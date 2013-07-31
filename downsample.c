#include <math.h>
#include <string.h>
#include "psrfits.h"

// TODO:  for these to work with OpenMP, we probably need
//        separate input and output arrays and then a copy.
//        Otherwise, the threads will step on each other.

void convert_4bit_to_8bit(unsigned char *indata, unsigned char *outdata, int N)
// This converts 4-bit indata to 8-bit outdata
// N is the total number of data points
{
    int ii;
    unsigned char uctmp;

    // Convert all the data from 4-bit to 8-bit
    for (ii = 0 ; ii < N / 2 ; ii++, indata++) {
        uctmp = *indata;
        *outdata++ = uctmp >> 4;   // 1st 4 bits (MSBs) are first nibble
        *outdata++ = uctmp & 0x0F; // 2nd 4 bits (LSBs) are second nibble
    }
}


void pf_4bit_to_8bit(struct psrfits *pf)
// This converts 4-bit pf->sub.rawdata to 8-bit pf->sub.data
{
    convert_4bit_to_8bit((unsigned char *)pf->sub.rawdata,
                         (unsigned char *)pf->sub.data,
                         pf->sub.bytes_per_subint * 2);
}


void convert_8bit_to_4bit(unsigned char *indata, unsigned char *outdata, int N)
// This converts 8-bit indata to 4-bit outdata
// N is the total number of data points
{
    int ii;

    // Convert all the data from 4-bit to 8-bit
    for (ii = 0 ; ii < N / 2 ; ii++, outdata++) {
        *outdata = *indata++ << 4;  // 1st 4 bits (MSBs) are first point
        *outdata += *indata++;      // 2nd 4 bits (LSBs) are second point
    }
}


void pf_8bit_to_4bit(struct psrfits *pf)
// This converts 8-bit pf->sub.data into 4-bit pf->sub.rawdata
{
    long long numoutsamp = pf->sub.bytes_per_subint * 2 / \
        (pf->hdr.ds_time_fact * pf->hdr.ds_freq_fact);
    convert_8bit_to_4bit((unsigned char *)pf->sub.data,
                         (unsigned char *)pf->sub.rawdata,
                         numoutsamp);
}


void get_stokes_I(struct psrfits *pf)
/* Move the Stokes I in place so that it is consecutive in the array */
{
    int ii;
    float *data;
    struct hdrinfo *hdr = &(pf->hdr);
    const int out_nchan = hdr->nchan / hdr->ds_freq_fact;

    // In this mode, average the polns first to make it like IQUV
    if (strncmp(hdr->poln_order, "AABBCRCI", 8)==0) {
        float *bbptr;
        int jj;
        for (ii = 0 ; ii < hdr->nsblk ; ii++) {
            data = pf->sub.fdata + ii * out_nchan * 4; // 4 polns
            bbptr = data + out_nchan;
            for (jj = 0 ; jj < out_nchan ; jj++, data++, bbptr++)
                *data = 0.5 * (*data + *bbptr); // Average AA and BB polns
        }
    }
    data = pf->sub.fdata;
    // Start from 1 since we don't need to move the 1st spectra
    for (ii = 1 ; ii < hdr->nsblk ; ii++) {
        memcpy(data + ii * out_nchan, 
               data + ii * 4 * out_nchan, 
               out_nchan * sizeof(float));
    }
}


void downsample_time(struct psrfits *pf)
/* Average adjacent time samples together in place */
/* This should be called _after_ make_subbands()   */
{
    int ii, jj, kk;
    struct hdrinfo *hdr = &(pf->hdr);
    float *data = pf->sub.fdata;
    float *indata, *outdata, *tmpspec;
    const int dsfact = hdr->ds_time_fact;
    // Treat the polns as being parts of the same spectrum
    int out_npol = hdr->npol;
    if (hdr->onlyI) out_npol = 1;
    const int in_nchan = hdr->nchan * out_npol;
    const int out_nchan = in_nchan / hdr->ds_freq_fact;
    const int out_nsblk = hdr->nsblk / dsfact;
    const float norm = 1.0 / dsfact;

    tmpspec = (float *)malloc(out_nchan * sizeof(float));
    indata = data;
    // Iterate over the output times
    for (ii = 0 ; ii < out_nsblk ; ii++) {
        // Initiaize the summation
        for (jj = 0 ; jj < out_nchan ; jj++) 
            tmpspec[jj] = 0.0;
        // Add up the samples in time in the tmp array
        for (jj = 0 ; jj < dsfact ; jj++) {
            outdata = tmpspec;
            for (kk = 0 ; kk < out_nchan ; kk++, indata++, outdata++)
                *outdata += *indata;
        }
        // Convert the sum to an average and put into the output array
        outdata = data + ii * out_nchan;
        for (jj = 0 ; jj < out_nchan ; jj++)
            outdata[jj] = tmpspec[jj] * norm;
    }
    free(tmpspec);
}


void guppi_update_ds_params(struct psrfits *pf)
/* Update the various output data arrays / values so that */
/* they are correct for the downsampled data.             */
{
    struct hdrinfo *hdr = &(pf->hdr);
    struct subint  *sub = &(pf->sub);

    int out_npol = hdr->npol;
    if (hdr->onlyI) out_npol = 1;
    int out_nchan = hdr->nchan / hdr->ds_freq_fact;
 
    if (hdr->ds_freq_fact > 1) {
        int ii;
        double dtmp;

        /* Note:  we don't need to malloc the subint arrays since */
        /*        their original values are longer by default.    */

        // The following correctly accounts for the middle-of-bin FFT offset
        dtmp = hdr->fctr - 0.5 * hdr->BW;
        dtmp += 0.5 * (hdr->ds_freq_fact - 1.0) * hdr->df;
        for (ii = 0 ; ii < out_nchan ; ii++)
            sub->dat_freqs[ii] = dtmp + ii * (hdr->df * hdr->ds_freq_fact);

        for (ii = 1 ; ii < out_npol ; ii++) {
            memcpy(sub->dat_offsets+ii*out_nchan,
                   sub->dat_offsets+ii*hdr->nchan,
                   sizeof(float)*out_nchan);
            memcpy(sub->dat_scales+ii*out_nchan,
                   sub->dat_scales+ii*hdr->nchan,
                   sizeof(float)*out_nchan);
        }
    }
}
