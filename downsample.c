#include <math.h>
#include <string.h>
#include "psrfits.h"

// TODO:  for these to work with OpenMP, we probably need
//        separate input and output arrays and then a copy.
//        Otherwise, the threads will step on each other.


void pack_8bit_to_2bit_unsigned(unsigned char *indata,
                                unsigned char *outdata, int N)
// packs (i.e. converts) 8-bit unsigned indata to 2-bit outdata
// N is total number of data points
{
    int ii;
    for (ii = 0; ii < N / 4; ii++, outdata++) {
        *outdata = *indata++ << 6;
        *outdata |= *indata++ << 4;
        *outdata |= *indata++ << 2;
        *outdata |= *indata++;
    }
}

void pack_8bit_to_2bit_signed(char *indata,
                              char *outdata, int N)
// packs (i.e. converts) 8-bit signed indata to 2-bit outdata
// N is total number of data points
{
    int ii;
    for (ii = 0; ii < N / 4; ii++, outdata++) {
        *outdata = (*indata++ & 0x03) << 6;
        *outdata |= (*indata++ & 0x03) << 4;
        *outdata |= (*indata++ & 0x03) << 2;
        *outdata |= *indata++ & 0x03;
    }
}

void pf_pack_8bit_to_2bit(struct psrfits *pf, int numunsigned)
// packs (i.e. converts) 8-bit indata to 2-bit outdata in psrfits struct
{
    int ii, poln;
    int nspec = pf->hdr.nsblk;
    int npol = pf->hdr.npol;
    int nchan = pf->hdr.nchan;
    for (ii = 0 ; ii < nspec ; ii++) {
        for (poln = 0 ; poln < npol ; poln++) {
            if (poln < numunsigned) { // unsigned
                unsigned char *indata = pf->sub.data +  \
                    nchan * (ii * npol + poln);
                unsigned char *outdata = pf->sub.rawdata + \
                    nchan * (ii * npol + poln) / 4;
                pack_8bit_to_2bit_unsigned(indata, outdata, nchan);
            } else { // signed
                char *indata = (char *) (pf->sub.data + \
                                         nchan * (ii * npol + poln));
                char *outdata = (char *) (pf->sub.rawdata + \
                                          nchan * (ii * npol + poln) / 4);
                pack_8bit_to_2bit_signed(indata, outdata, nchan);
            }
        }
    }
}

void pack_8bit_to_4bit_unsigned(unsigned char *indata,
                                unsigned char *outdata, int N)
// packs (i.e. converts) 8-bit unsigned indata to 4-bit outdata
// N is total number of data points
{
    int ii;
    for (ii = 0 ; ii < N / 2 ; ii++, outdata++) {
        *outdata = *indata++ << 4;
        *outdata |= *indata++;
    }
}

void pack_8bit_to_4bit_signed(char *indata,
                              char *outdata, int N)
// packs (i.e. converts) 8-bit signed indata to 4-bit outdata
// N is total number of data points
{
    int ii;
    for (ii = 0 ; ii < N / 2 ; ii++, outdata++) {
        *outdata = (*indata++ & 0x0F) << 4;
        *outdata |= *indata++ & 0x0F;
    }
}

void pf_pack_8bit_to_4bit(struct psrfits *pf, int numunsigned)
// packs (i.e. converts) 8-bit indata to 4-bit outdata in psrfits struct
{
    int ii, poln;
    int nspec = pf->hdr.nsblk;
    int npol = pf->hdr.npol;
    int nchan = pf->hdr.nchan;
    for (ii = 0 ; ii < nspec ; ii++) {
        for (poln = 0 ; poln < npol ; poln++) {
            if (poln < numunsigned) { // unsigned
                unsigned char *indata = pf->sub.data + \
                    nchan * (ii * npol + poln);
                unsigned char *outdata = pf->sub.rawdata + \
                    nchan * (ii * npol + poln) / 2;
                pack_8bit_to_4bit_unsigned(indata, outdata, nchan);
            } else { // signed
                char *indata = (char *) (pf->sub.data + \
                                         nchan * (ii * npol + poln));
                char *outdata = (char *) (pf->sub.rawdata + \
                                          nchan * (ii * npol + poln) / 2);
                pack_8bit_to_4bit_signed(indata, outdata, nchan);
            }
        }
    }
}

void unpack_2bit_to_8bit_unsigned(unsigned char *indata,
                                  unsigned char *outdata, int N)
// unpacks (i.e. converts) 2-bit unsigned indata to 8-bit outdata
// N is total number of data points
{
    int ii;
    unsigned char uctmp;
    for (ii = 0 ; ii < N / 4 ; ii++, indata++) {
        uctmp = *indata;
        *outdata++ = uctmp >> 6;
        *outdata++ = (uctmp >> 4) & 0x03;
        *outdata++ = (uctmp >> 2) & 0x03;
        *outdata++ = uctmp & 0x03;
    }
}

void unpack_2bit_to_8bit_signed(unsigned char *indata,
                                unsigned char *outdata, int N)
// unpacks (i.e. converts) 2-bit signed indata to 8-bit outdata
// N is total number of data points
{
    int ii;
    // This provides automatic sign extension (via a bitfield)
    // which is essential for twos complement signed numbers
    // https://graphics.stanford.edu/~seander/bithacks.html#FixedSignExtend
    struct {signed char x:2;} stmp;
    for (ii = 0 ; ii < N / 4 ; ii++, indata++) {
        stmp.x = *indata >> 6;
        *outdata++ = stmp.x;
        stmp.x = ((*indata >> 4) & 0x03);
        *outdata++ = stmp.x;
        stmp.x = ((*indata >> 2) & 0x03);
        *outdata++ = stmp.x;
        stmp.x = (*indata & 0x03);
        *outdata++ = stmp.x;
    }
}

void pf_unpack_2bit_to_8bit(struct psrfits *pf, int numunsigned)
// unpacks (i.e. converts) 2-bit indata to 8-bit outdata in psrfits struct
{
    int ii, poln;
    int nspec = pf->hdr.nsblk;
    int npol = pf->hdr.npol;
    int nchan = pf->hdr.nchan;
    for (ii = 0 ; ii < nspec ; ii++) {
        for (poln = 0 ; poln < npol ; poln++) {
            if (poln < numunsigned) { // unsigned
                unsigned char *indata = pf->sub.rawdata + \
                    ii * nspec * npol / 4 + poln * nchan / 4;
                unsigned char *outdata = pf->sub.data + ii * nspec * npol + \
                    poln * nchan;
                unpack_2bit_to_8bit_unsigned(indata, outdata, nchan);
            } else { // signed
                char *indata = (char *) (pf->sub.rawdata + \
                                         ii * nspec * npol / 4 + \
                                         poln * nchan / 4);
                char *outdata = (char *) (pf->sub.data + \
                                          ii * nspec * npol + \
                                          poln * nchan);
                unpack_2bit_to_8bit_signed(indata, outdata, nchan);
            }
        }
    }
}

void unpack_4bit_to_8bit_unsigned(unsigned char *indata,
                                  unsigned char *outdata, int N)
// unpacks (i.e. converts) 4-bit unsigned indata to 8-bit outdata
// N is total number of data points
{
    int ii;
    unsigned char uctmp;
    for (ii = 0 ; ii < N / 2 ; ii++, indata++) {
        uctmp = *indata;
        *outdata++ = uctmp >> 4;
        *outdata++ = uctmp & 0x0F;
    }
}

void unpack_4bit_to_8bit_signed(unsigned char *indata,
                                unsigned char *outdata, int N)
// unpacks (i.e. converts) 4-bit signed indata to 8-bit outdata
// N is total number of data points
{
    int ii;
    // This provides automatic sign extension (via a bitfield)
    // which is essential for twos complement signed numbers
    // https://graphics.stanford.edu/~seander/bithacks.html#FixedSignExtend
    struct {signed char x:4;} stmp;
    for (ii = 0 ; ii < N / 2 ; ii++, indata++) {
        stmp.x = *indata >> 4;
        *outdata++ = stmp.x;
        stmp.x = (*indata & 0x0F);
        *outdata++ = stmp.x;
    }
}

void pf_unpack_4bit_to_8bit(struct psrfits *pf, int numunsigned)
// unpacks (i.e. converts) 4-bit indata to 8-bit outdata in psrfits struct
{
    int ii, poln;
    int nspec = pf->hdr.nsblk;
    int npol = pf->hdr.npol;
    int nchan = pf->hdr.nchan;
    for (ii = 0 ; ii < nspec ; ii++) {
        for (poln = 0 ; poln < npol ; poln++) {
            if (poln < numunsigned) { // unsigned
                unsigned char *indata = pf->sub.rawdata + \
                    ii * nspec * npol / 2 + poln * nchan / 2;
                unsigned char *outdata = pf->sub.data + ii * nspec * npol + \
                    poln * nchan;
                unpack_4bit_to_8bit_unsigned(indata, outdata, nchan);
            } else { // signed
                char *indata = (char *) (pf->sub.rawdata + \
                                         ii * nspec * npol / 2 + \
                                         poln * nchan / 2);
                char *outdata = (char *) (pf->sub.data + \
                                          ii * nspec * npol + \
                                          poln * nchan);
                unpack_4bit_to_8bit_signed(indata, outdata, nchan);
            }
        }
    }
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
