#include <math.h>
#include <string.h>
#include "psrfits.h"

// TODO:  for these to work with OpenMP, we probably need
//        separate input and output arrays and then a copy.
//        Otherwise, the threads will step on each other.

void get_stokes_I(struct psrfits *pf)
/* Move the Stokes I in place so that it is consecutive in the array */
{
    int ii, inbytes, outbytes;
    struct hdrinfo *hdr = &(pf->hdr); // dereference the ptr to the header struct
    unsigned char *data = pf->sub.data;

    outbytes = hdr->nbits * hdr->nchan / 8;
    inbytes = outbytes * 4;  // 4 Stokes params
    // Start from 1 since we don't need to move the 1st spectra
    for (ii = 1 ; ii < hdr->nsblk ; ii++)
        memcpy(data + ii * outbytes, data + ii * inbytes, outbytes);
}

void downsample_freq(struct psrfits *pf)
/* Average adjacent frequency channels together in place    */
/* Note:  this only works properly for 8-bit data currently */
{
    int ii, jj, kk, itmp, iidx=0, oidx=0;
    struct hdrinfo *hdr = &(pf->hdr); // dereference the ptr to the header struct
    unsigned char *data = pf->sub.data;
    float norm = 1.0 / hdr->ds_freq_fact;

    // Treat the polns as being parts of the same spectrum
    int out_npol = hdr->npol;
    if (hdr->onlyI) out_npol = 1;
    int out_nchan = hdr->nchan * out_npol / hdr->ds_freq_fact;
    
    for (ii = 0 ; ii < hdr->nsblk ; ii++) { // Iterate over the time-slices
        for (jj = 0 ; jj < out_nchan ; jj++) { // ... output chans
            itmp = 0;
            for (kk = 0 ; kk < hdr->ds_freq_fact ; kk++) // ... adjacent chans
                itmp += data[oidx++]; 
            data[iidx++] = (int) rintf(((float) itmp) * norm);
        }
    }
}

void downsample_time(struct psrfits *pf)
/* Average adjacent time samples together in place */
/* This should be called _after_ downsample_freq() */
/* Note:  this only works properly for 8-bit data currently */
{
    int ii, jj, kk, itmp, iidx, oidx;
    struct hdrinfo *hdr = &(pf->hdr); // dereference the ptr to the header struct
    unsigned char *data = pf->sub.data;
    int dsfact = hdr->ds_time_fact;
    float norm = 1.0 / dsfact;

    // Treat the polns as being parts of the same spectrum
    int out_npol = hdr->npol;
    if (hdr->onlyI) out_npol = 1;
    int out_nchan = hdr->nchan * out_npol / hdr->ds_freq_fact;
    int out_nsblk = hdr->nsblk / dsfact;
    
    for (ii = 0 ; ii < out_nsblk ; ii++) { // Iterate over the output time-slices
        iidx = ii * out_nchan;
        for (jj = 0 ; jj < out_nchan ; jj++) { // ... channels
            itmp = 0;
            oidx = ii * dsfact * out_nchan + jj;
            for (kk = 0 ; kk < dsfact ; kk++) { // ... adjacent times
                itmp += data[oidx];
                oidx += out_nchan;
            }
            data[iidx+jj] = (int) rintf((float) itmp * norm);
        }
    }
}
