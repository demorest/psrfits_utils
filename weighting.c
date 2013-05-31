#include <math.h>
#include <string.h>
#include "psrfits.h"

float vector_sum(float *vector, int N)
{
    int ii;
    float result = 0.0;

    for (ii = 0 ; ii < N ; ii++)
        result += vector[ii];
    return result;
}


void combine_datamods(int N, float *inwgts, float *inscls, float *inoffs, 
                      float *outwgt, float *outscl, float *outoff)
{
    float inwgt, sumwgt = 0.0;
    int ii;

    *outwgt = *outscl = *outoff = 0.0;
    for (ii = 0 ; ii < N ; ii++) {
        inwgt = inwgts[ii];
        // weights are averaged
        sumwgt += inwgt;
        // scales and offsets are normal weighted average
        *outscl += inwgt * inscls[ii];
        *outoff += inwgt * inoffs[ii];
    }
    // avoid divide by zeros if weights are all zero
    *outscl = (sumwgt) ? *outscl / sumwgt : 1.0;
    *outoff = (sumwgt) ? *outoff / sumwgt : 0.0;
    // Since the stdev of the output data will be up-scaled
    // by sqrt(sumwgt), we will compensate for that using the
    // output weight.
    *outwgt = (sumwgt) ? 1.0 / sqrt(sumwgt): 0.0;
}


void make_weighted_datamods(struct psrfits *inpf, struct psrfits *outpf)
{
    int ii, jj, poln;
    const int out_nchan = inpf->hdr.nchan / outpf->hdr.ds_freq_fact;
    const int N = outpf->hdr.ds_freq_fact;
    float *inwgts, *inscls, *inoffs, *outwgts, *outscls, *outoffs;

    // Need to do this for each poln.  There are nchan * npol scales 
    // and weights, but only nchan weights.
    for (poln = 0 ; poln < inpf->hdr.npol ; poln++) {
        inwgts = inpf->sub.dat_weights;
        inscls = inpf->sub.dat_scales + poln * inpf->hdr.nchan;
        inoffs = inpf->sub.dat_offsets + poln * inpf->hdr.nchan;
        outwgts = outpf->sub.dat_weights;
        outscls = outpf->sub.dat_scales + poln * out_nchan;
        outoffs = outpf->sub.dat_offsets + poln * out_nchan;
        for (ii = 0, jj = 0 ; ii < out_nchan ; ii++, jj += N)
            combine_datamods(N, inwgts+jj, inscls+jj, inoffs+jj,
                             outwgts+ii, outscls+ii, outoffs+ii);
        //printf("%4d  %8.5f %8.5f %8.5f \n", ii, outpf->sub.dat_weights[ii], outpf->sub.dat_scales[ii], outpf->sub.dat_offsets[ii]);
        //printf("%4d  %8.5f %8.5f %8.5f \n", ii, *(outwgts+ii), *(outscls+ii), *(outoffs+ii));
    }
}

