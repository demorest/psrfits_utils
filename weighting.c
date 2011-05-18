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
    int ii;

    *outwgt = *outscl = *outoff = 0.0;
    for (ii = 0 ; ii < N ; ii++) {
        float inwgt = inwgts[ii];
        // weights are averaged
        *outwgt += inwgt;
        // scales are combined in quadrature
        *outscl += inwgt * inscls[ii] * inscls[ii];
        // offsets are weighted average
        *outoff += inwgt * inoffs[ii];
    }
    // avoid divide by zeros if weights are all zero
    *outscl = (*outwgt) ? sqrt(*outscl / *outwgt) : 1.0;
    *outoff = (*outwgt) ? *outoff / *outwgt : 0.0;
    *outwgt /= N;  // This is normal average of weights
}


void make_weighted_datamods(struct psrfits *inpf, struct psrfits *outpf)
{
    int ii, jj;
    const int out_nchan = inpf->hdr.nchan / outpf->hdr.ds_freq_fact;
    const int N = outpf->hdr.ds_freq_fact;
    float *inwgts = inpf->sub.dat_weights;
    float *inscls = inpf->sub.dat_scales;
    float *inoffs = inpf->sub.dat_offsets;
    float *outwgts = outpf->sub.dat_weights;
    float *outscls = outpf->sub.dat_scales;
    float *outoffs = outpf->sub.dat_offsets;

    for (ii = 0, jj = 0 ; ii < out_nchan ; ii++, jj += N) {
        combine_datamods(N, inwgts+jj, inscls+jj, inoffs+jj,
                         outwgts+ii, outscls+ii, outoffs+ii);
        // printf("%4d  %8.5f %8.5f %8.5f \n", ii, outpf->sub.dat_weights[ii], outpf->sub.dat_scales[ii], outpf->sub.dat_offsets[ii]);
    }
}

