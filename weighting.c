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
    *outscl = sqrt(*outscl / *outwgt);
    *outoff /= *outwgt;
    *outwgt /= N;  // This is normal average of weights
}


void make_weighted_datamods(struct psrfits *inpf, struct psrfits *outpf)
{
    int ii, jj;
    struct hdrinfo *inhdr = &(inpf->hdr);
    const int N = inhdr->ds_freq_fact;
    const int out_nchan = inhdr->nchan / N;
    float *inwgts = inpf->sub.dat_weights;
    float *inscls = inpf->sub.dat_scales;
    float *inoffs = inpf->sub.dat_offsets;
    float *outwgts = outpf->sub.dat_weights;
    float *outscls = outpf->sub.dat_scales;
    float *outoffs = outpf->sub.dat_offsets;

    for (ii = 0, jj = 0 ; ii < out_nchan ; ii++, jj += N) {
        combine_datamods(N, inwgts+jj, inscls+jj, inoffs+jj,
                         outwgts+ii, outscls+ii, outoffs+ii);
    }
}

