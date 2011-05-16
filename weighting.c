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
                      float *outwgt, float *outscl, float *outoff);
{
    int ii;
    float inwgt;
    
    *outwgt = *outscl = *outoff = 0.0;
    for (ii = 0 ; ii < N ; ii++) {
        inwgt = inwgts[ii];
        *outwgt += inwgt;
        *outscl += inwgt * inscls[ii] * inscls[ii];  // adding in quadrature
        *outoff += inwgt * inoffs[ii];
    }
    *outscl = sqrt(*outscl / *outwgt);
    *outoff /= *outwgt;
    *outwgt /= N;  // This is normal average of weights
}


void make_weighted_datamods(struct psrfits *inpf, struct psrfits *outpf)
{
    int ii, jj;
    struct hdrinfo *inhdr = &(pf->inhdr);
    const int N = inhdr->ds_freq_fact;
    const int out_nchan = inhdr->nchan / N;
    const float *inwgts = inpf->sub.dat_weights;
    const float *inscls = inpf->sub.dat_scales;
    const float *inoffs = inpf->sub.dat_offsets;
    const float *outwgts = outpf->sub.dat_weights;
    const float *outscls = outpf->sub.dat_scales;
    const float *outoffs = outpf->sub.dat_offsets;

    for (ii = 0, jj = 0 ; ii < out_nchan ; ii++, jj += N) {
        combine_datamods(N, inwgts+jj, inscls+jj, inoffs+jj,
                         outwgts+ii, outscls+ii, outoffs+ii);
    }
}

