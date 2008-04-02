#ifndef _FOLD_H
#define _FOLD_H
#include <fftw3.h>
int fold_8bit_power(struct polyco *pc, int npc, int imjd, double fmjd,
        char *data, int nsamp, int nchan, int npol, 
        double tsamp, double *foldbuf, long long *cntbuf, int nbin);
int normalize_folds(double *foldbuf, long long *cntbuf, 
        int nbin, int nchan, int npol);
#endif
