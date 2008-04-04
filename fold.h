#ifndef _FOLD_H
#define _FOLD_H
#include "polyco.h"

struct foldbuf {
    int nbin;
    int nchan;
    int npol;
    float *data;
    unsigned *count;
};

void malloc_foldbuf(struct foldbuf *f);

void free_foldbuf(struct foldbuf *f);

void clear_foldbuf(struct foldbuf *f);

int normalize_transpose_folds(float *out, struct foldbuf *f);

int fold_8bit_power(struct polyco *pc, int imjd, double fmjd, 
        char *data, int nsamp, double tsamp, struct foldbuf *f);

int accumulate_folds(struct foldbuf *ftot, struct foldbuf *f);

#endif
