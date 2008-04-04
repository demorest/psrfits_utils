/* Simple fold routines */
#include <math.h>
#include <string.h>
#include "fold.h"
#include "polyco.h"

void malloc_foldbuf(struct foldbuf *f) {
    // XXX align!
    f->data = (float *)malloc(sizeof(float) * f->nbin * f->npol * f->nchan);
    f->count = (unsigned *)malloc(sizeof(unsigned) * f->nbin);
}

void free_foldbuf(struct foldbuf *f) {
    if (f->data!=NULL) free(f->data);
    if (f->count!=NULL) free(f->count);
}

void clear_foldbuf(struct foldbuf *f) {
    memset(f->data, 0, sizeof(float) * f->nbin * f->npol * f->nchan);
    memset(f->count, 0, sizeof(float) * f->nbin);
}

static void vector_accumulate(float *out, float *in, int n) {
    int i;
    for (i=0; i<n; i++) { out[i] += in[i]; }
}

static void unpack_8bit(float *out, char *in, int n) {
    int i;
    for (i=0; i<n; i++) { out[i] = (float)in[i]; }
}

int fold_8bit_power(struct polyco *pc, int imjd, double fmjd, 
        char *data, int nsamp, double tsamp, struct foldbuf *f) {

    /* Find midtime */
    double fmjd_mid = fmjd + nsamp*tsamp/2.0/86400.0;

    /* Check polyco set */
    if (pc_out_of_range(pc, imjd, fmjd)) { return(-1); }

    /* Calc phase, phase step */
    double dphase=0.0;
    double phase = psr_phase(pc, imjd, fmjd, NULL);
    phase = fmod(phase, 1.0);
    if (phase<0.0) { phase += 1.0; }
    psr_phase(pc, imjd, fmjd_mid, &dphase);
    dphase *= tsamp;

    /* Fold em */
    int i, ibin;
    float *fptr, *dptr;
    dptr = (float *)malloc(sizeof(float)*f->nchan*f->npol); // XXX align!
    for (i=0; i<nsamp; i++) {
        ibin = (int)(phase * (double)f->nbin);
        if (ibin<0) { ibin+=f->nbin; }
        if (ibin>=f->nbin) { ibin-=f->nbin; }
        fptr = &f->data[ibin*f->nchan*f->npol];
        unpack_8bit(dptr, &data[i*f->nchan*f->npol], f->nchan * f->npol);
        vector_accumulate(fptr, dptr, f->npol * f->nchan);
        f->count[ibin]++;
        phase += dphase;
        if (phase>1.0) { phase -= 1.0; }
    }
    free(dptr);

    return(0);
}

int accumulate_folds(struct foldbuf *ftot, struct foldbuf *f) {
    if (ftot->nbin!=f->nbin || ftot->nchan!=f->nchan || ftot->npol!=f->npol) {
        return(-1);
    }
    int i;
    for (i=0; i<f->nbin; i++) { ftot->count[i] += f->count[i]; }
    vector_accumulate(ftot->data, f->data, f->nbin * f->nchan * f->npol);
    return(0);
}

/* normalize and transpose to psrfits order */
int normalize_transpose_folds(float *out, struct foldbuf *f) {
    int ibin, ii;
    for (ibin=0; ibin<f->nbin; ibin++) {
        for (ii=0; ii<f->nchan*f->npol; ii++) {
            out[ibin + ii*f->nchan*f->npol] =
                f->data[ii + ibin*f->nchan*f->npol] / (float)f->count[ibin];
        }
    }
    return(0);
}
