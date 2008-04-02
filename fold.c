/* Simple fold routines */
#include <math.h>
#include "polyco.h"

int fold_8bit_power(struct polyco *pc, int npc, int imjd, double fmjd,
        char *data, int nsamp, int nchan, int npol, 
        double tsamp, double *foldbuf, long long *cntbuf, int nbin) {

    /* Find midtime */
    double fmjd_mid = fmjd + nsamp*tsamp/2.0/86400.0;

    /* Find correct pc set */
    int ipc;
    for (ipc=0; ipc<npc; ipc++) {
        if (pc_out_of_range(&pc[ipc], imjd, fmjd_mid)==0) { break; }
    }
    if (ipc==npc) { return(-1); }

    /* Calc phase, phase step */
    double dphase=0.0;
    double phase = psr_phase(&pc[ipc], imjd, fmjd, NULL);
    phase = fmod(phase, 1.0);
    if (phase<0.0) { phase += 1.0; }
    psr_phase(&pc[ipc], imjd, fmjd_mid, &dphase);
    dphase *= tsamp;

    /* Fold em */
    int i, ipol, ichan, ibin;
    char *dptr;
    double *fptr;
    for (i=0; i<nsamp; i++) {
        ibin = (int)(phase * (double)nbin);
        if (ibin<0) { ibin+=nbin; }
        if (ibin>=nbin) { ibin-=nbin; }
        dptr = &data[i*nchan*npol];
        fptr = &foldbuf[ibin];
        for (ipol=0; ipol<npol; ipol++) {
            for (ichan=0; ichan<nchan; ichan++) {
                *fptr += (double)(*dptr);
                fptr+=nbin;
                dptr++;
            }
        }
        cntbuf[ibin]++;
        phase += dphase;
        if (phase>1.0) { phase -= 1.0; }
    }

    return(0);
}

int normalize_folds(double *foldbuf, long long *cntbuf, 
        int nbin, int nchan, int npol) {
    int ibin, ichan, ipol;
    for (ipol=0; ipol<npol; ipol++) {
        for (ichan=0; ichan<nchan; ichan++) {
            for (ibin=0; ibin<nbin; ibin++) {
                foldbuf[ibin+ichan*nbin+ipol*nbin*nchan] /=
                    (double)cntbuf[ibin];
            }
        }
    }
    return(0);
}
