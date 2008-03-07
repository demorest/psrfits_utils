/* psrfits.h */
#ifndef _PSRFITS_H_
#define _PSRFITS_H_
#include "fitsio.h"
#include "../polyco.h"
int psrfits_write_header(fitsfile *fptr, struct asp_params p, 
        int *status);
int psrfits_write_start_time(fitsfile *fptr, int imjd, double fmjd, 
        int *status);
int psrfits_create_polyco_tbl(fitsfile *fptr, int *status);
int psrfits_create_history_tbl(fitsfile *fptr, int *status);
int psrfits_create_subint_tbl(fitsfile *fptr, int nchan, int nbin, 
        int *status);
int psrfits_write_polyco(fitsfile *fptr, struct Polyco *pc, 
        struct asp_params p, int *status);
int psrfits_write_history(fitsfile *fptr, int nsub, int nbin, 
        double rf, int nchan, double ch_bw, int *status);
int psrfits_write_subint(fitsfile *fptr, int nchan, int nbin, double tint,
        double time, double per,
        float *rfs, float *scale, float *offset, float *data, 
        int *status);
#endif
