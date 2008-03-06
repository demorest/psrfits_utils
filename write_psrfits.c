/* psrfits.c
 *
 * Routines for writing psrfits files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "asp_params.h"
#include "psrfits.h"

/* Stupid helper routines */
void dec2hms(char *out, double in, int sflag) {
    int sign=1;
    char *ptr=out;
    int h, m;
    double s;
    if (in<0.0) { sign=-1; in=fabs(in); }
    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    s = in;
    if (sign==1 && sflag) { *ptr='+'; ptr++; }
    else if (sign==-1) { *ptr='-'; ptr++; }
    sprintf(ptr, "%2.2d:%2.2d:%06.3f", h, m, s);
}

/* Write the primary HDU keywords */
int psrfits_write_header(fitsfile *fptr, struct asp_params p, 
        int *status) {

    /* if fitsio status is non-zero, pass it thru */
    if (*status!=0) { return(*status); }

    /* Temp vars */
    int inttmp;
    float flttmp;
    double dbltmp;
    char strtmp[256];

    /* Point to primary HDU (should be already created) */
    fits_movabs_hdu(fptr, 1, NULL, status);

    /* Write em! */
    fits_write_key(fptr, TSTRING, "HDRVER", "2.12", NULL, status);
    fits_write_date(fptr, status); /* writes DATE keyword */
    fits_write_key(fptr, TSTRING, "OBSERVER", p.observer, NULL, status);
    fits_write_key(fptr, TSTRING, "PROJID", p.proj_id, NULL, status);
    fits_write_key(fptr, TSTRING, "TELESCOP", p.telescope, NULL, status);
    /* Skipped keywords: ANT_[XYZ] telescope coords */
    fits_write_key(fptr, TSTRING, "FRONTEND", p.front_end, NULL, status);
    switch (p.pol_mode[0]) {
        case 'l':
        case 'L':
            fits_write_key(fptr, TSTRING, "FD_POLN", "LIN", NULL, status);
            break;
        case 'c':
        case 'C':
            fits_write_key(fptr, TSTRING, "FD_POLN", "CIRC", NULL, status);
            break;
        default:
            fits_write_key(fptr, TSTRING, "FD_POLN", "UNKNOWN", NULL, status);
            break;
    }
    inttmp=1; fits_write_key(fptr, TINT, "FD_HAND", &inttmp,
            "Possibly not accurate", status);
    flttmp=0.0; fits_write_key(fptr, TFLOAT, "FD_SANG", &flttmp, 
            "[deg] Possibly not accurate", status);
    flttmp=0.0; fits_write_key(fptr, TFLOAT, "FD_XYPH", &flttmp,
            "[deg] Possibly not accurate", status);
    switch (p.telescope[0]) {
        case '1':
        case 'b':
            fits_write_key(fptr, TSTRING, "BACKEND", "GASP", NULL, status);
            break;
        case '3':
        default:
            fits_write_key(fptr, TSTRING, "BACKEND", "ASP", NULL, status);
            break;
    }
    fits_write_key(fptr, TSTRING, "BECONFIG", "config.xml", NULL, status);
    inttmp=0; fits_write_key(fptr, TINT, "BE_PHASE", &inttmp, 
            "Don't know what this means", status);
    inttmp=0; fits_write_key(fptr, TINT, "BE_DCC", &inttmp, 
            "Don't know what this means", status);
    dbltmp=0.0; fits_write_key(fptr, TDOUBLE, "TCYCLE", &dbltmp,
            "[s] Don't know what this means", status);
    inttmp=1; fits_write_key(fptr, TINT, "NRCVR", &inttmp, 
            "Probably correct", status);
    if (p.cal_scan) {
        fits_write_key(fptr, TSTRING, "OBS_MODE", "CAL", NULL, status);
    } else {
        fits_write_key(fptr, TSTRING, "OBS_MODE", "PSR", NULL, status);
    }
    /* TODO: convert p.imjd, p.fmjd to date, write to file as DATE-OBS */
    switch (p.n_ds) {
        case 1:
            flttmp = p.rf + (float)p.band_dir*p.n_chan*p.ch_bw/2.0; 
            flttmp += (float)p.band_dir*p.ch_bw/2.0;
            break;
        case 2:
            flttmp = p.rf + (float)p.band_dir*p.n_chan*p.ch_bw; 
            flttmp += (float)p.band_dir*p.ch_bw/2.0;
            break;
        case 3:
            flttmp = p.rf + (float)p.band_dir*p.n_chan*p.ch_bw; 
            flttmp += (float)p.band_dir*p.ch_bw/2.0;
            break;
        case 4:
            flttmp = p.rf + (float)p.band_dir*p.n_chan*p.ch_bw*2.0; 
            flttmp += (float)p.band_dir*p.ch_bw/2.0;
            break;
        default:
            flttmp = p.rf;
            break;
    }
    fits_write_key(fptr, TFLOAT, "OBSFREQ", &flttmp, "[MHz]", status);
    fits_write_key(fptr, TSTRING, "SRC_NAME", p.psr_name, NULL, status);
    fits_write_key(fptr, TSTRING, "COORD_MD", "J2000", NULL, status);
    sprintf(strtmp, "%.1f", p.epoch);
    fits_write_key(fptr, TSTRING, "EQUINOX", strtmp,
            "Why is this a string?", status);
    /* Assumes we're tracking */
    dec2hms(strtmp, p.ra, 0);
    fits_write_key(fptr, TSTRING, "RA", strtmp, NULL, status);
    fits_write_key(fptr, TSTRING, "STT_CRD1", strtmp, NULL, status);
    fits_write_key(fptr, TSTRING, "STP_CRD1", strtmp, NULL, status);
    dec2hms(strtmp, p.dec, 1);
    fits_write_key(fptr, TSTRING, "DEC", strtmp, NULL, status);
    fits_write_key(fptr, TSTRING, "STT_CRD2", strtmp, NULL, status);
    fits_write_key(fptr, TSTRING, "STP_CRD2", strtmp, NULL, status);
    flttmp=0.0; fits_write_key(fptr, TFLOAT, "BMAJ", &flttmp, 
            "[deg] Need?", status);
    flttmp=0.0; fits_write_key(fptr, TFLOAT, "BMIN", &flttmp, 
            "[deg] Need?", status);
    flttmp=0.0; fits_write_key(fptr, TFLOAT, "BPA", &flttmp, 
            "[deg] Need?", status);
    fits_write_key(fptr, TSTRING, "TRK_MODE", "TRACK", NULL, status);
    flttmp = p.t_dump;
    fits_write_key(fptr, TFLOAT, "SCANLEN", &flttmp, 
            "[s] Integration time or total time?", status);
    fits_write_key(fptr, TSTRING, "FD_MODE", "FA",
            "Don't know what this means", status);
    fits_write_key(fptr, TFLOAT, "FA_REQ", &p.feed_ang,
            "[deg]", status);
    if (p.cal_scan) {
        fits_write_key(fptr, TSTRING, "CAL_MODE", "EXT1",
                "Don't know what this means", status);
        flttmp=25.0; fits_write_key(fptr, TFLOAT, "CAL_FREQ", &flttmp,
                "[Hz]", status);
        flttmp=0.5; fits_write_key(fptr, TFLOAT, "CAL_DCYC", &flttmp,
                NULL, status);
        flttmp=0.0; fits_write_key(fptr, TFLOAT, "CAL_PHS", &flttmp,
                "[deg] Possibly not accurate", status);
    } else {
        fits_write_key(fptr, TSTRING, "CAL_MODE", "OFF", NULL, status);
        /* Need all the others in this case? */
    }

    /* Do these dates need to be accurate start times? 
     * Yes they do.. so we don't write them here!
     */
#if 0
    fits_write_key(fptr, TLONG, "STT_IMJD", &p.imjd, NULL, status);
    inttmp = (int)(p.fmjd * 86400.0);
    fits_write_key(fptr, TINT, "STT_SMJD", &inttmp, 
            "[s] NOTE: Not exact!", status);
    dbltmp = p.fmjd*86400.0 - (double)inttmp;
    fits_write_key(fptr, TDOUBLE, "STT_OFFS", &dbltmp, 
            "[s] NOTE: Not exact!", status);
#endif

    /* Skipped keyword STT_LST */

    /* All done */
    return(*status);

}

int psrfits_write_start_time(fitsfile *fptr, int imjd, double fmjd, 
        int *status) {

    int day_sec = (int)(fmjd * 86400.0);
    double frac = fmjd*86400.0 - day_sec;

    /* On nonzero status, skip out */
    if (*status) { return(*status); }

    /* Point to primary HDU (should be already created) */
    fits_movabs_hdu(fptr, 1, NULL, status);

    fits_write_key(fptr, TLONG, "STT_IMJD", &imjd, NULL, status);
    fits_write_key(fptr, TINT, "STT_SMJD", &day_sec, NULL, status);
    fits_write_key(fptr, TDOUBLE, "STT_OFFS", &frac, NULL, status);
    return(*status);
}

#define TSIZE 16

#define COL_NOUNIT(type,form) do { \
    ttype[i] = (char *)malloc(sizeof(char)*TSIZE); \
    tform[i] = (char *)malloc(sizeof(char)*TSIZE); \
    tunit[i] = NULL; \
    strncpy(ttype[i], type, TSIZE-1); \
    strncpy(tform[i], form, TSIZE-1); \
    ttype[i][TSIZE-1] = '\0'; \
    tform[i][TSIZE-1] = '\0'; \
    i++; } while (0)

#define COL_UNIT(type,form,unit) do { \
    ttype[i] = (char *)malloc(sizeof(char)*TSIZE); \
    tform[i] = (char *)malloc(sizeof(char)*TSIZE); \
    tunit[i] = (char *)malloc(sizeof(char)*TSIZE); \
    strncpy(ttype[i], type, TSIZE-1); \
    strncpy(tform[i], form, TSIZE-1); \
    strncpy(tunit[i], unit, TSIZE-1); \
    ttype[i][TSIZE-1] = '\0'; \
    tform[i][TSIZE-1] = '\0'; \
    tunit[i][TSIZE-1] = '\0'; \
    i++; } while (0)

    
    
int psrfits_create_polyco_tbl(fitsfile *fptr, int *status) { 
    int i=0,ncol;
    char *ttype[32], *tform[32], *tunit[32];

    /* On nonzero status, skip out */
    if (*status) { return(*status); }

    /* Set up table names, etc */
    COL_NOUNIT("DATE_PRO","24A");
    COL_NOUNIT("POLYVER","16A");
    COL_UNIT("NSPAN","1I","min");
    COL_NOUNIT("NCOEF","1I");
    COL_NOUNIT("NPBLK","1I");
    COL_NOUNIT("NSITE","8A");
    COL_UNIT("REF_FREQ","1D","MHz");
    COL_NOUNIT("PRED_PHS","1D");
    COL_NOUNIT("REF_MJD","1D");
    COL_NOUNIT("REF_PHS","1D");
    COL_UNIT("REF_F0","1D","Hz");
    COL_NOUNIT("LGFITERR","1D");
    COL_NOUNIT("COEFF","15D");
    ncol=i;

    /* Create table */
    fits_create_tbl(fptr, BINARY_TBL, 0, ncol, ttype, tform, tunit, 
            "POLYCO", status);

    /* Unalloc */
    for (i=0; i<ncol; i++) {
        free(ttype[i]);
        free(tform[i]);
        if (tunit[i]!=NULL) { free(tunit[i]); }
    }

    /* All done */
    return(*status);
}

int psrfits_create_history_tbl(fitsfile *fptr, int *status) {
    int i=0, ncol;
    char *ttype[32], *tform[32], *tunit[32];

    /* On nonzero status, skip out */
    if (*status) { return(*status); }

    /* Column names */
    COL_NOUNIT("DATE_PRO","24A");
    COL_NOUNIT("PROC_CMD","80A");
    COL_NOUNIT("POL_TYPE","8A");
    COL_NOUNIT("NSUB","1I");
    COL_NOUNIT("NPOL","1I");
    COL_NOUNIT("NBIN","1I");
    COL_UNIT("CTR_FREQ","1D","MHz");
    COL_NOUNIT("NCHAN","1I");
    COL_UNIT("CHAN_BW","1D","MHz");
    ncol=i;

    /* Create table */
    fits_create_tbl(fptr, BINARY_TBL, 0, ncol, ttype, tform, tunit, 
            "HISTORY", status);

    /* Unalloc */
    for (i=0; i<ncol; i++) {
        free(ttype[i]);
        free(tform[i]);
        if (tunit[i]!=NULL) { free(tunit[i]); }
    }

    /* All done */
    return(*status);
}

int psrfits_create_subint_tbl(fitsfile *fptr, int nchan, int nbin, 
        int *status) {
    int i=0,ncol,tmp;
    char stmp[32];
    char *ttype[32], *tform[32], *tunit[32];

    /* On nonzero status, skip out */
    if (*status) { return(*status); }

    /* Set up column names, types */
    COL_NOUNIT("ISUBINT","1J");
    // Skipped : INDEXVAL
    COL_UNIT("TSUBINT","1D","s");
    COL_UNIT("OFFS_SUB","1D","s");
    COL_UNIT("PERIOD","1D","s");
    // Skipped LST_SUB
    // Skipped RA_SUB
    // Skipped DEC_SUB
    // Skipped GLON_SUB
    // Skipped GLAT_SUB
    // Skipped FD_ANG
    // Skipped POS_ANG
    // Skipped PAR_ANG
    // Skipped TEL_AZ
    // Skipped TEL_ZEN
    sprintf(stmp, "%dE", nchan); 
    COL_UNIT("DAT_FREQ",stmp,"MHz");
    COL_NOUNIT("DAT_WTS",stmp);
    sprintf(stmp, "%dE", nchan*4); 
    COL_NOUNIT("DAT_OFFS",stmp);
    COL_NOUNIT("DAT_SCL",stmp);
    sprintf(stmp, "%dE", nbin*nchan*4); 
    COL_NOUNIT("DATA",stmp);
    ncol=i;

    /* Create table */
    fits_create_tbl(fptr, BINARY_TBL, 0, ncol, ttype, tform, tunit, 
            "SUBINT", status);

    /* Fill in keywords */
    fits_write_key(fptr, TSTRING, "INT_TYPE", "TIME", NULL, status);
    fits_write_key(fptr, TSTRING, "INT_UNIT", "PHS", NULL, status);
    fits_write_key(fptr, TINT, "NBIN", &nbin, NULL, status);
    tmp=1; fits_write_key(fptr, TINT, "NBITS", &tmp, NULL, status);
    fits_write_key(fptr, TINT, "NCH_FILE", &nchan, NULL, status);
    tmp=0; fits_write_key(fptr, TINT, "NCH_STRT", &tmp, NULL, status);
    tmp=4; fits_write_key(fptr, TINT, "NPOL", &tmp, NULL, status);
    tmp=1; fits_write_key(fptr, TINT, "NSBLK", &tmp, NULL, status);
    fits_write_key (fptr, TSTRING, "EPOCHS", "MIDTIME", NULL, status);

    /* Write TDIM for DATA column.
     * Assumes DATA is the last column in table.
     */
    long naxes[4];
    naxes[0] = nbin;
    naxes[1] = nchan;
    naxes[2] = 4;
    naxes[3] = 1;
    fits_write_tdim(fptr, ncol, 4, naxes, status);

    /* Unalloc */
    for (i=0; i<ncol; i++) {
        free(ttype[i]);
        free(tform[i]);
        if (tunit[i]!=NULL) { free(tunit[i]); }
    }

    /* Done */
    return(*status);
}

int psrfits_write_polyco(fitsfile *fptr, struct Polyco *pc, 
        struct asp_params p, int *status) {

    if (*status) { return(*status); }

    /* Point to correct hdu */
    fits_movnam_hdu(fptr, BINARY_TBL, "POLYCO", 0, status);

    /* Get date string */
    int tmp;
    double dtmp;
    char date[32];
    fits_get_system_time(date, &tmp, status);

    /* find out how many rows (so we know where to write) */
    long row;
    fits_get_num_rows(fptr, &row, status);
    row++;

    /* Go through columns:
     *  - retrieve colnum using the col name
     *  - write value
     */
    char *ptr;
    int i, col;

    ptr = date;
    fits_get_colnum(fptr, CASEINSEN, "DATE_PRO", &col, status);
    fits_write_col(fptr, TSTRING, col, row, 1, 1, &ptr, status);

    // skipped POLYVER .. what's it mean?

    fits_get_colnum(fptr, CASEINSEN, "NSPAN", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &pc->NMinutes, status);

    fits_get_colnum(fptr, CASEINSEN, "NCOEF", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &pc->NCoeff, status);

    /* Overwrite all NPBLKs so far */
    fits_get_colnum(fptr, CASEINSEN, "NPBLK", &col, status);
    for (i=1; i<=row; i++) {
        fits_write_col(fptr, TLONG, col, i, 1, 1, &row, status);
    }

    ptr = p.telescope;
    fits_get_colnum(fptr, CASEINSEN, "NSITE", &col, status);
    fits_write_col(fptr, TSTRING, col, row, 1, 1, &ptr, status);

    fits_get_colnum(fptr, CASEINSEN, "REF_FREQ", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &pc->RefFreq, status);

    // skipped PRED_PHS

    dtmp = pc->MjdMidInt + pc->MjdMidFrac;
    fits_get_colnum(fptr, CASEINSEN, "REF_MJD", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &dtmp, status);

    fits_get_colnum(fptr, CASEINSEN, "REF_PHS", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &pc->PhRotRef, status);

    fits_get_colnum(fptr, CASEINSEN, "REF_F0", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &pc->FRotRef, status);

    // skipped LGFITERR

    fits_get_colnum(fptr, CASEINSEN, "COEFF", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 15, pc->Coeff, status);

    /* All done */
    return(*status);
}

int psrfits_write_history(fitsfile *fptr, int nsub, int nbin, 
        double rf, int nchan, double ch_bw, int *status) {

    if (*status) { return(*status); }

    /* Point to correct hdu */
    fits_movnam_hdu(fptr, BINARY_TBL, "HISTORY", 0, status);

    /* Get date string */
    int tmp;
    double dtmp;
    char date[32];
    char strtmp[128];
    fits_get_system_time(date, &tmp, status);

    /* Go through columns:
     *  - retrieve colnum using the col name
     *  - write value
     * We're always writing row 1 here (new file).
     */
    char *ptr;
    int i, col, row=1;

    ptr = date;
    fits_get_colnum(fptr, CASEINSEN, "DATE_PRO", &col, status);
    fits_write_col(fptr, TSTRING, col, row, 1, 1, &ptr, status);

    sprintf(strtmp, "asp_result_psrfits");
    ptr = strtmp;
    fits_get_colnum(fptr, CASEINSEN, "PROC_CMD", &col, status);
    fits_write_col(fptr, TSTRING, col, row, 1, 1, &ptr, status);

    sprintf(strtmp, "XXYYCRCI");
    ptr = strtmp;
    fits_get_colnum(fptr, CASEINSEN, "POL_TYPE", &col, status);
    fits_write_col(fptr, TSTRING, col, row, 1, 1, &ptr, status);

    fits_get_colnum(fptr, CASEINSEN, "NSUB", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &nsub, status);

    tmp=4;
    fits_get_colnum(fptr, CASEINSEN, "NPOL", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &tmp, status);

    fits_get_colnum(fptr, CASEINSEN, "NBIN", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &nbin, status);

    fits_get_colnum(fptr, CASEINSEN, "CTR_FREQ", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &rf, status);

    fits_get_colnum(fptr, CASEINSEN, "NCHAN", &col, status);
    fits_write_col(fptr, TINT, col, row, 1, 1, &nchan, status);

    dtmp = -1.0*ch_bw;
    fits_get_colnum(fptr, CASEINSEN, "CHAN_BW", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &dtmp, status);

    /* All done */
    return(*status);
}

int psrfits_write_subint(fitsfile *fptr, int nchan, int nbin, double tint,
        double time, double per, 
        float *rfs, float *scale, float *offset, float *data, 
        int *status) {

    if (*status) { return(*status); }

    /* Go to correct hdu */
    fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, status);

    /* find out how many rows (so we know where to write) */
    long row;
    fits_get_num_rows(fptr, &row, status);
    row++;

    /* Go through columns */
    int col;

    fits_get_colnum(fptr, CASEINSEN, "ISUBINT", &col, status);
    fits_write_col(fptr, TLONG, col, row, 1, 1, &row, status);

    fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &tint, status);

    fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &time, status);

    fits_get_colnum(fptr, CASEINSEN, "PERIOD", &col, status);
    fits_write_col(fptr, TDOUBLE, col, row, 1, 1, &per, status);

    // skipped all position, angle, etc

    fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &col, status);
    fits_write_col(fptr, TFLOAT, col, row, 1, nchan, rfs, status);

    fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &col, status);
    fits_write_col(fptr, TFLOAT, col, row, 1, nchan, scale, status);

    fits_get_colnum(fptr, CASEINSEN, "DAT_OFFS", &col, status);
    fits_write_col(fptr, TFLOAT, col, row, 1, nchan*4, offset, status);

    fits_get_colnum(fptr, CASEINSEN, "DAT_SCL", &col, status);
    fits_write_col(fptr, TFLOAT, col, row, 1, nchan*4, scale, status);

    fits_get_colnum(fptr, CASEINSEN, "DATA", &col, status);
    fits_write_col(fptr, TFLOAT, col, row, 1, nbin*nchan*4, data, status);

    /* Done */
    return(*status);
}

