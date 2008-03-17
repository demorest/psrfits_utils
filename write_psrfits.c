#include <stdio.h>
#include <string.h>
#include <time.h>
#include "write_psrfits.h"

int psrfits_create_searchmode(fitsfile **fptr, 
                              struct psrfits_hdrinfo *pfh, 
                              int *status) {
    int itmp;
    long double ldtmp;
    double dtmp;
    char ctmp[40];
    time_t t;
    struct tm ts;

    // Create basic FITS file from our template
    fits_create_template(fptr, pfh->filename, PSRFITS_TEMPLATE, status);

    // Go to the primary HDU
    fits_movabs_hdu(*fptr, 1, NULL, status);

    // Update the keywords that need it
    sprintf(ctmp, "(1,%d,%d,%d)", pfh->nchan, pfh->npol, pfh->nsblk);
    // Note:  this is the date that the file was _written_, not the 
    // observation start date
    t = time(NULL);
    gmtime_r(&t, &ts);
    strftime(ctmp, 40, "%FT%T", &ts); 
    fits_update_key(*fptr, TSTRING, "DATE", ctmp, NULL, status);
    fits_update_key(*fptr, TSTRING, "OBSERVER", pfh->observer, NULL, status);
    fits_update_key(*fptr, TSTRING, "PROJID", pfh->project_id, NULL, status);
    fits_update_key(*fptr, TSTRING, "FRONTEND", pfh->frontend, NULL, status);
    if (pfh->summed_polns) {
        if (pfh->npol > 1) {
            printf("Warning!:  Can't have %d polarizations _and_ be summed!\n", 
                   pfh->npol);
        }
        itmp = 2;
        fits_update_key(*fptr, TLONG, "NRCVR", &itmp, NULL, status);
    } else {
        if (pfh->npol > 2) { // Can't have more than 2 polns
            itmp = 2;
            fits_update_key(*fptr, TLONG, "NRCVR", &itmp, NULL, status);
        } else {
            fits_update_key(*fptr, TLONG, "NRCVR", &(pfh->npol), NULL, status);
        }
    }
    fits_update_key(*fptr, TSTRING, "FD_POLN", pfh->poln_type, NULL, status);
    // Need to make entries here for the specific polarization 
    // settings PF_HAND< FD_SANG, FD_XYPH
    fits_update_key(*fptr, TSTRING, "DATE-OBS", pfh->date_obs, NULL, status);
    fits_update_key(*fptr, TDOUBLE, "OBSFREQ", &(pfh->fctr), NULL, status);
    fits_update_key(*fptr, TDOUBLE, "OBSBW", &(pfh->BW), NULL, status);
    fits_update_key(*fptr, TLONG, "OBSNCHAN", &(pfh->orig_nchan), NULL, status);
    fits_update_key(*fptr, TSTRING, "SRC_NAME", pfh->source, NULL, status);
    if (strcmp("TRACK", pfh->track_mode)) {
        printf("Warning!:  We don't currently handle non-tracking observations!\n");
        fits_update_key(*fptr, TSTRING, "TRK_MODE", pfh->track_mode, NULL, status);
    }
    // Note:  will need to change the following if we aren't tracking!
    fits_update_key(*fptr, TSTRING, "RA", pfh->ra_str, NULL, status);
    fits_update_key(*fptr, TSTRING, "DEC", pfh->dec_str, NULL, status);
    fits_update_key(*fptr, TSTRING, "STT_CRD1", pfh->ra_str, NULL, status);
    fits_update_key(*fptr, TSTRING, "STP_CRD1", pfh->ra_str, NULL, status);
    // Should update these at the end of the file (or obs?)
    fits_update_key(*fptr, TSTRING, "STT_CRD2", pfh->dec_str, NULL, status);
    fits_update_key(*fptr, TSTRING, "STP_CRD2", pfh->dec_str, NULL, status);
    fits_update_key(*fptr, TDOUBLE, "BMAJ", &(pfh->beam_FWHM), NULL, status);
    fits_update_key(*fptr, TDOUBLE, "BMIN", &(pfh->beam_FWHM), NULL, status);
    if (strcmp("OFF", pfh->cal_mode)) {
        printf("Yikes!\n");
        fits_update_key(*fptr, TDOUBLE, "CAL_FREQ", &(pfh->cal_freq), NULL, status);
        fits_update_key(*fptr, TDOUBLE, "CAL_DCYC", &(pfh->cal_dcyc), NULL, status);
        fits_update_key(*fptr, TDOUBLE, "CAL_PHS", &(pfh->cal_phs), NULL, status);
    }
    fits_update_key(*fptr, TDOUBLE, "SCANLEN", &(pfh->scanlen), NULL, status);
    itmp = (int) pfh->MJD_epoch;
    fits_update_key(*fptr, TLONG, "STT_IMJD", &itmp, NULL, status);
    ldtmp = (pfh->MJD_epoch - (long double) itmp) * 86400.0L;   // in sec
    itmp = (int) ldtmp;
    fits_update_key(*fptr, TLONG, "STT_SMJD", &itmp, NULL, status);
    ldtmp -= (long double) itmp;
    dtmp = (double) ldtmp;
    fits_update_key(*fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, status);
    // Note:  1 sidereal day = 86164.0905 seconds
    fits_update_key(*fptr, TDOUBLE, "STT_LST", &(pfh->start_lst), NULL, status);

    // Go to the SUBINT HDU
    fits_movabs_hdu(*fptr, 2, NULL, status);

    // Update the keywords that need it
    fits_update_key(*fptr, TLONG, "NPOL", &(pfh->npol), NULL, status);
    if (!pfh->summed_polns) {
        printf("Warning!: POL_TYPE might be incorrect!\n");
        if (pfh->npol==1)
            strcpy(ctmp, "AA");
        else if (pfh->npol==2)
            strcpy(ctmp, "AABB");
        else if (pfh->npol==4)
            strcpy(ctmp, "AABBCRCI");
        fits_update_key(*fptr, TSTRING, "POL_TYPE", ctmp, NULL, status);
    }
    fits_update_key(*fptr, TDOUBLE, "TBIN", &(pfh->dt), NULL, status);
    fits_update_key(*fptr, TLONG, "NBITS", &(pfh->nbits), NULL, status);
    fits_update_key(*fptr, TLONG, "NSUBOFFS", &(pfh->start_spec), NULL, status);
    fits_update_key(*fptr, TLONG, "NCHAN", &(pfh->nchan), NULL, status);
    fits_update_key(*fptr, TDOUBLE, "CHAN_BW", &(pfh->df), NULL, status);
    fits_update_key(*fptr, TLONG, "NSBLK", &(pfh->nsblk), NULL, status);

    // Update the column sizes for the colums containing arrays
    itmp = pfh->nchan;
    fits_modify_vector_len(*fptr, 13, itmp, status); // DAT_FREQ
    fits_modify_vector_len(*fptr, 14, itmp, status); // DAT_WTS
    itmp = pfh->nchan * pfh->npol;
    fits_modify_vector_len(*fptr, 15, itmp, status); // DAT_OFFS
    fits_modify_vector_len(*fptr, 16, itmp, status); // DAT_SCL
    itmp = (pfh->nbits * pfh->nchan * pfh->npol * pfh->nsblk) / 8;
    fits_modify_vector_len(*fptr, 17, itmp, status); // DATA

    // Update the TDIM field for the data column
    sprintf(ctmp, "(1,%d,%d,%d)", pfh->nchan, pfh->npol, pfh->nsblk);
    fits_update_key(*fptr, TSTRING, "TDIM17", ctmp, NULL, status);
    // The following is an alternate way to do that, but it copies
    // the TDIM17 field instead of updating it
    // long naxes[4] = {1, pfh->nchan, pfh->npol, pfh->nsblk};
    // fits_write_tdim(*fptr, 17, 4, naxes, status);
    
    return *status;
}


int psrfits_write_subint(fitsfile *fptr, 
                         struct psrfits_subint *pfs, 
                         int *status) {
    int row;
    float ftmp;

    row = (pfs->rownum + 1) % pfs->max_rows;
    fits_write_col(fptr, TDOUBLE, 1, row, 1, 1, &(pfs->tsubint), status);
    fits_write_col(fptr, TDOUBLE, 2, row, 1, 1, &(pfs->offs), status);
    fits_write_col(fptr, TDOUBLE, 3, row, 1, 1, &(pfs->lst), status);
    fits_write_col(fptr, TDOUBLE, 4, row, 1, 1, &(pfs->ra), status);
    fits_write_col(fptr, TDOUBLE, 5, row, 1, 1, &(pfs->dec), status);
    fits_write_col(fptr, TDOUBLE, 6, row, 1, 1, &(pfs->glon), status);
    fits_write_col(fptr, TDOUBLE, 7, row, 1, 1, &(pfs->glat), status);
    ftmp = (float) pfs->feed_ang;
    fits_write_col(fptr, TFLOAT, 8, row, 1, 1, &ftmp, status);
    ftmp = (float) pfs->pos_ang;
    fits_write_col(fptr, TFLOAT, 9, row, 1, 1, &ftmp, status);
    ftmp = (float) pfs->par_ang;
    fits_write_col(fptr, TFLOAT, 10, row, 1, 1, &ftmp, status);
    ftmp = (float) pfs->tel_az;
    fits_write_col(fptr, TFLOAT, 11, row, 1, 1, &ftmp, status);
    ftmp = (float) pfs->tel_zen;
    fits_write_col(fptr, TFLOAT, 12, row, 1, 1, &ftmp, status);
    fits_write_col(fptr, TFLOAT, 13, row, 1, pfs->nchan, 
                   pfs->dat_freqs, status);
    fits_write_col(fptr, TFLOAT, 14, row, 1, pfs->nchan, 
                   pfs->dat_weights, status);
    fits_write_col(fptr, TFLOAT, 15, row, 1, pfs->nchan * pfs->npol, 
                   pfs->dat_offsets, status);
    fits_write_col(fptr, TFLOAT, 16, row, 1, pfs->nchan * pfs->npol, 
                   pfs->dat_scales, status);
    fits_write_col(fptr, TBYTE, 17, row, 1, pfs->bytes_per_subint, 
                   pfs->data, status);

    return *status;
}
