#include <stdio.h>
#include <string.h>
#include <time.h>
#include "write_psrfits.h"

int psrfits_create_searchmode(struct psrfits *pf) {
    int itmp, *status;
    long double ldtmp;
    double dtmp;
    char ctmp[40];
    time_t t;
    struct tm ts;
    struct hdrinfo *hdr;

    hdr = &(pf->hdr);        // dereference the ptr to the header struct
    status = &(pf->status);  // dereference the ptr to the CFITSIO status

    // Initialize the key variables if needed
    if (pf->filenum == 0) {  // first time writing to the file
        pf->status = 0;
        pf->tot_rows = 0;
        pf->N = 0L;
        pf->T = 0.0;
        hdr->offset_subint = 0;
    }
    pf->filenum++;
    pf->rownum = 1;
    hdr->offset_subint = pf->tot_rows;

    // Update the filename
    sprintf(pf->filename, "%s_%04d.fits", pf->basefilename, pf->filenum);

    // Create basic FITS file from our template
    fits_create_template(&(pf->fptr), pf->filename, PSRFITS_TEMPLATE, status);

    // Go to the primary HDU
    fits_movabs_hdu(pf->fptr, 1, NULL, status);

    // Update the keywords that need it
    sprintf(ctmp, "(1,%d,%d,%d)", hdr->nchan, hdr->npol, hdr->nsblk);
    // Note:  this is the date that the file was _written_, not the 
    // observation start date
    t = time(NULL);
    gmtime_r(&t, &ts);
    strftime(ctmp, 40, "%FT%T", &ts); 
    fits_update_key(pf->fptr, TSTRING, "DATE", ctmp, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "OBSERVER", hdr->observer, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "PROJID", hdr->project_id, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "FRONTEND", hdr->frontend, NULL, status);
    if (hdr->summed_polns) {
        if (hdr->npol > 1) {
            printf("Warning!:  Can't have %d polarizations _and_ be summed!\n", 
                   hdr->npol);
        }
        itmp = 2;
        fits_update_key(pf->fptr, TLONG, "NRCVR", &itmp, NULL, status);
    } else {
        if (hdr->npol > 2) { // Can't have more than 2 polns
            itmp = 2;
            fits_update_key(pf->fptr, TLONG, "NRCVR", &itmp, NULL, status);
        } else {
            fits_update_key(pf->fptr, TLONG, "NRCVR", &(hdr->npol), NULL, status);
        }
    }
    fits_update_key(pf->fptr, TSTRING, "FD_POLN", hdr->poln_type, NULL, status);
    // Need to make entries here for the specific polarization 
    // settings PF_HAND< FD_SANG, FD_XYPH
    fits_update_key(pf->fptr, TSTRING, "DATE-OBS", hdr->date_obs, NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "OBSFREQ", &(hdr->fctr), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "OBSBW", &(hdr->BW), NULL, status);
    fits_update_key(pf->fptr, TLONG, "OBSNCHAN", &(hdr->orig_nchan), NULL, status);
    fits_update_key(pf->fptr, TSTRING, "SRC_NAME", hdr->source, NULL, status);
    if (strcmp("TRACK", hdr->track_mode)) {
        printf("Warning!:  We don't currently handle non-tracking observations!\n");
        fits_update_key(pf->fptr, TSTRING, "TRK_MODE", hdr->track_mode, NULL, status);
    }
    // Note:  will need to change the following if we aren't tracking!
    fits_update_key(pf->fptr, TSTRING, "RA", hdr->ra_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "DEC", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STT_CRD1", hdr->ra_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STP_CRD1", hdr->ra_str, NULL, status);
    // Should update these at the end of the file (or obs?)
    fits_update_key(pf->fptr, TSTRING, "STT_CRD2", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STP_CRD2", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "BMAJ", &(hdr->beam_FWHM), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "BMIN", &(hdr->beam_FWHM), NULL, status);
    if (strcmp("OFF", hdr->cal_mode)) {
        printf("Yikes!\n");
        fits_update_key(pf->fptr, TDOUBLE, "CAL_FREQ", &(hdr->cal_freq), NULL, status);
        fits_update_key(pf->fptr, TDOUBLE, "CAL_DCYC", &(hdr->cal_dcyc), NULL, status);
        fits_update_key(pf->fptr, TDOUBLE, "CAL_PHS", &(hdr->cal_phs), NULL, status);
    }
    fits_update_key(pf->fptr, TDOUBLE, "SCANLEN", &(hdr->scanlen), NULL, status);
    itmp = (int) hdr->MJD_epoch;
    fits_update_key(pf->fptr, TLONG, "STT_IMJD", &itmp, NULL, status);
    ldtmp = (hdr->MJD_epoch - (long double) itmp) * 86400.0L;   // in sec
    itmp = (int) ldtmp;
    fits_update_key(pf->fptr, TLONG, "STT_SMJD", &itmp, NULL, status);
    ldtmp -= (long double) itmp;
    dtmp = (double) ldtmp;
    fits_update_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, status);
    // Note:  1 sidereal day = 86164.0905 seconds
    // CALL sla_OBS (N, C, NAME, W, P, H)
    // sla_GMST (UT1)
    fits_update_key(pf->fptr, TDOUBLE, "STT_LST", &(hdr->start_lst), NULL, status);

    // Go to the SUBINT HDU
    fits_movabs_hdu(pf->fptr, 2, NULL, status);

    // Update the keywords that need it
    fits_update_key(pf->fptr, TLONG, "NPOL", &(hdr->npol), NULL, status);
    if (!hdr->summed_polns) {
        printf("Warning!: POL_TYPE might be incorrect!\n");
        if (hdr->npol==1)
            strcpy(ctmp, "AA");
        else if (hdr->npol==2)
            strcpy(ctmp, "AABB");
        else if (hdr->npol==4)
            strcpy(ctmp, "AABBCRCI");
        fits_update_key(pf->fptr, TSTRING, "POL_TYPE", ctmp, NULL, status);
    }
    fits_update_key(pf->fptr, TDOUBLE, "TBIN", &(hdr->dt), NULL, status);
    fits_update_key(pf->fptr, TLONG, "NBITS", &(hdr->nbits), NULL, status);
    fits_update_key(pf->fptr, TLONG, "NSUBOFFS", &(hdr->offset_subint), NULL, status);
    fits_update_key(pf->fptr, TLONG, "NCHAN", &(hdr->nchan), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "CHAN_BW", &(hdr->df), NULL, status);
    fits_update_key(pf->fptr, TLONG, "NSBLK", &(hdr->nsblk), NULL, status);

    // Update the column sizes for the colums containing arrays
    itmp = hdr->nchan;
    fits_modify_vector_len(pf->fptr, 13, itmp, status); // DAT_FREQ
    fits_modify_vector_len(pf->fptr, 14, itmp, status); // DAT_WTS
    itmp = hdr->nchan * hdr->npol;
    fits_modify_vector_len(pf->fptr, 15, itmp, status); // DAT_OFFS
    fits_modify_vector_len(pf->fptr, 16, itmp, status); // DAT_SCL
    itmp = (hdr->nbits * hdr->nchan * hdr->npol * hdr->nsblk) / 8;
    fits_modify_vector_len(pf->fptr, 17, itmp, status); // DATA

    // Update the TDIM field for the data column
    sprintf(ctmp, "(1,%d,%d,%d)", hdr->nchan, hdr->npol, hdr->nsblk);
    fits_update_key(pf->fptr, TSTRING, "TDIM17", ctmp, NULL, status);
    // The following is an alternate way to do that, but it copies
    // the TDIM17 field instead of updating it
    // long naxes[4] = {1, hdr->nchan, hdr->npol, hdr->nsblk};
    // fits_write_tdim(pf->fptr, 17, 4, naxes, status);
    
    return *status;
}


int psrfits_write_subint(struct psrfits *pf) {
    int row, *status, nchan, nivals;
    float ftmp;
    struct hdrinfo *hdr;
    struct subint *sub;

    hdr = &(pf->hdr);        // dereference the ptr to the header struct
    sub = &(pf->sub);        // dereference the ptr to the subint struct
    status = &(pf->status);  // dereference the ptr to the CFITSIO status
    nchan = hdr->nchan;
    nivals = hdr->nchan * hdr->npol;

    // Create the initial file or change to a new one if needed
    if (pf->filenum == 0 || pf->rownum > pf->rows_per_file) {
        if (pf->filenum) fits_close_file(pf->fptr, status);
        psrfits_create_searchmode(pf);
    }

    row = pf->rownum;
    fits_write_col(pf->fptr, TDOUBLE, 1, row, 1, 1, &(sub->tsubint), status);
    sub->offs = (pf->tot_rows + 0.5) * sub->tsubint;
    fits_write_col(pf->fptr, TDOUBLE, 2, row, 1, 1, &(sub->offs), status);
    fits_write_col(pf->fptr, TDOUBLE, 3, row, 1, 1, &(sub->lst), status);
    fits_write_col(pf->fptr, TDOUBLE, 4, row, 1, 1, &(sub->ra), status);
    fits_write_col(pf->fptr, TDOUBLE, 5, row, 1, 1, &(sub->dec), status);
    fits_write_col(pf->fptr, TDOUBLE, 6, row, 1, 1, &(sub->glon), status);
    fits_write_col(pf->fptr, TDOUBLE, 7, row, 1, 1, &(sub->glat), status);
    ftmp = (float) sub->feed_ang;
    fits_write_col(pf->fptr, TFLOAT, 8, row, 1, 1, &ftmp, status);
    ftmp = (float) sub->pos_ang;
    fits_write_col(pf->fptr, TFLOAT, 9, row, 1, 1, &ftmp, status);
    ftmp = (float) sub->par_ang;
    fits_write_col(pf->fptr, TFLOAT, 10, row, 1, 1, &ftmp, status);
    ftmp = (float) sub->tel_az;
    fits_write_col(pf->fptr, TFLOAT, 11, row, 1, 1, &ftmp, status);
    ftmp = (float) sub->tel_zen;
    fits_write_col(pf->fptr, TFLOAT, 12, row, 1, 1, &ftmp, status);
    fits_write_col(pf->fptr, TFLOAT, 13, row, 1, nchan, sub->dat_freqs, status);
    fits_write_col(pf->fptr, TFLOAT, 14, row, 1, nchan, sub->dat_weights, status);
    fits_write_col(pf->fptr, TFLOAT, 15, row, 1, nivals, sub->dat_offsets, status);
    fits_write_col(pf->fptr, TFLOAT, 16, row, 1, nivals, sub->dat_scales, status);
    // Need to change this for other data types...
    fits_write_col(pf->fptr, TBYTE, 17, row, 1, sub->bytes_per_subint, 
                   sub->data, status);

    // Flush the buffers
    if (pf->tot_rows==0) {
        fits_flush_file(pf->fptr, status);
    } else {
        fits_update_key(pf->fptr, TLONG, "NAXIS2", &(pf->rownum), NULL, status);
        fits_flush_buffer(pf->fptr, 0, status);
    }

    // Now update some key values
    pf->rownum++;
    pf->tot_rows++;
    pf->N += hdr->nsblk;
    pf->T = pf->N * hdr->dt;
    
    return *status;
}
