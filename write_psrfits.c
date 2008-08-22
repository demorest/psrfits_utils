/* write_psrfits.c */
#define _ISOC99_SOURCE  // For long double strtold
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "psrfits.h"
#include "polyco.h"

// Define different obs modes
static const int search=SEARCH_MODE, fold=FOLD_MODE;
int psrfits_obs_mode(const char *obs_mode) {
    if (strncmp("SEARCH", obs_mode, 6)==0) { return(search); }
    else if (strncmp("FOLD", obs_mode, 4)==0) { return(fold); }
    else if (strncmp("PSR", obs_mode, 3)==0) { return(fold); }
    else if (strncmp("CAL", obs_mode, 3)==0) { return(fold); }
    else {
        // TODO: what to do here? default to search for now
        printf("Warning: obs_mode '%s' not recognized, defaulting to SEARCH.\n",
                obs_mode);
        return(search);
    }
    return(search);
}

int psrfits_create(struct psrfits *pf) {
    int itmp, *status;
    long double ldtmp;
    double dtmp;
    char ctmp[40];
    struct hdrinfo *hdr;

    hdr = &(pf->hdr);        // dereference the ptr to the header struct
    status = &(pf->status);  // dereference the ptr to the CFITSIO status

    // Figure out what mode this is 
    int mode=0;
    mode = psrfits_obs_mode(hdr->obs_mode);

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
    // Fold mode template has additional tables (polyco, ephem)
    char *guppi_dir = getenv("GUPPI_DIR");
    char template_file[1024];
    if (guppi_dir==NULL) {
        fprintf(stderr, 
                "Error: GUPPI_DIR environment variable not set, exiting.\n");
        exit(1);
    }
    printf("Opening file '%s' ", pf->filename);
    if (mode==search) { 
        printf("in search mode.\n");
        sprintf(template_file, "%s/%s", guppi_dir, PSRFITS_SEARCH_TEMPLATE);
    } else if (mode==fold) { 
        printf("in fold mode.\n");
        sprintf(template_file, "%s/%s", guppi_dir, PSRFITS_FOLD_TEMPLATE);
    }
    fits_create_template(&(pf->fptr), pf->filename, template_file, status);

    // Go to the primary HDU
    fits_movabs_hdu(pf->fptr, 1, NULL, status);

    // Update the keywords that need it
    fits_get_system_time(ctmp, &itmp, status);
    // Note:  this is the date the file was _written_, not the obs start date
    fits_update_key(pf->fptr, TSTRING, "DATE", ctmp, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "TELESCOP", hdr->telescope,NULL, status);
    fits_update_key(pf->fptr, TSTRING, "OBSERVER", hdr->observer, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "PROJID", hdr->project_id, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "FRONTEND", hdr->frontend, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "BACKEND", hdr->backend, NULL, status);
    if (hdr->summed_polns) {
        if (hdr->npol > 1) {
            printf("Warning!:  Can't have %d polarizations _and_ be summed!\n", 
                   hdr->npol);
        }
        itmp = 2;
        fits_update_key(pf->fptr, TINT, "NRCVR", &itmp, NULL, status);
    } else {
        if (hdr->npol > 2) { // Can't have more than 2 polns
            itmp = 2;
            fits_update_key(pf->fptr, TINT, "NRCVR", &itmp, NULL, status);
        } else {
            fits_update_key(pf->fptr, TINT, "NRCVR", &(hdr->npol), NULL, status);
        }
    }
    fits_update_key(pf->fptr, TSTRING, "FD_POLN", hdr->poln_type, NULL, status);
    // TODO: Need to include specific poln settings PF_HAND< FD_SANG, FD_XYPH
    fits_update_key(pf->fptr, TSTRING, "DATE-OBS", hdr->date_obs, NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "OBSFREQ", &(hdr->fctr), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "OBSBW", &(hdr->BW), NULL, status);
    fits_update_key(pf->fptr, TINT, "OBSNCHAN", &(hdr->orig_nchan), NULL, status);
    fits_update_key(pf->fptr, TSTRING, "SRC_NAME", hdr->source, NULL, status);
    if (strcmp("TRACK", hdr->track_mode)) {
        printf("Warning!:  We don't currently handle non-tracking observations!\n");
        fits_update_key(pf->fptr, TSTRING, "TRK_MODE", hdr->track_mode, NULL, status);
    }
    // TODO: will need to change the following if we aren't tracking!
    fits_update_key(pf->fptr, TSTRING, "RA", hdr->ra_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "DEC", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STT_CRD1", hdr->ra_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STP_CRD1", hdr->ra_str, NULL, status);
    // TODO: update these at the end of the file or obs
    fits_update_key(pf->fptr, TSTRING, "STT_CRD2", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TSTRING, "STP_CRD2", hdr->dec_str, NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "BMAJ", &(hdr->beam_FWHM), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "BMIN", &(hdr->beam_FWHM), NULL, status);
    if (strcmp("OFF", hdr->cal_mode)) {
        fits_update_key(pf->fptr, TDOUBLE, "CAL_FREQ", &(hdr->cal_freq), NULL, status);
        fits_update_key(pf->fptr, TDOUBLE, "CAL_DCYC", &(hdr->cal_dcyc), NULL, status);
        fits_update_key(pf->fptr, TDOUBLE, "CAL_PHS", &(hdr->cal_phs), NULL, status);
    }
    fits_update_key(pf->fptr, TDOUBLE, "SCANLEN", &(hdr->scanlen), NULL, status);
    itmp = (int) hdr->MJD_epoch;
    fits_update_key(pf->fptr, TINT, "STT_IMJD", &itmp, NULL, status);
    ldtmp = (hdr->MJD_epoch - (long double) itmp) * 86400.0L;   // in sec
    itmp = (int) ldtmp;
    fits_update_key(pf->fptr, TINT, "STT_SMJD", &itmp, NULL, status);
    ldtmp -= (long double) itmp;
    dtmp = (double) ldtmp;
    fits_update_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "STT_LST", &(hdr->start_lst), NULL, status);

    // Go to the SUBINT HDU
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "SUBINT", 0, status);

    // Update the keywords that need it
    fits_update_key(pf->fptr, TINT, "NPOL", &(hdr->npol), NULL, status);
    if (!hdr->summed_polns) {
        // TODO:  These need to be updated for the real machine.
        if (hdr->npol==1)
            strcpy(ctmp, "AA");
        else if (hdr->npol==2)
            strcpy(ctmp, "AABB");
        else if (hdr->npol==4)
            strcpy(ctmp, "IQUV");
        fits_update_key(pf->fptr, TSTRING, "POL_TYPE", ctmp, NULL, status);
    } else {
        fits_update_key(pf->fptr, TSTRING, "POL_TYPE", "AA+BB", NULL, status);
    }
    // TODO what does TBIN mean in fold mode?
    fits_update_key(pf->fptr, TDOUBLE, "TBIN", &(hdr->dt), NULL, status);
    fits_update_key(pf->fptr, TINT, "NSUBOFFS", &(hdr->offset_subint), NULL, status);
    fits_update_key(pf->fptr, TINT, "NCHAN", &(hdr->nchan), NULL, status);
    fits_update_key(pf->fptr, TDOUBLE, "CHAN_BW", &(hdr->df), NULL, status);
    if (mode==search) {
        itmp = 1;
        fits_update_key(pf->fptr, TINT, "NSBLK", &(hdr->nsblk), NULL, status);
        fits_update_key(pf->fptr, TINT, "NBITS", &(hdr->nbits), NULL, status);
        fits_update_key(pf->fptr, TINT, "NBIN", &itmp, NULL, status);
    } else if (mode==fold) {
        itmp = 1;
        fits_update_key(pf->fptr, TINT, "NSBLK", &itmp, NULL, status);
        fits_update_key(pf->fptr, TINT, "NBITS", &itmp, NULL, status);
        fits_update_key(pf->fptr, TINT, "NBIN", &(hdr->nbin), NULL, status);
        fits_update_key(pf->fptr, TSTRING, "EPOCHS", "MIDTIME", NULL, status);
    }

    // Update the column sizes for the colums containing arrays
    itmp = hdr->nchan;
    fits_modify_vector_len(pf->fptr, 13, itmp, status); // DAT_FREQ
    fits_modify_vector_len(pf->fptr, 14, itmp, status); // DAT_WTS
    itmp = hdr->nchan * hdr->npol;
    fits_modify_vector_len(pf->fptr, 15, itmp, status); // DAT_OFFS
    fits_modify_vector_len(pf->fptr, 16, itmp, status); // DAT_SCL
    if (mode==search)  
        itmp = (hdr->nbits * hdr->nchan * hdr->npol * hdr->nsblk) / 8;
    else if (mode==fold) 
        itmp = (hdr->nbin * hdr->nchan * hdr->npol);
    fits_modify_vector_len(pf->fptr, 17, itmp, status); // DATA

    // Update the TDIM field for the data column
    if (mode==search)
        sprintf(ctmp, "(1,%d,%d,%d)", hdr->nchan, hdr->npol, hdr->nsblk);
    else if (mode==fold) 
        sprintf(ctmp, "(%d,%d,%d,1)", hdr->nbin, hdr->nchan, hdr->npol);
    fits_update_key(pf->fptr, TSTRING, "TDIM17", ctmp, NULL, status);

    fits_flush_file(pf->fptr, status);
    
    return *status;
}


int psrfits_write_subint(struct psrfits *pf) {
    int row, *status, nchan, nivals, mode;
    float ftmp;
    struct hdrinfo *hdr;
    struct subint *sub;

    hdr = &(pf->hdr);        // dereference the ptr to the header struct
    sub = &(pf->sub);        // dereference the ptr to the subint struct
    status = &(pf->status);  // dereference the ptr to the CFITSIO status
    nchan = hdr->nchan;
    nivals = hdr->nchan * hdr->npol;
    mode = psrfits_obs_mode(hdr->obs_mode);

    // Create the initial file or change to a new one if needed
    if (pf->filenum == 0 || pf->rownum > pf->rows_per_file) {
        if (pf->filenum) {
            printf("Closing file '%s'\n", pf->filename);
            fits_close_file(pf->fptr, status);
        }
        psrfits_create(pf);
    }

    row = pf->rownum;
    fits_write_col(pf->fptr, TDOUBLE, 1, row, 1, 1, &(sub->tsubint), status);
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
    if (mode==search) {
        // Need to change this for other data types...
        fits_write_col(pf->fptr, TBYTE, 17, row, 1, sub->bytes_per_subint, 
                   sub->data, status);
    } else if (mode==fold) { 
        // Fold mode writes floats for now..
        fits_write_col(pf->fptr, TFLOAT, 17, row, 1, sub->bytes_per_subint, 
                   sub->data, status);
    }

    // Flush the buffers if not finished with the file
    // Note:  this use is not entirely in keeping with the CFITSIO
    //        documentation recommendations.  However, manually 
    //        correcting NAXIS2 and using fits_flush_buffer()
    //        caused occasional hangs (and extrememly large
    //        files due to some infinite loop).
    fits_flush_file(pf->fptr, status);

    // Print status if bad
    fits_report_error(stderr, *status);

    // Now update some key values if no CFITSIO errors
    if (!(*status)) {
        pf->rownum++;
        pf->tot_rows++;
        pf->N += hdr->nsblk;
        pf->T = pf->N * hdr->dt;
    }
    
    return *status;
}

int psrfits_write_polycos(struct psrfits *pf, struct polyco *pc, int npc) {

    // Usual setup
    int *status = &(pf->status);

    // If mode!=fold, exit?

    // Save current HDU, move to polyco table
    int hdu;
    fits_get_hdu_num(pf->fptr, &hdu);
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "POLYCO", 0, status);

    int itmp;
    double dtmp;
    char datestr[32], ctmp[32];
    char *cptr;
    fits_get_system_time(datestr, &itmp, status);
    int i, row, col; 
    // XXX start at end of table?
    for (i=0; i<npc; i++) {

        row = i+1;

        cptr = datestr;
        fits_get_colnum(pf->fptr,CASEINSEN,"DATE_PRO",&col,status);
        fits_write_col(pf->fptr,TSTRING,col,row,1,1,&cptr,status);

        sprintf(ctmp, "11.005"); // Tempo version?
        cptr = ctmp;
        fits_get_colnum(pf->fptr,CASEINSEN,"POLYVER",&col,status);
        fits_write_col(pf->fptr,TSTRING,col,row,1,1,&cptr,status);

        fits_get_colnum(pf->fptr,CASEINSEN,"NSPAN",&col,status);
        fits_write_col(pf->fptr,TINT,col,row,1,1,&(pc[i].nmin),status);

        fits_get_colnum(pf->fptr,CASEINSEN,"NCOEF",&col,status);
        fits_write_col(pf->fptr,TINT,col,row,1,1,&(pc[i].nc),status);

        itmp = npc;
        fits_get_colnum(pf->fptr,CASEINSEN,"NPBLK",&col,status);
        fits_write_col(pf->fptr,TINT,col,row,1,1,&itmp,status);

        sprintf(ctmp,"%d", pc[i].nsite); // XXX convert to letter?
        cptr = ctmp;
        fits_get_colnum(pf->fptr,CASEINSEN,"NSITE",&col,status);
        fits_write_col(pf->fptr,TSTRING,col,row,1,1,&cptr,status);

        fits_get_colnum(pf->fptr,CASEINSEN,"REF_FREQ",&col,status);
        fits_write_col(pf->fptr,TFLOAT,col,row,1,1,&(pc[i].rf),status);

        // XXX needs to be accurate??
        dtmp=0.0;
        fits_get_colnum(pf->fptr,CASEINSEN,"PRED_PHS",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dtmp,status);

        dtmp = (double)pc[i].mjd + pc[i].fmjd;
        fits_get_colnum(pf->fptr,CASEINSEN,"REF_MJD",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dtmp,status);

        fits_get_colnum(pf->fptr,CASEINSEN,"REF_PHS",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&(pc[i].rphase),status);

        fits_get_colnum(pf->fptr,CASEINSEN,"REF_F0",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&(pc[i].f0),status);

        // XXX don't parse this yet
        dtmp=-6.0;
        fits_get_colnum(pf->fptr,CASEINSEN,"LGFITERR",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dtmp,status);

        fits_get_colnum(pf->fptr,CASEINSEN,"COEFF",&col,status);
        fits_write_col(pf->fptr,TDOUBLE,col,row,1,pc[i].nc,pc[i].c,status);
    }

    // Go back to orig HDU
    fits_movabs_hdu(pf->fptr, hdu, NULL, status);

    return *status;
}

int psrfits_write_ephem(struct psrfits *pf, FILE *parfile) {
    // Read a pulsar ephemeris (par file) and put it into
    // the psrfits PSREPHEM table.  Only minimal checking
    // is done.
   
    // Get status
    int *status = &(pf->status);

    // Save current HDU, move to psrephem table
    int hdu;
    fits_get_hdu_num(pf->fptr, &hdu);
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "PSREPHEM", 0, status);

    // Loop over lines in par file
    int row=1, col, dtype;
    double dval;
    int ival;
    long double ldval;
    char line[256], *ptr, *saveptr, *key, *val;
    while (fgets(line, 256, parfile)!=NULL) {

        // Convert tabs to spaces
        while ((ptr=strchr(line,'\t'))!=NULL) { *ptr=' '; }

        // strip leading whitespace
        ptr = line;
        while (*ptr==' ') { ptr++; }

        // Identify comments or blank lines
        if (line[0]=='\n' || line[0]=='#' || 
                (line[0]=='C' && line[1]==' '))
            continue;

        // Split into key/val (ignore fit flag and error)
        key = strtok_r(line,  " ", &saveptr);
        val = strtok_r(NULL, " ", &saveptr);
        if (key==NULL || val==NULL) continue; // TODO : complain?

        // Deal with any special cases here
        if (strncmp(key, "PSR", 3)==0)  {

            // PSR(J) -> PSR_NAME
            fits_get_colnum(pf->fptr,CASEINSEN,"PSR_NAME",&col,status);
            fits_write_col(pf->fptr,TSTRING,col,row,1,1,&val,status);

        } else if (strncmp(key, "F0", 2)==0) {

            // F is converted to mHz and split into int/frac
            ldval = strtold(val,NULL) * 1000.0; // Hz->mHz
            ival = (int)ldval;
            dval = ldval - (long double)ival;
            fits_get_colnum(pf->fptr,CASEINSEN,"IF0",&col,status);
            fits_write_col(pf->fptr,TINT,col,row,1,1,&ival,status);
            fits_get_colnum(pf->fptr,CASEINSEN,"FF0",&col,status);
            fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dval,status);

        } else if (strncmp(key, "TZRMJD", 6)==0) {

            // TZRMJD is split into int/frac
            ldval = strtold(val,NULL);
            ival = (int)ldval;
            dval = ldval - (long double)ival;
            fits_get_colnum(pf->fptr,CASEINSEN,"TZRIMJD",&col,status);
            fits_write_col(pf->fptr,TINT,col,row,1,1,&ival,status);
            fits_get_colnum(pf->fptr,CASEINSEN,"TZRFMJD",&col,status);
            fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dval,status);

        } else {

            // Find column, skip/warn if this one isn't present
            fits_get_colnum(pf->fptr,CASEINSEN,key,&col,status);
            if (*status==COL_NOT_FOUND) {
                fprintf(stderr, 
                        "psrfits_write_epherm warning: Couldn't find keyword %s "
                        "in ephemeris table.\n",
                        key);
                *status=0;
                continue;
            }

            // Need to convert string to appropriate column data type
            // and then write it to the column.  These should all be
            // either double int or string.
            fits_get_coltype(pf->fptr,col,&dtype,NULL,NULL,status);
            if (dtype==TDOUBLE || dtype==TFLOAT) { 
                dval = atof(val);
                fits_write_col(pf->fptr,TDOUBLE,col,row,1,1,&dval,status);
            } else if (dtype==TINT || dtype==TLONG || dtype==TSHORT) {
                ival = atoi(val);
                fits_write_col(pf->fptr,TINT,col,row,1,1,&ival,status);
            } else if (dtype==TSTRING) {
                fits_write_col(pf->fptr,TSTRING,col,row,1,1,&val,status);
            } else {
                fprintf(stderr, "psrfits_write_ephem warning: "
                        "Unhandled column datatype (key=%s)\n", key);
                continue;
            }
        }

        // sucess/failure
        if (*status) {
            fprintf(stderr, "psrfits_write_ephem failed: key=%s val=%s\n",
                    key, val);
            fits_report_error(stderr, *status);
            *status=0;
        } 
#if 0  // DEBUG
        else {
            fprintf(stderr, "psrfits_write_ephem success: key=%s val=%s\n",
                    key, val);
        }
#endif

    }

    // Go back to orig HDU
    fits_movabs_hdu(pf->fptr, hdu, NULL, status);

    return *status;
}

int psrfits_remove_polycos(struct psrfits *pf) {
    // Delete the polyco table
    
    int *status = &(pf->status);

    // Save current HDU, move to polyco table
    int hdu;
    fits_get_hdu_num(pf->fptr, &hdu);
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "POLYCO", 0, status);

    // Delete it
    fits_delete_hdu(pf->fptr, NULL, status);

    // Go to the SUBINT HDU
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "SUBINT", 0, status);

    return *status;
}

int psrfits_remove_ephem(struct psrfits *pf) {
    // Delete the polyco table
    
    int *status = &(pf->status);

    // Save current HDU, move to polyco table
    int hdu;
    fits_get_hdu_num(pf->fptr, &hdu);
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "PSREPHEM", 0, status);

    // Delete it
    fits_delete_hdu(pf->fptr, NULL, status);

    // Go to the SUBINT HDU
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "SUBINT", 0, status);

    return *status;
}

int psrfits_close(struct psrfits *pf) {
    if (!pf->status) {
        fits_close_file(pf->fptr, &(pf->status));
        printf("Closing file '%s'\n", pf->filename);
    }
    printf("Done.  Wrote %d subints (%f sec) in %d files (status = %d).\n",
           pf->tot_rows, pf->T, pf->filenum, pf->status);
    return pf->status;
}
