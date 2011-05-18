/* read_psrfits.c 
 * Paul Demorest, 05/2008
 */
#include <stdio.h>
#include <string.h>
#include "psrfits.h"

/* This function is similar to psrfits_create, except it
 * deals with reading existing files.  It is assumed that
 * basename and filenum are filled in correctly to point to 
 * the first file in the set OR that filename already contains
 * the correct file name.
 */
int psrfits_open(struct psrfits *pf) {

    int itmp;
    double dtmp;
    char ctmp[256];

    struct hdrinfo *hdr = &(pf->hdr);
    struct subint  *sub = &(pf->sub);
    struct foldinfo *fold = &(pf->fold);
    int *status = &(pf->status);

    sprintf(ctmp, "%s_%04d.fits", pf->basefilename, pf->filenum-1);
    if (pf->filename[0]=='\0' || 
        ((pf->filenum > 1) && (strcmp(ctmp, pf->filename)==0)))
        // The 2nd test checks to see if we are creating filenames ourselves
        sprintf(pf->filename, "%s_%04d.fits", pf->basefilename, pf->filenum);

    fits_open_file(&(pf->fptr), pf->filename, READONLY, status);
    pf->mode = 'r';

    // If file no exist, exit now
    if (*status) {
        return *status; 
    } else {
        printf("Opened file '%s'\n", pf->filename);
    }

    // Move to main HDU
    fits_movabs_hdu(pf->fptr, 1, NULL, status);

    // Figure out obs mode
    fits_read_key(pf->fptr, TSTRING, "OBS_MODE", hdr->obs_mode, NULL, status);
    int mode = psrfits_obs_mode(hdr->obs_mode);

    // Set the downsampling stuff to default values
    hdr->onlyI = 0;
    hdr->ds_time_fact = 1;
    hdr->ds_freq_fact = 1;

    // Blank parfile name, folding params
    fold->parfile[0] = '\0';
    fold->n_polyco_sets = 0;
    fold->pc = NULL;

    // Read some stuff
    fits_read_key(pf->fptr, TSTRING, "TELESCOP", hdr->telescope, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "OBSERVER", hdr->observer, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "PROJID", hdr->project_id, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "FRONTEND", hdr->frontend, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "BACKEND", hdr->backend, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "FD_POLN", hdr->poln_type, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "DATE-OBS", hdr->date_obs, NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "OBSFREQ", &(hdr->fctr), NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "OBSBW", &(hdr->BW), NULL, status);
    fits_read_key(pf->fptr, TINT, "OBSNCHAN", &(hdr->orig_nchan), NULL, status);
    hdr->orig_df = hdr->BW / hdr->orig_nchan;
    fits_read_key(pf->fptr, TDOUBLE, "CHAN_DM", &(hdr->chan_dm), NULL, status);
    if (*status==KEY_NO_EXIST) { hdr->chan_dm=0.0; *status=0; }
    fits_read_key(pf->fptr, TSTRING, "SRC_NAME", hdr->source, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "TRK_MODE", hdr->track_mode, NULL, status);
    // TODO warn if not TRACK?
    fits_read_key(pf->fptr, TSTRING, "RA", hdr->ra_str, NULL, status);
    fits_read_key(pf->fptr, TSTRING, "DEC", hdr->dec_str, NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "BMAJ", &(hdr->beam_FWHM), NULL, status);
    fits_read_key(pf->fptr, TSTRING, "CAL_MODE", hdr->cal_mode, NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "CAL_FREQ", &(hdr->cal_freq), NULL, 
            status);
    fits_read_key(pf->fptr, TDOUBLE, "CAL_DCYC", &(hdr->cal_dcyc), NULL, 
            status);
    fits_read_key(pf->fptr, TDOUBLE, "CAL_PHS", &(hdr->cal_phs), NULL, status);
    fits_read_key(pf->fptr, TSTRING, "FD_MODE", hdr->feed_mode, NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "FA_REQ", &(hdr->feed_angle), NULL, 
            status);
    fits_read_key(pf->fptr, TDOUBLE, "SCANLEN", &(hdr->scanlen), NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "FD_SANG", &(hdr->fd_sang), NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "FD_XYPH", &(hdr->fd_xyph), NULL, status);
    fits_read_key(pf->fptr, TINT, "FD_HAND", &(hdr->fd_hand), NULL, status);
    fits_read_key(pf->fptr, TINT, "BE_PHASE", &(hdr->be_phase), NULL, status);

    fits_read_key(pf->fptr, TINT, "STT_IMJD", &itmp, NULL, status);
    hdr->MJD_epoch = (long double)itmp;
    hdr->start_day = itmp;
    fits_read_key(pf->fptr, TDOUBLE, "STT_SMJD", &dtmp, NULL, status);
    hdr->MJD_epoch += dtmp/86400.0L;
    hdr->start_sec = dtmp;
    fits_read_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, status);
    hdr->MJD_epoch += dtmp/86400.0L;
    hdr->start_sec += dtmp;

    fits_read_key(pf->fptr, TDOUBLE, "STT_LST", &(hdr->start_lst), NULL, 
            status);

    // Move to first subint
    fits_movnam_hdu(pf->fptr, BINARY_TBL, "SUBINT", 0, status);

    // Read some more stuff
    fits_read_key(pf->fptr, TINT, "NPOL", &(hdr->npol), NULL, status);
    fits_read_key(pf->fptr, TSTRING, "POL_TYPE", &(hdr->poln_order), NULL, status);
    if (strncmp(hdr->poln_order, "AA+BB", 6)==0) hdr->summed_polns=1;
    else hdr->summed_polns=0;
    fits_read_key(pf->fptr, TDOUBLE, "TBIN", &(hdr->dt), NULL, status);
    fits_read_key(pf->fptr, TINT, "NBIN", &(hdr->nbin), NULL, status);
    fits_read_key(pf->fptr, TINT, "NSUBOFFS", &(hdr->offset_subint), NULL, 
            status);
    fits_read_key(pf->fptr, TINT, "NCHAN", &(hdr->nchan), NULL, status);
    fits_read_key(pf->fptr, TDOUBLE, "CHAN_BW", &(hdr->df), NULL, status);
    fits_read_key(pf->fptr, TINT, "NSBLK", &(hdr->nsblk), NULL, status);
    fits_read_key(pf->fptr, TINT, "NBITS", &(hdr->nbits), NULL, status);

    if (mode==SEARCH_MODE) {
        long long lltmp = hdr->nsblk;  // Prevents a possible overflow in numerator below
        lltmp = (lltmp * hdr->nbits * hdr->nchan * hdr->npol) / 8L;
        sub->bytes_per_subint = (int) lltmp;
    } else if (mode==FOLD_MODE) {
        sub->bytes_per_subint = 
            (hdr->nbin * hdr->nchan * hdr->npol); // XXX data type??
    }

    // Init counters
    pf->rownum = 1;
    fits_read_key(pf->fptr, TINT, "NAXIS2", &(pf->rows_per_file), NULL, status);

    return *status;
}


void apply_scales_and_offsets(int numchan, int numpol, int numspect,
                              float *scales, float *offsets,
                              unsigned char *inbuf, float *outbuf)
{
    int ii, jj;
    unsigned char *inptr = inbuf;
    float *outptr = outbuf;
    const int N = numchan * numpol;

    for (ii = 0 ; ii < numspect ; ii++) {
        float *sptr = scales;
        float *optr = offsets;
        for (jj = 0 ; jj < N ; jj++, sptr++, optr++, inptr++, outptr++) {
            *outptr = *sptr * (float)(*inptr) + *optr;
        }
    }
}


void scale_and_offset_data(struct psrfits *pf) {
    // Make sure that pf->sub.fdata has been allocated!
    apply_scales_and_offsets(pf->hdr.nchan, pf->hdr.npol, pf->hdr.nsblk,
                             pf->sub.dat_scales, pf->sub.dat_offsets,
                             pf->sub.data, pf->sub.fdata);
}


/* Read next subint from the set of files described
 * by the psrfits struct.  It is assumed that all files
 * form a consistent set.  Read automatically goes to the
 * next file when one ends.  Arrays should be allocated
 * outside this routine.
 */
int psrfits_read_subint(struct psrfits *pf) {

    struct hdrinfo *hdr = &(pf->hdr);
    struct subint  *sub = &(pf->sub);
    int colnum = 0, *status = &(pf->status);

    // See if we need to move to next file
    if (pf->rownum > pf->rows_per_file) {
        printf("Closing file '%s'\n", pf->filename);
        fits_close_file(pf->fptr, status);
        pf->filenum++;
        psrfits_open(pf);
        if (*status==104) {
            printf("Finished with all input files.\n");
            pf->filenum--;
            *status = 1;
            return *status;
        }
    }

    int mode = psrfits_obs_mode(hdr->obs_mode);
    int nchan = hdr->nchan;
    int nivals = hdr->nchan * hdr->npol;
    int row = pf->rownum;

    // TODO: bad! really need to base this on column names
    fits_get_colnum(pf->fptr, 0, "TSUBINT", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->tsubint),
            NULL, status);
    double last_offs = sub->offs;
    fits_get_colnum(pf->fptr, 0, "OFFS_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->offs),
            NULL, status);
    // Hack to fix wrapping in coherent data
    if (pf->tot_rows > 0) {
        double delta_offs = sub->offs - last_offs;
	double wrap_offs = 4294967296L * hdr->dt;
        if (delta_offs < -0.5*wrap_offs) {
            sub->offs += wrap_offs;
	    fprintf(stderr, "Warning: detected likely counter wrap, attempting to fix it.\n");
        }
    }
    fits_get_colnum(pf->fptr, 0, "LST_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->lst),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "RA_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->ra),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DEC_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->dec),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "GLON_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->glon),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "GLAT_SUB", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->glat),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "FD_ANG", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->feed_ang),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "POS_ANG", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->pos_ang),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "PAR_ANG", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->par_ang),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "TEL_AZ", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->tel_az),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "TEL_ZEN", &colnum, status);
    fits_read_col(pf->fptr, TDOUBLE, colnum, row, 1, 1, NULL, &(sub->tel_zen),
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DAT_FREQ", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, row, 1, nchan, NULL, sub->dat_freqs,
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DAT_WTS", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, row, 1, nchan, NULL, sub->dat_weights,
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DAT_OFFS", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, row, 1, nivals, NULL, sub->dat_offsets,
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DAT_SCL", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, row, 1, nivals, NULL, sub->dat_scales,
            NULL, status);
    fits_get_colnum(pf->fptr, 0, "DATA", &colnum, status);
    if (mode==SEARCH_MODE) {
        fits_read_col(pf->fptr, TBYTE, colnum, row, 1, sub->bytes_per_subint,
                      NULL, sub->rawdata, NULL, status);
        if (hdr->nbits==4) pf_4bit_to_8bit(pf);
    } else if (mode==FOLD_MODE) {
        fits_read_col(pf->fptr, TFLOAT, colnum, row, 1, sub->bytes_per_subint,
                      NULL, sub->data, NULL, status);
    }

    // Complain on error
    fits_report_error(stderr, *status);

    // Update counters
    if (!(*status)) {
        pf->rownum++;
        pf->tot_rows++;
        pf->N += hdr->nsblk;
        pf->T = pf->N * hdr->dt;
    }

    return *status;
}

/* Read only part (N "spectra" in time) of the DATA column from the
 * current subint from the set of files described by the psrfits
 * struct.  Put the data into the buffer "buffer".  It is assumed that
 * all files form a consistent set.  Read automatically goes to the
 * next file when one ends.  Arrays should be allocated outside this
 * routine.  Counters are _not_ updated as they are in
 * psrfits_read_subint().
 */
int psrfits_read_part_DATA(struct psrfits *pf, int N, float *fbuffer) {

    struct hdrinfo *hdr = &(pf->hdr);
    int colnum = 0, *status = &(pf->status);

    // See if we need to move to next file
    if (pf->rownum > pf->rows_per_file) {
        printf("Closing file '%s'\n", pf->filename);
        fits_close_file(pf->fptr, status);
        pf->filenum++;
        psrfits_open(pf);
        if (*status==104) {
            printf("Finished with all input files.\n");
            pf->filenum--;
            *status = 1;
            return *status;
        }
    }

    int nivals = hdr->nchan * hdr->npol;
    long long numdata = (long long) nivals * (long long) N;
    long long bytes_to_read = (hdr->nbits * numdata) / 8L;
    float *offsets = (float *)malloc(sizeof(float) * nivals);
    float *scales = (float *)malloc(sizeof(float) * nivals);
    unsigned char *buffer = (unsigned char *)malloc(numdata);
    unsigned char *rawbuffer = buffer;
    if (hdr->nbits==4) {
        rawbuffer = (unsigned char *)malloc(bytes_to_read);
    }

    // Now read the data
    fits_get_colnum(pf->fptr, 0, "DAT_OFFS", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, pf->rownum, 1, nivals,
                  NULL, offsets, NULL, status);
    fits_get_colnum(pf->fptr, 0, "DAT_SCL", &colnum, status);
    fits_read_col(pf->fptr, TFLOAT, colnum, pf->rownum, 1, nivals,
                  NULL, scales, NULL, status);
    fits_get_colnum(pf->fptr, 0, "DATA", &colnum, status);
    fits_read_col(pf->fptr, TBYTE, colnum, pf->rownum, 1, bytes_to_read,
                  NULL, rawbuffer, NULL, status);
    if (hdr->nbits==4) {
        convert_4bit_to_8bit(rawbuffer, (unsigned char *)buffer, numdata);
        free(rawbuffer);
    }
    
    // Now convert the 8-bit data to floats using the scales and offsets
    apply_scales_and_offsets(hdr->nchan, hdr->npol, N,
                             scales, offsets, buffer, fbuffer);
    free(offsets);
    free(scales);
    free(buffer);

    // Complain on error
    fits_report_error(stderr, *status);
    
    return *status;
}
