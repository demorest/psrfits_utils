/* read_psrfits.c */
#include <stdio.h>
#include <string.h>
#include "psrfits.h"

extern void pf_unpack_2bit_to_8bit(struct psrfits *pf, int numunsigned);
extern void pf_unpack_4bit_to_8bit(struct psrfits *pf, int numunsigned);
extern void unpack_2bit_to_8bit_unsigned(unsigned char *indata,
                                         unsigned char *outdata, int N);
extern void unpack_4bit_to_8bit_unsigned(unsigned char *indata,
                                         unsigned char *outdata, int N);

int is_search_PSRFITS(char *filename)
// Return 1 if the file described by filename is a PSRFITS file
// Return 0 otherwise.
{
    fitsfile *fptr;
    int status=0;
    char ctmp[80], comment[120];

    // Read the primary HDU
    fits_open_file(&fptr, filename, READONLY, &status);
    if (status) return 0;

    // Make the easy check first
    fits_read_key(fptr, TSTRING, "FITSTYPE", ctmp, comment, &status);
    if (status || strcmp(ctmp, "PSRFITS")) return 0;

    // See if the data are search-mode
    fits_read_key(fptr, TSTRING, "OBS_MODE", ctmp, comment, &status);
    if (status || (strcmp(ctmp, "SEARCH") &&
                   strcmp(ctmp, "SRCH"))) return 0;

    fits_close_file(fptr, &status);
    return 1;  // it is search-mode PSRFITS
}


void psrfits_set_files(struct psrfits *pf, int numfiles, char *filenames[])
/* Determine whether we are using explicit filenames (and */
/* how many there are, or whether we are building our own */
/* using sequence numbers and the basefilename.           */
{
    int ii;

    // Check if we were passed a basefilename
    if (numfiles==1 && (!is_search_PSRFITS(filenames[0]))) {
        strncpy(pf->basefilename, filenames[0], 200);
        printf("Using '%s' as our base PSRFITS file name.\n", pf->basefilename);
        pf->numfiles = 0;  // dynamically determined
        pf->filenum = 1;
        pf->filename[0] = '\0';
        pf->filenames = NULL;
        return;
    }
    // Using real filenames instead...
    pf->basefilename[0]='\0';
    pf->numfiles = numfiles;
    pf->filenum = 0;
    pf->filenames = filenames;
    // Test that all of the input files are valid PSRFITS
    for (ii = 0; ii < numfiles; ii++) {
        // printf("Checking '%s'...\n", filenames[ii]);
        if (!is_search_PSRFITS(filenames[ii]))
            fprintf(stderr, "Error:  '%s' is not PSRFITS.  Exiting.\n", filenames[ii]);
    }
    printf("Found %d valid PSRFITS files for input.\n", numfiles);
    return;
}


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

    if (pf->numfiles==0) {
        // Dynamically generated file names
        sprintf(pf->filename, "%s_%04d.fits", pf->basefilename, pf->filenum);
    } else {
        // Using explicit filenames
        if (pf->filenum < pf->numfiles) {
            strncpy(pf->filename, pf->filenames[pf->filenum], 200);
        } else {
            *status = FILE_NOT_OPENED;
            return *status;
        }
    }

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
                              int numunsigned, 
                              float *scales, float *offsets,
                              unsigned char *inbuf, float *outbuf)
{
    int ii, jj, poln, N;
    float *outptr = outbuf;

    N = numchan * numpol;
    for (ii = 0 ; ii < numspect ; ii++) {
        for (poln = 0 ; poln < numpol ; poln++) {
            float *sptr = scales + poln * numchan;
            float *optr = offsets + poln * numchan;
            if (poln < numunsigned) {
                unsigned char *inptr = inbuf + ii * N + poln * numchan;
                for (jj = 0 ; jj < numchan ; jj++, sptr++, optr++, inptr++, outptr++)
                    *outptr = *sptr * (float)(*inptr) + *optr;
            } else {
                signed char *inptr = (signed char *)(inbuf + ii * N + poln * numchan);
                for (jj = 0 ; jj < numchan ; jj++, sptr++, optr++, inptr++, outptr++)
                    *outptr = *sptr * (float)(*inptr) + *optr;
            }
        }
    }
}


void scale_and_offset_data(struct psrfits *pf, int numunsigned) {
    // Make sure that pf->sub.fdata has been allocated!
    apply_scales_and_offsets(pf->hdr.nchan, pf->hdr.npol, pf->hdr.nsblk,
                             numunsigned,
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
        if (*status==FILE_NOT_OPENED) {
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
    int numunsigned = hdr->npol;
    if (hdr->npol==4) {
        if (strncmp(hdr->poln_order, "AABBCRCI", 8)==0)
            numunsigned = 2;
        if (strncmp(hdr->poln_order, "IQUV", 4)==0)
            numunsigned = 1;
    }

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
        if (hdr->nbits==2) pf_unpack_2bit_to_8bit(pf, numunsigned);
        else if (hdr->nbits==4) pf_unpack_4bit_to_8bit(pf, numunsigned);
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
int psrfits_read_part_DATA(struct psrfits *pf, int N, int numunsigned, 
                           float *fbuffer) {

    struct hdrinfo *hdr = &(pf->hdr);
    int colnum = 0, *status = &(pf->status);

    // See if we need to move to next file
    if (pf->rownum > pf->rows_per_file) {
        printf("Closing file '%s'\n", pf->filename);
        fits_close_file(pf->fptr, status);
        pf->filenum++;
        psrfits_open(pf);
        if (*status==FILE_NOT_OPENED) {
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
    if (hdr->nbits==4 || hdr->nbits==2) {
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
        unpack_4bit_to_8bit_unsigned(rawbuffer,
                                     (unsigned char *)buffer, numdata);
        free(rawbuffer);
    } else if (hdr->nbits==2) {
        unpack_2bit_to_8bit_unsigned(rawbuffer,
                                     (unsigned char *)buffer, numdata);
        free(rawbuffer);
    }
    
    // Now convert the 8-bit data to floats using the scales and offsets
    apply_scales_and_offsets(hdr->nchan, hdr->npol, N, numunsigned,
                             scales, offsets, buffer, fbuffer);
    free(offsets);
    free(scales);
    free(buffer);

    // Complain on error
    fits_report_error(stderr, *status);
    
    return *status;
}
