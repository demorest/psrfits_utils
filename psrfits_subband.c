// This code is to partially de-disperse and subband
// PSRFITS search-mode data.  Currently it is specifically
// for GUPPI data, however, I intend to make it more general
// eventually.   S. Ransom  Oct 2008
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include "psrfits.h"
#include "psrfits_subband_cmd.h"

// This tests to see if 2 times are within 100ns of each other
#define TEST_CLOSE(a, b) (fabs((a)-(b)) <= 1e-7 ? 1 : 0)

#define DEBUG 0

extern double delay_from_dm(double dm, double freq_emitted);
extern int split_root_suffix(char *input, char **root, char **suffix);
extern void avg_std(char *x, int n, double *mean, double *std);
extern short transpose_bytes(unsigned char *a, int nx, int ny, 
                             unsigned char *move, int move_size);
extern void downsample_freq(struct psrfits *pf);

struct subband_info {
    int nsub;
    int nchan;
    int chan_per_sub;
    int npol;
    int max_early;
    int max_late;
    int max_overlap;
    int buflen;  // Number of samples (in time) in a block
    int bufwid;  // Number of channels times number of polns
    int Tbufsize;
    double dm;
    double sub_df;
    double *sub_delays;
    double *chan_delays;
    int *idelays;
    float *sub_freqs;
    float *weights;
    float *offsets;
    float *scales;
    float *chan_avgs;
    float *chan_stds;
    unsigned char *buffer;
    unsigned char *tmpbuf;
    unsigned char *outbuffer;
    unsigned char *Tbuf;
};

void get_chan_stats(struct psrfits *pfi, 
                    struct subband_info *si) {
    int ii, offset=0;
    double avg, std;

    for (ii = 0 ; ii < si->bufwid ; ii++) {
        avg_std((char *)pfi->sub.data+offset, si->buflen/8, &avg, &std);
        si->chan_avgs[ii] = avg;
        si->chan_stds[ii] = std;
        offset += si->buflen;
#ifdef DEBUG1
        printf("xxx %d %f %f\n", ii, avg, std);
#endif        
    }
}

void get_sub_stats(struct psrfits *pfi, 
                   struct subband_info *si) {
    int ii, offset=0;
    double avg, std;

    for (ii = 0 ; ii < si->bufwid/si->chan_per_sub ; ii++) {
        avg_std((char *)pfi->sub.data+offset, si->buflen, &avg, &std);
        offset += si->buflen;
#ifdef DEBUG
        printf("xxx %d %f %f\n", ii, avg, std);
#endif        
    }
}

void get_next_row(struct psrfits *pfi, 
                  struct subband_info *si) {
    static int firsttime = 1, num_pad_blocks = 0;
    static double last_offs, row_duration;
    double diff_offs, dnum_blocks;
    int ii;
    
    if (firsttime) {
        row_duration = pfi->sub.tsubint;
        last_offs = pfi->sub.offs-row_duration;
        firsttime = 0;
    }

    if (num_pad_blocks==0) {  // Try to read the PSRFITS file

        // Read the first row of data
        psrfits_read_subint(pfi);
        diff_offs = pfi->sub.offs - last_offs;

        if (!TEST_CLOSE(diff_offs, row_duration) || 
            pfi->status) { // Missing rows
            if (!pfi->status) {
                dnum_blocks = diff_offs/row_duration - 1.0;
                num_pad_blocks = (int)(dnum_blocks + 1e-12);
#ifdef DEBUG
                printf("At row %d, found %d dropped rows.\n", 
                       pfi->rownum, num_pad_blocks);
#endif        
                pfi->rownum--;   // Will re-read when no more padding
                pfi->tot_rows--; // Only count "real" rows towards tot_rows
                pfi->N -= pfi->hdr.nsblk;  // Will be re-added below for padding
            } else { // End of the files
                num_pad_blocks = 1;
            }
            // Fill the data buffer with with the chan_avgs so that it
            // acts like a previously read block (or row)
            for (ii = 0 ; ii < si->bufwid ; ii++) {
                unsigned char *tind = pfi->sub.data + ii * si->buflen;
                memset((char *)tind, (int)rint(si->chan_avgs[ii]), si->buflen);
            }

        } else { // Return the row from the file

            // Transpose the data so that each channel is together in RAM
            transpose_bytes(pfi->sub.data, si->buflen, si->bufwid, 
                            si->Tbuf, si->Tbufsize);
            // Determine channel statistics
            get_chan_stats(pfi, si);

            last_offs = pfi->sub.offs;
            return;
        }
    }

    // Return the padding from before
    last_offs += row_duration;
    pfi->N += pfi->hdr.nsblk;
    pfi->T = pfi->N * pfi->hdr.dt;
    num_pad_blocks--;
}


void convert_next_to_current(struct psrfits *pfi, 
                             struct subband_info *si) {
    int ii;

    // Prep and shift the data in the output buffer
    for (ii = 0 ; ii < si->bufwid ; ii++) {
        int absdelay = abs(si->idelays[ii]);
        int num = si->buflen - absdelay;
        int ioffset = 0;
        int ooffset = 0;
        unsigned char *iind = pfi->sub.data + ii * si->buflen;
        unsigned char *oind = si->buffer + ii * si->buflen;
        unsigned char *tind = si->tmpbuf + ii * si->max_overlap;
        if (si->idelays[ii] < 0) {
            // negative delays arrived before the subband center, 
            // so we get them from tmpbuf
            ooffset = absdelay;
            memcpy(oind, tind, absdelay);
            // Put the stuff we need to buffer into the tmpbuf
            memcpy(tind, iind+num, absdelay);
        } else {
            // positive delays arrive later than the subband center,
            // so some of them come from the next block
            ioffset = absdelay;
            // Put the stuff from the new block into buffer
            memcpy(oind+num, iind, absdelay);
        }
        // Now place the data we have access to into buffer
        memcpy(oind+ooffset, iind+ioffset, num);
    }
}


void complete_current_from_next(struct psrfits *pfi, 
                                struct subband_info *si) {
    int ii;
    
    // Prep and shift the data in the output buffer
    for (ii = 0 ; ii < si->bufwid ; ii++) {
        int absdelay = abs(si->idelays[ii]);
        int ooffset = 0;
        unsigned char *iind = pfi->sub.data + ii * si->buflen;
        unsigned char *oind = si->buffer + ii * si->buflen;
        if (si->idelays[ii] > 0) {
            // positive delays arrive later than the subband center,
            // so some of them come from the next block
            ooffset = si->buflen - absdelay;
            // Put the stuff from the new block into buffer
            memcpy(oind+ooffset, iind, absdelay);
        }
    }
}


void init_subbanding(int nsub, double dm, 
                     struct psrfits *pfi, 
                     struct subband_info *si) {
    int ii, jj, kk, cindex;
    double lofreq, dtmp;
    
    si->nsub = nsub;
    si->nchan = pfi->hdr.nchan;
    si->npol = pfi->hdr.npol;
    si->chan_per_sub = si->nchan / si->nsub;
    si->bufwid = si->nchan * si->npol; // Freq * polns
    si->buflen = pfi->hdr.nsblk;  // Time
    if (si->nchan % si->nsub) {
        fprintf(stderr, 
                "Error!  %d channels is not evenly divisible by %d subbands!\n", 
                si->nchan, si->nsub);
        exit(1);
    }
    si->dm = dm;
    si->sub_df = pfi->hdr.orig_df * si->chan_per_sub;
    si->sub_freqs = (float *)malloc(sizeof(float) * si->nsub);
    si->chan_delays = (double *)malloc(sizeof(double) * si->nchan);
    si->sub_delays = (double *)malloc(sizeof(double) * si->nsub);
    si->idelays = (int *)malloc(sizeof(int) * si->nchan * si->npol);
    si->weights = (float *)malloc(sizeof(float) * si->nsub);
    si->offsets = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->scales = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->chan_avgs = (float *)malloc(sizeof(float) * si->bufwid);
    si->chan_stds = (float *)malloc(sizeof(float) * si->bufwid);

    // Following should be an upper limit (if storing data as floats)
    //si->buffer = (unsigned char *)malloc(sizeof(float) * 
    //                                     si->buflen * si->bufwid);
    //si->outbuffer = (unsigned char *)malloc(sizeof(float) * 
    //                                        si->nsub * si->npol * si->buflen);
    si->buffer = (unsigned char *)calloc(si->buflen * si->bufwid, 
                                         sizeof(unsigned char));
    // Will need the following for out-of-place subbanding
    //si->outbuffer = (unsigned char *)calloc(si->nsub * si->npol * si->buflen, 
    //                                        sizeof(unsigned char));
    // For in-place subbanding...
    si->outbuffer = si->buffer;
    si->Tbufsize = si->bufwid * si->buflen / 2;
    si->Tbuf = (unsigned char *)malloc(sizeof(unsigned char) * si->Tbufsize);

    /* Alloc data buffers for the input PSRFITS file */
    pfi->sub.dat_freqs = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_weights = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_offsets = (float *)malloc(sizeof(float)
                                          * pfi->hdr.nchan * pfi->hdr.npol);
    pfi->sub.dat_scales  = (float *)malloc(sizeof(float)
                                          * pfi->hdr.nchan * pfi->hdr.npol);
    pfi->sub.data = (unsigned char *)malloc(pfi->sub.bytes_per_subint);
        
    // Read the first row of data
    psrfits_read_subint(pfi);
    // Reset the read counters since we'll re-read
    pfi->rownum--;
    pfi->tot_rows--;
    pfi->N -= pfi->hdr.nsblk;

    // Compute the subband properties, DM delays and offsets
    lofreq = pfi->sub.dat_freqs[0] - pfi->hdr.orig_df * 0.5;
    for (ii = 0, cindex = 0 ; ii < si->nsub ; ii++) {
        dtmp = lofreq + ((double)ii + 0.5) * si->sub_df;
        si->sub_freqs[ii] = dtmp;
        si->sub_delays[ii] = delay_from_dm(si->dm, dtmp);
#ifdef DEBUG1
        printf("%4d  %.8f  %.4f\n", ii, dtmp, si->sub_delays[ii]*1e6);
#endif
// These need fixed based on input data
        si->weights[ii] = 1.0;
        for (jj = 0 ; jj < si->npol ; jj++) {
            si->offsets[jj*si->nsub+ii] = 0.0;
            si->scales[jj*si->nsub+ii] = 1.0;
        }
        // Determine the dispersion delays and convert them
        // to offsets in units of sample times
        for (jj = 0 ; jj < si->chan_per_sub ; jj++, cindex++) {
            si->chan_delays[cindex] = delay_from_dm(si->dm, 
                                                    pfi->sub.dat_freqs[cindex]);
            si->chan_delays[cindex] -= si->sub_delays[ii];
            si->idelays[cindex] = (int)rint(si->chan_delays[cindex] / pfi->hdr.dt);
            // Copy the delays if we have more than 1 poln
            for (kk = 1 ; kk > si->npol ; kk++)
                si->idelays[kk*si->nchan+cindex] = si->idelays[cindex];
                
        }
    }

    // Now determine the earliest and latest delays
    si->max_early = si->max_late = 0;
    for (ii = 0 ; ii < si->nchan ; ii++) {
        if (si->idelays[ii] < si->max_early)
            si->max_early = si->idelays[ii];
        if (si->idelays[ii] > si->max_late)
            si->max_late = si->idelays[ii];
    }
    si->max_overlap = abs(si->max_early) + si->max_late;
#ifdef DEBUG1
    printf("max_early   = %d\n", si->max_early);
    printf("max_late    = %d\n", si->max_late);
    printf("max_overlap = %d\n", si->max_overlap);
#endif
    si->tmpbuf = (unsigned char *)malloc(sizeof(float) * 
                                         si->max_overlap * si->bufwid);
    
    // re-read the first row (i.e. for "real" this time)
    get_next_row(pfi, si);
    
    // Now fill tmpbuf with the chan_avgs so that it acts like
    // a previously read block (or row)
    for (ii = 0 ; ii < si->bufwid ; ii++) {
        unsigned char *tind = si->tmpbuf + ii * si->max_overlap;
        memset((char *)tind, (int)rint(si->chan_avgs[ii]), si->max_overlap);
    }

    // Now shift most of the data into the output buffer
    convert_next_to_current(pfi, si);
}

void set_output_vals(struct psrfits *pfi, 
                     struct psrfits *pfo, 
                     struct subband_info *si) {
    // Copy everything
    *pfo = *pfi;
    pfo->filenum = 0; // This causes the output files to be created
    pfo->filename[0] = '\0';
    pfo->rownum = 1;
    pfo->tot_rows = 0;
    pfo->N = 0;
    sprintf(pfo->basefilename, "%s_subs", pfi->basefilename);
    // Reset different params
    pfo->sub.dat_freqs = si->sub_freqs;
    pfo->sub.dat_weights = si->weights;
    pfo->sub.dat_offsets = si->offsets;
    pfo->sub.dat_scales  = si->scales;
    //pfo->hdr.nbits = sizeof(short)*8;
    //pfo->sub.FITS_typecode = TSHORT;
    //pfo->hdr.nbits = 8;
    //pfo->sub.FITS_typecode = TBYTE;
    // For in-place subbanding using the code in downsample.c
    pfo->hdr.ds_freq_fact = si->chan_per_sub;
    pfo->hdr.chan_dm = si->dm;
    //pfo->hdr.nchan = si->nsub;
    //pfo->hdr.npol = si->npol;
    //pfo->sub.bytes_per_subint = (pfo->hdr.nbits * pfo->hdr.nchan * 
    //                             pfo->hdr.npol * pfo->hdr.nsblk) / 8;
    pfo->sub.data = si->outbuffer;
}


int main(int argc, char *argv[]) {
    Cmdline *cmd;
    struct psrfits pfi, pfo;
    struct subband_info si;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);

    // Open the input PSRFITs files
    pfi.tot_rows = pfi.N = pfi.T = pfi.status = 0;
    pfi.filenum = 1;
    pfi.filename[0] = '\0';
    sprintf(pfi.basefilename, cmd->argv[0]);
    int rv = psrfits_open(&pfi);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    // Initialize the subbanding
    // (including reading the first row of data and
    //  putting it in si->buffer and si->tmpbuf)
    init_subbanding(cmd->nsub, cmd->dm, &pfi, &si);

    // Update the output PSRFITS structure
    set_output_vals(&pfi, &pfo, &si);

    // Loop through the data
    do {
        // Read the next row (or padding)
        get_next_row(&pfi, &si);
        
        // Put the overlapping parts from it into si->buffer
        complete_current_from_next(&pfi, &si);
        
        // Sum si->buffer to create subbands
        // must un-transpose and then sum
        transpose_bytes(si.buffer, si.bufwid, si.buflen, 
                        si.Tbuf, si.Tbufsize);
        downsample_freq(&pfo);
        //printf("Input row: %d   Output row: %d\n", pfi.rownum, pfo.rownum);
        //transpose_bytes(si.buffer, si.buflen, si.bufwid/si.chan_per_sub, 
        //                si.Tbuf, si.Tbufsize);
        //get_sub_stats(&pfo, &si);

        // Write the new row to the output file
        // Need to update the other params as well
        pfo.sub.offs = ((pfo.rownum-1)+0.5) * pfo.sub.tsubint;
        psrfits_write_subint(&pfo);

        // Shift the psr->sub.data to si->buffer and si->tmpbuf
        convert_next_to_current(&pfi, &si);

    } while (pfi.status == 0);

    psrfits_close(&pfi);
    if (rv) { fits_report_error(stderr, rv); }
    psrfits_close(&pfo);
    if (rv) { fits_report_error(stderr, rv); }
    exit(0);
}
