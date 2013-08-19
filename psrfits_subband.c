// This code is to partially de-disperse and subband
// PSRFITS search-mode data.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include "psrfits.h"
#include "psrfits_subband_cmd.h"

// This tests to see if 2 times are within 100ns of each other
#define TEST_CLOSE(a, b) (fabs((a)-(b)) <= 1e-7 ? 1 : 0)

extern double delay_from_dm(double dm, double freq_emitted);
extern int split_root_suffix(char *input, char **root, char **suffix);
extern void avg_std(float *x, int n, double *mean, double *std, int stride);
extern void split_path_file(char *input, char **path, char **file);
extern void get_stokes_I(struct psrfits *pf);
extern void downsample_time(struct psrfits *pf);

struct subband_info {
    int nsub;
    int nchan;
    int numunsigned;  // Number of polns using unsigned, latter ones are signed
    int chan_per_sub;
    int npol;
    int max_early;
    int max_late;
    int max_overlap;
    int buflen;  // Number of spectra in time in a block / row
    int bufwid;  // Number of channels times number of polns
    double dm;
    double sub_df;
    double *sub_delays;
    double *chan_delays;
    int *idelays;
    float *sub_freqs;
    float *userwgts;
    float *weights;
    float *offsets;
    float *scales;
    float *chan_avgs;
    float *chan_stds;
    float *fbuffer;
    float *outfbuffer;
};


static void print_percent_complete(int current, int number, int reset)
{
   static int newper = 0, oldper = -1;

   if (reset) {
      oldper = -1;
      newper = 0;
   } else {
      newper = (int) (current / (float) (number) * 100.0);
      if (newper < 0)
         newper = 0;
      if (newper > 100)
         newper = 100;
      if (newper > oldper) {
         printf("\r%3d%% ", newper);
         fflush(stdout);
         oldper = newper;
      }
   }
}

void print_raw_chan_stats(unsigned char *data, int nspect, int nchan, int npol){
    int ii, jj;
    double avg, std;
    float *tmparr;

    tmparr = (float *)malloc(sizeof(float) * nspect);
    for (ii = 0 ; ii < nchan ; ii++) {
        // Grab the raw data
        for (jj = 0 ; jj < nspect ; jj++)
            tmparr[jj] = (float)data[npol * nchan * jj + ii];
        avg_std(tmparr, nspect, &avg, &std, 1);
        printf("%4d: %10.3f %10.3f\n", ii, avg, std);
    }
    free(tmparr);
}


void get_chan_stats(struct psrfits *pfi, struct subband_info *si){   
    int ii;
    double avg, std;

    for (ii = 0 ; ii < si->bufwid ; ii++) {
        // Only use 1/8 of the total length in order to speed things up
        avg_std(pfi->sub.fdata + ii, si->buflen / 8, &avg, &std, si->bufwid);
        //printf("%d %f %f\n", ii, avg, std);
        si->chan_avgs[ii] = avg;
        si->chan_stds[ii] = std;
    }
}


void new_scales_and_offsets(struct psrfits *pfo, int numunsigned) {
    int ii, poln;
    double avg, std;
    const int nspec = pfo->hdr.nsblk / pfo->hdr.ds_time_fact;
    const int nchan = pfo->hdr.nchan / pfo->hdr.ds_freq_fact;
    const int npoln = (pfo->hdr.onlyI) ? 1 : pfo->hdr.npol;
    const int bufwid = npoln * nchan;
    const float target_std = 20.0;  // This is reasonable for 8-bit data

    for (poln = 0 ; poln < npoln ; poln++) {
        float target_avg = (poln < numunsigned) ? 128.0 : 0.0;
        float *out_scls = pfo->sub.dat_scales + poln * nchan;
        float *out_offs = pfo->sub.dat_offsets + poln * nchan;
        for (ii = 0 ; ii < nchan ; ii++) {
            float *fptr = pfo->sub.fdata + poln * nchan + ii;
            avg_std(fptr, nspec, &avg, &std, bufwid);
            out_scls[ii] = std / target_std;
            out_offs[ii] = avg - (target_avg * out_scls[ii]);
        }
    }
}


void new_weights(struct psrfits *inpf, struct psrfits *outpf)
{
    int ii, jj;
    const int out_nchan = inpf->hdr.nchan / outpf->hdr.ds_freq_fact;
    const int N = outpf->hdr.ds_freq_fact;
    float *outwgts = outpf->sub.dat_weights;

    for (ii = 0 ; ii < out_nchan ; ii++) {
        float *inwgts = inpf->sub.dat_weights + ii * N;
        float sumwgt = 0.0;
        for (jj = 0 ; jj < N ; jj++, inwgts++)
            sumwgt += *inwgts;
        // avoid divide by zeros if weights are all zero
        outwgts[ii] = (sumwgt) ? sumwgt / N : 0.0;
    }
}


void fill_chans_with_avgs(int N, int samp_per_spect, float *buffer, float *avgs)
{
    int ii, jj;
    float *fdataptr = buffer;
    for (ii = 0 ; ii < N ; ii++) {
        float *avgptr = avgs;
        for (jj = 0 ; jj < samp_per_spect ; jj++, fdataptr++, avgptr++) {
            *fdataptr = *avgptr;
        }
    }
}

// This routine removes the offsets from the floating point
// data, divides by the scales, and converts (with clipping)
// the resulting value into unsigned chars in the sub.data buffer
void un_scale_and_offset_data(struct psrfits *pf, int numunsigned)
{
    int ii, jj, poln;
    float *inptr = pf->sub.fdata;
    unsigned char *outptr = pf->sub.data;
    const int nspec = pf->hdr.nsblk / pf->hdr.ds_time_fact;
    const int nchan = pf->hdr.nchan / pf->hdr.ds_freq_fact;
    const int npoln = (pf->hdr.onlyI) ? 1 : pf->hdr.npol;

    for (ii = 0 ; ii < nspec ; ii++) {
        for (poln = 0 ; poln < npoln ; poln++) {
            float *sptr = pf->sub.dat_scales + poln * nchan;
            float *optr = pf->sub.dat_offsets + poln * nchan;
            if (poln < numunsigned) {
                for (jj = 0 ; jj < nchan ; jj++, sptr++, optr++, inptr++, outptr++) {
                    //printf("%d  %f  %f  %f\n", jj, *inptr, *sptr, *optr);
                    float ftmp = (*inptr - *optr) / *sptr + 0.5;
                    ftmp = (ftmp >= 255.0) ? 255.0 : ftmp;
                    ftmp = (ftmp < 0.0) ? 0.0 : ftmp;
                    *outptr = (unsigned char) ftmp;
                }
            } else {
                for (jj = 0 ; jj < nchan ; jj++, sptr++, optr++, inptr++, outptr++) {
                    //printf("%d  %f  %f  %f\n", jj, *inptr, *sptr, *optr);
                    float ftmp = (*inptr - *optr) / *sptr + 0.5;
                    ftmp = (ftmp >= 128.0) ? 128.0 : ftmp;
                    ftmp = (ftmp < -127.0) ? -127 : ftmp;
                    *outptr = (signed char) ftmp;
                }
            }
        }
    }
}


int get_current_row(struct psrfits *pfi, struct subband_info *si) {
    static int firsttime = 1, num_pad_blocks = 0;
    static double last_offs, row_duration;
    double diff_offs, dnum_blocks;
    
    if (firsttime) {
        row_duration = pfi->sub.tsubint;
        last_offs = pfi->sub.offs-row_duration;
        firsttime = 0;
    }

    print_percent_complete(pfi->rownum, pfi->rows_per_file, 
                           pfi->rownum==1 ? 1 : 0);

#if 0
    printf("row %d\n", pfi->rownum);
#endif

    if (num_pad_blocks==0) {  // Try to read the PSRFITS file

        // Read the current row of data
        psrfits_read_subint(pfi);
        diff_offs = pfi->sub.offs - last_offs;
        if (si->userwgts) // Always overwrite if using user weights
            memcpy(pfi->sub.dat_weights, si->userwgts,
                   pfi->hdr.nchan * sizeof(float));

        if (!TEST_CLOSE(diff_offs, row_duration) || pfi->status) {
            if (pfi->status) { // End of the files
                num_pad_blocks = 1;
            } else { // Missing row(s)
                dnum_blocks = diff_offs/row_duration - 1.0;
                num_pad_blocks = (int)(dnum_blocks + 1e-7);
                pfi->rownum--;   // Will re-read when no more padding
                pfi->tot_rows--; // Only count "real" rows towards tot_rows
#if 1
                printf("At row %d, found %d dropped rows.\n", 
                       pfi->rownum, num_pad_blocks);
                printf("Adding a missing row (#%d) of padding to the subbands.\n", 
                       pfi->tot_rows);
#endif        
                pfi->N -= pfi->hdr.nsblk;  // Will be re-added below for padding
            }
            // Now fill the main part of si->fbuffer with the chan_avgs so that
            // it acts like a correctly read block (or row)
            fill_chans_with_avgs(si->buflen, si->bufwid,
                                 pfi->sub.fdata, si->chan_avgs);
        } else { // Return the row from the file
            // Compute the float representations of the data
            scale_and_offset_data(pfi, si->numunsigned);
            // Determine channel statistics
            get_chan_stats(pfi, si);
            last_offs = pfi->sub.offs;
            return 0;
        }
    }

    // Return the same padding as before
    last_offs += row_duration;
    pfi->N += pfi->hdr.nsblk;
    pfi->T = pfi->N * pfi->hdr.dt;
    num_pad_blocks--;
    return num_pad_blocks;
}


/* Average adjacent frequency channels, including dispersion, together
 * to make de-dispersed subbands.  It works for all weights, scales,
 * and offsets.
 */
void make_subbands(struct psrfits *pfi, struct subband_info *si) {
    int ii, jj, kk, ll;
    float *weights = pfi->sub.dat_weights;
    float *indata = pfi->sub.fdata;
    float *outdata = si->outfbuffer;
    const int dsfact = si->chan_per_sub;
    const int in_bufwid = si->bufwid;

    // Compute the inv_sumwgts vector for normalizing
    float *inv_sumwgts = (float *)malloc(sizeof(float) * si->nsub);
    for (ii = 0 ; ii < si->nsub ; ii++) {
        float sumwgts = 0.0;
        for (jj = 0 ; jj < dsfact ; jj++, weights++)
            sumwgts += *weights;
        // Take reciprocal if sum is not zero
        inv_sumwgts[ii] = (sumwgts) ? 1.0 / sumwgts : 0.0;
    }

    // Iterate over the times 
    for (ii = 0 ; ii < si->buflen ; ii++) {
        // Iterate over the polarizations
        for (jj = 0 ; jj < si->npol ; jj++) {
            weights = pfi->sub.dat_weights;
            int *idelays = si->idelays;
            float *norm = inv_sumwgts;
            // Iterate over the output channels
            for (kk = 0 ; kk < si->nsub ; kk++, norm++, outdata++) {
                // Iterate over the input channels
                float ftmp = 0.0;
                for (ll = 0 ; ll < dsfact ; ll++, idelays++, indata++, weights++) {
                    ftmp += *weights * *(indata + *idelays * in_bufwid);
                }
                // Now convert the sum to a weighted average
                *outdata = ftmp * *norm;
            }
        }
    }
    free(inv_sumwgts);
}


void init_subbanding(struct psrfits *pfi,
                     struct psrfits *pfo,
                     struct subband_info *si,
                     Cmdline *cmd) {
    int ii, jj, kk, cindex;
    double lofreq, dtmp;
    
    si->nsub = cmd->nsub;
    si->nchan = pfi->hdr.nchan;
    si->npol = pfi->hdr.npol;
    si->numunsigned = si->npol;
    if (si->npol==4) {
        if (strncmp(pfi->hdr.poln_order, "AABBCRCI", 8)==0)
            si->numunsigned = 2;
        if (strncmp(pfi->hdr.poln_order, "IQUV", 4)==0)
            si->numunsigned = 1;
    }
    si->chan_per_sub = si->nchan / si->nsub;
    si->bufwid = si->nchan * si->npol; // Freq * polns
    si->buflen = pfi->hdr.nsblk;  // Number of spectra in each row
    // Check the downsampling factor in time
    if (si->buflen % cmd->dstime) {
        fprintf(stderr,
                "Error!:  %d spectra per row is not evenly divisible by -dstime of %d!\n",
                si->buflen, cmd->dstime);
        exit(1);
    }
    // Check the downsampling factor in frequency
    if (si->nchan % si->nsub) {
        fprintf(stderr, 
                "Error!  %d channels is not evenly divisible by %d subbands!\n",
                si->nchan, si->nsub);
        exit(1);
    }
    si->dm = cmd->dm;
    si->sub_df = pfi->hdr.df * si->chan_per_sub;
    si->sub_freqs = (float *)malloc(sizeof(float) * si->nsub);
    si->chan_delays = (double *)malloc(sizeof(double) * si->nchan);
    si->sub_delays = (double *)malloc(sizeof(double) * si->nsub);
    si->idelays = (int *)malloc(sizeof(int) * si->nchan);
    si->weights = (float *)malloc(sizeof(float) * si->nsub);
    si->offsets = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->scales = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->chan_avgs = (float *)malloc(sizeof(float) * si->bufwid);
    si->chan_stds = (float *)malloc(sizeof(float) * si->bufwid);

    /* Alloc data buffers for the input PSRFITS file */
    pfi->sub.dat_freqs = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_weights = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_offsets = (float *)malloc(sizeof(float)
                                           * pfi->hdr.nchan * pfi->hdr.npol);
    pfi->sub.dat_scales  = (float *)malloc(sizeof(float)
                                           * pfi->hdr.nchan * pfi->hdr.npol);
    pfi->sub.rawdata = (unsigned char *)malloc(pfi->sub.bytes_per_subint);
    if (pfi->hdr.nbits!=8) {
        pfi->sub.data = (unsigned char *)malloc(pfi->sub.bytes_per_subint *
                                                (8 / pfi->hdr.nbits));
    } else {
        pfi->sub.data = pfi->sub.rawdata;
    }

    // Read the first row of data
    psrfits_read_subint(pfi);
    if (si->userwgts) // Always overwrite if using user weights
        memcpy(pfi->sub.dat_weights, si->userwgts, pfi->hdr.nchan * sizeof(float));

    // Reset the read counters since we'll re-read
    pfi->rownum--;
    pfi->tot_rows--;
    pfi->N -= pfi->hdr.nsblk;

    // Compute the subband properties, DM delays and offsets
    lofreq = pfi->sub.dat_freqs[0] - pfi->hdr.df * 0.5;
    for (ii = 0, cindex = 0 ; ii < si->nsub ; ii++) {
        dtmp = lofreq + ((double)ii + 0.5) * si->sub_df;
        si->sub_freqs[ii] = dtmp;
        si->sub_delays[ii] = delay_from_dm(si->dm, dtmp);
        // Determine the dispersion delays and convert them
        // to offsets in units of sample times
        for (jj = 0 ; jj < si->chan_per_sub ; jj++, cindex++) {
            si->chan_delays[cindex] = delay_from_dm(si->dm, 
                                                    pfi->sub.dat_freqs[cindex]);
            si->chan_delays[cindex] -= si->sub_delays[ii];
            si->idelays[cindex] = (int)rint(si->chan_delays[cindex] / pfi->hdr.dt);
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

    // This buffer will hold the float-converted input data, plus the bits 
    // of data from the previous and next blocks
    si->fbuffer = (float *)calloc((si->buflen + 2 * si->max_overlap) *
                                  si->bufwid, sizeof(float));
    // The input data will be stored directly in the buffer space
    // So the following is really just an offset into the bigger buffer
    pfi->sub.fdata = si->fbuffer + si->max_overlap * si->bufwid;
    
    // Now start setting values for the output arrays
    *pfo = *pfi;
    // Determine the length of the outputfiles to use
    if (cmd->filetimeP) {
        pfo->rows_per_file = 10 * \
            (int) rint(0.1 * (cmd->filetime / pfi->sub.tsubint));
    } else if (cmd->filelenP) {        
        long long filelen;
        int bytes_per_subint;
        filelen = cmd->filelen * (1L<<30);  // In GB
        bytes_per_subint = (pfo->hdr.nbits * pfo->hdr.nchan * 
                            pfo->hdr.npol * pfo->hdr.nsblk) / \
            (8 * si->chan_per_sub * cmd->dstime * (cmd->onlyIP ? 4 : 1));
        pfo->rows_per_file = filelen / bytes_per_subint;
    } else {  // By default, keep the filesize roughly constant
        pfo->rows_per_file = pfi->rows_per_file * si->chan_per_sub * 
            cmd->dstime * (cmd->onlyIP ? 4 : 1);
    }
    pfo->filenum = 0; // This causes the output files to be created
    pfo->filename[0] = '\0';
    pfo->rownum = 1;
    pfo->tot_rows = 0;
    pfo->N = 0;
    // Set the "orig" values to those of the input file
    pfo->hdr.orig_nchan = pfi->hdr.nchan;
    pfo->hdr.orig_df = pfi->hdr.df;
    {
        char *inpath, *infile;
        split_path_file(pfi->basefilename, &inpath, &infile);
        sprintf(pfo->basefilename, "%s_subs", infile);
        free(inpath);
        free(infile);
    }
    // Reset different params
    pfo->sub.dat_freqs = si->sub_freqs;
    pfo->sub.dat_weights = si->weights;
    pfo->sub.dat_offsets = si->offsets;
    pfo->sub.dat_scales  = si->scales;
    pfo->hdr.ds_freq_fact = si->chan_per_sub;
    pfo->hdr.ds_time_fact = cmd->dstime;
    pfo->hdr.onlyI = cmd->onlyIP;
    pfo->hdr.chan_dm = si->dm;
    pfo->sub.rawdata = (unsigned char *)malloc(si->nsub * si->npol * si->buflen);
    if (pfo->hdr.nbits!=8) {
        pfo->sub.data = (unsigned char *)malloc(si->nsub * si->npol * si->buflen *
                                                (8 / pfi->hdr.nbits));
    } else {
        pfo->sub.data = pfo->sub.rawdata;
    }
    si->outfbuffer = (float *)calloc(si->nsub * si->npol * si->buflen,
                                     sizeof(float));
    pfo->sub.fdata = si->outfbuffer;

    // Now re-read the first row (i.e. for "real" this time)
    get_current_row(pfi, si);

    // Set the new weights properly
    new_weights(pfi, pfo);

    // Now fill the first part of si->fbuffer with the chan_avgs so that
    // it acts like a previously read block (or row)
    fill_chans_with_avgs(si->max_overlap, si->bufwid,
                         si->fbuffer, si->chan_avgs);
}


void read_weights(char *filenm, int *numchan, float **weights)
{
    FILE *infile;
    int ii, N, chan;
    float wgt;
    char line[80];

    infile = fopen(filenm, "r");

    // Read the input file once to count the lines
    N = 0;
    while (!feof(infile)){
        fgets(line, 80, infile);
        if (line[0]!='#') {
            sscanf(line, "%d %f\n", &chan, &wgt);
            N++;
        }
    }
    N--;
    *numchan = N;

    // Allocate the output arrays
    *weights = (float *)malloc(N * sizeof(float));

    // Rewind and read the EVENTs for real
    rewind(infile);
    ii = 0;
    while (ii < N) {
        fgets(line, 80, infile);
        if (line[0]!='#') {
            sscanf(line, "%d %f\n", &chan, (*weights)+ii);
            ii++;
        }
    }
    fclose(infile);
}


int main(int argc, char *argv[]) {
    Cmdline *cmd;
    struct psrfits pfi, pfo;
    struct subband_info si;
    int stat=0, padding=0, userN=0;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);

    // Open the input PSRFITs files
    psrfits_set_files(&pfi, cmd->argc, cmd->argv);
    // Use the dynamic filename allocation
    if (pfi.numfiles==0) pfi.filenum = cmd->startfile;
    pfi.tot_rows = pfi.N = pfi.T = pfi.status = 0;
    int rv = psrfits_open(&pfi);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    // Read the user weights if requested
    si.userwgts = NULL;
    if (cmd->wgtsfileP) {
        read_weights(cmd->wgtsfile, &userN, &si.userwgts);
        if (userN != pfi.hdr.nchan) {
            printf("Error!:  Input data has %d channels, but '%s' contains only %d weights!\n",
                   pfi.hdr.nchan, cmd->wgtsfile, userN);
            exit(0);
        }
        printf("Overriding input channel weights with those in '%s'\n",
               cmd->wgtsfile);
    }

    // Initialize the subbanding
    // (including reading the first row of data and
    //  putting it in si->fbuffer)
    init_subbanding(&pfi, &pfo, &si, cmd);

    if (cmd->outputbasenameP)
      strcpy(pfo.basefilename, cmd->outputbasename);

    // Loop through the data
    do {
        // Put the overlapping parts from the next block into si->buffer
        float *ptr = pfi.sub.fdata + si.buflen * si.bufwid;
        if (padding==0)
            stat = psrfits_read_part_DATA(&pfi, si.max_overlap, si.numunsigned, ptr);
        if (stat || padding) { // Need to use padding since we ran out of data
            printf("Adding a missing row (#%d) of padding to the subbands.\n",
                   pfi.tot_rows);
            // Now fill the last part of si->fbuffer with the chan_avgs so that
            // it acts like a correctly read block (or row)
            fill_chans_with_avgs(si.max_overlap, si.bufwid,
                                 ptr, si.chan_avgs);
        }
        //print_raw_chan_stats(pfi.sub.data, pfi.hdr.nsblk, 
        //                     pfi.hdr.nchan, pfi.hdr.npol);

        // Now create the subbanded row in the output buffer
        make_subbands(&pfi, &si);

        // Output only Stokes I (in place via floats)
        if (pfo.hdr.onlyI && pfo.hdr.npol==4)
            get_stokes_I(&pfo);

        // Downsample in time (in place via floats)
        if (pfo.hdr.ds_time_fact > 1)
            downsample_time(&pfo);

        // Compute new scales and offsets so that we can pack
        // into 8-bits reliably
        new_scales_and_offsets(&pfo, si.numunsigned);

        // Convert the floats back to bytes in the output array
        un_scale_and_offset_data(&pfo, si.numunsigned);
        //print_raw_chan_stats(pfo.sub.data, pfo.hdr.nsblk / pfo.hdr.ds_time_fact,  
        //                     pfo.hdr.nchan / pfo.hdr.ds_freq_fact, pfo.hdr.npol);

        // Write the new row to the output file
        pfo.sub.offs = (pfo.tot_rows+0.5) * pfo.sub.tsubint;
        psrfits_write_subint(&pfo);

        // Break out of the loop here if stat is set
        if (stat) break;

        // shift the last part of the current row into the "last-row" 
        // part of the data buffer
        memcpy(si.fbuffer, si.fbuffer + si.buflen * si.bufwid,
               si.max_overlap * si.bufwid * sizeof(float));

        // Read the next row (or padding)
        padding = get_current_row(&pfi, &si);

        // Set the new weights properly
        new_weights(&pfi, &pfo);

    } while (pfi.status == 0);

    rv = psrfits_close(&pfi);
    if (rv>100) { fits_report_error(stderr, rv); }
    rv = psrfits_close(&pfo);
    if (rv>100) { fits_report_error(stderr, rv); }
    exit(0);
}
