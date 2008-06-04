/* fold_psrfits.c
 *
 * Fold PSRFITS search data into PSRFITS folded format.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <getopt.h>
#include <signal.h>
#include <pthread.h>
#include "polyco.h"
#include "fold.h"
#include "psrfits.h"

/* Signal handler */
int run=1;
void cc(int sig) { run=0; }

void usage() {
    printf(
            "Usage: fold_psrfits [options] input_filename_base\n"
            "Options:\n"
            "  -h, --help               Print this\n"
            "  -o name, --output=name   Output base filename\n"
            "  -b nn, --nbin=nn         Number of profile bins (256)\n"
            "  -t nn, --tsub=n          Folded subintegration time, sec (60)\n"
            "  -j nn, --nthread=nn      Max number of threads (4)\n"
            "  -i nn, --initial=nn      Starting input file number\n"
            "  -f nn, --final=nn        Ending input file number\n"
            "  -s src, --src=src        Override source name from file\n"
            "  -q, --quiet              No progress indicator\n"
          );
}

int main(int argc, char *argv[]) {

    /* Cmd line */
    static struct option long_opts[] = {
        {"output",  1, NULL, 'o'},
        {"nbin",    1, NULL, 'b'},
        {"tsub",    1, NULL, 't'},
        {"nthread", 1, NULL, 'j'},
        {"initial", 1, NULL, 'i'},
        {"final",   1, NULL, 'f'},
        {"source",  1, NULL, 's'},
        {"quiet",   0, NULL, 'q'},
        {"help",    0, NULL, 'h'},
        {0,0,0,0}
    };
    int opt, opti;
    int nbin=256, nthread=4, fnum_start=1, fnum_end=0;
    int quiet=0;
    double tfold = 60.0; 
    char output_base[256] = "fold_out";
    char source[24];  source[0]='\0';
    while ((opt=getopt_long(argc,argv,"o:b:t:j:i:f:s:qh",long_opts,&opti))!=-1) {
        switch (opt) {
            case 'o':
                strncpy(output_base, optarg, 255);
                output_base[255]='\0';
                break;
            case 'b':
                nbin = atoi(optarg);
                break;
            case 't':
                tfold = atof(optarg);
                break;
            case 'j':
                nthread = atoi(optarg);
                break;
            case 'i':
                fnum_start = atoi(optarg);
                break;
            case 'f':
                fnum_end = atoi(optarg);
                break;
            case 's':
                strncpy(source, optarg, 24);
                source[23]='\0';
                break;
            case 'q':
                quiet=1;
                break;
            case 'h':
            default:
                usage();
                exit(0);
                break;
        }

    }
    if (optind==argc) { 
        usage();
        exit(1);
    }

    /* Open first file */
    struct psrfits pf;
    sprintf(pf.basefilename, argv[optind]);
    pf.filenum = fnum_start;
    pf.status = 0;
    int rv = psrfits_open(&pf);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    /* Check any constraints */
    if (pf.hdr.nbits!=8) { 
        fprintf(stderr, "Only implemented for 8-bit data (read nbits=%d).\n",
                pf.hdr.nbits);
        exit(1);
    }

    /* Set up output file */
    struct psrfits pf_out;
    memcpy(&pf_out, &pf, sizeof(struct psrfits));
    sprintf(pf_out.basefilename, output_base);
    sprintf(pf_out.hdr.obs_mode, "PSR");
    if (source[0]!='\0') { strncpy(pf_out.hdr.source, source, 24); }
    pf_out.fptr = NULL;
    pf_out.filenum=0;
    pf_out.status=0;
    pf_out.hdr.nbin=nbin;
    pf_out.sub.FITS_typecode = TFLOAT;
    pf_out.sub.bytes_per_subint = 
        pf_out.hdr.nchan * pf_out.hdr.npol * pf_out.hdr.nbin;
    rv = psrfits_create(&pf_out);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    /* Alloc data buffers */
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf_out.sub.dat_freqs = pf.sub.dat_freqs;
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf_out.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) 
            * pf.hdr.nchan * pf.hdr.npol);
    pf_out.sub.dat_offsets = (float *)malloc(sizeof(float) 
            * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales  = (float *)malloc(sizeof(float) 
            * pf.hdr.nchan * pf.hdr.npol);
    pf_out.sub.dat_scales  = (float *)malloc(sizeof(float) 
            * pf.hdr.nchan * pf.hdr.npol);
    /* Now see later thread management section for input data buffer */
    //pf.sub.data  = (unsigned char *)malloc(pf.sub.bytes_per_subint);
    pf_out.sub.data  = (unsigned char *)malloc(sizeof(float) 
            * pf_out.sub.bytes_per_subint);

    /* Output scale/offset */
    int i;
    for (i=0; i<pf.hdr.nchan*pf.hdr.npol; i++) {
        pf_out.sub.dat_scales[i] = 1.0;
        pf_out.sub.dat_offsets[i] = 0.0;
    }
    for (i=0; i<pf.hdr.nchan; i++) { pf_out.sub.dat_weights[i]=1.0; }

    /* Read polycos */
    FILE *pcfile = fopen("polyco.dat", "r");
    if (pcfile==NULL) { 
        fprintf(stderr, "Couldn't open polyco file.\n");
        exit(1);
    }
    int npc=0, ipc=0;
    struct polyco *pc = NULL;
    do { 
        pc = (struct polyco *)realloc(pc, sizeof(struct polyco) * (npc+1));
        rv = read_one_pc(pcfile, &pc[npc]);
        npc++;
    } while (rv==0); 
    npc--; // Final "read" is really a error or EOF.
    if (npc==0) {
        fprintf(stderr, "Error parsing polyco file.\n");
        exit(1);
    }
    fclose(pcfile);

    /* For now, just write all polycos */
    rv = psrfits_write_polycos(&pf_out, pc, npc);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    /* Alloc total fold buf */
    struct foldbuf fb;
    fb.nchan = pf.hdr.nchan;
    fb.npol = pf.hdr.npol;
    fb.nbin = pf_out.hdr.nbin;
    malloc_foldbuf(&fb);
    clear_foldbuf(&fb);

    /* Set up thread management */
    pthread_t *thread_id;
    struct fold_args *fargs;
    thread_id = (pthread_t *)malloc(sizeof(pthread_t) * nthread);
    fargs = (struct fold_args *)malloc(sizeof(struct fold_args) * nthread);
    for (i=0; i<nthread; i++) { 
        thread_id[i] = 0; 
        fargs[i].data = (char *)malloc(sizeof(char)*pf.sub.bytes_per_subint);
        fargs[i].fb = (struct foldbuf *)malloc(sizeof(struct foldbuf));
        fargs[i].fb->nbin = pf_out.hdr.nbin;
        fargs[i].fb->nchan = pf.hdr.nchan;
        fargs[i].fb->npol = pf.hdr.npol;
        fargs[i].nsamp = pf.hdr.nsblk;
        fargs[i].tsamp = pf.hdr.dt;
        malloc_foldbuf(fargs[i].fb);
        clear_foldbuf(fargs[i].fb);
    }

    /* Main loop */
    rv=0;
    int imjd;
    double fmjd, fmjd0=0, fmjd_next=0, fmjd_epoch;
    double offs0=0, offs1=0;
    //double phase=0.0, freq=1.0;
    int first=1, subcount=0;
    int cur_thread = 0;
    signal(SIGINT, cc);
    while (run) { 

        /* Read data block */
        pf.sub.data = (unsigned char *)fargs[cur_thread].data;
        rv = psrfits_read_subint(&pf);
        if (rv) { run=0; break; }

        /* If we've passed final file, exit */
        if (fnum_end && pf.filenum>fnum_end) { run=0; break; }

        /* Get start date, etc */
        imjd = (int)pf.hdr.MJD_epoch;
        fmjd = (double)(pf.hdr.MJD_epoch - (long double)imjd);
        fmjd += (pf.sub.offs-0.5*pf.sub.tsubint)/86400.0;

        /* First time stuff */
        if (first) {
            fmjd0 = fmjd;
            fmjd_next = fmjd + tfold/86400.0;
            pf_out.sub.offs=0.0;
            offs0 = pf.sub.offs - 0.5*pf.sub.tsubint;
            offs1 = pf.sub.offs + 0.5*pf.sub.tsubint;
            first=0;
        }

        /* Keep track of timestamp */
        // TODO also pointing stuff.
        pf_out.sub.offs += pf.sub.offs;
        subcount++;

        /* Update output block end time */
        offs1 = pf.sub.offs + 0.5*pf.sub.tsubint;

        /* Select polyco set */
        ipc = select_pc(pc, npc, pf_out.hdr.source, imjd, fmjd);
        if (ipc<0) { 
            fprintf(stderr, "No matching polycos (src=%s, imjd=%d, fmjd=%f)\n",
                    pf_out.hdr.source, imjd, fmjd);
            break;
        }

        /* TODO: deal with scale/offset */

        /* Fold this subint */
        fargs[cur_thread].pc = &pc[ipc];
        fargs[cur_thread].imjd = imjd;
        fargs[cur_thread].fmjd = fmjd;
        rv = pthread_create(&thread_id[cur_thread], NULL, 
                fold_8bit_power_thread, &fargs[cur_thread]);
        if (rv) {
            fprintf(stderr, "Thread creation error.\n");
            exit(1);
        }
        cur_thread++;

        /* Combine thread results if needed */
        if (cur_thread==nthread || fmjd>fmjd_next) {

            /* Loop over active threads */
            for (i=0; i<cur_thread; i++) {

                /* Wait for thread to finish */
                rv = pthread_join(thread_id[i], NULL);
                if (rv) { 
                    fprintf(stderr, "Thread join error.\n");
                    exit(1);
                }

                /* Combine its result into total fold */
                accumulate_folds(&fb, fargs[i].fb);

                /* Reset thread info */
                clear_foldbuf(fargs[i].fb);
                thread_id[i] = 0;

            }

            /* Reset active thread count */
            cur_thread = 0;
        }

        /* See if integration needs to be written, etc */
        if (fmjd > fmjd_next) {

            /* Figure out timestamp */
            pf_out.sub.offs /= (double)subcount;
            pf_out.sub.tsubint = offs1 - offs0;
            fmjd_epoch = fmjd0 + pf_out.sub.offs/86400.0;
            /*
            // Don't need this stuff if we set EPOCHS=MIDTIME
            ipc = select_pc(pc, npc, pf.hdr.source, imjd, fmjd_epoch); 
            if (ipc<0) { 
                fprintf(stderr, "Polyco error, exiting.\n");
                exit(1);
            }
            phase = psr_phase(&pc[ipc], imjd, fmjd_epoch, &freq);
            phase = fmod(phase, 1.0);
            pf_out.sub.offs -= phase/freq; // ref epoch needs 0 phase
            */

            /* Transpose, output subint */
            normalize_transpose_folds((float *)pf_out.sub.data, &fb);
            psrfits_write_subint(&pf_out);

            /* Clear counters, avgs */
            clear_foldbuf(&fb);
            pf_out.sub.offs = 0.0;
            offs0 = pf.sub.offs - 0.5*pf.sub.tsubint;
            subcount=0;

            /* Set next output time */
            fmjd_next = fmjd + tfold/86400.0;
            if (!quiet) printf("\rWrote subint   \n");
        }


        /* Progress report */
        if (!quiet) {
            printf("\rFile %d %5.1f%%", pf.filenum, 
                    100.0 * (float)(pf.rownum-1)/(float)pf.rows_per_file);
            fflush(stdout);
        }
    }

    /* Join any running threads */
    for (i=0; i<cur_thread; i++)  
        if (thread_id[i]) pthread_join(thread_id[i], NULL);

    psrfits_close(&pf_out);
    psrfits_close(&pf);

    if (rv) { fits_report_error(stderr, rv); }
    exit(0);
}
