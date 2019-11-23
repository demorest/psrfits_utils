#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "vectors.h"
#include "psrfits2psrfits_cmd.h"
#include "psrfits.h"
#include "rescale.h"

//If this flag is 1, and 4-bit conversion is requested, the scaling will be
//done for 4-bit conversion, but scaled data will be written out as 8-bit
#define DEBUG 0

int main(int argc, char *argv[])
{
    int numfiles, ii, numrows, rownum, ichan, itsamp, datidx;
    int spec_per_row, status;
    float offset, scale, datum, packdatum, maxval, fulltsubint;
    float *datachunk;
    FILE **infiles;
    struct psrfits pfin, pfout;
    Cmdline *cmd;
    fitsfile *infits, *outfits;
    char outfilename[128], tform[8];
    char *pc1, *pc2;
    int first = 1, dummy = 0, nclipped;
    short int *inrowdata;
    unsigned char *outrowdata;

    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(1);
    }
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);
    numfiles = cmd->argc;
    infiles = (FILE **) malloc(numfiles * sizeof(FILE *));

#ifdef DEBUG
    showOptionValues();
#endif

    printf("\n         PSRFITS 16-bit to 4-bit Conversion Code\n");
    printf("         by J. Deneva, S. Ransom, & S. Chatterjee\n\n");

    // Open the input files
    status = 0;                 //fits_close segfaults if this is not initialized
    printf("Reading input data from:\n");
    for (ii = 0; ii < numfiles; ii++) {
        printf("  '%s'\n", cmd->argv[ii]);

        //Get the file basename and number from command-line argument
        //(code taken from psrfits2fil)
        pc2 = strrchr(cmd->argv[ii], '.');      // at .fits
        *pc2 = 0;               // terminate string
        pc1 = pc2 - 1;
        while ((pc1 >= cmd->argv[ii]) && isdigit(*pc1))
            pc1--;
        if (pc1 <= cmd->argv[ii]) {     // need at least 1 char before filenum
            puts("Illegal input filename. must have chars before the filenumber");
            exit(1);
        }
        pc1++;                  // we were sitting on "." move to first digit
        pfin.filenum = atoi(pc1);
        pfin.fnamedigits = pc2 - pc1;   // how many digits in filenumbering scheme.
        *pc1 = 0;               // null terminate the basefilename
        strcpy(pfin.basefilename, cmd->argv[ii]);
        pfin.initialized = 0;   // set to 1 in  psrfits_open()
        pfin.status = 0;
        //(end of code taken from psrfits2fil)

        //Open the existing psrfits file
        if (psrfits_open(&pfin, READONLY) != 0) {
            fprintf(stderr, "error opening file\n");
            fits_report_error(stderr, pfin.status);
            exit(1);
        }
        // Create the subint arrays
        if (first) {
            pfin.sub.dat_freqs = (float *) malloc(sizeof(float) * pfin.hdr.nchan);
            pfin.sub.dat_weights = (float *) malloc(sizeof(float) * pfin.hdr.nchan);
            pfin.sub.dat_offsets =
                (float *) malloc(sizeof(float) * pfin.hdr.nchan * pfin.hdr.npol);
            pfin.sub.dat_scales =
                (float *) malloc(sizeof(float) * pfin.hdr.nchan * pfin.hdr.npol);
            //first is set to 0 after data buffer allocation further below
        }

        infits = pfin.fptr;
        spec_per_row = pfin.hdr.nsblk;
        numrows = pfin.rows_per_file;

        if (ii == 0) {

            // Create the output PSRFITS file
            sprintf(outfilename, "%s.%0*d.fits", cmd->outfile, pfin.fnamedigits,
                    pfin.filenum);
            fits_create_file(&outfits, outfilename, &status);

            //Instead of copying HDUs one by one, can move to the SUBINT
            //HDU, and copy all the HDUs preceding it
            fits_movnam_hdu(infits, BINARY_TBL, "SUBINT", 0, &status);
            fits_copy_file(infits, outfits, 1, 0, 0, &status);

            //Copy the SUBINT table header
            fits_copy_header(infits, outfits, &status);
            fits_flush_buffer(outfits, 0, &status);

            //Set NAXIS2 in the output SUBINT table to 0 b/c we haven't 
            //written any rows yet
            dummy = 0;
            fits_update_key(outfits, TINT, "NAXIS2", &dummy, NULL, &status);

            //Edit the NBITS key
            if (DEBUG) {
                dummy = 8;
                fits_update_key(outfits, TINT, "NBITS", &dummy, NULL, &status);
            } else {
                fits_update_key(outfits, TINT, "NBITS", &(cmd->numbits), NULL,
                                &status);
            }

            //Edit the TFORM17 column: # of data bytes per row 
            //fits_get_colnum(outfits,1,"DATA",&dummy,&status);
            if (DEBUG)
                sprintf(tform, "%dB",
                        pfin.hdr.nsblk * pfin.hdr.nchan * pfin.hdr.npol);
            else
                sprintf(tform, "%dB", pfin.hdr.nsblk * pfin.hdr.nchan *
                        pfin.hdr.npol * cmd->numbits / 8);

            fits_update_key(outfits, TSTRING, "TTYPE17", "DATA", NULL, &status);
            fits_update_key(outfits, TSTRING, "TFORM17", tform, NULL, &status);

            //Edit NAXIS1: row width in bytes
            fits_read_key(outfits, TINT, "NAXIS1", &dummy, NULL, &status);
            if (DEBUG) {
                dummy = dummy - pfin.hdr.nsblk * pfin.hdr.nchan *
                    pfin.hdr.npol * (pfin.hdr.nbits - 8) / 8;
            } else {
                dummy = dummy - pfin.hdr.nsblk * pfin.hdr.nchan *
                    pfin.hdr.npol * (pfin.hdr.nbits - cmd->numbits) / 8;
            }
            fits_update_key(outfits, TINT, "NAXIS1", &dummy, NULL, &status);
            fits_close_file(outfits, &status);

            //Now reopen the file so that the pfout structure is initialized
            pfout.status = 0;
            pfout.initialized = 0;
            pfout.fnamedigits = pfin.fnamedigits;
            pfout.filenum = pfin.filenum;
            sprintf(pfout.basefilename, "%s.", cmd->outfile);

            if (psrfits_open(&pfout, READWRITE) != 0) {
                fprintf(stderr, "error opening file\n");
                fits_report_error(stderr, pfout.status);
                exit(1);
            }
            outfits = pfout.fptr;
            maxval = pow(2, cmd->numbits) - 1;
            fprintf(stderr, "maxval: %f\n", maxval);

            //These are not initialized in psrfits_open but are needed 
            //in psrfits_write_subint (not obvious what are the corresponding 
            //fields in any of the psrfits table headers)
            pfout.hdr.ds_freq_fact = 1;
            pfout.hdr.ds_time_fact = 1;

            rownum = 0;
        }


        while (psrfits_read_subint(&pfin, first) == 0) {
            fprintf(stderr, "Working on row %d\n", ++rownum);

            //If this is the first row, store the length of a full subint
            if (rownum == 1)
                fulltsubint = pfin.sub.tsubint;

            //If this is the last row for the observation, and it's partial,
            //drop it.
            if (rownum == numrows && pfin.sub.tsubint != fulltsubint) {
                fprintf(stderr,
                        "Dropping partial row of length %f s (full row is %f s)\n",
                        pfin.sub.tsubint, fulltsubint);
                break;
            }
            //Copy the subint struct from pfin to pfout, but correct 
            //elements that are not the same 
            pfout.sub = pfin.sub;       //this copies array pointers too
            pfout.sub.bytes_per_subint =
                pfin.sub.bytes_per_subint * pfout.hdr.nbits / pfin.hdr.nbits;
            pfout.sub.dataBytesAlloced = pfout.sub.bytes_per_subint;
            pfout.sub.FITS_typecode = TBYTE;

            if (first) {
                //Allocate scaling buffer and output buffer
                datachunk = gen_fvect(spec_per_row);
                outrowdata = gen_bvect(pfout.sub.bytes_per_subint);

                first = 0;
            }
            pfout.sub.data = outrowdata;

            inrowdata = (short int *) pfin.sub.data;
            nclipped = 0;

            // Loop over all the channels:
            for (ichan = 0; ichan < pfout.hdr.nchan * pfout.hdr.npol; ichan++) {
                // Populate datachunk[] by picking out all time samples for ichan
                for (itsamp = 0; itsamp < spec_per_row; itsamp++)
                    datachunk[itsamp] = (float) (inrowdata[ichan + itsamp *
                                                           pfout.hdr.nchan *
                                                           pfout.hdr.npol]);

                // Compute the statistics here, and put the offsets and scales in
                // pf.sub.dat_offsets[] and pf.sub.dat_scales[]

                if (rescale(datachunk, spec_per_row, cmd->numbits, &offset, &scale)
                    != 0) {
                    printf("Rescale routine failed!\n");
                    return (-1);
                }
                pfout.sub.dat_offsets[ichan] = offset;
                pfout.sub.dat_scales[ichan] = scale;

                // Since we have the offset and scale ready, rescale the data:
                for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                    datum = (scale == 0.0) ? 0.0 :
                        roundf((datachunk[itsamp] - offset) / scale);
                    if (datum < 0.0) {
                        datum = 0;
                        nclipped++;
                    } else if (datum > maxval) {
                        datum = maxval;
                        nclipped++;
                    }

                    inrowdata[ichan + itsamp * pfout.hdr.nchan * pfout.hdr.npol] =
                        (short int) datum;

                }
                // Now inrowdata[ichan] contains rescaled ints.
            }

            // Then do the conversion and store the
            // results in pf.sub.data[] 
            if (cmd->numbits == 8 || DEBUG) {
                for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                    datidx = itsamp * pfout.hdr.nchan * pfout.hdr.npol;
                    for (ichan = 0; ichan < pfout.hdr.nchan * pfout.hdr.npol;
                         ichan++, datidx++) {
                        pfout.sub.data[datidx] = (unsigned char) inrowdata[datidx];
                    }
                }
            } else if (cmd->numbits == 4) {
                for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                    datidx = itsamp * pfout.hdr.nchan * pfout.hdr.npol;
                    for (ichan = 0; ichan < pfout.hdr.nchan * pfout.hdr.npol;
                         ichan += 2, datidx += 2) {

                        packdatum = inrowdata[datidx] * 16 + inrowdata[datidx + 1];
                        pfout.sub.data[datidx / 2] = (unsigned char) packdatum;
                    }
                }
            } else if (cmd->numbits == 2) {
                for (itsamp = 0; itsamp < spec_per_row; itsamp++) {
                    datidx = itsamp * pfout.hdr.nchan * pfout.hdr.npol;
                    for (ichan = 0; ichan < pfout.hdr.nchan * pfout.hdr.npol;
                         ichan += 4, datidx += 4) {

                        packdatum = inrowdata[datidx] * 256 + inrowdata[datidx + 1];
                        pfout.sub.data[datidx / 4] = (unsigned char) packdatum;
                    }
                }
            } else {
                fprintf(stderr, "Only 2, 4, or 8-bit output formats supported.\n");
                fprintf(stderr, "Bits per sample requested: %d\n", cmd->numbits);
                exit(1);
            }


            //pfout.sub.offs = (pfout.tot_rows+0.5) * pfout.sub.tsubint;
            fprintf(stderr, "nclipped: %d fraction clipped: %f\n", nclipped,
                    (float) nclipped / (pfout.hdr.nchan * pfout.hdr.npol *
                                        pfout.hdr.nsblk));

            // Now write the row. 
            status = psrfits_write_subint(&pfout);
            if (status) {
                printf("\nError (%d) writing PSRFITS...\n\n", status);
                break;
            }
        }

        //Close the files 
        fits_close_file(infits, &status);
    }

    fits_close_file(outfits, &status);

    // Free the structure arrays too...
    free(datachunk);
    free(infiles);

    free(pfin.sub.dat_freqs);
    free(pfin.sub.dat_weights);
    free(pfin.sub.dat_offsets);
    free(pfin.sub.dat_scales);

    free(pfin.sub.data);
    free(pfout.sub.data);
    free(pfin.sub.stat);

    return 0;
}
