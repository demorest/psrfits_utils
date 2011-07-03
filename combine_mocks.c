//This code combines the two frequency bands from the Mock 
//spectrometers at Arecibo. K. Stovall  Oct 2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <libgen.h>
#include "psrfits.h"
#include "combine_mocks_cmd.h"

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

//Routine taken from PRESTO
void avg_var(float *x, int n, double *mean, double *var)
/* For a float vector, *x, of length n, this routine  */
/* returns the mean and variance of *x.               */
{
   long i;
   double an = 0.0, an1 = 0.0, dx;

   /*  Modified (29 June 98) C version of the following:        */
   /*  ALGORITHM AS 52  APPL. STATIST. (1972) VOL.21, P.226     */
   /*  Returned values were checked with Mathematica 3.01       */

   if (n < 1) {
      printf("\vVector length must be > 0 in avg_var().  Exiting\n");
      exit(1);
   } else {
      *mean = (double) x[0];
      *var = 0.0;
   }

   for (i = 1; i < n; i++) {
      an = (double) (i + 1);
      an1 = (double) (i);
      dx = (x[i] - *mean) / an;
      *var += an * an1 * dx * dx;
      *mean += dx;
   }

   if (n > 1)
      *var /= an1;

   return;
}

//End of routine taken from PRESTO

int main(int argc, char *argv[])
{
   Cmdline *cmd;
   struct psrfits pfupper, pflower, pfo;
   fitsfile *infits, *outfits;
   char *pc1, *pc2;
   char outfilename[200];       //Name of outfile if not specified on command line
   int stat = 0, padding = 0, userN = 0, status;

   // Call usage() if we have no command line arguments
   if (argc == 1) {
      Program = argv[0];
      usage();
      exit(0);
   }
   // Parse the command line using the excellent program Clig
   cmd = parseCmdline(argc, argv);
   pfupper.tot_rows = pfupper.N = pfupper.T = pfupper.status = 0;       //Initialize upper band
   pflower.tot_rows = pflower.N = pflower.T = pflower.status = 0;       //Initialize lower band
   pfupper.filenum = pflower.filenum = 1;
   pfo.tot_rows = pfo.N = pfo.T = pfo.status = pfo.multifile = 0;       //Initialize output
   sprintf(pfupper.filename, cmd->argv[0]);     //Copy filename specified on command line to
   sprintf(pflower.filename, cmd->argv[0]);     //upper and lower bands, will correct filenames shortly
   if ((pc2 = strstr(pfupper.filename, "s1")) != NULL)  //Upper contains s1, change to s0
      strncpy(pc2, "s0", 2);
   else if ((pc2 = strstr(pflower.filename, "s0")) != NULL)     //Lower contains s0, change to s1
      strncpy(pc2, "s1", 2);
   else {
      printf("Unable to determine which sideband is which\n");
      exit(EXIT_FAILURE);
   }
   //Setting the name of the output file, setting as same name as input file, but removing s0/s1. 
   pc1 = strstr(pflower.filename, "s1");
   pc2 = strrchr(pflower.filename, '.');        //At '.fits'
   pc2--;
   while ((pc2 >= pflower.filename) && isdigit(*pc2))   //Move through the digits to the separation char.
      pc2--;
   strncpy(outfilename, pflower.filename, pc1 - pflower.filename);      //Copy everything up to s1 into outfilename
   strncpy(outfilename + (pc1 - pflower.filename), pc1 + 2, pc2 - pc1 - 2);     //Concatenate from after s1 to char before the separation char.
   pc1 = outfilename + (pc2 - pflower.filename - 2);
   *pc1 = 0;
   int rv = psrfits_open(&pfupper);   //Open upper band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   rv = psrfits_open(&pflower);       //Open lower band
   if (rv) {
      fits_report_error(stderr, rv);
      exit(1);
   }
   pfo = pflower;               //Copy all lower band variables into the output struct
   if (!cmd->outputbasenameP)
      sprintf(pfo.basefilename, basename(outfilename));
   else
      sprintf(pfo.basefilename, cmd->outputbasename);
   pfo.filenum = 0;
   sprintf(pfo.filename, "\0"); //Set filename to null so psrfits_open will create the filename for me
   pfo.rownum = 1;
   pfo.tot_rows = 0;
   pfo.N = 0;
   if (pfupper.rows_per_file != pflower.rows_per_file) {        //Sanity check for the two input frequency bands
      fprintf(stderr, "rows_per_file in input files do not match!\n");
      exit(1);
   }

   double upperfreqoflower, nextfromlower, lowerfreqofupper, numchandiff;       //Used to find which frequencies to take from each band
   double offsetfactor, scalefactor;    //Factors which will be applied to offsets and scales
   int upchanskip, lowchanskip; //Number of channels to skip in each banda

   //Variables used to make code cleaner
   int extrachanoffset, outoffset, upperoffset, numtocopyupper, loweroffset_skip,
       loweroffset, numtocopylower, newuppernchan, newlowernchan;
   double df = pflower.hdr.df;
   int nchan = pflower.hdr.nchan;
   int outnchan;
   int npol = pflower.hdr.npol;
   int nbits = pflower.hdr.nbits;
   int nsblk = pflower.hdr.nsblk;
   //Allocate memory for all upper and lower data
   pflower.sub.dat_freqs = (float *) malloc(sizeof(float) * nchan);
   pflower.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pflower.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pflower.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pflower.sub.rawdata = (unsigned char *) malloc(pflower.sub.bytes_per_subint);
   pflower.sub.data = (unsigned char *) malloc(pflower.sub.bytes_per_subint*2);

   pfupper.sub.dat_freqs = (float *) malloc(sizeof(float) * nchan);
   pfupper.sub.dat_weights = (float *) malloc(sizeof(float) * nchan);
   pfupper.sub.dat_offsets = (float *) malloc(sizeof(float) * nchan * npol);
   pfupper.sub.dat_scales = (float *) malloc(sizeof(float) * nchan * npol);
   pfupper.sub.rawdata = (unsigned char *) malloc(pfupper.sub.bytes_per_subint);
   pfupper.sub.data = (unsigned char *) malloc(pfupper.sub.bytes_per_subint*2);

   int firsttime = 1;           //First time through do while loop
   do {
      print_percent_complete(pflower.rownum, pflower.rows_per_file,
                             pflower.rownum == 1 ? 1 : 0);
      psrfits_read_subint(&pflower);
      psrfits_read_subint(&pfupper);
      if (firsttime) {          //First time through loop, calculate factors for scales and offsets and number of channels to skip
         firsttime = 0;
         //Find the number of channels in the upper band which will be skipped
         if (df < 0) {          //Find channel order, low to high or high to low
            upperfreqoflower = pflower.sub.dat_freqs[0];        //Highest frequency channel in lower band
            lowerfreqofupper = pfupper.sub.dat_freqs[nchan - 1];        //Lowest frequency channel in upper band
         } else {
            upperfreqoflower = pflower.sub.dat_freqs[nchan - 1];        //Highest frequency channel in lower band
            lowerfreqofupper = pfupper.sub.dat_freqs[0];        //Lowest frequency channel in upper band
         }
         nextfromlower = upperfreqoflower + fabs(df);   //Second highest channel in lower band
         numchandiff = (nextfromlower - lowerfreqofupper) / fabs(df);   //Number of channels to skip in float form
         int chanskip;
         if (numchandiff > 0) { //Make sure there are channels which need to be skipped
            if (numchandiff - (double) ((int) numchandiff) > .5)        // See whether we need to round up integer channels to skip
               chanskip = (int) numchandiff + 1;
            else
               chanskip = (int) numchandiff;
         } else
            chanskip = 0;       //No need to skip any channels
         if (chanskip % 2 == 1) {       //Odd number of channels, give lower band the extra channel
            upchanskip = chanskip / 2;
            lowchanskip = chanskip / 2 + 1;
         } else                 //Even number of channels to skip
            upchanskip = lowchanskip = chanskip / 2;
         if (upchanskip % 2 == 1) {     //We want an even number of channels in upper band for 4-bit data to get copied correctly
            ++lowchanskip;
            --upchanskip;
         }
         //Find new values given the number of channels skipped
         pfo.hdr.nchan = outnchan = nchan + nchan - chanskip + 2;       //New number of channels, plus 2 to make nchan=960 (many factors of 2)
         pfo.hdr.BW = (double) outnchan *fabs(df);      //New bandwidth
         pfo.hdr.fctr =         //New center frequency
             (pflower.hdr.fctr - (double) (nchan / 2) * fabs(df)) + pfo.hdr.BW / 2.0;
         pfo.sub.bytes_per_subint =     //Calculate new number of bytes in each subint
             outnchan * nsblk * nbits / 8 * npol;
         //Allocate space for output data now that we know the new number of channels
         pfo.sub.dat_freqs = (float *) malloc(sizeof(float) * outnchan);
         pfo.sub.dat_weights = (float *) malloc(sizeof(float) * outnchan);
         pfo.sub.dat_offsets = (float *) malloc(sizeof(float) * outnchan * npol);
         pfo.sub.dat_scales = (float *) malloc(sizeof(float) * outnchan * npol);
         pfo.sub.rawdata = (unsigned char *) malloc(pfo.sub.bytes_per_subint);
         pfo.sub.data = (unsigned char *) malloc(pfo.sub.bytes_per_subint*2);
         newuppernchan = nchan - upchanskip;    //The number of channels to copy from the upper sideband.
         newlowernchan = nchan - lowchanskip;   //The number of channels to copy from the lower sideband.

         extrachanoffset = 2;     //Offset for 2 extra freq channels making nchan 960 in bytes
         outoffset = (outnchan * npol);     //Offset in each loop due to previously written data
         upperoffset = (nchan * npol);      //Offset in loop for upper band
         numtocopyupper = (newuppernchan * npol);   //Number of bytes to copy from upper band
         loweroffset_skip = (lowchanskip * npol);   //Number of bytes to skip when copying lower band due to 
         //having written upper band
         loweroffset =          //Number of bytes to skip due to having written previous lower band data
             (nchan * npol);
         numtocopylower = (newlowernchan * npol);   //Number of bytes to copy from lower band
         double upmean, upvar, lowmean, lowvar;
         avg_var(pfupper.sub.dat_offsets + (nchan - upchanskip),        //Find the mean and variance of the upper band's offsets
                 upchanskip, &upmean, &upvar);
         printf("Upper offset stats: mean=%f variance=%f\n", upmean, upvar);
         avg_var(pflower.sub.dat_offsets, lowchanskip, &lowmean, &lowvar);      //Find the mean and variance of the lower band's offsets
         printf("Lower offset stats: mean=%f variance=%f\n", lowmean, lowvar);
         printf("Applying factor of %f to upper offsets\n", (lowmean / upmean));
         offsetfactor = lowmean / upmean;       //Set offset factor used to correct variance differences in the two bands
         avg_var(pfupper.sub.dat_scales + (nchan - upchanskip), //Find the mean and var. of the upper band's scales
                 upchanskip, &upmean, &upvar);
         printf("Upper scales stats: mean=%f variance=%f\n", upmean, upvar);
         avg_var(pflower.sub.dat_scales, lowchanskip, &lowmean, &lowvar);       //Find the mean and var. of the lower band's scales
         printf("Lower scales stats: mean=%f variance=%f\n", lowmean, lowvar);
         printf("Applying factor of %f to upper scales\n", (lowmean / upmean));
         scalefactor = lowmean / upmean;        //Set scale factor used to correct variance differences in the two bands
      }
      if (pflower.status == 0 && pfupper.status == 0) {
         //Copy info from the lower band subint struct to the output file's subint struct
         pfo.sub.tsubint = pflower.sub.tsubint;
         pfo.sub.offs = pflower.sub.offs;
         pfo.sub.lst = pflower.sub.lst;
         pfo.sub.ra = pflower.sub.ra;
         pfo.sub.dec = pflower.sub.dec;
         pfo.sub.glon = pflower.sub.glon;
         pfo.sub.glat = pflower.sub.glat;
         pfo.sub.feed_ang = pflower.sub.feed_ang;
         pfo.sub.pos_ang = pflower.sub.pos_ang;
         pfo.sub.par_ang = pflower.sub.par_ang;
         pfo.sub.tel_az = pflower.sub.tel_az;
         pfo.sub.tel_zen = pflower.sub.tel_zen;
         pfo.sub.FITS_typecode = pflower.sub.FITS_typecode;

         //Create variables to reduce column width of lines below
         float *dat_freqs = pfo.sub.dat_freqs;
         float *udat_freqs = pfupper.sub.dat_freqs;
         float *ldat_freqs = pflower.sub.dat_freqs;
         float *dat_weights = pfo.sub.dat_weights;
         float *udat_weights = pfupper.sub.dat_weights;
         float *ldat_weights = pflower.sub.dat_weights;
         float *dat_offsets = pfo.sub.dat_offsets;
         float *udat_offsets = pfupper.sub.dat_offsets;
         float *ldat_offsets = pflower.sub.dat_offsets;
         float *dat_scales = pfo.sub.dat_scales;
         float *udat_scales = pfupper.sub.dat_scales;
         float *ldat_scales = pflower.sub.dat_scales;
         unsigned char *data = pfo.sub.data;
         unsigned char *udata = pfupper.sub.data;
         unsigned char *ldata = pflower.sub.data;

         if (df < 0) {
            //Copy frequency labels
            dat_freqs[1] = udat_freqs[0] + fabs(df);    //Calculate the frequency labels
            dat_freqs[0] = dat_freqs[1] + fabs(df);     //for our two empty frequency channels
            int newuppernchan = nchan - upchanskip;     //The number of channels to copy from the upper band
            int newlowernchan = nchan - lowchanskip;    //The number of channels to copy from the lower band
            memcpy(dat_freqs + 2, udat_freqs, sizeof(float) * newuppernchan);   //Copy from the upper band, skipping first two chans.
            memcpy(dat_freqs + newuppernchan + 2,       //Copy from the lower band
                   ldat_freqs + lowchanskip, sizeof(float) * newlowernchan);
            //Copy weights
            dat_weights[0] = dat_weights[1] = 0;        //Set the weights of first two channels to 0, so they shouldn't be used in calculations
            memcpy(dat_weights + 2, udat_weights,       //Copy weights from the upper band
                   sizeof(float) * newuppernchan);
            memcpy(dat_weights + 2 + newuppernchan,     //Copy weights from the lower band
                   ldat_weights + lowchanskip, sizeof(float) * newlowernchan);
            //Copy offsets
            dat_offsets[0] = dat_offsets[1] =   //Set offsets of first two channels to the same as upper's first channel
                udat_offsets[0];        //(shouldn't matter since they should be ignored)
            int ii;
            for (ii = 0; ii < newuppernchan; ++ii)      //Apply offset factor to upper band
               udat_offsets[ii] = udat_offsets[ii] * (offsetfactor);
            memcpy(dat_offsets + 2 * npol, udat_offsets,        //Copy upper offsets
                   sizeof(float) * newuppernchan * npol);
            memcpy(dat_offsets + (newuppernchan + 2) * npol,    //Copy lower offsets
                   ldat_offsets + lowchanskip, sizeof(float) * newlowernchan * npol);
            //Copy scales
            for (ii = 0; ii < newuppernchan; ++ii)      //Apply scale factor to upper band
               udat_scales[ii] = udat_scales[ii] * (scalefactor);
            dat_scales[0] = dat_scales[1] = udat_scales[0];
            memcpy(dat_scales + 2 * npol, udat_scales,  //Copy upper scales
                   sizeof(float) * newuppernchan * npol);
            memcpy(dat_scales + (newuppernchan + 2) * npol,     //Copy lower scales
                   ldat_scales + lowchanskip, sizeof(float) * newlowernchan * npol);
            //Copy the data
            for (ii = 0; ii < nsblk; ++ii) {    //Loop through data copying into place
               memcpy(data + ii * outoffset + extrachanoffset,
                      udata + ii * upperoffset, numtocopyupper);
               memcpy(data + ii * outoffset + extrachanoffset +
                      numtocopyupper,
                      ldata + ii * loweroffset + loweroffset_skip, numtocopylower);
            }
            psrfits_write_subint(&pfo);
         } else {
         }
      }
   } while (pfo.rownum <= pfo.rows_per_file);
   printf("Closing file '%s'\n", pflower.filename);
   fits_close_file(pfupper.fptr, &status);
   printf("Closing file '%s'\n", pfupper.filename);
   fits_close_file(pflower.fptr, &status);
   exit(0);
}
