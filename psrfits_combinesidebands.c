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
#include "psrfits_combinesidebands_cmd.h"

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

int main(int argc, char *argv[]) {
    Cmdline *cmd;
    struct psrfits pfupper,pflower, pfo;
    int stat=0, padding=0, userN=0;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);
    pfupper.tot_rows = pfupper.N = pfupper.T = pfupper.status = 0;
    pflower.tot_rows = pflower.N = pflower.T = pflower.status = 0;
    pfupper.filenum = pflower.filenum = 1;
    pfo.tot_rows = pfo.N = pfo.T = pfo.status = 0;
    sprintf(pfupper.filename, cmd->argv[0]);
    sprintf(pflower.filename, cmd->argv[0]);
    char *loc;
    if((loc=strstr(pfupper.filename,"s1"))!=NULL)
      strncpy(loc,"s0",2);
    else if((loc=strstr(pflower.filename,"s0"))!=NULL)
      strncpy(loc,"s1",2);
    else
    {
      printf("Unable to determine which sideband is which\n");
      exit(EXIT_FAILURE);
    }
    char outfilename[200]="\0";
    loc=strstr(pflower.filename,"s1");
    strncpy(outfilename,pflower.filename,loc-pflower.filename);
    strncpy(outfilename+(loc-pflower.filename),loc+2,4);
    int rv=psrfits_open(&pfupper);
    if (rv) { fits_report_error(stderr, rv); exit(1); }
    rv=psrfits_open(&pflower);
    if (rv) { fits_report_error(stderr, rv); exit(1); }
    pfo=pflower;
    printf("filenum=%d\n",pflower.filenum);
    printf("filenum=%d\n",pfo.filenum);
    sprintf(pfo.basefilename,outfilename);
    printf("pfo.basefilename=%s\n",pfo.basefilename);
    pfo.filenum = 0;
    pfo.filename[0] = '\0';
    pfo.rownum = 1;
    pfo.tot_rows = 0;
    pfo.N = 0;
    if(pfupper.rows_per_file!=pflower.rows_per_file)
    {
      fprintf(stderr,"rows_per_file in input files do not match!\n");
      exit(1);
    }
    
    //Find the number of channels in the upper sideband which will be skipped
/*    double upperfreqoflower=pflower.hdr.fctr+(double)(pflower.hdr.nchan/2)*fabs(pflower.hdr.df);
    double nextfromlower=upperfreqoflower+fabs(pflower.hdr.df);
    double lowerfreqofupper=pfupper.hdr.fctr-(double)(pfupper.hdr.nchan/2)*fabs(pfupper.hdr.df);
    double numchandiff=(nextfromlower-lowerfreqofupper)/fabs(pfupper.hdr.df);
    int chanskip;
    if(numchandiff-(double)((int)numchandiff)>.5)
      chanskip=(int)numchandiff+1;
    else
      chanskip=(int)numchandiff;*/

    double upperfreqoflower,nextfromlower,lowerfreqofupper,numchandiff;
    int chanskip;

/*    //Using the number of skipped channels, find new values for nchan,BW, and fctr
    pfo.hdr.nchan=pfupper.hdr.nchan+pflower.hdr.nchan-chanskip;
    pfo.hdr.BW=(double)pfo.hdr.nchan*fabs(pflower.hdr.df);
    pfo.hdr.fctr=(pflower.hdr.fctr-(double)(pflower.hdr.nchan/2)*fabs(pflower.hdr.df))+pfo.hdr.BW/2.0;
    pfo.sub.bytes_per_subint=pfo.hdr.nchan*pfo.hdr.nsblk*pfo.hdr.nbits/8*pfo.hdr.npol;*/

    pflower.sub.dat_freqs = (float *)malloc(sizeof(float) * pflower.hdr.nchan);
    pflower.sub.dat_weights = (float *)malloc(sizeof(float) * pflower.hdr.nchan);
    pflower.sub.dat_offsets = (float *)malloc(sizeof(float) * pflower.hdr.nchan * pflower.hdr.npol);
    pflower.sub.dat_scales  = (float *)malloc(sizeof(float) * pflower.hdr.nchan * pflower.hdr.npol);
    pflower.sub.data = (unsigned char *)malloc(pflower.sub.bytes_per_subint);

    pfupper.sub.dat_freqs = (float *)malloc(sizeof(float) * pfupper.hdr.nchan);
    pfupper.sub.dat_weights = (float *)malloc(sizeof(float) * pfupper.hdr.nchan);
    pfupper.sub.dat_offsets = (float *)malloc(sizeof(float) * pfupper.hdr.nchan * pfupper.hdr.npol);
    pfupper.sub.dat_scales  = (float *)malloc(sizeof(float) * pfupper.hdr.nchan * pfupper.hdr.npol);
    pfupper.sub.data = (unsigned char *)malloc(pfupper.sub.bytes_per_subint);

/*    pfo.sub.dat_freqs = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
    pfo.sub.dat_weights = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
    pfo.sub.dat_offsets = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
    pfo.sub.dat_scales = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
    pfo.sub.data = (unsigned char *)malloc(pfo.sub.bytes_per_subint);*/

    int firsttime=1;
    int loops=0;
    do {
      print_percent_complete(pflower.rownum, pflower.rows_per_file, pflower.rownum==1 ? 1 : 0);
      psrfits_read_subint(&pflower);
      psrfits_read_subint(&pfupper);
      if(firsttime)
      {
        firsttime=0;
        //Find the number of channels in the upper sideband which will be skipped
        if(pflower.hdr.df<0)
        {
          upperfreqoflower=pflower.sub.dat_freqs[0];
          lowerfreqofupper=pfupper.sub.dat_freqs[pfupper.hdr.nchan-1];
        }
        else
        {
          upperfreqoflower=pflower.sub.dat_freqs[pflower.hdr.nchan-1];
          lowerfreqofupper=pfupper.sub.dat_freqs[0];
        }
        nextfromlower=upperfreqoflower+fabs(pflower.hdr.df);
        numchandiff=(nextfromlower-lowerfreqofupper)/fabs(pflower.hdr.df);
        chanskip;
        if(numchandiff>0)
        {
          if(numchandiff-(double)((int)numchandiff)>.5)
            chanskip=(int)numchandiff+1;
          else
            chanskip=(int)numchandiff;
        }
        else
          chanskip=0;
        printf("upperfreqoflower=%f\n",upperfreqoflower);
        printf("lowerfreqofupper=%f\n",lowerfreqofupper);
        printf("nextfromlower=%f\n",nextfromlower);
        printf("numchandiff=%f\n",numchandiff);
        printf("chanskip=%d\n",chanskip);
        printf("nbits=%d\n",pfo.hdr.nbits);
        //Using the number of skipped channels, find new values for nchan,BW, and fctr
        pfo.hdr.nchan=pfupper.hdr.nchan+pflower.hdr.nchan-chanskip+2;
        pfo.hdr.BW=(double)pfo.hdr.nchan*fabs(pflower.hdr.df);
        pfo.hdr.fctr=(pflower.hdr.fctr-(double)(pflower.hdr.nchan/2)*fabs(pflower.hdr.df))+pfo.hdr.BW/2.0;
        pfo.sub.bytes_per_subint=pfo.hdr.nchan*pfo.hdr.nsblk*pfo.hdr.nbits/8*pfo.hdr.npol;
        printf("pfo.sub.bytes_per_subint=%d\n",pfo.sub.bytes_per_subint);
        pfo.sub.dat_freqs = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
        pfo.sub.dat_weights = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
        pfo.sub.dat_offsets = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
        pfo.sub.dat_scales = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
        pfo.sub.data = (unsigned char *)malloc(pfo.sub.bytes_per_subint);
      }
      if(pflower.status == 0 && pfupper.status == 0)
      {
        pfo.sub.tsubint=pflower.sub.tsubint;
        pfo.sub.offs=pflower.sub.offs;
        pfo.sub.lst=pflower.sub.lst;
        pfo.sub.ra=pflower.sub.ra;
        pfo.sub.dec=pflower.sub.dec;
        pfo.sub.glon=pflower.sub.glon;
        pfo.sub.glat=pflower.sub.glat;
        pfo.sub.feed_ang=pflower.sub.feed_ang;
        pfo.sub.pos_ang=pflower.sub.pos_ang;
        pfo.sub.par_ang=pflower.sub.par_ang;
        pfo.sub.tel_az=pflower.sub.tel_az;
        pfo.sub.tel_zen=pflower.sub.tel_zen;
        pfo.sub.FITS_typecode=pflower.sub.FITS_typecode;

        if(pflower.hdr.df<0)
        {
          //Copy frequency labels
          pfo.sub.dat_freqs[1]=pfupper.sub.dat_freqs[0]+fabs(pflower.hdr.df);
          pfo.sub.dat_freqs[0]=pfo.sub.dat_freqs[1]+fabs(pflower.hdr.df);
          memcpy(pfo.sub.dat_freqs+2,pfupper.sub.dat_freqs,sizeof(float)*(pfupper.hdr.nchan-chanskip));
          memcpy((pfo.sub.dat_freqs)+pfupper.hdr.nchan-chanskip+2,pflower.sub.dat_freqs,sizeof(float)*pflower.hdr.nchan);
          //Copy weights
          pfo.sub.dat_weights[0]=pfo.sub.dat_weights[1]=0;
          memcpy(pfo.sub.dat_weights+2,pfupper.sub.dat_weights,sizeof(float)*(pfupper.hdr.nchan-chanskip));
          memcpy((pfo.sub.dat_weights+2)+pfupper.hdr.nchan-chanskip,pflower.sub.dat_weights,sizeof(float)*pflower.hdr.nchan);
          //Copy offsets
          pfo.sub.dat_offsets[0]=pfo.sub.dat_offsets[1]=0;
          memcpy(pfo.sub.dat_offsets+(2*pfo.hdr.npol),pfupper.sub.dat_offsets,sizeof(float)*(pfupper.hdr.nchan-chanskip)*pfupper.hdr.npol);
          memcpy((pfo.sub.dat_offsets)+(pfupper.hdr.nchan-chanskip)*pfo.hdr.npol+(2*pfo.hdr.npol),pflower.sub.dat_offsets,sizeof(float)*(pflower.hdr.nchan)*pflower.hdr.npol);
          //Copy scales
          pfo.sub.dat_scales[0]=pfo.sub.dat_scales[1]=0;
          memcpy(pfo.sub.dat_scales+(2*pfo.hdr.npol),pfupper.sub.dat_scales,sizeof(float)*(pfupper.hdr.nchan-chanskip)*pfupper.hdr.npol);
          memcpy((pfo.sub.dat_scales)+(pfupper.hdr.nchan-chanskip)*pfo.hdr.npol+(2*pfo.hdr.npol),pflower.sub.dat_scales,sizeof(float)*(pflower.hdr.nchan)*pflower.hdr.npol);
          //Copy the data
          int i=0;
          //int k=0;
          for(i=0;i<pfo.hdr.nsblk;++i)
          {
            int j=0;
//            printf("upper %d %d %d\n",(i*(int)((double)pfo.hdr.nchan*(double)pfo.hdr.npol*((double)pfo.hdr.nbits/8.0))+(int)(2.0*((double)pfo.hdr.nbits/8.0))),i*(int)((double)pfupper.hdr.nchan*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)),(int)((double)(pfupper.hdr.nchan-chanskip)*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)));
            memcpy(pfo.sub.data+(i*(int)((double)pfo.hdr.nchan*(double)pfo.hdr.npol*((double)pfo.hdr.nbits/8.0))+(int)(2.0*((double)pfo.hdr.nbits/8.0))),pfupper.sub.data+i*(int)((double)pfupper.hdr.nchan*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)),(int)((double)(pfupper.hdr.nchan-chanskip)*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)));
//            printf("lower %d %d %d\n",(i*(int)((double)pfo.hdr.nchan*(double)pfo.hdr.npol*((double)pfo.hdr.nbits/8.0))+(int)(2.0*((double)pfo.hdr.nbits/8.0)))+(int)((double)pfupper.hdr.nchan-chanskip*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)),i*(int)((double)pflower.hdr.nchan*(double)pflower.hdr.npol*((double)pflower.hdr.nbits/8.0)),(int)((double)pflower.hdr.nchan*(double)pflower.hdr.npol*((double)pflower.hdr.nbits/8.0)));
            memcpy(pfo.sub.data+(i*(int)((double)pfo.hdr.nchan*(double)pfo.hdr.npol*((double)pfo.hdr.nbits/8.0))+(int)(2.0*((double)pfo.hdr.nbits/8.0)))+(int)((double)(pfupper.hdr.nchan-chanskip)*(double)pfupper.hdr.npol*((double)pfupper.hdr.nbits/8.0)),pflower.sub.data+i*(int)((double)pflower.hdr.nchan*(double)pflower.hdr.npol*((double)pflower.hdr.nbits/8.0)),(int)((double)pflower.hdr.nchan*(double)pflower.hdr.npol*((double)pflower.hdr.nbits/8.0)));
            //for(j=0;j<pfupper.hdr.nchan-chanskip;++j,++k)
            //  printf("%d %d %d %d\n",i*pfupper.hdr.nchan+j,pfupper.sub.data[i*pfupper.hdr.nchan+j],k,pfo.sub.data[k]);
            //for(j=0;j<pflower.hdr.nchan;++j,++k)
            //  printf("%d %d %d %d\n",i*pflower.hdr.nchan+j,pflower.sub.data[i*pflower.hdr.nchan+j],k,pfo.sub.data[k]);
          }
          psrfits_write_subint(&pfo);
          //++loops;
          //if(loops==2)
          //  exit(0);
        }
        else
        {
        }
      }
    } while(pflower.status == 0 && pfupper.status == 0);
    
    exit(0);
}
