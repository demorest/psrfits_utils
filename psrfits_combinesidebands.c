//This code combines the two sidebands from the Mock 
//spectrometers at Arecibo. K. Stovall  Oct 2010

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
    fitsfile *infits,*outfits;
    char outfilename[128];
    int stat=0, padding=0, userN=0,status;

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
    pfo.tot_rows = pfo.N = pfo.T = pfo.status = pfo.multifile = 0;
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
    loc=strstr(pflower.filename,"s1");
    strncpy(outfilename,pflower.filename,loc-pflower.filename);
    strncpy(outfilename+(loc-pflower.filename),loc+2,2);
    int rv=psrfits_open(&pfupper,READONLY);
    if (rv) { fits_report_error(stderr, rv); exit(1); }
    rv=psrfits_open(&pflower,READONLY);
    if (rv) { fits_report_error(stderr, rv); exit(1); }
    pfo=pflower;
    if(!cmd->outputbasenameP)
      sprintf(pfo.basefilename,outfilename);
    else
      sprintf(pfo.basefilename,cmd->outputbasename);

    pfo.filenum = 0;
    sprintf(pfo.filename,"\0");
    pfo.rownum = 1;
    pfo.tot_rows = 0;
    pfo.N = 0;
    if(pfupper.rows_per_file!=pflower.rows_per_file)
    {
      fprintf(stderr,"rows_per_file in input files do not match!\n");
      exit(1);
    }
    
    double upperfreqoflower,nextfromlower,lowerfreqofupper,numchandiff;
    int upchanskip,lowchanskip;
    int extrachanoffset,outoffset,upperoffset,numtocopyupper,loweroffset_skip,loweroffset,numtocopylower,newuppernchan,newlowernchan;

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

    int firsttime=1;
    int loops=0;
    do
    {
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
        int chanskip;
        if(numchandiff>0)
        {
          if(numchandiff-(double)((int)numchandiff)>.5)
            chanskip=(int)numchandiff+1;
          else
            chanskip=(int)numchandiff;
        }
        else
          chanskip=0;
        if(chanskip%2==1)
        {
          upchanskip=chanskip/2;
          lowchanskip=chanskip/2+1;
        }
        else
          upchanskip=lowchanskip=chanskip/2;
        if(upchanskip%2==1)
        {
          ++lowchanskip;
          --upchanskip;
        }
        //Using the number of skipped channels, find new values for nchan,BW, and fctr
        pfo.hdr.nchan=pfupper.hdr.nchan+pflower.hdr.nchan-chanskip+2;
        pfo.hdr.BW=(double)pfo.hdr.nchan*fabs(pflower.hdr.df);
        pfo.hdr.fctr=(pflower.hdr.fctr-(double)(pflower.hdr.nchan/2)*fabs(pflower.hdr.df))+pfo.hdr.BW/2.0;
        pfo.sub.bytes_per_subint=pfo.hdr.nchan*pfo.hdr.nsblk*pfo.hdr.nbits/8*pfo.hdr.npol;
        pfo.sub.dat_freqs = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
        pfo.sub.dat_weights = (float *)malloc(sizeof(float) * pfo.hdr.nchan);
        pfo.sub.dat_offsets = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
        pfo.sub.dat_scales = (float *)malloc(sizeof(float) * pfo.hdr.nchan * pfo.hdr.npol);
        pfo.sub.data = (unsigned char *)malloc(pfo.sub.bytes_per_subint);
        newuppernchan=pfupper.hdr.nchan-upchanskip; //The number of channels to copy from the upper sideband.
        newlowernchan=pflower.hdr.nchan-lowchanskip; //The number of channels to copy from the lower sideband.

        extrachanoffset=(2*pfo.hdr.nbits)/8; //Offset for 2 extra freq channels
        outoffset=(pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nbits)/8;
        upperoffset=(pfupper.hdr.nchan*pfupper.hdr.npol*pfo.hdr.nbits)/8;
        numtocopyupper=(newuppernchan*pfupper.hdr.npol*pfupper.hdr.nbits)/8;
        loweroffset_skip=(lowchanskip*pflower.hdr.npol*pflower.hdr.nbits)/8;
        loweroffset=(pflower.hdr.nchan*pflower.hdr.npol*pflower.hdr.nbits)/8;
        numtocopylower=(newlowernchan*pflower.hdr.npol*pflower.hdr.nbits)/8;
        char fitsoutfilename[128];
        sprintf(fitsoutfilename,"%s_%04d.fits",pfo.basefilename,pflower.filenum);
        fits_create_file(&outfits,fitsoutfilename,&status);
        psrfits_close(&pflower);
        psrfits_open(&pflower,READONLY);
        infits=pflower.fptr;
        fits_movnam_hdu(infits,BINARY_TBL,"SUBINT",0,&status);
        fits_copy_file(infits,outfits,1,0,0,&status);
        fits_copy_header(infits, outfits, &status);
        fits_flush_buffer(outfits,0,&status);
        int dummy;
        char tform[8];
        char tdim[10];
        sprintf(tform,"%dE",pfo.hdr.nchan);
        fits_update_key(outfits, TINT, "NCHAN", &pfo.hdr.nchan, NULL, &status);
        fits_update_key(outfits,TSTRING,"TFORM13",tform,NULL,&status);
        fits_update_key(outfits,TSTRING,"TFORM14",tform,NULL,&status);
        fits_update_key(outfits,TSTRING,"TFORM15",tform,NULL,&status);
        fits_update_key(outfits,TSTRING,"TFORM16",tform,NULL,&status);
        sprintf(tform,"%dB",(pfo.hdr.nsblk*pfo.hdr.nchan*pfo.hdr.npol*pfo.hdr.nbits)/8);
        fits_update_key(outfits, TSTRING, "TFORM17", tform, NULL, &status);
        fits_read_key(outfits, TINT, "NAXIS1", &dummy, NULL, &status);
        dummy=dummy+(pfo.hdr.nsblk*(pfo.hdr.nchan-pflower.hdr.nchan)*pfo.hdr.npol*pfo.hdr.nbits)/8+(pfo.hdr.nchan-pflower.hdr.nchan)*16;
        sprintf(tdim,"(1, %d, 1, %d)",pfo.hdr.nchan,pfo.hdr.nsblk);
        fits_update_key(outfits, TSTRING, "TDIM17",tdim,NULL,&status);
        fits_update_key(outfits, TINT, "NAXIS1", &dummy, NULL, &status);
        dummy=0;
        fits_update_key(outfits, TINT, "NAXIS2", &dummy, NULL, &status);
        fits_movabs_hdu(outfits,1,NULL,&status);
        fits_update_key(outfits,TDOUBLE,"OBSFREQ",&pfo.hdr.fctr,NULL,&status);
        fits_update_key(outfits,TINT,"OBSNCHAN",&pfo.hdr.nchan,NULL,&status);
        fits_update_key(outfits,TDOUBLE,"OBSBW",&pfo.hdr.BW,NULL,&status);
        fits_movnam_hdu(outfits,BINARY_TBL,"SUBINT",0,&status);
        fits_close_file(outfits,&status);
        pfo.filenum++;
        if(psrfits_open(&pfo, READWRITE)!=0)
        {
          fprintf(stderr, "error opening file\n");
          fits_report_error(stderr, pfo.status);
          exit(1);
        } 

        pflower.status=0;
        psrfits_read_subint(&pflower);
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
          int newuppernchan=pfupper.hdr.nchan-upchanskip; //The number of channels to copy from the upper sideband.
          int newlowernchan=pflower.hdr.nchan-lowchanskip; //The number of channels to copy from the lower sideband.
          memcpy(pfo.sub.dat_freqs+2,pfupper.sub.dat_freqs,sizeof(float)*newuppernchan); //The first two frequency channels are skipped
          memcpy(pfo.sub.dat_freqs+newuppernchan+2,pflower.sub.dat_freqs+lowchanskip,sizeof(float)*newlowernchan);
          //Copy weights
          pfo.sub.dat_weights[0]=pfo.sub.dat_weights[1]=0;
          memcpy(pfo.sub.dat_weights+2,pfupper.sub.dat_weights,sizeof(float)*newuppernchan);
          memcpy(pfo.sub.dat_weights+2+newuppernchan,pflower.sub.dat_weights+lowchanskip,sizeof(float)*newlowernchan);
          //Copy offsets
          pfo.sub.dat_offsets[0]=pfo.sub.dat_offsets[1]=pfupper.sub.dat_offsets[0];
          int a,b,ii;
          float upmean=0,lowmean=0;
          for(ii=0;ii<upchanskip;++ii)
            upmean=upmean+*(pfupper.sub.dat_offsets+(pfupper.hdr.nchan-upchanskip+ii));
          upmean=upmean/(float)upchanskip;
          for(ii=0;ii<lowchanskip;++ii)
            lowmean=lowmean+*(pflower.sub.dat_offsets+ii);
          lowmean=lowmean/(float)lowchanskip;
          for(ii=0;ii<pfupper.hdr.nchan;++ii)
            pfupper.sub.dat_offsets[ii]=pfupper.sub.dat_offsets[ii]*(lowmean/upmean);
          memcpy(pfo.sub.dat_offsets+2*pfo.hdr.npol,pfupper.sub.dat_offsets,sizeof(float)*newuppernchan*pfupper.hdr.npol);
          memcpy(pfo.sub.dat_offsets+(newuppernchan+2)*pfo.hdr.npol,pflower.sub.dat_offsets+lowchanskip,sizeof(float)*newlowernchan*pflower.hdr.npol);
          //Copy scales
          upmean=lowmean=0;
          for(ii=0;ii<upchanskip;++ii)
            upmean=upmean+*(pfupper.sub.dat_scales+(pfupper.hdr.nchan-upchanskip+ii));
          upmean=upmean/(float)upchanskip;
          for(ii=0;ii<lowchanskip;++ii)
            lowmean=lowmean+*(pflower.sub.dat_scales+ii);
          lowmean=lowmean/(float)lowchanskip;
          for(ii=0;ii<pfupper.hdr.nchan;++ii)
            pfupper.sub.dat_scales[ii]=pfupper.sub.dat_scales[ii]*(lowmean/upmean);
          pfo.sub.dat_scales[0]=pfo.sub.dat_scales[1]=pfupper.sub.dat_scales[0];
          memcpy(pfo.sub.dat_scales+2*pfo.hdr.npol,pfupper.sub.dat_scales,sizeof(float)*newuppernchan*pfupper.hdr.npol);
          memcpy(pfo.sub.dat_scales+(newuppernchan+2)*pfo.hdr.npol,pflower.sub.dat_scales+lowchanskip,sizeof(float)*newlowernchan*pflower.hdr.npol);
          //Copy the data
          int i=0;
          for(i=0;i<pfo.hdr.nsblk;++i)
          {
            memcpy(pfo.sub.data+i*outoffset+extrachanoffset,pfupper.sub.data+i*upperoffset,numtocopyupper);
            memcpy(pfo.sub.data+i*outoffset+extrachanoffset+numtocopyupper,pflower.sub.data+i*loweroffset+loweroffset_skip,numtocopylower);
          }
          psrfits_write_subint(&pfo);
        }
        else
        {
        }
      }
    } while(pfo.rownum<=pfo.rows_per_file);
    printf("Closing file '%s'\n", pflower.filename);
    fits_close_file(pfupper.fptr,&status);
    printf("Closing file '%s'\n", pfupper.filename);
    fits_close_file(pflower.fptr,&status);
    exit(0);
}
