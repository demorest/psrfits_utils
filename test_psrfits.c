#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "write_psrfits.h"

void dec2hms(char *out, double in, int sflag) {
    int sign = 1;
    char *ptr = out;
    int h, m;
    double s;
    if (in<0.0) { sign = -1; in = fabs(in); }
    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    s = in;
    if (sign==1 && sflag) { *ptr='+'; ptr++; }
    else if (sign==-1) { *ptr='-'; ptr++; }
    sprintf(ptr, "%2.2d:%2.2d:%07.4f", h, m, s);
}

int main(int argc, char *argv[]) {
    int status = 0, ii;
    double dtmp;
    fitsfile *fptr;
    struct psrfits_hdrinfo pfh;
    struct psrfits_subint pfs;
    
    // First, set values for our hdrinfo structure
    strcpy(pfh.filename, "test_psrfits.fits");
    strcpy(pfh.observer, "John Doe");
    strcpy(pfh.source, "Cool PSR A");
    strcpy(pfh.frontend, "L-band");
    strcpy(pfh.project_id, "GBT09A-001");
    strcpy(pfh.date_obs, "2010-01-01T05:15:30.000");
    strcpy(pfh.poln_type, "LIN");
    strcpy(pfh.track_mode, "TRACK");
    strcpy(pfh.cal_mode, "OFF");
    strcpy(pfh.feed_mode, "FA");
    pfh.dt = 0.000050;
    pfh.fctr = 1400.0;
    pfh.nchan = 256;
    pfh.BW = 800.0;
    pfh.orig_df = pfh.df = pfh.BW / pfh.nchan;
    pfh.ra2000 = 302.0876876;
    dec2hms(pfh.ra_str, pfh.ra2000/15.0, 0);
    pfh.dec2000 = -3.456987698;
    dec2hms(pfh.dec_str, pfh.dec2000, 1);
    pfh.azimuth = 123.123;
    pfh.zenith_ang = 23.0;
    pfh.beam_FWHM = 0.25;
    pfh.start_lst = 10000.0;
    pfh.start_sec = 25000.82736876;
    pfh.scanlen = 200;
    pfh.start_day = 55000;
    pfh.scan_number = 3;
    pfh.orig_nchan = pfh.nchan;
    pfh.rcvr_polns = 2;
    pfh.npol = 1;
    pfh.summed_polns = 1;
    pfh.start_spec = 0;
    pfh.nbits = 8;
    pfh.nsblk = 200;
    pfh.MJD_epoch = 55555.123123123123123123L;  // Note the "L" for long double

    // Now set values for our subint structure
    pfs.tsubint = pfh.nsblk * pfh.dt;
    pfs.rownum = 0;
    pfs.offs = (pfs.rownum + 0.5) * pfs.tsubint;
    pfs.lst = pfh.start_lst;
    pfs.ra = pfh.ra2000;
    pfs.dec = pfh.dec2000;
    // Need to fix these.  Link with SLALIB?
    pfs.glon = 0.0;
    pfs.glat = 0.0;
    pfs.feed_ang = 0.0;
    pfs.pos_ang = 0.0;
    pfs.par_ang = 0.0;
    pfs.tel_az = pfh.azimuth;
    pfs.tel_zen = pfh.zenith_ang;
    pfs.nbits = pfh.nbits;
    pfs.nchan = pfh.nchan;
    pfs.npol = pfh.npol;
    pfs.nsblk = pfh.nsblk;
    pfs.bytes_per_subint = (pfs.nbits * pfs.nchan * pfs.npol * pfs.nsblk) / 8;
    pfs.FITS_typecode = TBYTE;  // 11 = byte
    pfs.max_rows = 100;  // Need to set this based on PSRFITS_MAXFILELEN

    // Create and initialize the subint arrays
    pfs.dat_freqs = (float *)malloc(sizeof(float) * pfs.nchan);
    pfs.dat_weights = (float *)malloc(sizeof(float) * pfs.nchan);
    dtmp = pfh.fctr - 0.5 * pfh.BW + 0.5 * pfh.df;
    for (ii = 0 ; ii < pfs.nchan ; ii++) {
        pfs.dat_freqs[ii] = dtmp + ii * pfh.df;
        pfs.dat_weights[ii] = 1.0;
    }
    pfs.dat_offsets = (float *)malloc(sizeof(float) * pfs.nchan * pfs.npol);
    pfs.dat_scales = (float *)malloc(sizeof(float) * pfs.nchan * pfs.npol);
    for (ii = 0 ; ii < pfs.nchan * pfs.npol ; ii++) {
        pfs.dat_offsets[ii] = 0.0;
        pfs.dat_scales[ii] = 1.0;
    }
 
    pfs.data = (unsigned char *)malloc(pfs.bytes_per_subint);
    for (ii = 0 ; ii < pfs.bytes_per_subint ; ii++) {
        pfs.data[ii] = ii % 256;
    }

    // Create the PSRFITS file
    psrfits_create_searchmode(&fptr, &pfh, &status);

    // Now write several subints worth of data
    psrfits_write_subint(fptr, &pfs, &status);
    pfs.rownum++;
    pfs.offs = (pfs.rownum + 0.5) * pfs.tsubint;

    psrfits_write_subint(fptr, &pfs, &status);
    pfs.rownum++;
    pfs.offs = (pfs.rownum + 0.5) * pfs.tsubint;

    psrfits_write_subint(fptr, &pfs, &status);

    // Close the file
    fits_close_file(fptr, &status);
    printf("status = %d\n", status);

    exit(0);
}
