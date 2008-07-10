/* polyco.c
 * routines to read/use polyco.dat
 */

#include "polyco.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int read_one_pc(FILE *f, struct polyco *pc) {

    int i, j;
    char *rv;
    int ret;
    char buf[90];
    /* Read in polyco chunk */
    rv = fgets(buf, 90, f);
    if (rv==NULL) { return(-1); }
    strncpy(pc->psr, &buf[0], 10);  pc->psr[10] = '\0';
    pc->mjd = atoi(&buf[31]);
    pc->fmjd = atof(&buf[39]);
    if ((rv=strchr(pc->psr, ' '))!=NULL) { *rv='\0'; }
    rv = fgets(buf,90,f);
    if (rv==NULL) { return(-1); }
    pc->rphase = fmod(atof(&buf[0]),1.0);
    pc->f0 = atof(&buf[20]);
    pc->nsite = atoi(&buf[42]);
    pc->nmin = atoi(&buf[43]);
    pc->nc = atoi(&buf[50]);
    pc->rf = atof(&buf[55]);
    for (i=0; i<pc->nc/3 + (pc->nc%3)?1:0; i++) {
        rv=fgets(buf, 90, f);
        if (rv==NULL) { return(-1); }
        for (j=0; j<90; j++) { if (buf[j]=='D' || buf[j]=='d') buf[j]='e'; }
        ret=sscanf(buf, "%lf %lf %lf", 
                &(pc->c[3*i]), &(pc->c[3*i+1]), &(pc->c[3*i+2]));
        if (ret!=3) { return(-1); }
    }

    return(0);

}

int read_pc(FILE *f, struct polyco *pc, const char *psr, int mjd, double fmjd) {

    /* Read through until we get to right psr, mjd */
    int done=0, nomatch=0;
    int i, j;
    char *rv;
    int ret;
    char buf[90];
    float tdiff;
    while (!done) {
        /* Read in polyco chunk */
        rv = fgets(buf, 90, f);
        if (rv==NULL) { done=1; nomatch=1; continue; }
        strncpy(pc->psr, &buf[0], 10);  pc->psr[10] = '\0';
        pc->mjd = atoi(&buf[31]);
        pc->fmjd = atof(&buf[39]);
        if ((rv=strchr(pc->psr, ' '))!=NULL) { *rv='\0'; }
        rv = fgets(buf,90,f);
        pc->rphase = fmod(atof(&buf[0]),1.0);
        pc->f0 = atof(&buf[20]);
        pc->nsite = atoi(&buf[42]);
        pc->nmin = atoi(&buf[43]);
        pc->nc = atoi(&buf[50]);
        pc->rf = atof(&buf[55]);
        for (i=0; i<pc->nc/3 + (pc->nc%3)?1:0; i++) {
            rv=fgets(buf, 90, f);
            if (rv==NULL) { return(-1); }
            for (j=0; j<90; j++) { if (buf[j]=='D' || buf[j]=='d') buf[j]='e'; }
            ret=sscanf(buf, "%lf %lf %lf", 
                    &(pc->c[3*i]), &(pc->c[3*i+1]), &(pc->c[3*i+2]));
            if (ret!=3) { return(-1); }
        }
        /* check for correct psr - null psrname matches any */
        if (psr!=NULL) { if (strcmp(pc->psr, psr)!=0) { continue; } }
        tdiff = 1440.0*((double)(mjd-pc->mjd) + (fmjd-pc->fmjd));
        if (fabs(tdiff) > (float)pc->nmin/2.0) { continue; }
        done=1;
    }

    return(-1*nomatch);

}

/* Select appropriate polyco set */
int select_pc(const struct polyco *pc, int npc, const char *psr,
        int imjd, double fmjd) {
    int ipc;
    char *tmp = psr;
    if (tmp[0]=='J' || tmp[0]=='B') tmp++;
    for (ipc=0; ipc<npc; ipc++) {
        if (psr!=NULL) { if (strcmp(pc[ipc].psr,psr)!=0) { continue; } }
        if (pc_out_of_range(&pc[ipc],imjd,fmjd)==0) { break; }
    }
    if (ipc<npc) { return(ipc); }
    return(-1);
}

/* Compute pulsar phase given polyco struct and mjd */
double psr_phase(const struct polyco *pc, int mjd, double fmjd, double *freq) {
    double dt = 1440.0*((double)(mjd-pc->mjd)+(fmjd-pc->fmjd));
    int i;
    double phase = pc->c[pc->nc-1];
    double f = 0.0;
    if (fabs(dt)>(double)pc->nmin/2.0) { return(-1.0); }
    for (i=pc->nc-1; i>0; i--) {
        phase = dt*(phase) + pc->c[i-1];
        f = dt*(f) + (double)i*pc->c[i];
    }
    f = pc->f0 + (1.0/60.0)*f;
    phase += pc->rphase + dt*60.0*pc->f0;
    if (freq!=NULL) { *freq = f; }
    return(phase);
}

double psr_fdot(const struct polyco *pc, int mjd, double fmjd, double *fdot) {
    double dt = 1440.0*((double)(mjd-pc->mjd)+(fmjd-pc->fmjd));
    if (fabs(dt)>(double)pc->nmin/2.0) { return(-1.0); }
    double fd=0.0;
    int i;
    for (i=pc->nc-1; i>1; i--) {
        fd = dt*fd + ((double)i)*((double)i-1.0)*pc->c[i];
    }
    fd /= 60.0;
    if (fdot!=NULL) { *fdot=fd; }
    return(fd);
}

double psr_phase_avg(const struct polyco *pc, int mjd, 
        double fmjd1, double fmjd2) {
    double dt1 = 1440.0*((double)(mjd-pc->mjd)+(fmjd1-pc->fmjd));
    double dt2 = 1440.0*((double)(mjd-pc->mjd)+(fmjd2-pc->fmjd));
    if (fabs(dt1)>(double)pc->nmin/2.0) { return(-1.0); }
    if (fabs(dt2)>(double)pc->nmin/2.0) { return(-1.0); }
    double pavg;
    int i;
    double tmp1=0.0, tmp2=0.0;
    for (i=pc->nc-1; i>=0; i--) {
        tmp1 = dt1*tmp1 + pc->c[i]/((double)i+1.0);
        tmp2 = dt2*tmp2 + pc->c[i]/((double)i+1.0);
    }
    tmp1 *= dt1; tmp2 *= dt2;
    pavg = (tmp2-tmp1)/(dt2-dt1) + pc->rphase + (dt1+dt2)*60.0*pc->f0/2.0;
    return(pavg);
}

int pc_range_check(const struct polyco *pc, int mjd, double fmjd) {
    double dt;
    dt = (double)(mjd - pc->mjd) + (fmjd - pc->fmjd);
    dt *= 1440.0;
    if (dt < -1.0*(double)pc->nmin/2.0) { return(-1); }
    else if (dt > (double)pc->nmin/2.0) { return(1); }
    else { return(0); }
}

int pc_out_of_range(const struct polyco *pc, int mjd, double fmjd) {
    double dt;
    dt = (double)(mjd - pc->mjd) + (fmjd - pc->fmjd);
    dt *= 1440.0;
    if (fabs(dt)>(double)pc->nmin/2.0) { return(1); }
    return(0);
}
