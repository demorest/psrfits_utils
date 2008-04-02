
#ifndef _POLYCO_H
#define _POLYCO_H

#include <stdio.h>
#include <stdlib.h>

struct polyco {
    char psr[15];
    int mjd;
    double fmjd;
    double rphase;
    double f0;
    int nsite;
    int nmin;
    int nc;
    float rf;
    double c[15];
};

int read_one_pc(FILE *f, struct polyco *pc);
int read_pc(FILE *f, struct polyco *pc, const char *psr, int mjd, double fmjd);
double psr_phase(struct polyco *pc, int mjd, double fmjd, double *freq);
double psr_fdot(struct polyco *pc, int mjd, double fmjd, double *fdot);
double psr_phase_avg(struct polyco *pc, int mjd, double fmjd1, double fmjd2);
int pc_range_check(struct polyco *pc, int mjd, double fmjd);
int pc_out_of_range(struct polyco *pc, int mjd, double fmjd);
#endif
