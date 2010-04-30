#include <stdlib.h>
#include <stdio.h>
#include "psrfits.h"
#include "mpi.h"

#define HDRLEN 14400

int main(int argc, char *argv[])
{
    int ii, jj, status;
    int numprocs, myid;
    struct psrfits pf, pftemp;
    char hostname[20];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    {
        FILE *hostfile;
        char tmpname[100];
        
        hostfile = chkfopen("/etc/hostname", "r");
        fscanf(hostfile, "%s\n", tmpname);
        hostname = (char *) calloc(strlen(tmpname) + 1, 1);
        memcpy(hostname, tmpname, strlen(tmpname));
        fclose(hostfile);
        if (hostname != NULL) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0) printf("\n");
            fflush(NULL);
            for (ii=0 ; ii < numprocs ; ii++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (myid == ii)
                    printf("Process %3d is on machine %s\n", myid, hostname);
                fflush(NULL);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            fflush(NULL);
        }
    }

    // Call usage() if we have no command line arguments
    
    if (argc == 1) {
        if (myid == 0) {
            printf("usage:  mpimerge_psrfits basefilename\n\n");
        }
        MPI_Finalize();
        exit(1);
    }
    
    if (myid == 0) { // Master proc only
        printf("\n\n");
        printf("      MPI Search-mode PSRFITs Combiner\n");
        printf("              by Scott M. Ransom\n\n");
    } else { // all other procs
        sprintf(pf.basefilename, "/data/gpu/%s/partial/%s", 
                hostname, argv[1]);
        
    }

    // Initialize some key parts of the PSRFITS structure
    pf.tot_rows = pf.N = pf.T = pf.status = 0;
    pf.filenum = 1;
    pf.filename[0] = '\0';

    if (myid == 1) {
        FILE *psrfitsfile;
        char hdr[HDRLEN], filenm[200];

        // Read the header info
        sprintf(filenm, "%s_0001.fits", pf.basefilename);
        psrfitsfile = chkfopen(filenm, 'r');
        chkfread(&hdr, 1, HDRLEN, psrfitsfile);
        fclose(psrfitsfile);

        // Send the heafer to the master proc
        MPI_Send(hdr, HDRLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

    } else if (myid == 0) {
        FILE *psrfitsfile;
        char hdr[HDRLEN], tmpfilenm[80];

        // Receive the header info from proc 1
        MPI_Recv(hdr, HDRLEN, MPI_CHAR, 1, 0, MPI_COMM_WORLD);

        // Now write that header to a temp file
        strcpy(tmpfilenm, "mpi_merge_psrfits.XXXXXX");
        mkstemp(tmpfilenm); 
        psrfitsfile = chkfopen(tmpfilenm, 'w');
        chkfrite(&hdr, 1, HDRLEN, psrfitsfile);
        fclose(psrfitsfile);

        // And read the key information into a PSRFITS struct
        status = psrfits_open(&pf);
        status = psrfits_close(&pf)
        
        // Now create the output PSTFITS file
        strcpy(pf.basefilename, argv[1]);
        pf.filenum = 0;
        pf.multifile = 1;
        status = psrfits_create(&pf);
        status = psrfits_close(&pf)
    }

    if (myid > 0) {
        status = psrfits_open(&pf);
        status = psrfits_close(&pf)
    }

    MPI_Finalize();
    exit(0);
}
