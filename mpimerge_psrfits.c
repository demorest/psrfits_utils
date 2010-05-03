#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "psrfits.h"
#include "mpi.h"

#define HDRLEN 14400

int main(int argc, char *argv[])
{
    int ii, jj, status, statsum;
    int numprocs, numbands, myid, gpuid, *gpuids;
    struct psrfits pf;
    char hostname[100];
    MPI_Status mpistat;
    MPI_Group all_group, gpu_group; 
    MPI_Comm gpu_comm; 

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    numbands = numprocs - 1;
    {
        FILE *hostfile;
        
        hostfile = fopen("/etc/hostname", "r");
        fscanf(hostfile, "%s\n", hostname);
        fclose(hostfile);
        if (hostname != NULL) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0) printf("\n");
            fflush(NULL);
            for (ii = 0 ; ii < numprocs ; ii++) {
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

    gpuids = (int *)malloc(sizeof(int), numbands);
    for (ii = 1 ; ii <= numbands ; ii++) gpuids[ii] = ii;

    // Extract the original group handle
    MPI_Comm_group(MPI_COMM_WORLD, &all_group); 
    // Add the GPU nodes to the new group
    MPI_Group_incl(all_group, numbands, gpuids, &gpu_group);
    free(gpuids);

    // Create new communicator

    MPI_Comm_create(MPI_COMM_WORLD, gpu_group, &gpu_comm); 
    MPI_Group_rank(gpu_group, &gpuid); 


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
        sprintf(pf.basefilename, "/data/gpu/partial/%s/%s", 
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
        psrfitsfile = fopen(filenm, "r");
        fread(&hdr, 1, HDRLEN, psrfitsfile);
        fclose(psrfitsfile);

        // Send the heafer to the master proc
        MPI_Send(hdr, HDRLEN, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

    } else if (myid == 0) {
        FILE *psrfitsfile;
        char hdr[HDRLEN], tmpfilenm[80];

        // Receive the header info from proc 1
        MPI_Recv(hdr, HDRLEN, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &mpistat);

        // Now write that header to a temp file
        strcpy(tmpfilenm, "mpi_merge_psrfits.XXXXXX");
        mkstemp(tmpfilenm); 
        psrfitsfile = fopen(tmpfilenm, "w");
        fwrite(&hdr, 1, HDRLEN, psrfitsfile);
        fclose(psrfitsfile); 
        sprintf(pf.filename, "%s", tmpfilenm);

        // And read the key information into a PSRFITS struct
        status = psrfits_open(&pf);
        status = psrfits_close(&pf);
        remove(tmpfilenm);

        // Now create the output PSTFITS file
        strcpy(pf.basefilename, argv[1]);
        pf.filenum = 0;
        pf.multifile = 1;
        pf.filename[0] = '\0';
        pf.hdr.orig_nchan *= numbands;
        pf.hdr.nchan *= numbands;
        pf.hdr.fctr = pf.hdr.fctr - 0.5 * pf.hdr.BW + numbands/2.0 * pf.hdr.BW;
        pf.hdr.BW *= numbands;
        pf.sub.bytes_per_subint *= numbands;
        status = psrfits_create(&pf);
    }

    // Open the input PSRFITs files for real
    if (myid > 0) status = psrfits_open(&pf);

    // Alloc data buffers for the PSRFITS files
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_offsets = (float *)malloc(sizeof(float)
                                         * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales  = (float *)malloc(sizeof(float)
                                         * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);

    // Now loop over the rows...
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid > 0) {
        do {
            status = psrfits_read_subint(&pf);
            // Combine statuses of GPU nodes by summing....
            MPI_Allreduce(&status, &statsum, 1, MPI_INT, MPI_SUM, gpu_comm); 
            if (statsum) {
                // send statsum to master
                break;
            }
            
            // Step over all the spectra * polns
            for (jj = 0 ; jj < pf.hdr.nsblk * pf.hdr.npol ; jj++) {
            }
            
        } while (statsum == 0);
    } else { // master
        do {
            // Write the new row to the output file
            status = psrfits_write_subint(&pf);
            if (status) break;
        } while (statsum == 0);
    }


    // Loop over each spectra and re-arrange them so they are a single band

    status = psrfits_close(&pf);
    MPI_Finalize();
    exit(0);
}
