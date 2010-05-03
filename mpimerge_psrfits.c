#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "psrfits.h"
#include "mpi.h"

#define HDRLEN 14400

int main(int argc, char *argv[])
{
    int ii, nc, ncnp, status=0, statsum=0;
    int numprocs, numbands, myid;
    int *counts, *offsets;
    struct psrfits pf;
    char hostname[100];
    MPI_Status mpistat;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    numbands = numprocs - 1;

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
    }

    // Determine the hostnames of the processes
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
    
    // Basefilenames for the GPU nodes
    if (myid > 0) {
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

        // Send the header to the master proc
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
        nc = pf.hdr.nchan;
        ncnp = pf.hdr.nchan * pf.hdr.npol;
        pf.hdr.orig_nchan *= numbands;
        pf.hdr.nchan *= numbands;
        pf.hdr.fctr = pf.hdr.fctr - 0.5 * pf.hdr.BW + numbands/2.0 * pf.hdr.BW;
        pf.hdr.BW *= numbands;
        pf.sub.bytes_per_subint *= numbands;
        status = psrfits_create(&pf);
    }

    // Open the input PSRFITs files for real
    if (myid > 0) {
        status = psrfits_open(&pf);
        nc = pf.hdr.nchan;
        ncnp = pf.hdr.nchan * pf.hdr.npol;
    }

    // Alloc data buffers for the PSRFITS files
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) * 
                                         pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales  = (float *)malloc(sizeof(float) * 
                                         pf.hdr.nchan * pf.hdr.npol);
    pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);

    // Counts and offsets for MPI_Gatherv
    counts = (int *)malloc(sizeof(int) * numprocs);
    offsets = (int *)malloc(sizeof(int) * numprocs);
    counts[0] = offsets[0] = 0;  //  master sends nothing

    // Now loop over the rows (i.e. subints)...
    do {
        MPI_Barrier(MPI_COMM_WORLD);
        // Read the current subint from each of the "slave" nodes
        if (myid > 0) status = psrfits_read_subint(&pf);
        // Combine statuses of all nodes by summing....
        MPI_Allreduce(&status, &statsum, 1, 
                      MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (statsum) break;
            
        if (myid == 1) { // Send all of the non-band-specific parts to master
            MPI_Send(&pf.sub.tsubint, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.offs, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.lst, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.ra, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.dec, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.glon, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.glat, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.feed_ang, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.pos_ang, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.par_ang, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.tel_az, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&pf.sub.tel_zen, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else if (myid == 0) { // Receive all of the non-data parts
            MPI_Recv(&pf.sub.tsubint, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.offs, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.lst, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.ra, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.dec, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.glon, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.glat, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.feed_ang, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.pos_ang, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.par_ang, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.tel_az, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
            MPI_Recv(&pf.sub.tel_zen, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &mpistat);
        }

        // Now gather the vector quantities...

        // Vectors of length nchan
        for (ii = 1 ; ii < numprocs ; ii++) {
            counts[ii] = nc;
            offsets[ii] = (ii - 1) * nc;
        }
        status = MPI_Gatherv(pf.sub.dat_freqs, nc, MPI_FLOAT, 
                             pf.sub.dat_freqs, counts, offsets, MPI_FLOAT, 
                             0, MPI_COMM_WORLD);
        status = MPI_Gatherv(pf.sub.dat_weights, nc, MPI_FLOAT, 
                             pf.sub.dat_weights, counts, offsets, MPI_FLOAT, 
                             0, MPI_COMM_WORLD);
        
        // Vectors of length nchan * npol
        for (ii = 1 ; ii < numprocs ; ii++) {
            counts[ii] = ncnp;
            offsets[ii] = (ii - 1) * ncnp;
        }
        status = MPI_Gatherv(pf.sub.dat_offsets, ncnp, MPI_FLOAT, 
                             pf.sub.dat_offsets, counts, offsets, MPI_FLOAT, 
                             0, MPI_COMM_WORLD);
        status = MPI_Gatherv(pf.sub.dat_scales, ncnp, MPI_FLOAT, 
                             pf.sub.dat_scales, counts, offsets, MPI_FLOAT, 
                             0, MPI_COMM_WORLD);

        // Vectors of length nchan for the raw data
        for (ii = 1 ; ii < numprocs ; ii++) {
            counts[ii] = nc;
            offsets[ii] = (ii - 1) * nc;
        }
        // Step over all the spectra * polns (Note: this does a "corner turn)
        for (ii = 0 ; ii < pf.hdr.nsblk * pf.hdr.npol ; ii++) {
            status = MPI_Gatherv(pf.sub.data + ii * nc, nc, MPI_UNSIGNED_CHAR, 
                                 pf.sub.data + ii * nc * numbands, 
                                 counts, offsets, MPI_UNSIGNED_CHAR, 
                                 0, MPI_COMM_WORLD);
        }
            
        // Write the new row to the output file
        if (myid == 0)
            status = psrfits_write_subint(&pf);

    } while (statsum == 0);

    // Free the arrays
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);
    free(counts);
    free(offsets);
    
    // Close the files and finalize things
    status = psrfits_close(&pf);
    MPI_Finalize();
    exit(0);
}
