CFLAGS = -g -O3 -Wall
PROGS = check_guppi_status clean_guppi_shmem test_udp_recv \
	check_guppi_databuf test_net_thread
OBJS = guppi_net_thread.o guppi_rawdisk_thread.o\
       guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
       guppi_params.o \
       hget.o hput.o
LIBS = -lcfitsio -lm -lpthread
all: $(PROGS) 
clean:
	rm -f $(PROGS) *~ *.o
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS)
