CFLAGS = -g -O3 -Wall
PROGS = check_guppi_databuf check_guppi_status clean_guppi_shmem test_udp_recv 
OBJS  = guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
        guppi_params.o \
        hget.o hput.o
THREAD_PROGS = test_net_thread
THREAD_OBJS  = guppi_net_thread.o guppi_rawdisk_thread.o 
LIBS = -lcfitsio -lm -lpthread
all: $(PROGS) $(THREAD_PROGS)
clean:
	rm -f $(PROGS) $(THREAD_PROGS) *~ *.o
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS)
$(THREAD_PROGS): $$@.c $(THREAD_OBJS) $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(THREAD_OBJS) $(OBJS) $(LIBS)
