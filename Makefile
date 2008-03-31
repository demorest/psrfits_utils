CFLAGS = -g -O3 -Wall -I/users/sransom/64bit/include
PY_INCLUDE = /usr/include/python2.4
PROGS = check_guppi_databuf check_guppi_status clean_guppi_shmem test_udp_recv test_psrfits
OBJS  = guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
        guppi_params.o write_psrfits.o hget.o hput.o sla.o
THREAD_PROGS = test_net_thread
THREAD_OBJS  = guppi_net_thread.o guppi_rawdisk_thread.o 
LIBS = -L/users/sransom/64bit/lib -lcfitsio -L/home/sransom/presto/lib -lsla -lm -lpthread
all: $(PROGS) $(THREAD_PROGS) _possem.so
clean:
	rm -f $(PROGS) $(THREAD_PROGS) *~ *.o test_psrfits_0*.fits _possem.so
INSTALL_DIR = ../bin
install: $(PROGS) $(THREAD_PROGS)
	mkdir -p $(INSTALL_DIR) && \
	cp -f -p $(PROGS) $(THREAD_PROGS) $(INSTALL_DIR)
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS)
$(THREAD_PROGS): $$@.c $(THREAD_OBJS) $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(THREAD_OBJS) $(OBJS) $(LIBS)
_possem.so:
	$(CC) -I$(PY_INCLUDE) -o $@ -shared -fPIC -O2 possem_wrap.c
