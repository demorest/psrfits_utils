OPT64 = /users/sransom/64bit
ifeq ($(shell hostname),beef) 
    OPT64 = /opt/64bit
endif
ifeq ($(shell hostname),tofu) 
    OPT64 = /opt/64bit
endif
CFLAGS = -O3 -Wall -DFOLD_USE_INTRINSICS -I$(OPT64)/include
#CFLAGS = -g -Wall -I$(OPT64)/include
PROGS = check_guppi_databuf check_guppi_status clean_guppi_shmem \
	test_udp_recv test_psrfits test_psrfits_read fold_psrfits \
	fix_psrfits_polyco
OBJS  = guppi_status.o guppi_databuf.o guppi_udp.o guppi_error.o \
	guppi_params.o guppi_time.o guppi_thread_args.o \
	write_psrfits.o read_psrfits.o misc_utils.o \
	fold.o polyco.o hget.o hput.o sla.o downsample.o
THREAD_PROGS = test_net_thread guppi_daq_fold
THREAD_OBJS  = guppi_net_thread.o guppi_rawdisk_thread.o \
	       guppi_psrfits_thread.o guppi_fold_thread.o \
	       guppi_null_thread.o
LIBS = -L$(OPT64)/lib -lcfitsio -L$(PRESTO)/lib -lsla -lm -lpthread
all: $(PROGS) $(THREAD_PROGS) guppi_daq psrfits_subband psrfits_subband_orig
clean:
	rm -f $(PROGS) $(THREAD_PROGS) guppi_daq psrfits.tgz *~ *.o test_psrfits_0*.fits
INSTALL_DIR = ../bin
install: $(PROGS) $(THREAD_PROGS) guppi_daq psrfits_subband
	mkdir -p $(INSTALL_DIR) && \
	cp -f -p $(PROGS) $(THREAD_PROGS) guppi_daq psrfits_subband $(INSTALL_DIR)
psrfits.tgz: psrfits.h read_psrfits.c write_psrfits.c polyco.c polyco.h \
	guppi_PSRFITS_v3.4_fold_template.txt \
	guppi_PSRFITS_v3.4_search_template.txt
	tar cvzf $@ $^
guppi_daq: guppi_daq.c guppi_daq_cmd.o $(THREAD_OBJS) $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ guppi_daq_cmd.o $(THREAD_OBJS) $(OBJS) $(LIBS)
find_dropped_blocks: find_dropped_blocks.o 
	$(CC) $(CFLAGS) $< -o $@ -L$(OPT64)/lib -lcfitsio -lm
psrfits_subband: psrfits_subband.c psrfits_subband_cmd.o $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ psrfits_subband_cmd.o $(OBJS) $(LIBS)
psrfits_subband_orig: psrfits_subband_orig.c psrfits_subband_cmd.o $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ psrfits_subband_cmd.o $(OBJS) $(LIBS)
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS) $(THREAD_LIBS)
$(THREAD_PROGS): $$@.c $(THREAD_OBJS) $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(THREAD_OBJS) $(OBJS) $(LIBS)
