CC = gcc
CFLAGS = -g -O2 -Wall
PROGS = check_guppi_status
OBJS = guppi_status.o guppi_error.o
LIBS = -lcfitsio -lm -lpthread
all: $(PROGS) 
clean:
	rm -f $(PROGS) *~ *.o
.SECONDEXPANSION:
$(PROGS): $$@.c $(OBJS)
	$(CC) $(CFLAGS) $< -o $@ $(OBJS) $(LIBS)
