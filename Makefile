CC = gcc
CFLAGS = -g -O2 -Wall
PROGS = 
all: $(PROGS) guppi_status.o
clean:
	rm -f $(PROG) *~ *.o
