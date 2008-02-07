CC = gcc
CFLAGS = -g -O2 -WALL
PROGS = 
all: $(PROGS)
clean:
	rm -f $(PROG) *~ *.o
