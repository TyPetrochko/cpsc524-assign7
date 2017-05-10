CC = icc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99
FC = ifort
FFLAGS = -g -O3 -xHost -fno-alias
MPICC = mpicc

all: parallel

parallel: parallel.o timing.o util.o
	$(MPICC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f mvcp-ga *.o

