CC = icpc

CFLAGS= -g -O3 -xHost -mkl -fno-alias -qopenmp
#CFLAGS= -g -O0 -xHost -mkl -fno-alias -qopenmp

OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt 

TARGETS = mvcp-ga
TARGETOBJECTS = mvcp-ga.o

.SUFFIXES: .o .c

all: $(TARGETS)


$(TARGETS): $(TARGETOBJECTS)
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $(OPTFLAGS) -o $@ $<

clean: 
	rm -f $(TARGETOBJECTS) $(TARGETS) *.optrpt


