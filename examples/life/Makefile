all: all-redirect

# ?= is a GNU make extension that means "set only if not set".  
MPI_DIR ?= /usr/local/mpich-1.2.5.2/
# MPI_DIR=/home/rross/cvs/build-mpich2
MPICC=$(MPI_DIR)/bin/mpicc
MPIEXEC = $(MPI_DIR)/bin/mpiexec
PNCDIR ?= /usr/local/parallel-netcdf-0.9.4/

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE

STDOBJS= mlife.o mlife-io-stdout.o mlife-pt2pt.o
MPIOBJS= mlife.o mlife-io-mpiio.o mlife-pt2pt.o
PNCOBJS= mlife.o mlife-io-pnetcdf.o mlife-pt2pt.o

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)

EXECS = mlife mlife-mpiio mlife2d

default: mlife

all-redirect: $(EXECS)

mlife-io-pnetcdf.o: mlife-io-pnetcdf.c
	$(MPICC) $(CFLAGS) -I$(PNCDIR)/include -c $(?) -o $(@)

mlife: $(STDOBJS)
	$(MPICC) -o $@ $(LDFLAGS) $(STDOBJS) $(LIBS)

mlife-mpiio: $(MPIOBJS)
	$(MPICC) -o $@ $(LDFLAGS) $(MPIOBJS) $(LIBS)

mlife-pnetcdf: $(PNCOBJS)
	$(MPICC) -o $@ $(LDFLAGS) -L$(PNCDIR)/lib $(PNCOBJS) $(LIBS) -lpnetcdf

mlife2d: mlife2d.o mlife2d-io-stdout.o
	$(MPICC) -o $@ $(LDFLAGS) mlife2d.o mlife2d-io-stdout.o $(LIBS)

clean::
	rm -f *.o mlife *~

run: $(PROG)
	$(MPIEXEC) -n $(NP) $(PROG) $(ARGS)
