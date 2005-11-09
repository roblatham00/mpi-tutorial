all: all-redirect

# ?= is a GNU make extension that means "set only if not set".  
MPI_DIR ?= /usr/local/mpich2-1.0.3/
# Some systems may want to add -static to the MPICC line
MPICC   = $(MPI_DIR)/bin/mpicc 
MPIF77  = $(MPI_DIR)/bin/mpif77
MPIEXEC = $(MPI_DIR)/bin/mpiexec
PNCDIR ?= /usr/local/parallel-netcdf-1.0.0/

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
# _POSIX_C_SOURCE=199506L and HAVE_NANOSLEEP should be set if
# nanosleep is available
CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=199506L -DHAVE_NANOSLEEP

STDOBJS= mlife.o mlife-io-stdout.o mlife-pt2pt.o
MPIOBJS= mlife.o mlife-io-mpiio.o mlife-pt2pt.o
PNCOBJS= mlife.o mlife-io-pnetcdf.o mlife-pt2pt.o

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)
.f.o:   $?
	$(MPIF77) $(FFLAGS) -c $(?) -o $(@)

EXECS = mlife mlife-mpiio mlife-fence mlife-pscw mlife2d 

default: mlife

all-redirect: $(EXECS)

mlife-io-pnetcdf.o: mlife-io-pnetcdf.c
	$(MPICC) $(CFLAGS) -I$(PNCDIR)/include -c $(?) -o $(@)

mlife: $(STDOBJS)
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) $(STDOBJS) $(LIBS)

mlife-mpiio: $(MPIOBJS)
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) $(MPIOBJS) $(LIBS)

mlife-pnetcdf: $(PNCOBJS)
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) -L$(PNCDIR)/lib $(PNCOBJS) $(LIBS) -lpnetcdf

mlife2d: mlife2d.o mlife2d-io-stdout.o mlife2d-pt2pt.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife2d.o mlife2d-io-stdout.o mlife2d-pt2pt.o $(LIBS)

mlife-fence: mlife.o mlife-io-stdout.o mlife-fence.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife.o mlife-io-stdout.o mlife-fence.o $(LIBS)

mlife-pscw: mlife.o mlife-io-stdout.o mlife-pscw.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife.o mlife-io-stdout.o mlife-pscw.o $(LIBS)

mlife2d-fence: mlife2d.o mlife2d-io-stdout.o mlife2d-fence.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife2d.o mlife2d-io-stdout.o mlife2d-fence.o $(LIBS)

mlife2d1-fence: mlife2d.o mlife2d-io-stdout1.o mlife2d-fence.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife2d.o mlife2d-io-stdout1.o mlife2d-fence.o $(LIBS)

mlife2d1: mlife2d.o mlife2d-io-stdout1.o mlife2d-pt2pt.o
	$(MPICC) $(CFLAGS) -o $@ $(LDFLAGS) mlife2d.o mlife2d-io-stdout1.o mlife2d-pt2pt.o $(LIBS)

#
# Fortran versions (not yet complete)
mlife2df: mlife2df.o mlife2df-io-stdout.o mlife2df-pt2pt.o
	$(MPIF77) $(FFLAGS) -o $@ $(LDFLAGS) mlife2df.o mlife2df-io-stdout.o mlife2df-pt2pt.o $(LIBS)
mlife2df-fence: mlife2df.o mlife2df-io-stdout.o mlife2df-fence.o
	$(MPIF77) $(FFLAGS) -o $@ $(LDFLAGS) mlife2df.o mlife2df-io-stdout.o mlife2df-fence.o $(LIBS)

clean::
	rm -f *.o mlife *~ *.ps *.pdf mlife-mpiio mlife-pnetcdf mlife mlife2d mlife2d-fence mlife-pscw *.chkpt *.exe *.stackdump

run: $(PROG)
	$(MPIEXEC) -n $(NP) $(PROG) $(ARGS)
