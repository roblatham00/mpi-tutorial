MPI_DIR=/usr/local/mpich-1.2.5.2/
# MPI_DIR=/home/rross/cvs/build-mpich2
MPICC=$(MPI_DIR)/bin/mpicc
PNCDIR=/usr/local/parallel-netcdf-0.9.4/

CFLAGS+= -g -O0 -Wall

STDOBJS= mlife.o mlife-io-stdout.o
MPIOBJS= mlife.o mlife-io-mpiio.o
PNCOBJS= mlife.o mlife-io-pnetcdf.o

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)

default: mlife

all: mlife

mlife-io-pnetcdf.o: mlife-io-pnetcdf.c
	$(MPICC) $(CFLAGS) -I$(PNCDIR)/include -c $(?) -o $(@)

mlife: $(STDOBJS)
	$(MPICC) -o $@ $(LDFLAGS) $(STDOBJS) $(LIBS)

mlife-mpiio: $(MPIOBJS)
	$(MPICC) -o $@ $(LDFLAGS) $(MPIOBJS) $(LIBS)

mlife-pnetcdf: $(PNCOBJS)
	$(MPICC) -o $@ $(LDFLAGS) -L$(PNCDIR)/lib $(PNCOBJS) $(LIBS) -lpnetcdf

clean::
	rm -f *.o mlife *~
