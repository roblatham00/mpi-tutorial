MPI_DIR ?= /usr/local/mpich-1.2.5.2/
# MPI_DIR=/home/rross/cvs/build-mpich2
MPICC=$(MPI_DIR)/bin/mpicc
PNCDIR=/usr/local/parallel-netcdf-0.9.4/

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)

default: csrio-mpiio.o

all: mlife

test1: test1.o csrio-stdout.o
	$(MPICC) -o $@ $(LDFLAGS) test1.o csrio-stdout.o $(LIBS)

test1-mpiio: test1.o csrio-mpiio.o
	$(MPICC) -o $@ $(LDFLAGS) test1.o csrio-mpiio.o $(LIBS)

clean::
	rm -f *.o *~ test1 test1-mpiio

# Enscript command:
# enscript -1rC -E --color --style emacs -f "Helvetica16" \
# -F "Helvetica16" --header='File: $n||Page $% of $=' $(?) \
# -o $(@)

