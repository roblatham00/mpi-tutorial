MPI_DIR = /Users/gropp/installs/mpich2-44/bin
MPICC   = /Users/gropp/installs/mpich2-44/bin/mpicc
PNCDIR  = 

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
#CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE
CFLAGS  = -g -O2

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

