all: all-redirect

MPI_DIR = 
# Some systems may want to add -static to the MPICC line
MPICC   = /home/gropp/tmp/mpi2-inst/bin/mpicc
MPIF77  = /home/gropp/tmp/mpi2-inst/bin/mpif77
MPIEXEC = /home/gropp/tmp/mpi2-inst/bin/mpiexec
MPIEXEC_NP = -n

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
#CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=199506L
CFLAGS = -g -O2

STDOBJS= 
MPIOBJS= 

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)
.f.o:   $?
	$(MPIF77) $(FFLAGS) -c $(?) -o $(@)

EXECS = 

default:

all-redirect: $(EXECS)

clean::
	rm -f *.o *.ps *.pdf *.exe *.stackdump

run: $(PROG)
	$(MPIEXEC) -n $(NP) $(PROG) $(ARGS)
