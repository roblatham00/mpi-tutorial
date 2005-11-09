all: all-redirect

# ?= is a GNU make extension that means "set only if not set".  
MPI_DIR ?= /usr/local/mpich2-1.0.3/
# Some systems may want to add -static to the MPICC line
MPICC   = $(MPI_DIR)/bin/mpicc 
MPIF77  = $(MPI_DIR)/bin/mpif77
MPIEXEC = $(MPI_DIR)/bin/mpiexec

# _XOPEN_SOURCE ensures that drand48 and friends are defined in stdlib.h
CFLAGS+= -g -O0 -Wall -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=199506L

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
