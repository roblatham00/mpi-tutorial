all: all-redirect

# ?= is a GNU make extension that means "set only if not set".  
MPI_DIR ?= /sandbox/lusk/install-mpd-oct25
# MPI_DIR= /usr/local/mpich-1.2.5.2/
MPICC=$(MPI_DIR)/bin/mpicc
MPIEXEC = $(MPI_DIR)/bin/mpiexec

EXECS = pneo_fence pneo_pscw
all-redirect: $(EXECS)

default: $(EXECS) 

pneo_fence: pneo_fence.o
	$(MPICC) -o $@ pneo_fence.o

pneo_pscw: pneo_pscw.o
	$(MPICC) -o $@ pneo_pscw.o

.c.o:	$?
	$(MPICC) -c $(?) -o $(@)

clean::
	rm -f *.o *~ pneo_fence pneo_pscw

run: $(PROG)
	$(MPIEXEC) -n $(NP) $(PROG) pneofile

# Enscript command:
# enscript -1rC -E --color --style emacs -f "Helvetica16" \
# -F "Helvetica16" --header='File: $n||Page $% of $=' $(?) \
# -o $(@)

