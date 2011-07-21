all: all-redirect

# ?= is a GNU make extension that means "set only if not set".  
MPI_DIR = /Users/gropp/installs/mpich2-44/bin
MPICC   = /Users/gropp/installs/mpich2-44/bin/mpicc
MPIEXEC = /Users/gropp/installs/mpich2-44/bin/mpiexec
MPIEXEC_NP = -n

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
	$(MPIEXEC) $(MPIEXEC_NP) $(NP) $(PROG) pneofile

# Enscript command:
# enscript -1rC -E --color --style emacs -f "Helvetica16" \
# -F "Helvetica16" --header='File: $n||Page $% of $=' $(?) \
# -o $(@)

