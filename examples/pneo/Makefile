MPI_DIR ?= /sandbox/lusk/install-mpd-oct25
# MPI_DIR= /usr/local/mpich-1.2.5.2/
MPICC=$(MPI_DIR)/bin/mpicc

default: pneo_fence pneo_pscw

pneo_fence: pneo_fence.o
	$(MPICC) -o $@ pneo_fence.o

pneo_pscw: pneo_pscw.o
	$(MPICC) -o $@ pneo_pscw.o

.c.o:	$?
	$(MPICC) -c $(?) -o $(@)

clean::
	rm -f *.o *~ pneo_fence pneo_pscw

# Enscript command:
# enscript -1rC -E --color --style emacs -f "Helvetica16" \
# -F "Helvetica16" --header='File: $n||Page $% of $=' $(?) \
# -o $(@)

