MPI_DIR=/usr/local/mpich-1.2.5.2/
# MPI_DIR=/home/rross/cvs/build-mpich2
MPICC=$(MPI_DIR)/bin/mpicc

CFLAGS+= -g -O0 -Wall

OBJS= mlife.o mlife-io-stdout.o

.c.o:	$?
	$(MPICC) $(CFLAGS) -c $(?) -o $(@)

default: mlife

all: mlife

mlife: $(OBJS)
	$(MPICC) -o $@ $(LDFLAGS) $(OBJS) $(LIBS)

clean::
	rm -f $(OBJS) mlife *~
