/*  This is an abstraction of the pNeo brain simulation program,
    written to illustrate the MPI one-sided operations

    Run with

      pneo <filename>

    where <filename> is a file of connections, one on each line consisting
    of a blank-separated pair of integers in character format.

*/

#include <stdio.h>
#include "mpi.h"
#define MAX_CONNECTIONS 10000
#define MAX(x,y) ((x) > (y) ? x : y)

typedef struct {
    int source;
    int dest;
} connection;

connection connarray[MAX_CONNECTIONS];
int connarray_count;

typedef struct {
    int source;
    int inspike;
} inconnection;

inconnection *inarray;
int inarray_count;

typedef struct {
    int dest;
    int disp;
    int outspike;
} outconnection;

outconnection *outarray;
int outarray_count;

int state;

int numprocs, myrank, maxcell;  /* why global? */
int itercount;
int max_steps = 100;
MPI_Win win;

int  setup_connarray( char *);
void dump_connarray( void );
void setup_local_arrays( void );
void dump_local_arrays( void );
void init_state( void );
void compute_state( void );
void output_spikes( void );

int main(int argc, char *argv[])
{

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    maxcell = setup_connarray(argv[1]);
    /* dump_connarray(); */
    if (maxcell > numprocs-1) {
	printf("This needs to be run with at least %d processes\n", maxcell+1);
	MPI_Finalize();
	exit(-1);
    }

    init_state();
    setup_local_arrays();
    dump_local_arrays();
    printf( "state for rank %d = %d\n", myrank, state);

    /* make input arrays the windows */
    MPI_Win_create(inarray, inarray_count * sizeof(int), sizeof(int),
		   MPI_INFO_NULL, MPI_COMM_WORLD, &win); 


    for (itercount = 0; itercount < max_steps; itercount++ ) {
	if (myrank == 0)
	    printf("iteration %d:\n", itercount);
	compute_state();
	if (state > 4) {
	    output_spikes();
	    state = 0;
	}
    }

    MPI_Win_free(&win);
    MPI_Finalize();
    return 0;
}

int setup_connarray(char *filename)
{
    FILE *confile;
    int i, cellcount, maxpair, maxcell;

    if ((confile = fopen(filename, "r")) == NULL) {
	printf("could not open connection file %s\n", filename);
	return -1;
    }

    maxcell = connarray_count = i = 0;
    inarray_count = outarray_count = 0;

    while ((fscanf(confile, "%d %d", &connarray[i].source, &connarray[i].dest)) == 2) {
	connarray_count++; 
	maxpair = MAX(connarray[i].source, connarray[i].dest);
	if (maxpair > maxcell)
	    maxcell = maxpair;
	if (connarray[i].dest == myrank)
	    inarray_count++;
	if (connarray[i].source == myrank)
	    outarray_count++;
	i++;
    }
    return maxcell;
}

void dump_connarray()
{
    int i;

    for (i = 0; i < connarray_count; i++) {
	printf("%d %d\n", connarray[i].source, connarray[i].dest);
    }
}    

void setup_local_arrays()
{
    int i, j, k, n, dispcnt;

    inarray  = (inconnection *) malloc(inarray_count * sizeof(inconnection));
    outarray = (outconnection *) malloc(outarray_count * sizeof(outconnection));

    for (i = j = k = 0 ; i < connarray_count; i++) {
	if (connarray[i].dest == myrank) {
	    inarray[j].source = connarray[i].source;
	    inarray[j].inspike = 0;
	    j++;
	}
	inarray_count = j;

        if (connarray[i].source == myrank) {
	    dispcnt = 0;
	    for (n = 0; n < connarray_count; n++) {
		if (connarray[n].dest == connarray[i].dest) {
		    if (connarray[n].source == myrank) 
			break;
		    else
			dispcnt++; /* increment counter on destination process i */
		}
	    }
	    outarray[k].dest = connarray[i].dest;
	    outarray[k].disp = dispcnt;
	    outarray[k].outspike = 0;
	    k++;
	}
	outarray_count = k;
    }
}

void dump_local_arrays()
{
    int i;
     
    printf("inarray for process %d:\n", myrank);
    for (i = 0; i < inarray_count; i++) 
	printf("%d %d\n", inarray[i].source, inarray[i].inspike);
    printf("outarray for process %d:\n", myrank);
    for (i = 0; i < outarray_count; i++) 
	printf("%d %d %d\n", outarray[i].dest, outarray[i].disp, outarray[i].outspike);
}

void init_state()
{
    state = (myrank + 5) % numprocs; /* random */
}

void compute_state()
{
    int i;
    int num_incoming = 0;
    
    MPI_Win_fence(0, win);
    for (i = 0; i < inarray_count; i++) {
	num_incoming += inarray[i].inspike;
    }

    state = state + num_incoming + 1;
}

void output_spikes()
{
    int i;
    int spike = 1;

    for (i = 0; i < outarray_count; i++) {
	printf("putting spike from %d to %d in interation %d\n",
	       myrank, outarray[i].dest, itercount);
	MPI_Put(&spike, 1, MPI_INT, outarray[i].dest, outarray[i].disp, 1, MPI_INT, win);
    }
}
