/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "mlife.h"
#include "mlife-io.h"

static int MLIFE_nextstate(int **matrix, int y, int x);
static int MLIFE_exchange(int **matrix, int LRows, int dimsz, MPI_Comm comm,
			  int prev, int next, int left, int right );
static int MLIFE_parse_args(int argc, char **argv);

/* options */
static int opt_dimsz = 50, opt_iter = 10, opt_restart_iter = -1;
static char opt_prefix[64] = "mlife";

/* decomposition of the domain in 2-d */
static int opt_pcols = 0, opt_prows = 0;

/* The Life function */
double life(int matrix_size, int ntimes, MPI_Comm comm)
{
    int      rank, nprocs;
    int      next, prev, left, right;
    int      i, j, k;
    int    **matrix, **temp, **addr;
    int      *matrixData, *tempData;
    double   slavetime, totaltime, starttime;
    int      LRows, LCols;     /* Number of rows and columns in local mesh
				  (without halo) */
    int      GFirstRow, GFirstCol; /* Index of the first row and column
				      of this mesh in the global mesh, 
				      relative to the [1][1] entry of the 
				      local mesh */

    /* Determine size and my rank in communicator */
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Set neighbors and determine my part of the matrix, 
       row-block distribution */
    MLIFE_MeshDecomp( rank, nprocs, 
		      matrix_size, matrix_size,
		      &left, &right, &prev, &next, 
		      &LRows, &LCols, &GFirstRow, &GFirstCol );

    /* allocate the memory dynamically for the matrix */

    /* Allocate a single block for the entire (local) matrix.
       This simplifies the process of moving data between processes */
    
    matrixData = (int *)malloc( sizeof(int) * (LRows+2)*(LCols+2) );
    tempData   = (int *)malloc( sizeof(int) * (LRows+2)*(LCols+2) );
       
    matrix = (int **)malloc(sizeof(int *)*(LRows+2)) ;
    temp   = (int **)malloc(sizeof(int *)*(LRows+2)) ;
    matrix[0] = matrixData;
    temp[0]   = tempData;
    for (i = 1; i < LRows+2; i++) {
	matrix[i] = matrix[i-1] + LCols + 2;
	temp[i]   = temp[i-1]   + LCols + 2;
    }

    /* Initialize the boundaries of the life matrix */
    for (j = 0; j < LCols+2; j++)
	matrix[0][j] = matrix[LRows+1][j] = temp[0][j] = temp[LRows+1][j] = DIES ;
    for (i = 0; i < LRows+2; i++)
	matrix[i][0] = matrix[i][LCols+1] = temp[i][0] = temp[i][LCols+1] = DIES ;

    /* Initialize the life matrix */
    for (i = 1; i <= LRows; i++)  {
	srand48((long)(1000^(i + GFirstRow-1)));
	/* advance to the random number generator to the first *owned* cell */
	for (j=1; j<GFirstCol; j++) {    
	    (void)drand48();
	}

	for (j = 1; j<= LCols; j++)
	    if (drand48() > 0.5)  
		matrix[i][j] = BORN ;
	    else
		matrix[i][j] = DIES ;
    }

    /* Play the game of life for given number of iterations */
    starttime = MPI_Wtime() ;
    for (k = 0; k < ntimes; k++)
    {
	MLIFE_exchange(matrix, LRows, LCols, comm, prev, next, left, right );

	/* Calculate new state */
	for (i = 1; i <= LRows; i++) {
	    for (j = 1; j < LCols+1; j++) {
		temp[i][j] = MLIFE_nextstate(matrix, i, j);
	    }
	}

	/* Swap the matrices */
	addr = matrix ;
	matrix = temp ;
	temp = addr ;

	MLIFEIO_Checkpoint(opt_prefix, matrix, matrix_size, matrix_size, 
			   k, MPI_INFO_NULL);
    }

    /* Return the average time taken/processor */
    slavetime = MPI_Wtime() - starttime;
    MPI_Reduce(&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    return(totaltime/(double)nprocs);
}

int main(int argc, char *argv[])
{
    int myargs[6];
    int rank, N, iters;
    double time ;
  
    MPI_Init (&argc, &argv);
    MLIFEIO_Init(MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

    /* Parse command line arguments */
    /* TODO: PUT THIS ALL IN ONE FUNCTION */
    if (rank == 0) {
	MLIFE_parse_args(argc, argv);
	myargs[0] = opt_dimsz;
	myargs[1] = opt_iter;
	myargs[2] = opt_restart_iter;
	myargs[3] = strlen(opt_prefix) + 1;
	myargs[4] = opt_pcols;
	myargs[5] = opt_prows;
    }
    MPI_Bcast(myargs, 6, MPI_INT, 0, MPI_COMM_WORLD);
    N     = myargs[0];
    iters = myargs[1];

    MPI_Bcast(opt_prefix, myargs[3], MPI_CHAR, 0, MPI_COMM_WORLD);

    /* Call the life routine */
    time = life(N, iters, MPI_COMM_WORLD);

    /* Print the total time taken */
    if (rank == 0)
	printf("[%d] Life finished in %lf seconds\n",rank,time/100);

    MLIFEIO_Finalize();
    MPI_Finalize();

    return 0;
}

static int MLIFE_nextstate(int **matrix,
			   int y /* row */,
			   int x /* col */)
{
    int sum;

    sum = matrix[y-1][x-1] + matrix[y-1][x] + matrix[y-1][x+1] +
	matrix[y][x-1] + matrix[y][x+1] +
	matrix[y+1][x-1] + matrix[y+1][x] + matrix[y+1][x+1];

    if (sum < 2 || sum > 3) {
	return DIES;
    }
    else if (sum == 3) {
	return BORN;
    }
    else {
	return matrix[y][x]; /* no change */
    }
}
	
/* Compute coordinates of this patch, given the rank and size of the
   global mesh.  Also return the neighbors */

void MLIFE_MeshDecomp( int rank, int nprocs, 
		       int GRows, int GCols,
		       int *leftP, int *rightP, int *topP, int *bottomP, 
		       int *LRowsP, int *LColsP, 
		       int *GFirstRowP, int *GFirstColP )
{
    int dims[2];
    int top, bottom, left, right;
    int prow, pcol;
    int firstrow, lastrow, firstcol, lastcol;

    /* Form the decomposition */
    dims[0] = opt_prows;
    dims[1] = opt_pcols;
    MPI_Dims_create( nprocs, 2, dims );

    /* Compute the cartesian coords of this process; number across
       rows changing column by 1 changes rank by 1) */
    prow = rank / dims[1];
    pcol = rank % dims[1];
    
    /* Compute the neighbors */
    left = right = top = bottom = MPI_PROC_NULL;
    if (prow > 0) {
	top = rank - dims[1];
    }
    if (pcol > 0) {
	left = rank - 1;
    }
    if (prow < dims[0]-1) {
	bottom = rank + dims[1];
    }
    if (pcol < dims[1]-1) {
	right = rank + 1;
    }
    if (leftP) {
	/* Allow a NULL for leftP to suppress all of the neighbor 
	   information */
	*leftP = left; *rightP = right; *topP = top; *bottomP = bottom;
    }
    
    /* Compute the decomposition of the global mesh.
       These are for the "active" part of the mesh, and range from 1 to GRows
       by 1 to GCols */
    firstcol = 1 + pcol * (GCols / dims[1]);
    firstrow = 1 + prow * (GRows / dims[0]);
    if (pcol == dims[1] - 1) {
	lastcol = GCols;
    }
    else {
	lastcol  = 1 + (pcol + 1) * (GCols / dims[1]) - 1;
    }
    if (prow == dims[0] - 1) {
	lastrow = GRows;
    }
    else {
	lastrow = 1 + (prow + 1) * (GRows / dims[0]) - 1;
    }

    *LRowsP     = lastrow - firstrow + 1;
    *LColsP     = lastcol - firstcol + 1;
    *GFirstRowP = firstrow;
    *GFirstColP = firstcol;
}

static int MLIFE_parse_args(int argc, char **argv)
{
    int ret;

    while ((ret = getopt(argc, argv, "s:i:p:r:x:y:")) >= 0)
    {
	switch(ret) {
  	    case 's':
		opt_dimsz = atoi(optarg);
		break;
	    case 'i':
		opt_iter = atoi(optarg);
		break;
	    case 'r':
		opt_restart_iter = atoi(optarg);
	    case 'p':
		strncpy(opt_prefix, optarg, 63);
		break;
	    case 'x':
	        opt_pcols = atoi(optarg);
	        break;
	    case 'y':
		opt_prows = atoi(optarg);
		break;
	    default:
		break;
	}
    }

    return 0;
}

static int MLIFE_exchange(int **matrix,
			  int LRows,
			  int LCols,
			  MPI_Comm comm,
			  int prev /* rank */,
			  int next /* rank */,
			  int left /* rank */, 
			  int right /* rank */ )
{
    int err;
    MPI_Request reqs[4];
    MPI_Status  statuses[4];
    static MPI_Datatype vectype = MPI_DATATYPE_NULL;

    /* Send and receive boundary information */
    /* TODO: POST IRECVS BEFORE ISENDS? */
    /* TODO: ERROR CHECKING? */

    if (vectype == MPI_DATATYPE_NULL) {
	MPI_Type_vector( LRows, 1, LCols+2, MPI_INT, &vectype );
    }
    /* First, move the left, right edges */
    MPI_Isend( &matrix[1][1], 1, vectype, left, 0, comm, reqs );
    MPI_Irecv( &matrix[1][0], 1, vectype, left, 0, comm, reqs+1 );
    MPI_Isend( &matrix[1][LCols], 1, vectype, right, 0, comm, reqs+2 );
    MPI_Irecv( &matrix[1][LCols+1], 1, vectype, right, 0, comm, reqs+3 );
    err = MPI_Waitall( 4, reqs, statuses );

    /* Now move the top, bottom edges (including the diagonal points) */
    MPI_Isend(&matrix[1][0], LCols + 2, MPI_INT, prev, 0, comm, reqs);
    MPI_Irecv(&matrix[0][0], LCols + 2, MPI_INT, prev, 0, comm, reqs+1);
    MPI_Isend(&matrix[LRows][0], LCols + 2, MPI_INT, next, 0, comm, reqs+2);
    MPI_Irecv(&matrix[LRows+1][0], LCols + 2, MPI_INT, next, 0, comm, reqs+3);

    err = MPI_Waitall(4, reqs, statuses);


    return err;
}
