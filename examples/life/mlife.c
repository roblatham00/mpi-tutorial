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
static int MLIFE_parse_args(int argc, char **argv);

/* options */
static int opt_rows = 50, opt_cols = 70, opt_iter = 10, opt_restart_iter = -1;
static char opt_prefix[64] = "mlife";

/* The Life function */
double life(int rows, int cols, int ntimes, MPI_Comm comm)
{
    int      rank, nprocs;
    int      next, prev;
    int      i, j, k;
    int      mysize;
    int    **matrix, **temp, **addr;
    int     *matrixdata, *tempdata;
    double   slavetime, totaltime, starttime;
    int      myoffset;

    /* Determine size and my rank in communicator */
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    /* Set neighbors */
    if (rank == 0) 
	prev = MPI_PROC_NULL;
    else
	prev = rank-1;
    if (rank == nprocs - 1)
	next = MPI_PROC_NULL;
    else
	next = rank+1;

    /* Determine my part of the matrix, row-block distribution */
    mysize   = MLIFE_myrows(rows, rank, nprocs);
    myoffset = MLIFE_myrowoffset(rows, rank, nprocs);

    /* allocate the memory dynamically for the matrix */
    matrix = (int **)malloc(sizeof(int *)*(mysize+2)) ;
    temp = (int **)malloc(sizeof(int *)*(mysize+2)) ;
    matrixdata = (int *) malloc((mysize + 2)*(cols + 2) * sizeof(int));
    tempdata = (int *) malloc((mysize + 2)*(cols + 2) * sizeof(int));

    /* set up pointers for convenience */
    matrix[0] = matrixdata;
    temp[0]   = tempdata;
    for (i = 1; i < mysize+2; i++) {
	matrix[i] = matrix[i-1] + cols + 2;
	temp[i]   = temp[i-1] + cols + 2;
    }

    /* Initialize the boundaries of the life matrix */
    for (j = 0; j < cols+2; j++)
	matrix[0][j] = matrix[mysize+1][j] = temp[0][j] = temp[mysize+1][j] = DIES ;
    for (i = 0; i < mysize+2; i++)
	matrix[i][0] = matrix[i][cols+1] = temp[i][0] = temp[i][cols+1] = DIES ;

    /* Initialize the life matrix */
    for (i = 1; i <= mysize; i++)  {
	srand48((long)(1000^(i + myoffset)));
	for (j = 1; j<= cols; j++)
	    if (drand48() > 0.5)  
		matrix[i][j] = BORN ;
	    else
		matrix[i][j] = DIES ;
    }

    MLIFE_exchange_init(comm, &matrix[0][0], &temp[0][0], mysize, cols, prev,
			next);

    /* Play the game of life for given number of iterations */
    starttime = MPI_Wtime() ;
    for (k = 0; k < ntimes; k++)
    {
	MLIFE_exchange(matrix, mysize, cols);

	/* Calculate new state */
	for (i = 1; i <= mysize; i++) {
	    for (j = 1; j < cols+1; j++) {
		temp[i][j] = MLIFE_nextstate(matrix, i, j);
	    }
	}

	/* Swap the matrices */
	addr = matrix ;
	matrix = temp ;
	temp = addr ;

	MLIFEIO_Checkpoint(opt_prefix, matrix, rows, cols, 
			   k, MPI_INFO_NULL);
    }

    /* Return the average time taken/processor */
    slavetime = MPI_Wtime() - starttime;
    MPI_Reduce(&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    MLIFE_exchange_finalize();
    free(matrix);
    free(temp);
    free(matrixdata);
    free(tempdata);

    return(totaltime/(double)nprocs);
}

int main(int argc, char *argv[])
{
    int rank;
    double time;
  
    MPI_Init (&argc, &argv);
    MLIFEIO_Init(MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MLIFE_parse_args(argc, argv);

    /* Call the life routine */
    time = life(opt_rows, opt_cols, opt_iter, MPI_COMM_WORLD);

    /* Print the total time taken */
    if (rank == 0)
	printf("[%d] Life finished in %lf seconds of calculation\n",
	       rank, time/100.0);

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
	

int MLIFE_myrows(int rows, int rank, int nprocs)
{
    int mysize;

    mysize = rows / nprocs;

    if (rank < rows % nprocs) {
	mysize += 1;
    }

    return mysize;
}

int MLIFE_myrowoffset(int rows, int rank, int nprocs)
{
    int myoffset;

    myoffset = rank * (rows / nprocs);

    if (rank > (rows % nprocs)) {
	myoffset += (rows % nprocs);
    }
    else {
	myoffset += rank;
    }

    return myoffset;
}

static int MLIFE_parse_args(int argc, char **argv)
{
    int ret;
    int rank;
    int myargs[5];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
	while ((ret = getopt(argc, argv, "x:y:i:p:r:")) >= 0)
	{
	    switch(ret) {
		case 'x':
		    opt_cols = atoi(optarg);
		    break;
		case 'y':
		    opt_rows = atoi(optarg);
		    break;
		case 'i':
		    opt_iter = atoi(optarg);
		    break;
		case 'r':
		    opt_restart_iter = atoi(optarg);
		case 'p':
		    strncpy(opt_prefix, optarg, 63);
		    break;
		default:
		    break;
	    }
	}

	myargs[0] = opt_rows;
	myargs[1] = opt_cols;
	myargs[2] = opt_iter;
	myargs[3] = opt_restart_iter;
	myargs[4] = strlen(opt_prefix) + 1;
    }

    MPI_Bcast(myargs, 5, MPI_INT, 0, MPI_COMM_WORLD);
    opt_rows = myargs[0];
    opt_cols = myargs[1];
    opt_iter = myargs[2];

    MPI_Bcast(opt_prefix, myargs[4], MPI_CHAR, 0, MPI_COMM_WORLD);

    return 0;
}
