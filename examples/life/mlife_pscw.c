/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "mlife.h"
#include "mlife-io.h"

static int MLIFE_nextstate(int **matrix, int y, int x);
static int MLIFE_exchange(int **matrix, int mysize, int dimsz, 
                          MPI_Win top_row_win,
                          MPI_Win bottom_row_win,
			  int prev, int next,
                          MPI_Group prev_grp,
                          MPI_Group next_grp);
static int MLIFE_parse_args(int argc, char **argv);

/* options */
static int opt_dimsz = 50, opt_iter = 10, opt_restart_iter = -1;
static char opt_prefix[64] = "mlife";

/* The Life function */
double life(int matrix_size, int ntimes, MPI_Comm comm)
{
    int      rank, nprocs;
    int      next, prev;
    int      i, j, k;
    int      mysize;
    int    **matrix, **temp, **addr;
    double   slavetime, totaltime, starttime;
    int      myoffset;
    MPI_Win top_row_win, bottom_row_win;
    MPI_Win matrix_top_row_win, matrix_bottom_row_win;
    MPI_Win temp_top_row_win, temp_bottom_row_win;
    MPI_Group prev_grp=MPI_GROUP_EMPTY, next_grp=MPI_GROUP_EMPTY, comm_world_grp;

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
    mysize   = MLIFE_myrows(matrix_size, rank, nprocs);
    myoffset = MLIFE_myrowoffset(matrix_size, rank, nprocs);

    /* allocate the memory dynamically for the matrix */
    matrix = (int **)malloc(sizeof(int *)*(mysize+2)) ;
    temp = (int **)malloc(sizeof(int *)*(mysize+2)) ;
    for (i = 0; i < mysize+2; i++) {
	matrix[i] = (int *)malloc(sizeof(int)*(matrix_size+2)) ;
	temp[i] = (int *)malloc(sizeof(int)*(matrix_size+2)) ;
    }

    /* Initialize the boundaries of the life matrix */
    for (j = 0; j < matrix_size+2; j++)
	matrix[0][j] = matrix[mysize+1][j] = temp[0][j] = temp[mysize+1][j] = DIES ;
    for (i = 0; i < mysize+2; i++)
	matrix[i][0] = matrix[i][matrix_size+1] = temp[i][0] = temp[i][matrix_size+1] = DIES ;

    /* Initialize the life matrix */
    for (i = 1; i <= mysize; i++)  {
	srand48((long)(1000^(i + myoffset)));
	for (j = 1; j<= matrix_size; j++)
	    if (drand48() > 0.5)  
		matrix[i][j] = BORN ;
	    else
		matrix[i][j] = DIES ;
    }

    MPI_Win_create(matrix[0], (matrix_size+2)*sizeof(int), 1, MPI_INFO_NULL, 
                   comm, &matrix_top_row_win);
    MPI_Win_create(matrix[mysize+1], (matrix_size+2)*sizeof(int), 1, 
                   MPI_INFO_NULL, comm, &matrix_bottom_row_win);

    MPI_Win_create(&temp[0][0], (matrix_size+2)*sizeof(int), 1, 
                   MPI_INFO_NULL, comm, &temp_top_row_win);
    MPI_Win_create(&temp[mysize+1][0], (matrix_size+2)*sizeof(int), 1, 
                   MPI_INFO_NULL, comm, &temp_bottom_row_win);

    MPI_Comm_group(MPI_COMM_WORLD, &comm_world_grp);

    if (prev != MPI_PROC_NULL)
        MPI_Group_incl(comm_world_grp, 1, &prev, &prev_grp);

    if (next != MPI_PROC_NULL)
        MPI_Group_incl(comm_world_grp, 1, &next, &next_grp);

    MPI_Group_free(&comm_world_grp);

    top_row_win = matrix_top_row_win;
    bottom_row_win = matrix_bottom_row_win;

    /* Play the game of life for given number of iterations */
    starttime = MPI_Wtime() ;
    for (k = 0; k < ntimes; k++)
    {
	MLIFE_exchange(matrix, mysize, matrix_size, top_row_win, 
                       bottom_row_win, prev, next, prev_grp, next_grp);

	/* Calculate new state */
	for (i = 1; i <= mysize; i++) {
	    for (j = 1; j < matrix_size+1; j++) {
		temp[i][j] = MLIFE_nextstate(matrix, i, j);
	    }
	}

	/* Swap the matrices */
	addr = matrix ;
	matrix = temp ;
	temp = addr ;

        /* Swap the win ptrs */
        top_row_win = (top_row_win == matrix_top_row_win) ? 
            temp_top_row_win : matrix_top_row_win;
        bottom_row_win = (bottom_row_win == matrix_bottom_row_win) ? 
            temp_bottom_row_win : matrix_bottom_row_win;

	MLIFEIO_Checkpoint(opt_prefix, matrix, matrix_size, k, MPI_INFO_NULL);
    }

    /* Return the average time taken/processor */
    slavetime = MPI_Wtime() - starttime;
    MPI_Reduce(&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

    MPI_Win_free(&matrix_top_row_win);
    MPI_Win_free(&matrix_bottom_row_win);
    MPI_Win_free(&temp_top_row_win);
    MPI_Win_free(&temp_bottom_row_win);

    return(totaltime/(double)nprocs);
}

int main(int argc, char *argv[])
{
    int myargs[4];
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
    }
    MPI_Bcast(myargs, 4, MPI_INT, 0, MPI_COMM_WORLD);
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
	

int MLIFE_myrows(int dimsz, int rank, int nprocs)
{
    int mysize;

    mysize = dimsz / nprocs;

    if (rank < dimsz % nprocs) {
	mysize += 1;
    }

    return mysize;
}

int MLIFE_myrowoffset(int dimsz, int rank, int nprocs)
{
    int myoffset;

    myoffset = rank * (dimsz / nprocs);

    if (rank > (dimsz % nprocs)) {
	myoffset += (dimsz % nprocs);
    }
    else {
	myoffset += rank;
    }

    return myoffset;
}

static int MLIFE_parse_args(int argc, char **argv)
{
    int ret;

    while ((ret = getopt(argc, argv, "s:i:p:r:")) >= 0)
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
	    default:
		break;
	}
    }

    return 0;
}

static int MLIFE_exchange(int **matrix,
			  int mysize,
			  int dimsz,
			  MPI_Win top_row_win,
                          MPI_Win bottom_row_win,
			  int prev /* rank */,
			  int next /* rank */,
                          MPI_Group prev_grp,
                          MPI_Group next_grp)
{
    int err=MPI_SUCCESS;

    /* Send and receive boundary information */

    MPI_Win_post( prev_grp, 0, top_row_win);
    MPI_Win_post( next_grp, 0, bottom_row_win);
    
    MPI_Win_start( prev_grp, 0, bottom_row_win);
    MPI_Win_start( next_grp, 0, top_row_win);
    
    MPI_Put(&matrix[1][0], dimsz + 2, MPI_INT, prev, 0, dimsz + 2, 
            MPI_INT, bottom_row_win);
    MPI_Put(&matrix[mysize][0], dimsz + 2, MPI_INT, next, 0, dimsz + 2, 
            MPI_INT, top_row_win);
    
    MPI_Win_complete(top_row_win);
    MPI_Win_complete(bottom_row_win);
    
    MPI_Win_wait(top_row_win);
    MPI_Win_wait(bottom_row_win);

    return err;
}
