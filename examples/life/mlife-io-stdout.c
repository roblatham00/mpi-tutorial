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

/* stdout implementation of checkpoint (no restart) for MPI Life
 *
 * Data output in matrix order: spaces represent dead cells,
 * '*'s represent live ones.
 *
 */

static int MLIFEIO_Type_create_rowblk(int **matrix, int myrows,
				      int dimsz, MPI_Datatype *newtype);
static void MLIFEIO_Row_print(int *data, int dimsz, int rownr);

static MPI_Comm mlifeio_comm = MPI_COMM_NULL;

int MLIFEIO_Init(MPI_Comm comm)
{
    int err;

    err = MPI_Comm_dup(comm, &mlifeio_comm);

    return err;
}

int MLIFEIO_Finalize(void)
{
    int err;

    err = MPI_Comm_free(&mlifeio_comm);

    return err;
}

int MLIFEIO_Checkpoint(char    *prefix,
		       int    **matrix,
		       int      dimsz,
		       int      iter,
		       MPI_Info info)
{
    int err;
    int rank, nprocs;
    int myrows, myoffset;

    MPI_Datatype type;

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    myrows   = MLIFE_myrows(dimsz, rank, nprocs);
    myoffset = MLIFE_myrowoffset(dimsz, rank, nprocs);

    if (rank == 0) {
	int i, procrows, totrows;

	printf("%c# Iteration %d\n### Begin Data ###\n", 12, iter);

	/* print rank 0 data first */
	for (i=0; i < myrows; i++) {
	    MLIFEIO_Row_print(&matrix[i][1], dimsz, i);
	}
	totrows = myrows;

	/* receive and print others' data */
	for (i=1; i < nprocs; i++) {
	    MPI_Status status;
	    int j, *data;

	    procrows = MLIFE_myrows(dimsz, i, nprocs);
	    data = (int *) malloc(procrows * dimsz * sizeof(int));

	    MPI_Recv(data, procrows * dimsz, MPI_INT, i, 1, mlifeio_comm,
		     &status);

	    for (j=0; i < procrows; j++) {
		MLIFEIO_Row_print(&data[i * dimsz], dimsz, totrows + j);
	    }
	    totrows += procrows;

	    free(data);
	}

	printf("### End Data ###\n");
    }
    else {
	/* send all data to rank 0 */

	MLIFEIO_Type_create_rowblk(matrix, myrows, dimsz, &type);
	MPI_Type_commit(&type);
	err = MPI_Send(MPI_BOTTOM, 1, type, 0, 1, mlifeio_comm);
	MPI_Type_free(&type);
    }

    return err;
}

static void MLIFEIO_Row_print(int *data,
			      int dimsz,
			      int rownr)
{
    int i;

    printf("%3d: ", rownr);
    for (i=0; i < dimsz; i++) {
	printf("%c", (data[i] == BORN) ? '*' : ' ');
    }
    printf("\n");
}

int MLIFEIO_Restart(char    *prefix,
		    int    **matrix,
		    int      dimsz,
		    int      iter,
		    MPI_Info info)
{
    return MPI_ERR_IO;
}

int MLIFEIO_Can_restart(void)
{
    return 0;
}

static int MLIFEIO_Type_create_rowblk(int         **matrix,
				      int           myrows,
				      int           dimsz,
				      MPI_Datatype *newtype)
{
    int i, err;
    MPI_Aint *rowptrs;

    rowptrs = (MPI_Aint *) malloc(myrows * sizeof(MPI_Aint));
    /* error handling */

    /* TODO: USE AN HVECTOR PLUS AN HINDEXED INSTEAD */

    /* create type describing rows, skipping boundary cells */
    for (i = 0; i < myrows; i++) {
	/* TODO: IS THIS A BAD HABIT :)? */
	rowptrs[i] = ((MPI_Aint) &matrix[i+1][1]) / sizeof(int);
    }
    err = MPI_Type_create_indexed_block(myrows, dimsz, rowptrs,
					MPI_INTEGER, newtype);

    free(rowptrs);

    return err;
}
