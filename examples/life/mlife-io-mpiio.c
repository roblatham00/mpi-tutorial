/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "mlife-io.h"

/* MPI-IO implementation of checkpoint and restart for MPI Life
 *
 * Data stored in matrix order, with a header consisting of two
 * integers: matrix size in one dimension, and iteration number.
 *
 * Each checkpoint is stored in its own file.
 */

static int MLIFEIO_Type_create_rowblk(int **matrix, int myrows,
				      int dimsz, MPI_Datatype *newtype);
static int MLIFEIO_Type_create_header_and_rowblk(int **matrix,
						 int myrows,
						 int *dimsz_p,
						 int *iter_p,
						 MPI_Datatype *newtype);

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
    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;
    int rank, nprocs;
    int myrows, myoffset;

    MPI_File fh;
    MPI_Status status;

    MPI_Datatype type;
    MPI_Offset myfileoffset;

    char filename[64];

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    myrows   = MLIFE_myrows(dimsz, rank, nprocs);
    myoffset = MLIFE_myrowoffset(dimsz, rank, nprocs);

    snprintf(filename, 63, "%s-%d.chkpt", prefix, iter);

    err = MPI_File_open(mlifeio_comm, filename, amode, info, &fh);
    if (err != MPI_SUCCESS) {
	fprintf(stderr, "Error opening %s.\n", filename);
	return err;
    }

    if (rank == 0) {
	MLIFEIO_Type_create_header_and_rowblk(matrix, myrows, &dimsz,
					      &iter, &type);
	myfileoffset = 0;
    }
    else {
	MLIFEIO_Type_create_rowblk(matrix, myrows, dimsz, &type);
	myfileoffset = ((myoffset * dimsz) + 2) * sizeof(int);
    }

    MPI_Type_commit(&type);
    err = MPI_File_write_at_all(fh, myfileoffset, MPI_BOTTOM, 1,
				    type, &status);
    MPI_Type_free(&type);

    err = MPI_File_close(&fh);
    return err;
}

int MLIFEIO_Restart(char    *prefix,
		    int    **matrix,
		    int      dimsz,
		    int      iter,
		    MPI_Info info)
{
    int err;
    int amode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;
    int rank, nprocs;
    int myrows, myoffset;
    int buf[2]; /* dimsz, iteration */

    MPI_File fh;
    MPI_Status status;

    MPI_Datatype type;
    MPI_Offset myfileoffset;

    char filename[64];

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    myrows   = MLIFE_myrows(dimsz, rank, nprocs);
    myoffset = MLIFE_myrowoffset(dimsz, rank, nprocs);

    snprintf(filename, 63, "%s-%d.chkpt", prefix, iter);

    err = MPI_File_open(mlifeio_comm, filename, amode, info, &fh);
    if (err != MPI_SUCCESS) return err;

    /* check that dimsz matches */
    err = MPI_File_read_at_all(fh, 0, buf, 2, MPI_INT, &status);
    /* TODO: err handling */

    if (buf[0] != dimsz) {
	printf("restart failed.\n");
	/* TODO: FIX THIS ERROR */
	return -1;
    }

    MLIFEIO_Type_create_rowblk(matrix, myrows, dimsz, &type);
    myfileoffset = ((myoffset * dimsz) + 2) * sizeof(int);

    MPI_Type_commit(&type);
    err = MPI_File_read_at_all(fh, myfileoffset, MPI_BOTTOM, 1,
			       type, &status);
    MPI_Type_free(&type);

    err = MPI_File_close(&fh);
    return err;
}

int MLIFEIO_Can_restart(void)
{
    return 1;
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

static int MLIFEIO_Type_create_header_and_rowblk(int **matrix,
						 int myrows,
						 int *dimsz_p,
						 int *iter_p,
						 MPI_Datatype *newtype)
{
    int err;
    int lens[3] = { 1, 1, 1 };
    MPI_Aint disps[3];
    MPI_Datatype types[3];
    MPI_Datatype rowblk;

    MLIFEIO_Type_create_rowblk(matrix, myrows, *dimsz_p, &rowblk);
    
    disps[0] = (MPI_Aint) dimsz_p;
    disps[1] = (MPI_Aint) iter_p;
    disps[2] = (MPI_Aint) MPI_BOTTOM;
    types[0] = MPI_INTEGER;
    types[1] = MPI_INTEGER;
    types[2] = rowblk;

#if 0    
    MPI_Type_create_struct(3, lens, disps, types, newtype);
#endif
    err = MPI_Type_struct(3, lens, disps, types, newtype);

    MPI_Type_free(&rowblk);

    return err;
}
