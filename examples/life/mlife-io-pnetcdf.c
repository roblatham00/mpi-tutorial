/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pnetcdf.h>

#include "mlife-io.h"

/* Parallel netCDF implementation of checkpoint and restart for MPI Life
 *
 * Data stored in a 2D variable in matrix order.
 *
 * Each checkpoint is stored in its own file.
 *
 * TODO: ADD ERROR CHECKING.
 */

static int MLIFEIO_Type_create_rowblk(int **matrix, int myrows,
				      int dimsz, MPI_Datatype *newtype);

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
		       int      rows,
		       int      cols,
		       int      iter,
		       MPI_Info info)
{
    int err;
    int cmode = 0;
    int rank, nprocs;
    int myrows, myoffset;

    int ncid, varid, coldim, rowdim, dims[2];
    MPI_Offset start[2];
    MPI_Offset count[2];
    int i, j, *buf;

    MPI_Datatype type;
    MPI_Offset myfileoffset;

    char filename[64];

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    myrows   = MLIFE_myrows(rows, rank, nprocs);
    myoffset = MLIFE_myrowoffset(rows, rank, nprocs);

    snprintf(filename, 63, "%s-%d.nc", prefix, iter);

    err = ncmpi_create(mlifeio_comm, filename, cmode, info, &ncid);
    if (err != 0) {
	fprintf(stderr, "Error opening %s.\n", filename);
	return err;
    }

    ncmpi_def_dim(ncid, "col", cols, &coldim);
    ncmpi_def_dim(ncid, "row", rows, &rowdim);
    dims[0] = coldim; /* TODO: WHICH ONE CHANGES MOST QUICKLY? */
    dims[1] = rowdim;
    ncmpi_def_var(ncid, "matrix", NC_INT, 2, dims, &varid);

    /* TODO: store iteration and dimsz as attributes? */

    ncmpi_enddef(ncid);

    start[0] = 0; /* col start */
    start[1] = myoffset; /* row start */
    count[0] = cols;
    count[1] = myrows;

#if 0
    /* if flexible data mode interface were finished... */
    MLIFEIO_Type_create_rowblk(matrix, myrows, cols, &type);
    MPI_Type_commit(&type);

    ncmpi_put_vara_all(ncid, varid, start, count, MPI_BOTTOM, 1, type);

    MPI_Type_free(&type);
#else
    /* need to pack and then write, or write in rows. */
    buf = (int *) malloc(cols * myrows * sizeof(int));

    for (i=0; i < myrows; i++) {
	for (j=0; j < cols; j++) {
	    buf[(i*cols) + j] = matrix[i+1][j];
	}
    }

    ncmpi_put_vara_int_all(ncid, varid, start, count, buf);
    free(buf);
#endif

    ncmpi_close(ncid);
    return err;
}

int MLIFEIO_Restart(char    *prefix,
		    int    **matrix,
		    int      rows,
		    int      cols,
		    int      iter,
		    MPI_Info info)
{
    int err = MPI_SUCCESS;
    int rank, nprocs;
    int myrows, myoffset;
    int buf[2]; /* dimsz, iteration */

    MPI_Datatype type;
    MPI_Offset myfileoffset;

    char filename[64];

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    myrows   = MLIFE_myrows(rows, rank, nprocs);
    myoffset = MLIFE_myrowoffset(rows, rank, nprocs);

    snprintf(filename, 63, "%s-%d.chkpt", prefix, iter);

#if 0
    err = MPI_File_open(mlifeio_comm, filename, amode, info, &fh);
    if (err != MPI_SUCCESS) return err;

    /* check that dimsz matches */
    err = MPI_File_read_at_all(fh, 0, buf, 2, MPI_INT, &status);
    /* TODO: err handling */

    if (buf[0] != rows) {
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
#endif
    return err;
}

int MLIFEIO_Can_restart(void)
{
    return 1;
}

static int MLIFEIO_Type_create_rowblk(int         **matrix,
				      int           myrows,
				      int           cols,
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
    err = MPI_Type_create_indexed_block(myrows, cols, rowptrs,
					MPI_INTEGER, newtype);

    free(rowptrs);

    return err;
}
