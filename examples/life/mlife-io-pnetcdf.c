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

    myrows   = MLIFE_myrows(dimsz, rank, nprocs);
    myoffset = MLIFE_myrowoffset(dimsz, rank, nprocs);

    snprintf(filename, 63, "%s-%d.nc", prefix, iter);

    err = ncmpi_create(mlifeio_comm, filename, cmode, info, &ncid);
    if (err != 0) {
	fprintf(stderr, "Error opening %s.\n", filename);
	return err;
    }

    ncmpi_def_dim(ncid, "col", dimsz, &coldim);
    ncmpi_def_dim(ncid, "row", dimsz, &rowdim);
    dims[0] = coldim; /* TODO: WHICH ONE CHANGES MOST QUICKLY? */
    dims[1] = rowdim;
    ncmpi_def_var(ncid, "matrix", NC_INT, 2, dims, &varid);

    /* TODO: store iteration and dimsz as attributes? */

    ncmpi_enddef(ncid);

    start[0] = 0; /* col start */
    start[1] = myoffset; /* row start */
    count[0] = dimsz;
    count[1] = myrows;

#if 0
    /* if flexible data mode interface were done... */
    MLIFEIO_Type_create_rowblk(matrix, myrows, dimsz, &type);
    MPI_Type_commit(&type);

    ncmpi_put_vara_all(ncid, varid, start, count, MPI_BOTTOM, 1, type);

    MPI_Type_free(&type);
#else
    /* need to pack and then write, or write in rows. */
    buf = (int *) malloc(dimsz * myrows * sizeof(int));

    for (i=0; i < myrows; i++) {
	for (j=0; j < dimsz; j++) {
	    buf[(i*dimsz) + j] = matrix[i+1][j];
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
