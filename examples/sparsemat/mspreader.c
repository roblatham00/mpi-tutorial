/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

/* 
 * A simple program to read a sparse matrix, stored in 
 * compressed-sparse-row (CSR) format.
 *
 * CSR format is defined by 3 arrays and one integer:
 * n - number of rows in the matrix (integer)
 * ia[n+1] - index of data for each row.  Last (n+1st) element gives the 
 * total length of the arrays a and ja (integer)
 * a[nz]   - values of matrix entries (double)
 * ja[nz]  - column number for corresponding a entry
 * All index values are 1-origin (this format was developed for Fortran)
 *
 * For example, the sparse matrix
 *      ( 1  2  0
 *        3  0  0
 *        0  4  5 )
 * is stored as:
 * n = 3
 * ia = ( 1, 3, 4, 6 )
 * a  = ( 1, 2, 3, 4, 5 )
 * ja = ( 1, 2, 1, 2, 3 )
 *
 * The number of elements on the ith row is ia[i+1] - ia[i] (which is why ia
 * has n+1, not just n entries).
 *
 * A file format such as the Boeing-Harwell exchange format will often 
 * include the number of non-zeros separately from the ia[n+1] value.
 *
 * For the purposes of our example, we will store the data in binary (native)
 * representation (later, we can look at external32)
 *
 * title (char, 80 bytes, fixed size)
 * n (int)
 * nz (int)
 * ia[i], i=1,...,n+1 (int array)
 * ja[i], i=1,...,nz  (int array)
 * a[i],  i=1,...,nz  (double array)
 *
 * A simple parallel reader does the following
 *  read_all( title, n, nz )   - All processes read the header
 *  compute location of entries to read:
 *    row range = (n/comm_size) * comm_rank to (n/comm_size) * (comm_rank + 1)-1 *  read row range of ia into a local ia array (ialocal)
 *  compute elements of file to read based on the indices in ialocal:
 *    JAOffset = sizeofHeader + (n+1) * sizeof(int) + (ialocal[o]-1) * sizeof(int)
 *    AOffset = sizeofHeader + (n+1) * sizeof(int) + nz * sizeof(int) + 
 *                (ialocal[o]-1) * sizeof(double)
 *    nelements = ialocal[lastLocalRow+1] - ialocal[0]
 *    read nelements ints at JAOffset
 *    read nelements doubles at AOffset
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

static MPI_Comm csrio_comm = MPI_COMM_NULL;

int CSRIO_Init(MPI_Comm comm)
{
    int err;

    err = MPI_Comm_dup(comm, &csrio_comm);

    return err;
}

int CSRIO_Finalize(void)
{
    int err;

    err = MPI_Comm_free(&csrio_comm);

    return err;
}


/* TODO: move info to init? */
int CSRIO_Read_header(char *filename,
                      char **title_p, /* TODO: pass in 80 char buf instead? */
                      int *nr_p,
                      int *nz_p,
                      MPI_Info info)
{
    int err = 0, ioerr;
    char *buf;
    int nrnz[2];

    int rank;
    int amode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

    MPI_Datatype type;
    int lens[3];
    MPI_Aint disps[3];
    MPI_Datatypes types[3];

    buf = (char *) malloc(80 * sizeof(char));
    if (buf != NULL) return MPI_ERR_UNKNOWN;

    MPI_Comm_rank(csrio_comm, &rank);

    /* often it is faster for only one process to open/access/close
     * the file when only a small amount of data is going to be
     * accessed.
     */
    if (rank == 0) {
        MPI_File fh;
        MPI_Status status;

        err = MPI_File_open(MPI_COMM_SELF, filename, amode, info, &fh);
        if (err == MPI_SUCCESS) {
            err = MPI_File_read_at(fh, 0, buf, 80, MPI_CHAR, &status);
        }
        if (err == MPI_SUCCESS) {
            err = MPI_File_read_at(fh, 80, nrnz, 2, MPI_INT, &status);
        }

        MPI_File_close(&fh);
    }
    ioerr = err;
    
    /* define a struct that describes all our data */
    lens[0] = 80;
    lens[1] = 2;
    lens[2] = 1;
    disps[0] = buf; /* TODO: USE MPI_ADDRESS? */
    disps[1] = nrnz;
    disps[2] = &ioerr;
    types[0] = MPI_CHAR;
    types[1] = MPI_INT;
    types[2] = MPI_INT;

    MPI_Type_struct(3, lens, disps, types, &type);

    /* broadcast the header data to everyone */
    err = MPI_Bcast(MPI_BOTTOM, 1, type, 0, csrio_comm);
    if (err == MPI_SUCCESS && ioerr == MPI_SUCCESS) {
        *nr_p = nrnz[0];
        *nz_p = nrnz[1];
        *title_p = buf; /* app frees */

        return MPI_SUCCESS;
    }
    else {
        free(buf);

        return (err != MPI_SUCCESS) ? err : ioerr;
    }
}

/* CSRIO_Read_rows
 *
 * Parameters:
 * n         - number of rows in matrix
 * nz        - number of nonzero values in matrix
 * row_start - first row to read (0-origin)
 * row_end   - last row to read
 * my_ia     - pointer to local memory for storing row start indices
 * my_ja_p   - pointer to local memory for storing column indices
 * my_a_p    - data values
 */
int CSRIO_Read_rows(char *filename,
		    int n,
		    int nz,
		    int row_start,
		    int row_end,
		    int *my_ia,
		    int **my_ja_p,
		    double **my_a_p,
		    MPI_Info info)
{
    int err;

    MPI_Aint my_ia_off, my_ja_off, my_a_off;

    int next_row_ia;
    int lens[2];
    MPI_Aint disps[2];

    int amode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;
    MPI_File fh;
    MPI_Status status;
    MPI_Datatype type;

    err = MPI_File_open(csrio_comm, filename, amode, info, &fh);
    if (err != MPI_SUCCESS) {
	return err;
    }

    my_ia_off = 80 * sizeof(char) + 2 * sizeof(int) + row_start * sizeof(int);

    /* must read one more row start to calculate number of elements */
    lens[0] = row_end - row_start;
    lens[1] = 1;
    disps[0] = my_ia;
    disps[1] = &next_row_ia;

    MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
    MPI_Type_commit(&type);

    err = MPI_File_read_at_all(fh, my_ia_off, MPI_BOTTOM, 1, type, &status);
    if (err != MPI_SUCCESS) {
	return err;
    }
    MPI_Type_free(&type);

    count = next_row_ia - my_ia[0];

    *my_ja_p = (int *) malloc(count * sizeof(int));
    if (*my_ja_p == NULL) {
	return MPI_ERR_IO;
    }

    *my_a_p = (double *) malloc(count * sizeof(double));
    if (*my_a_p == NULL) {
	return MPI_ERR_IO;
    }

    /* read local portion of ja */
    my_ja_off = 80 * sizeof(char) + 2 * sizeof(int) + n * sizeof(int) +
	my_ia[0] * sizeof(int);

    err = MPI_File_read_at_all(fh, my_ja_off, *my_ja_p, count, MPI_INT,
			       &status);
    if (err != MPI_SUCCESS) {
	return err;
    }

    /* read local portion of a */
    my_a_off = 80 * sizeof(char) + (2 + n + nz) * sizeof(int) +
	my_ia[0] * sizeof(double);

    err = MPI_File_read_at_all(fh, my_a_off, *my_a_p, count, MPI_DOUBLE,
			       &status);
    if (err != MPI_SUCCESS) {
	return err;
    }

    return MPI_SUCCESS;
}
		    


/* assumption: processes have all elements from a given row */
int CSRIO_Write(char *filename,
                char *title,
                int dimsz, /* total # of rows */
                int mystartrow,
                int myrows, /* local rows */
                int mynonzero,
                double *data,  /* a */
                int *rowstart, /* ia */
                int *col,      /* ja */
                MPI_Info info)
{
    int err;

    int prevnonzero, totalnonzero;

    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;
    int rank, nprocs;

    MPI_File fh;
    MPI_Status status;

    MPI_Datatype type;
    MPI_Offset myfilerowoffset, myfilecoloffset, myfiledataoffset;

    MPI_Comm_size(mlifeio_comm, &nprocs);
    MPI_Comm_rank(mlifeio_comm, &rank);

    
    /* TODO: Use exscan */
    err = MPI_Scan(&mynonzero, &prevnonzero, 1, MPI_INT, MPI_SUM, csrio_comm);
    prevnonzero -= mynonzero; /* MPI_Scan is inclusive */

    err = MPI_Allgather(&mynonzero, 1, MPI_INT,
                        &totalnonzero, 1, MPI_INT, 0, csrio_comm);

    printf("rank %d has %d elements, will start at %d\n", mynonzero,
           prevnonzero);

    err = MPI_File_open(csrio_comm, filename, amode, info, &fh);
    if (err != MPI_SUCCESS) return err;

    if (rank == 0) {
        /* rank 0 writes out the title, # of rows, and count of nonzeros */
        char titlebuf[80];
        int intbuf[2];

        memset(titlebuf, 0, 80);
        strncpy(titlebuf, title, 79);
        err = MPI_File_write_at(fh, 0, titlebuf, 80, MPI_CHAR, &status);

        intbuf[0] = dimsz;
        intbuf[1] = totalnonzero;
        err = MPI_File_write_at(fh, 80, intbuf, 2, MPI_INT, &status);
    }

    /* everyone writes their row offsets, columns, and data into the
     * correct location
     */
    myfilerowoffset = 80 * sizeof(char) + 2 * sizeof(int) + 
        mystartrow * sizeof(int);

    myfilecoloffset = 80 * sizeof(char) + 2 * sizeof(int) +
        dimsz * sizeof(int) +
        prevnonzero * sizeof(int);
        
    myfiledataoffset = 80 * sizeof(char) + 2 * sizeof(int) +
        dimsz * sizeof(int) +
        totalnonzero * sizeof(int) +
        prevnonzero * sizeof(double);

    /* TODO: combine the first two steps? */
    err = MPI_File_write_at_all(fh, myfilerowoffset, rowstart, myrows,
                                MPI_INT, &status);

    err = MPI_File_write_at_all(fh, myfilecoloffset, col, mynonzero,
                                MPI_INT, &status);

    err = MPI_File_write_at_all(fh, mydataoffset, data, mynonzero,
                                MPI_DOUBLE, &status);

    return err; /* TODO: better error handling */
}
                

/* One elaboration is to compute an "optimal" decomposition by trying to give
 * each process almost the same number of non-zeros, which may change the
 * number of consecutive rows of the matrix given to each process.  This can
 * be done using scan on the ia array (either all processes read all of ia and
 * do duplicate computation or read part and do MPI_Exscan), where each
 * process first starts with the number of nonzeros provided by the uniform
 * block decomposition of the ia array.  My looking at the scan values and the
 * total number of non zeros, you can figure out if there are too many non
 * zeros in the processes to your left or right, and then move some in the
 * correct direction.  The details aren't as important as allowing the number
 * of consecutive blocks to be computed or provided by some outside
 * function. -- Bill Gropp
 */
