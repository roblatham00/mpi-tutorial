/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

/* 
 * Our storage format for CSR will use native byte format.  The file will
 * contain:
 *
 * title (char, 80 bytes, fixed size)
 * n (int)
 * nz (int)
 * ia[i], i=1,...,n+1 (int array)
 * ja[i], i=1,...,nz  (int array)
 * a[i],  i=1,...,nz  (double array)
 *
 * See README.txt for more information on CSR.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

static MPI_Comm csrio_comm = MPI_COMM_NULL;
static MPI_Info csrio_info = MPI_INFO_NULL;

/* CSRIO_Init
 *
 * Parameters:
 * comm - communicator describing group of processes that will perform I/O
 * info - set of hints passed to CSRIO calls
 */
int CSRIO_Init(MPI_Comm comm, MPI_Info info)
{
    int err;

    err = MPI_Comm_dup(comm, &csrio_comm);
    if (err == MPI_SUCCESS) {
	err = MPI_Info_dup(info, &csrio_info);
    }

    return err;
}

/* CSRIO_Finalize
 */
int CSRIO_Finalize(void)
{
    MPI_Comm_free(&csrio_comm);
    MPI_Info_free(&csrio_info);

    return MPI_SUCCESS;
}


/* CSRIO_Read_header
 *
 * Parameters:
 * filename - name of file from which header will be read
 * title    - pointer to buffer of at least 80 characters which will hold
 *            title from file if call completes successfully
 * n_p      - address of integer in which number of rows is stored on success
 * nz_p     - address of integer in which number of nonzeros is stored
 *            on success
 *
 * Returns MPI_SUCCESS on success, MPI error code on error.
 */
int CSRIO_Read_header(char *filename, char *title, int  *n_p, int  *nz_p)
{
    int err = 0, ioerr;
    int nrnz[2];

    int rank;
    int amode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

    MPI_Datatype type;
    int lens[3];
    MPI_Aint disps[3];
    MPI_Datatype types[3];

    MPI_Comm_rank(csrio_comm, &rank);

    /* often it is faster for only one process to open/access/close
     * the file when only a small amount of data is going to be
     * accessed.
     */
    if (rank == 0) {
        MPI_File fh;
        MPI_Status status;

        err = MPI_File_open(MPI_COMM_SELF, filename, amode, csrio_info, &fh);
        if (err == MPI_SUCCESS) {
            err = MPI_File_read_at(fh, 0, title, 80, MPI_CHAR, &status);
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
    MPI_Address(title, &disps[0]);
    MPI_Address(nrnz, &disps[1]);
    MPI_Address(&ioerr, &disps[2]);
    types[0] = MPI_CHAR;
    types[1] = MPI_INT;
    types[2] = MPI_INT;

    MPI_Type_struct(3, lens, disps, types, &type);

    /* broadcast the header data to everyone */
    err = MPI_Bcast(MPI_BOTTOM, 1, type, 0, csrio_comm);
    if (err == MPI_SUCCESS && ioerr == MPI_SUCCESS) {
        *n_p = nrnz[0];
        *nz_p = nrnz[1];

        return MPI_SUCCESS;
    }
    else {
        return (err != MPI_SUCCESS) ? err : ioerr;
    }
}


/* CSRIO_Read_rows
 *
 * Parameters:
 * n         - number of rows in matrix
 * nz        - number of nonzero values in matrix
 * my_nz     - maximum number of nonzero values to read into local buffer
 *             (ignored if equal to 0); holds actual number of nonzero
 *             values on success
 * row_start - first row to read (0-origin)
 * row_end   - last row to read
 * my_ia     - pointer to local memory for storing row start indices
 * my_ja_p   - address of pointer to local memory for storing column indices
 *             (if not NULL, then region must be large enough for my_nz values)
 * my_a_p    - address of pointer to local memory for storing data values
 *             (if not NULL, then region must be large enough for my_nz values)
 *
 * Notes:
 * If (*my_ja_p == NULL) then memory will be allocated.  Likewise, if
 * (*my_a_p == NULL) then memory will be allocated.
 *
 * Returns MPI_SUCCESS on success, MPI error code on error.
 */
int CSRIO_Read_rows(char *filename, int n, int nz, int *my_nz_p, int row_start,
		    int row_end, int *my_ia, int **my_ja_p, double **my_a_p)
{
    int i, count, err;

    MPI_Aint my_ia_off, my_ja_off, my_a_off;

    int next_row_ia, my_nz_ok, my_ja_ok = 1, my_a_ok = 1, my_mem_ok,
	all_mem_ok;

    int lens[2];
    MPI_Aint disps[2];

    int amode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;
    MPI_File fh;
    MPI_Status status;
    MPI_Datatype type;

    err = MPI_File_open(csrio_comm, filename, amode, csrio_info, &fh);
    if (err != MPI_SUCCESS) {
	return err;
    }

    my_ia_off = 80 * sizeof(char) + 2 * sizeof(int) + row_start * sizeof(int);

    /* must read one more row start to calculate number of elements */
    lens[0] = row_end - row_start + 1;
    lens[1] = 1;
    MPI_Address(my_ia, &disps[0]);
    MPI_Address(&next_row_ia, &disps[1]);

    MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
    MPI_Type_commit(&type);

    err = MPI_File_read_at_all(fh, my_ia_off, MPI_BOTTOM, 1, type, &status);
    if (err != MPI_SUCCESS) {
	return err;
    }
    MPI_Type_free(&type);

    count = next_row_ia - my_ia[0];

    /* verify local nz value, allocate memory as necessary */
    my_nz_ok = (*my_nz_p == 0 || *my_nz_p >= count) ? 1 : 0;

    if (*my_ja_p == NULL) {
	*my_ja_p = (int *) malloc(count * sizeof(int));
	if (*my_ja_p == NULL) my_ja_ok = 0;
    }

    if (*my_a_p == NULL) {
	*my_a_p = (double *) malloc(count * sizeof(double));
	if (*my_a_p == NULL) my_a_ok = 0;
    }

    /* verify everyone has adequate memory regions and abort now if not */
    my_mem_ok = (my_nz_ok && my_a_ok && my_ja_ok) ? 1 : 0;

    MPI_Allreduce(&my_mem_ok, &all_mem_ok, 1, MPI_INT, MPI_MIN, csrio_comm);
    if (!all_mem_ok) {
	return MPI_ERR_IO;
    }

    /* save actual number of local nonzeros */
    *my_nz_p = count;

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

    /* convert ia values to local references */
    for (i=1; i < row_end - row_start + 1; i++) {
	my_ia[i] -= my_ia[0];
    }
    my_ia[0] = 0;

    return MPI_SUCCESS;
}


/* CSRIO_Write
 *
 * Parameters:
 * filename  - name of file to hold data
 * n         - number of rows in matrix
 * my_nz     - number of nonzero values to read into local buffer
 * row_start - first row to write
 * row_end   - last row to write
 * my_ia     - local row start indices
 * my_ja     - column indices for local array values
 * my_a      - data values
 * 
 * Returns MPI_SUCCESS on success, MPI error code on error.
 */
int CSRIO_Write(char *filename, char *title, int n, int my_nz, int row_start,
                int row_end, const int my_ia[], const int my_ja[],
		const double my_a[])
{
    int i, err;
    int *tmp_ia;

    int prev_nz, tot_nz;

    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE |
	MPI_MODE_UNIQUE_OPEN;
    int rank, nprocs;

    MPI_File fh;
    MPI_Status status;

    MPI_Offset myfilerowoffset, myfilecoloffset, myfiledataoffset;

    MPI_Comm_size(csrio_comm, &nprocs);
    MPI_Comm_rank(csrio_comm, &rank);
    
    /* TODO: Use exscan */
    err = MPI_Scan(&my_nz, &prev_nz, 1, MPI_INT, MPI_SUM, csrio_comm);
    prev_nz -= my_nz; /* MPI_Scan is inclusive */

    err = MPI_Allgather(&my_nz, 1, MPI_INT, &tot_nz, 1, MPI_INT, csrio_comm);

    printf("rank %d has %d elements, will start at %d\n", rank, my_nz,
           prev_nz);

    err = MPI_File_open(csrio_comm, filename, amode, csrio_info, &fh);
    if (err != MPI_SUCCESS) return err;

    /* rank 0 writes out the title, # of rows, and count of nonzeros */
    if (rank == 0) {
        char titlebuf[80];
        int intbuf[2];

        memset(titlebuf, 0, 80);
        strncpy(titlebuf, title, 79);
        err = MPI_File_write_at(fh, 0, titlebuf, 80, MPI_CHAR, &status);

        intbuf[0] = n;
        intbuf[1] = tot_nz;
        err = MPI_File_write_at(fh, 80, intbuf, 2, MPI_INT, &status);
    }

    /* everyone writes their row offsets, columns, and data into the
     * correct location
     */
    myfilerowoffset = 80 * sizeof(char) + 2 * sizeof(int) + 
        row_start * sizeof(int);

    myfilecoloffset = 80 * sizeof(char) + 2 * sizeof(int) +
        n * sizeof(int) + prev_nz * sizeof(int);
        
    myfiledataoffset = 80 * sizeof(char) + 2 * sizeof(int) +
        n * sizeof(int) + tot_nz * sizeof(int) + prev_nz * sizeof(double);

    /* TODO: combine the first two steps? */

    /* copy ia; adjust to be relative to global data */
    tmp_ia = (int *) malloc((row_end - row_start + 1) * sizeof(int));
    if (tmp_ia == NULL) return MPI_ERR_IO; /* TODO: BETTER ERRS */

    for (i=0; i < row_end - row_start + 1; i++) {
	tmp_ia[i] = my_ia[i] + prev_nz;
    }
    err = MPI_File_write_at_all(fh, myfilerowoffset, tmp_ia,
				row_end - row_start + 1,
                                MPI_INT, &status);
    free(tmp_ia);

    err = MPI_File_write_at_all(fh, myfilecoloffset, (void *) my_ja, my_nz,
                                MPI_INT, &status);

    err = MPI_File_write_at_all(fh, myfiledataoffset, (void *) my_a, my_nz,
                                MPI_DOUBLE, &status);

    return err; /* TODO: better error handling */
}
