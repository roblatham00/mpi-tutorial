/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef MLIFE_H
#define MLIFE_H

#if 0
extern void   srand48();
extern double drand48();
extern char * malloc();
#endif

extern char *optarg;

int MLIFE_myrows(int dimsz, int rank, int nprocs);
int MLIFE_myrowoffset(int dimsz, int rank, int nprocs);

double life(int matrix_size, int ntimes, MPI_Comm comm);

#define BORN 1
#define DIES 0

#endif
