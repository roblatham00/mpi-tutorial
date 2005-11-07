C
C
C  (C) 2004 by University of Chicago.
C      See COPYRIGHT in top-level directory.
C

C stdout implementation of checkpoint (no restart) for MPI Life
C
C Data output in matrix order: spaces represent dead cells,
C '*'s represent live ones.
C


       subroutine MLIFEIO_Init(comm)
       integer comm
       integer ierr
       integer mlifeio_comm
       common /mlifeio/ mlifeio_comm
       save /mlifeio/
       call mpi_comm_dup( comm, mlifeio_comm, ierr )
       return
       end
C
       subroutine MLIFEIO_Finalize()
       integer ierr
       
       if (mlifeio_comm .ne. MPI_COMM_NULL) then
          call mpi_comm_free( mlifeio_comm, ierr )
       endif
       end

       subroutine MLIFEIO_Checkpoint(prefix, matrix, lRows, lCols,        &
     &                               GRows, GCols, iter, info )
       character*(*) prefix
       integer lRows, lCols, GRows, GCols, iter, info
       integer matrix(lRows,lCols)

       integer err, rank, nprocs, r, i
       integer GFirstRow, GFirstCol

       err = 0

       call mpi_comm_size( mlifeio_comm, nprocs, ierr )
       call mpi_comm_rank( mlifeio_comm, rank, ierr )

       call MLIFE_MeshDecomp(rank, nprocs, GRows, GCols,                  &
     &                NULL, NULL, NULL, NULL,                             &
     &                LRows, LCols, GFirstRow, GFirstCol)

C   each proc. writes its part of the display, in rank order
       if (rank .eq. 0) then 
          print *, "[H[2J# Iteration ", iter
       endif
       

C     Slow but simple ...
      do r=0, nprocs-1
         if (rank .eq. r) then
C     print rank 0 data first 
            do i=1, LRows
               if (GFirstCol .eq. 1) then
                  print *, "[%03d;%03dH", 1+(i-1+GFirstRow+1),          &
     &                 GFirstCol
               else 
                  print *, "[%03d;%03dH", 1+(i-1+GFirstRow+1),          &
     &                 GFirstCol+5
               endif
               
               call MLIFEIO_Row_print(matrix(i,1), LCols,                 &
     &              i+GFirstRow-1, GFirstCol .eq. 1)
            enddo
         endif
         call mpi_barrier(mlifeio_comm, ierr )
      enddo
      
      if (rank .eq. nprocs-1) then
         print *, "[%03d;%03dH", GRows+3, 1
      endif

C give time to see the results 
      call MLIFEIO_msleep(250)
      
      end

      subroutine MLIFEIO_Row_print(data, cols, rownr, labelrow )
      integer data(*), cols, rownr
      logical labelrow
      integer i
      character*160 line

      if (labelrow) then
         print *, "%3d: ", rownr
      endif
      do i=1, cols
         if (data(i) .eq. 1) then
            line(i:i) = "*";
         endif
      enddo
      print *, line(1:cols)
      return
      end

      integer function MLIFEIO_Restart(prefix, matrix, GRows, GCols,     &
     &     iter, info)
      character*(*) prefix
      integer matrix(*), GRows, GCols, iter, info
      MLIFEIO_Restart = MPI_ERR_IO
      return
      end
C
      integer function MLIFEIO_Can_restart()
      MLIFEIO_Can_restart = 0
      return 
      end
C
      subroutine MLIFEIO_msleep( msec )
      integer msec
      double precision t0
      include 'mpif.h'
C
      t0 = mpi_wtime()
      do while (mpi_wtime() - t0 < 0.001 * msec ) 
      enddo
      return
      end
