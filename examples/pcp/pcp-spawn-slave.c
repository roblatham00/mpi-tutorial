/* pcp from SUT, in MPI */
#include "mpi.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#define BUFSIZE    256*1024
#define CMDSIZE    80
int main( int argc, char *argv[] )
{
    int mystatus, allstatus, done, numread;
    char outfilename[128], controlmsg[80];
    int outfd;
    char buf[BUFSIZE];
    MPI_Comm slavecomm;

    MPI_Init( &argc, &argv );

    MPI_Comm_get_parent( &slavecomm );
    MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, slavecomm );
    if ( strcmp( controlmsg, "exit" ) == 0 ) {
        MPI_Finalize();
	return( -1 );
    }

    MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, slavecomm );
    if ( (outfd = open( controlmsg, O_CREAT|O_WRONLY|O_TRUNC,
			S_IRWXU ) ) == -1 ) 
        mystatus = -1;
    else
        mystatus = 0;
    MPI_Allreduce( &mystatus, &allstatus, 1, MPI_INT, MPI_MIN,
		   slavecomm );
    if ( allstatus == -1 ) {
	MPI_Finalize();
	return( -1 );
    }

    /* at this point all files have been successfully opened */
    
    done = 0;
    while ( !done ) {
	MPI_Bcast( &numread, 1, MPI_INT, 0, slavecomm );
	if ( numread > 0 ) {
	    MPI_Bcast( buf, numread, MPI_BYTE, 0, slavecomm );
	    write( outfd, buf, numread ); 
	}
	else {	  
	    close( outfd );
	    done = 1;
	}
    }
    MPI_Comm_free( &slavecomm );
    MPI_Finalize();
    return 0;
}
