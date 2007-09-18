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
    int infd, outfd;
    char buf[BUFSIZE];
    MPI_Info hostinfo;
    MPI_Comm pcpslaves;

    MPI_Init( &argc, &argv );

    makehostlist( argv[1], "targets" );
    MPI_Info_create( &hostinfo );
    MPI_Info_set( hostinfo, "file", "targets" );
    MPI_Info_set( hostinfo, "soft", "0-1000" );
    MPI_Comm_spawn( "pcp_slave", MPI_ARGV_NULL, 1000, hostinfo,
		    0, MPI_COMM_SELF, &pcpslaves,
		    MPI_ERRCODES_IGNORE );
    MPI_Info_free( &hostinfo );

    strcpy( outfilename, argv[3] );
    if ( (infd = open( argv[2], O_RDONLY ) ) == -1 ) {
        fprintf( stderr, "input %s does not exist\n", argv[2] );
	sprintf( controlmsg, "exit" );
	MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, pcpslaves );
	MPI_Finalize();
	return( -1 );
    }
    else {
        sprintf( controlmsg, "ready" );
        MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, pcpslaves );
    }

    sprintf( controlmsg, outfilename );
    MPI_Bcast( controlmsg, CMDSIZE, MPI_CHAR, 0, pcpslaves );
    if ( (outfd = open( outfilename, O_CREAT|O_WRONLY|O_TRUNC,
			S_IRWXU ) ) == -1 )
        mystatus = -1;
    else
        mystatus = 0;
    MPI_Allreduce( &mystatus, &allstatus, 1, MPI_INT, MPI_MIN,
		   pcpslaves );
    if ( allstatus == -1 ) {
        fprintf( stderr, "output file %s could not be opened\n",
		 outfilename );
	MPI_Finalize();
	return( -1 );
    }

    /* at this point all files have been successfully opened */
    
    done = 0;
    while ( !done ) {
        numread = read( infd, buf, BUFSIZE );
	MPI_Bcast( &numread, 1, MPI_INT, 0, pcpslaves );
	if ( numread > 0 ) {
	    MPI_Bcast( buf, numread, MPI_BYTE, 0, pcpslaves );
	    write( outfd, buf, numread );
	}
	else {	  
	    close( outfd );
	    done = 1;
	}
    }
    MPI_Comm_free( &pcpslaves );
    MPI_Finalize();
}

int makehostlist( char spec[80], char filename[80] )
{

}

