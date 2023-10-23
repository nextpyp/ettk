/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include "mpi.h"

int main( int argc, char *argv[] )
{
    int rank;
    int size;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Init( 0, 0 );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name,&namelen);

    printf( "Hello world from process %d of %d, running on %s\n", rank, size, processor_name );
    MPI_Finalize();
    return 0;
}
