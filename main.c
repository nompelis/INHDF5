/******************************************************************************
  *** The driver code for the examples of how to use the INHDF5 API
  *** 
  *** Copyright (c) 2020-2024 Ioannis Nompelis
  *** Ioannis Nompelis <nompelis@nobelware.com>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <mpi.h>
#include <hdf5.h>

#include "inUtils_HDF5.h"

// a demo that uses 2 ranks (processes) found in "test.c"
int oldtest();

// a deme that handles a 2D dataset
int test( MPI_Comm* );

// some testing...
int test2( MPI_Comm* );

int main( int argc, char *argv[] )
{
   int ierr=0;
   int irank,nrank;
   MPI_Comm comm;


   //-- start the day with MPI
   MPI_Init( &argc, &argv );
   MPI_Comm_dup( MPI_COMM_WORLD, &comm );
   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );


// (void) oldtest();
// (void) test( &comm );
   (void) test2( &comm );


   MPI_Finalize();

   return 0;
}

