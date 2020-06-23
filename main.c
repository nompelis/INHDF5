/******************************************************************************
  *** The driver code of an example of how to use the INHDF5 API
  *** 
  *** Copyright (c) 2020 Ioannis Nompelis
  *** Ioannis Nompelis <nompelis@nobelware.com>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <mpi.h>
#include <hdf5.h>

#include "inUtils_HDF5.h"


int main( int argc, char *argv[] )
{
   int ierr=0;
   int irank,nrank;
   MPI_Comm comm;

   hid_t id_alist, id_file, id_grp, id_data;
   hid_t id_type, id_space;
   hsize_t idims[1],ioff[1],icnt[1];
   int* idum;
   int n;


   //-- start the day with MPI
   MPI_Init( &argc, &argv );
   MPI_Comm_dup( MPI_COMM_WORLD, &comm );
   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   if( nrank > 2 ) {
      printf("This demo is meant to run with 2 processes \n");
      MPI_Finalize();
      return 1;
   }

   //-- open the HDF5 API
   H5open();

   // create a file
   ierr = inUtils_HDF_CreateFile( "test_file", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("File creation success\n");
   // close the file to subsequently test the opening of a file
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("File closing success\n");

   // open a file
   ierr = inUtils_HDF_OpenFile( "test_file", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("File opening success\n");

   // create a group
   ierr = inUtils_HDF_CreateGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf("Group creation success\n");
   // close this group the test subsequent opening of the group
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf("Group closing success\n");

   // open a group
   ierr = inUtils_HDF_OpenGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf("Group opening success\n");

   // create a dataspace (to use in creating a dataset)
   idims[0] = 20;
   id_space = H5Screate_simple( 1, idims, idims );

   // create a dataset
   id_type = H5T_NATIVE_INT;
   ierr = inUtils_HDF_CreateDataset( "my_dataset", comm, id_grp,
                                     id_type, id_space,
                                     &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset creation success\n");
   // close the dataset so that we can test the re-opening
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset closing success\n");

   // close the space that is no longer needed
   H5Sclose( id_space );

   // open a dataset
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset opening success\n");

   // create an array to write to the dataset
   idum = (int*) malloc( 10 * sizeof(int) );    // no error-trapping
   for(n=0;n<10;++n) idum[n] = irank + 1;

   // write to the dataset
   if( irank == 0 ) ioff[0] = 0;
   if( irank == 1 ) ioff[0] = 10;
   icnt[0] = 10;
   ierr = inUtils_HDF_DatasetWrite( comm, id_data, id_type, id_space,
                                    1, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf("Write to dataset successful\n");

   // close the dataset so that we can test the re-opening and re-reading
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset closing success\n");

   // close the space
   H5Sclose( id_space );

   // change the stored data to something other than what we will read
   for(n=0;n<10;++n) idum[n] = irank + 1;

   // open the dataset so that we can test reading
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset opening success\n");

   // read the dataset
   if( irank == 0 ) ioff[0] = 0;
   if( irank == 1 ) ioff[0] = 10;
   icnt[0] = 10;
   ierr = inUtils_HDF_DatasetRead( comm, id_data, id_type, id_space,
                                   1, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf("Reading from dataset successful\n");

   // close the dataset for doog
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("Dataset closing success\n");

   // close the space
   H5Sclose( id_space );

   // show the contents of the buffer as we read it
   if( irank == 0 ) {
      printf(" Rank 0 data: ");
      for(n=0;n<10;++n) printf(" %d",idum[n]);
      printf(" \n");
      sleep(1);
      MPI_Barrier( comm );
   } else {
      MPI_Barrier( comm );
      printf(" Rank 1 data: ");
      for(n=0;n<10;++n) printf(" %d",idum[n]);
      printf(" \n");
      sleep(1);
   }

   // drop buffer
   free( idum );

   printf("Final termination \n");
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );

   //-- close the HDF5 and MPI APIs
   H5close();
   MPI_Finalize();

   return 0;
}

