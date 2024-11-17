/******************************************************************************
  *** Example code of how to use the INHDF5 API
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


// The original demonstration that uses 2 ranks (processes)

int oldtest()
{
   int ierr=0;
   int irank,nrank;
   MPI_Comm comm;

   hid_t id_alist, id_file, id_grp, id_data;
   hid_t id_type, id_space;
   hsize_t idims[1],ioff[1],icnt[1];
   int* idum;
   int n;


   MPI_Comm_dup( MPI_COMM_WORLD, &comm );
   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   if( nrank > 2 ) {
      printf("This demo is meant to run with 2 processes \n");
      return 1;
   }

   //-- open the HDF5 API
   H5open();

   // create a file
   ierr = inUtils_HDF_CreateFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("File creation success\n");
   // close the file to subsequently test the opening of a file
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("File closing success\n");

   // open a file
   ierr = inUtils_HDF_OpenFile( "TEST_FILE", comm, &id_alist, &id_file );
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

   return 0;
}


// A demonstration involving more complex dataset

int test( MPI_Comm* comm_ )
{
   int ierr=0;
   int irank,nrank;
   MPI_Comm comm;

   MPI_Comm_dup( *comm_, &comm );
   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   hid_t id_alist, id_file, id_grp, id_data;
   hid_t id_type, id_space;
   hsize_t idims[2],ioff[2],icnt[2];
   long int* idum;
   long int n,nseg,nc,nmax,noff;
   double t0 = MPI_Wtime(), t1;

   // build a large array on this process (with rectangular data)
   nseg = 1000*1000;       // number of "rows" in every rank's segment
   nc = 10;                // number of "columns" for the data table
   idum = (long int*) malloc( (size_t) (nseg*nc) * sizeof(long int) );
   if( idum == NULL ) {
      printf("Allocation error\n");
      exit(1);   // abnormal exit; the MPI runtime env will terminate them all
   }
   nmax = ((long int) nrank) * nseg;
   noff = ((long int) irank) * nseg;

   //-- open the HDF5 API
   H5open();

   // create a file
   ierr = inUtils_HDF_CreateFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File creation success\n");
   // close the file (and its list) to subsequently test the opening of a file
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("< File closing success\n");

   // open a file
   ierr = inUtils_HDF_OpenFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File opening success\n");

   // create a group
   ierr = inUtils_HDF_CreateGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group creation success\n");
   // close this group the test subsequent opening of the group
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf("<< Group closing success\n");

   // open a group
   ierr = inUtils_HDF_OpenGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group opening success\n");

   // create a dataspace (to use in creating a dataset)
   idims[0] = (hsize_t) nmax;
   idims[1] = (hsize_t) nc;
   id_space = H5Screate_simple( 2, idims, idims );

   // create a dataset
   if( irank == 0 ) printf(">>> Creating dataset\n");
   id_type = H5T_NATIVE_LONG;
   ierr = inUtils_HDF_CreateDataset( "my_dataset", comm, id_grp,
                                     id_type, id_space,
                                     &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset creation success\n");
   // close the dataset so that we can test the re-opening
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("<<< Dataset closing success\n");
   // close the space that is no longer needed
   H5Sclose( id_space );

   // open a dataset
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset opening success\n");

   if( irank == 0 ) printf(" Constructing data\n");
   // put data in the array, of the form "row * 100 + column"
   for(n=0;n<nseg;++n) {
      for(long int k=0;k<nc;++k) {
         idum[ n*nc + k ] = (noff + n+1)*100 + k;
      }
   }

   // write to the dataset
   if( irank == 0 ) printf(">>> Writing to dataset\n");
   ioff[0] = (hsize_t) noff;
   ioff[1] = (hsize_t) 0;
   icnt[0] = (hsize_t) nseg;
   icnt[1] = (hsize_t) nc;
   ierr = inUtils_HDF_DatasetWrite( comm, id_data, id_type, id_space,
                                    2, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Write to dataset successful\n");

   // close the dataset so that we can test the re-opening and re-reading
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("<<< Dataset closing success\n");
   // close the space
   H5Sclose( id_space );

   // change the stored data to something other than what we will read
   if( irank == 0 ) printf(" Corrupting data in memory\n");
   for(n=0;n<nseg;++n) {
      for(long int k=0;k<nc;++k) {
         idum[ n*nc + k ] = ~( (noff + n)*100 + k );  // flip all bits
      }
   }
   if( irank == 0 ) {
      printf(" Rank 0 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
      sleep(1);
      MPI_Barrier( comm );
   } else if( irank == nrank-1 ) {
      MPI_Barrier( comm );
      printf(" Rank 1 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
   } else {
      MPI_Barrier( comm );
   }

   // open the dataset so that we can test reading
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset opening success\n");

   // read the dataset
   ioff[0] = (hsize_t) noff;
   ioff[1] = (hsize_t) 0;
   icnt[0] = (hsize_t) nseg;
   icnt[1] = (hsize_t) nc;
   ierr = inUtils_HDF_DatasetRead( comm, id_data, id_type, id_space,
                                   2, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Reading dataset successful\n");

   // close the dataset
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf("<<< Dataset closing success\n");
   // close the space
   H5Sclose( id_space );

   // show the contents of the buffer as we read it
   if( irank == 0 ) {
      printf(" Rank 0 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
      sleep(1);
      MPI_Barrier( comm );
   } else if( irank == nrank-1 ) {
      MPI_Barrier( comm );
      printf(" Rank 1 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
   } else {
      MPI_Barrier( comm );
   }

   // drop buffer
   free( idum );

   if( irank == 0 ) printf("<< Closing group\n");
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( irank == 0 ) printf("< Closing file\n");
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   t1 = MPI_Wtime();
   if( irank == 0 ) printf(" Final termination after %lf sec\n",t1-t0);

   //-- close the HDF5 and MPI APIs
   H5close();

   return 0;
}


// Some testing of mine...

int test2( MPI_Comm* comm_ )
{
   int ierr=0;
   int irank,nrank;
   MPI_Comm comm;

   MPI_Comm_dup( *comm_, &comm );
   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   hid_t id_alist, id_file, id_grp, id_data;
   hid_t id_type, id_space;
   int nd=2;   // change to "1" to use a 1D layout (dimension 2 will be ignored)
   hsize_t idims[2],ioff[2],icnt[2];
   long int* idum;
   long int n,nseg,nc,nmax,noff;
   double t0,t1;

   // build a large array on this process (with rectangular data)
   nseg = 1000*1000;       // number of "rows" in every rank's segment
   nc = 10;                // number of "columns" for the data table
   idum = (long int*) malloc( (size_t) (nseg*nc) * sizeof(long int) );
   if( idum == NULL ) {
      printf("Allocation error\n");
      exit(1);   // abnormal exit; the MPI runtime env will terminate them all
   }
   nmax = ((long int) nrank) * nseg;
   noff = ((long int) irank) * nseg;


   //-- open the HDF5 API
   H5open();

   t0 = MPI_Wtime();

   // create a file IN PARALLEL
   ierr = inUtils_HDF_CreateFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File creation success\n");
   // close the file to subsequently test the opening of a file
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File closing success\n");

   // open a file IN PARALLEL
   ierr = inUtils_HDF_OpenFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File opening success\n");

   // create a group IN PARALLEL
   ierr = inUtils_HDF_CreateGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group creation success\n");
   // close this group the test subsequent opening of the group
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group closing success\n");

   // open a group IN PARALLEL
   ierr = inUtils_HDF_OpenGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group opening success\n");

   // create a dataspace (to use in creating a dataset)
   idims[0] = (hsize_t) nmax;
   if( nd == 1 ) idims[0] = (hsize_t) (nc*nmax);
   idims[1] = (hsize_t) nc;
   id_space = H5Screate_simple( nd, idims, idims );

   // create a dataset IN PARALLEL
   if( irank == 0 ) printf(">>> Creating dataset...\n");
   id_type = H5T_NATIVE_LONG;
   ierr = inUtils_HDF_CreateDataset( "my_dataset", comm, id_grp,
                                     id_type, id_space,
                                     &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset creation success\n");
   // close the dataset so that we can test the re-opening
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset closing success\n");
   // close the space that is no longer needed
   H5Sclose( id_space );

   // open a dataset IN PARALLEL (not requesting back type and space)
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, NULL, NULL );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset opening success\n");

   // close everything
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset closing success\n");
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group closing success\n");
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File closing success\n");

   t1 = MPI_Wtime();
   if( irank == 0 ) printf("Time elapsed: %lf sec\n",t1-t0);


   // put data in the array, of the form "row * 100 + column"
   for(n=0;n<nseg;++n) {
      for(long int k=0;k<nc;++k) {
         idum[ n*nc + k ] = (noff + n+1)*100 + k;
      }
   }


   if( irank == 0 ) printf("Opening and writing in PARALLEL \n");
   sleep(1);

   t0 = MPI_Wtime();

   ierr = inUtils_HDF_OpenFile( "TEST_FILE", comm, &id_alist, &id_file );
   if( ierr == 0 ) if( irank == 0 ) printf("> File opening success\n");
   ierr = inUtils_HDF_OpenGroup( "my_group", comm, id_file, &id_grp );
   if( ierr == 0 ) if( irank == 0 ) printf(">> Group opening success\n");
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset opening success\n");
   if( irank == 0 ) printf("[Space number of dims: %d]\n",
                            H5Sget_simple_extent_ndims( id_space ) );

/////// cleanup if we do not continue ////////
// ierr = inUtils_HDF_CloseDataset( comm, &id_data );
// H5Sclose( id_space );
// ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
// ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );

   // write to the dataset
   if( irank == 0 ) printf(">>> Writing to dataset...\n");
   ioff[0] = (hsize_t) noff;
   ioff[1] = (hsize_t) 0;
   icnt[0] = (hsize_t) nseg;
   icnt[1] = (hsize_t) nc;
   if( nd == 1 ) {
      ioff[0] = (hsize_t) (nc*noff);
      icnt[0] = (hsize_t) (nc*nseg);
   }
   ierr = inUtils_HDF_DatasetWrite( comm, id_data, id_type, id_space,
                                    nd, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Write to dataset successful\n");

   // close the dataset so that we can test the re-opening and re-reading
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset closing success\n");
   H5Sclose( id_space );

/////// cleanup if we do not continue ////////
// ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
// if( ierr == 0 ) if( irank == 0 ) printf(">> Group closing success\n");
// ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );
// if( ierr == 0 ) if( irank == 0 ) printf("> File closing success\n");

   t1 = MPI_Wtime();
   if( irank == 0 ) printf("Time elapsed: %lf sec\n",t1-t0);


   // change the stored data to something other than what we will read
   for(n=0;n<nseg;++n) {
      for(long int k=0;k<nc;++k) {
         idum[ n*nc + k ] = ~( (noff + n)*100 + k );  // flip all bits
      }
   }
   if( irank == 0 ) {
      printf(" Rank 0 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
      sleep(1);
      MPI_Barrier( comm );
   } else if( irank == nrank-1 ) {
      MPI_Barrier( comm );
      printf(" Rank 1 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
   } else {
      MPI_Barrier( comm );
   }

   t0 = MPI_Wtime();
   if( irank == 0 ) printf("About to read back the data from the open file\n");

   // open the dataset so that we can test reading
   if( irank == 0 ) printf(">>> Opening dataset\n");
   ierr = inUtils_HDF_OpenDataset( "my_dataset", comm, id_grp,
                                   &id_data, &id_type, &id_space );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset opening success\n");

   // read the dataset
   if( irank == 0 ) printf(">>> Reading from dataset...\n");
   ioff[0] = (hsize_t) noff;
   ioff[1] = (hsize_t) 0;
   icnt[0] = (hsize_t) nseg;
   icnt[1] = (hsize_t) nc;
   if( nd == 1 ) {
      ioff[0] = (hsize_t) (nc*noff);
      icnt[0] = (hsize_t) (nc*nseg);
   }
   ierr = inUtils_HDF_DatasetRead( comm, id_data, id_type, id_space,
                                   nd, ioff, icnt, idum );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Reading successful\n");

   // close the dataset
   ierr = inUtils_HDF_CloseDataset( comm, &id_data );
   if( ierr == 0 ) if( irank == 0 ) printf(">>> Dataset closing success\n");
   // close the space
   H5Sclose( id_space );

   // show the contents of the buffer as we read it
   if( irank == 0 ) {
      printf(" Rank 0 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
      sleep(1);
      MPI_Barrier( comm );
   } else if( irank == nrank-1 ) {
      MPI_Barrier( comm );
      printf(" Rank 1 data: ");
      for(n=0;n<10;++n) printf(" %ld",idum[n]);
      printf(" \n");
   } else {
      MPI_Barrier( comm );
   }

   if( irank == 0 ) printf(">> Closing group\n");
   ierr = inUtils_HDF_CloseGroup( comm, &id_grp );
   if( irank == 0 ) printf("> Closing file\n");
   ierr = inUtils_HDF_CloseFile( comm, &id_alist, &id_file );

   t1 = MPI_Wtime();
   if( irank == 0 ) printf("Final termination after %lf sec\n",t1-t0);

   // drop buffer
   free( idum );

   //-- close the HDF5 and MPI APIs
   H5close();

   return 0;
}

