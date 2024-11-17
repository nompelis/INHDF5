/******************************************************************************
  *** An example of how to use the INHDF5 API
  *** 
  *** Copyright (c) 2020-2024 Ioannis Nompelis
  *** Ioannis Nompelis <nompelis@nobelware.com>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <hdf5.h>

#include "inUtils_HDF5.h"


int inUtils_HDF_CreateFile( const char* filename, MPI_Comm comm,
                            hid_t* id_alist, 
                            hid_t* id_file )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_alist_tmp, id_clist_tmp, id_file_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // create creation list and file access property list with parallel I/O
   id_clist_tmp = H5Pcreate( H5P_FILE_CREATE );
   id_alist_tmp = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio( id_alist_tmp, comm, MPI_INFO_NULL );

   id_file_tmp = H5Fcreate( filename, H5F_ACC_TRUNC,
                            id_clist_tmp, id_alist_tmp );
   if( id_file_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) printf("Could not create file: \"%s\"\n",filename);
      if( id_file_tmp >= 0 ) H5Fclose( id_file_tmp );
      H5Pclose( id_alist_tmp );
      H5Pclose( id_clist_tmp );
      ierr = 1;
      goto return_point;
   }

   // copy structures to returned pointers
   memcpy( id_alist, &id_alist_tmp, sizeof(hid_t) );
   memcpy( id_file, &id_file_tmp, sizeof(hid_t) );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to open a file and return handle and access list
// (closes the access property list as well)
//

int inUtils_HDF_OpenFile( const char* filename, MPI_Comm comm,
                          hid_t* id_alist, 
                          hid_t* id_file )
{
   return inUtils_HDF_OpenFileOpt( filename, comm, id_alist, id_file, 1 );
}

int inUtils_HDF_OpenFileOpt( const char* filename, MPI_Comm comm,
                             hid_t* id_alist, 
                             hid_t* id_file,
                             int iop )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_alist_tmp, id_file_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );


   // create file access property list with parallel I/O
   id_alist_tmp = H5Pcreate( H5P_FILE_ACCESS );
   H5Pset_fapl_mpio( id_alist_tmp, comm, MPI_INFO_NULL );

   // open the file
   if( iop == 0 ) {
      id_file_tmp = H5Fopen( filename, H5F_ACC_RDONLY, id_alist_tmp );
   } else {
      id_file_tmp = H5Fopen( filename, H5F_ACC_RDWR, id_alist_tmp );
   }
   if( id_file_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error opening file: \"%s\" \n",filename);
      if( id_file_tmp >= 0 ) H5Fclose( id_file_tmp );
      H5Pclose( id_alist_tmp );
      ierr = 1;
      goto return_point;
   }

   // copy the handles to return
   memcpy( id_alist, &id_alist_tmp, sizeof(hid_t) );
   memcpy( id_file, &id_file_tmp, sizeof(hid_t) );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to close a file that is open in parallel
// (closes the access property list as well)
//

int inUtils_HDF_CloseFile( MPI_Comm comm,
                           hid_t* id_alist, 
                           hid_t* id_file )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // close the file
   ierr = H5Fclose( *id_file );
   if( ierr != 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error closeing file. (Access list stays open)\n");
      ierr = 1;
      goto return_point;
   }

   // close the access list
   H5Pclose( *id_alist );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//--- Groups ---

//
// function to create a group within a handle
//

int inUtils_HDF_CreateGroup( const char* grpname, MPI_Comm comm,
                             hid_t id_tree,
                             hid_t* id_grp )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_grp_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // create a group in parallel
   id_grp_tmp = H5Gcreate( id_tree, grpname,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   if( id_grp_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error creating group: \"%s\" \n",grpname);
      if( id_grp_tmp >= 0 ) H5Gclose( id_grp_tmp );
      ierr = 1;
      goto return_point;
   }

   // copy the handle to return
   memcpy( id_grp, &id_grp_tmp, sizeof(hid_t) );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to close an open group that has been opened in parallel
//

int inUtils_HDF_CloseGroup( MPI_Comm comm, hid_t* id_grp )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // close the group in parallel
   ierr = H5Gclose( *id_grp );
   if( ierr != 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error closing group \n");
      ierr = 1;
      goto return_point;
   }

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to open a group in a handle and return handle
//

int inUtils_HDF_OpenGroup( const char* grpname, MPI_Comm comm,
                           hid_t id_tree,
                           hid_t* id_grp )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_grp_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // open the group
   id_grp_tmp = H5Gopen( id_tree, grpname );
   if( id_grp_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error opening group: \"%s\" \n",grpname);
      if( id_grp_tmp >= 0 ) H5Gclose( id_grp_tmp );
      ierr = 1;
      goto return_point;
   }

   // copy the handle to return
   memcpy( id_grp, &id_grp_tmp, sizeof(hid_t) );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//--- Datasets ---


//
// function to create a dataset of the incoming data type
//

int inUtils_HDF_CreateDataset( const char* dataname, MPI_Comm comm,
                               hid_t id_tree,
                               hid_t id_type,
                               hid_t id_space,
                               hid_t* id_data )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_data_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // create teh dataset
   id_data_tmp = H5Dcreate( id_tree, dataname, id_type, id_space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   if( id_data_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error creating dataset: \"%s\" \n",dataname);
      if( id_data_tmp >= 0 ) H5Dclose( id_data_tmp );
      ierr = 1;
      goto return_point;
   }

   // copy the handle to return
   memcpy( id_data, &id_data_tmp, sizeof(hid_t) );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to close an open dataset that has been opened in parallel
//

int inUtils_HDF_CloseDataset( MPI_Comm comm, hid_t* id_data )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // close the dataset in parallel
   ierr = H5Dclose( *id_data );
   if( ierr != 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error closing dataset \n");
      ierr = 1;
      goto return_point;
   }

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// function to open a dataset in parallel
// (returns space and type of the data)
//

int inUtils_HDF_OpenDataset( const char* dataname, MPI_Comm comm, hid_t id_tree,
                             hid_t* id_data,
                             hid_t* id_type,
                             hid_t* id_space )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_data_tmp, id_type_tmp, id_space_tmp;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );

   // open the dataset in parallel
   id_data_tmp = H5Dopen( id_tree, dataname );
   if( id_data_tmp < 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      fprintf( stdout, " e  Error opening dataset: \"%s\" \n",dataname);
      if( id_data_tmp >= 0 ) H5Dclose( id_data_tmp );
      ierr = 1;
      goto return_point;
   }

   /// retrieve the space of the dataset
   id_space_tmp = H5Dget_space( id_data_tmp );

   // retrieve the datatype in the dataset
   id_type_tmp = H5Dget_type( id_data_tmp );

   // Potentially provide info about the dataset using these handles...

   // copy the handles to return
   memcpy( id_data, &id_data_tmp, sizeof(hid_t) );
   if( id_type != NULL ) {
      memcpy( id_type, &id_type_tmp, sizeof(hid_t) );
   } else {
      H5Tclose( id_type_tmp );
   }
   if( id_space != NULL ) {
      memcpy( id_space, &id_space_tmp, sizeof(hid_t) );
   } else {
      H5Sclose( id_space_tmp );
   }

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}


//
// functions to write to a dataset that has been opened in parallel
//

int inUtils_HDF_DatasetWriteCall( int iop,
                              MPI_Comm comm, hid_t id_data,
                              hid_t id_type, hid_t id_space,
                              int nd,
                              hsize_t* ioff,
                              hsize_t* icnt,
                              void* data )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_type_tmp, id_space_tmp, id_mem, id_plist;
   htri_t ht;
   int nds,ndd, n;
   hsize_t *idims,*idim1,*idimd,*idim2,*istd,*iblk;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );


   // retrieve the datatype in the dataset
   id_type_tmp = H5Dget_type( id_data );

   // check for concistency of data types
   ht = H5Tequal( id_type, id_type_tmp );

   // immediatelly release the datatype
   H5Tclose( id_type_tmp );
   if( ht <= 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Datatype is inconcistent \n");
      ierr = 1;
      goto return_point;
   }

   // retrieve the space of the dataset
   id_space_tmp = H5Dget_space( id_data );
   // check for concistency of data space rank
   nds = H5Sget_simple_extent_ndims( id_space );
   ndd = H5Sget_simple_extent_ndims( id_space_tmp );
   if( nds != ndd || nds != nd || ndd != nd ) {
      if( irank == 0 ) fprintf( stdout, " e  Sizes are inconcistent \n");
      H5Sclose( id_space_tmp );
      ierr = 2;
      goto return_point;
   }

   idims = (hsize_t *) malloc(6*((size_t) nd) * sizeof(hsize_t) );
   if( idims == NULL ) ierr = -1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Internal memory error\n");
      if( idims != NULL ) free( idims );
      H5Sclose( id_space_tmp );
      ierr = -1;
      goto return_point;
   }
   idim1 = &( idims[ nd*1 ] );
   idimd = &( idims[ nd*2 ] );
   idim2 = &( idims[ nd*3 ] );
   istd = &( idims[ nd*4 ] );
   iblk = &( idims[ nd*5 ] );

   // check for concistency of data space dimensions
   H5Sget_simple_extent_dims( id_space, idims, idim1 );
   H5Sget_simple_extent_dims( id_space_tmp, idimd, idim2 );

   // immediatelly drop the space handle
   H5Sclose( id_space_tmp );
   ierr = 0;
   for(n=0;n<nds;++n) {
      if( idims[n] != idimd[n] ) ierr = ierr + 1;
   }
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Dataspace dimension different \n");
      free( idims );
      ierr = 3;
      goto return_point;
   }

   // create a transfer property list with or without parallel I/O
   if( iop == 0 ) {
      id_plist = H5Pcreate( H5P_DATASET_XFER );
      H5Pset_dxpl_mpio( id_plist, H5FD_MPIO_COLLECTIVE );
   } else if( iop == 1 ) {
      id_plist = H5Pcreate( H5P_DATASET_XFER );
      H5Pset_dxpl_mpio( id_plist, H5FD_MPIO_INDEPENDENT );
   } else {
      id_plist = H5P_DEFAULT;
   }

   // create a memory dataspace
   id_mem = H5Screate_simple( nd, icnt, icnt );

   // select the area where the data will be written
   for(n=0;n<nds;++n) {
      istd[n] = 1;
      iblk[n] = 1;
   }
   H5Sselect_hyperslab( id_space, H5S_SELECT_SET, ioff, istd, icnt, iblk );

   // perform the parallel write operations
   ierr = H5Dwrite( id_data, id_type, id_mem, id_space, id_plist, data );
   if( ierr != 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  There was a problem writing \n");
      H5Pclose( id_plist );
      ierr = 100;
      goto return_point;
   }

   // drop handles no longer needed
   H5Sclose( id_mem );
   if( iop == 0 || iop == 1 ) H5Pclose( id_plist );

   //drop the sizes array
   free( idims );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}

// the high level call
int inUtils_HDF_DatasetWrite( MPI_Comm comm, hid_t id_data,
                              hid_t id_type, hid_t id_space,
                              int nd,
                              hsize_t* ioff,
                              hsize_t* icnt,
                              void* data )
{
   return inUtils_HDF_DatasetWriteCall( 0, comm, id_data, id_type, id_space,
                                        nd, ioff, icnt, data );
}

//
// functions to read from a dataset that has been opened in parallel
//

int inUtils_HDF_DatasetReadCall( int iop,
                             MPI_Comm comm, hid_t id_data,
                             hid_t id_type, hid_t id_space,
                             int nd,
                             hsize_t* ioff,
                             hsize_t* icnt,
                             void* data )
{
   int ierr=0;
   int irank,nrank;
   void *old_err_client_data;   // to store the old error handler
   H5E_auto_t old_err_func;     // to store the old error handler
   hid_t id_estack;             // new error stack
   hid_t id_type_tmp, id_space_tmp, id_mem, id_plist;
   htri_t ht;
   int nds,ndd, n;
   hsize_t *idims,*idim1,*idimd,*idim2,*istd,*iblk;


   MPI_Comm_rank( comm, &irank );
   MPI_Comm_size( comm, &nrank );

   // turn off HDF error-printing; store old handler to restore later
   id_estack = H5Ecreate_stack();
   H5Eget_auto( id_estack, &old_err_func, &old_err_client_data );
   H5Eset_auto( id_estack, NULL, NULL );


   // retrieve the datatype in the dataset
   id_type_tmp = H5Dget_type( id_data );

   // check for concistency of data types
   ht = H5Tequal( id_type, id_type_tmp );

   // immediatelly release the datatype
   H5Tclose( id_type_tmp );
   if( ht <= 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Datatype is inconcistent \n");
      ierr = 1;
      goto return_point;
   }

   // retrieve the space of the dataset
   id_space_tmp = H5Dget_space( id_data );
   // check for concistency of data space rank
   nds = H5Sget_simple_extent_ndims( id_space );
   ndd = H5Sget_simple_extent_ndims( id_space_tmp );
   if( nds != ndd || nds != nd || ndd != nd ) {
      if( irank == 0 ) fprintf( stdout, " e  Sizes are inconcistent \n");
      H5Sclose( id_space_tmp );
      ierr = 2;
      goto return_point;
   }

   idims = (hsize_t *) malloc(6*((size_t) nd) * sizeof(hsize_t) );
   if( idims == NULL ) ierr = -1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Internal memory error\n");
      if( idims != NULL ) free( idims );
      H5Sclose( id_space_tmp );
      ierr = -1;
      goto return_point;
   }
   idim1 = &( idims[ nd*1 ] );
   idimd = &( idims[ nd*2 ] );
   idim2 = &( idims[ nd*3 ] );
   istd = &( idims[ nd*4 ] );
   iblk = &( idims[ nd*5 ] );

   // check for concistency of data space dimensions
   H5Sget_simple_extent_dims( id_space, idims, idim1 );
   H5Sget_simple_extent_dims( id_space_tmp, idimd, idim2 );

   // immediatelly drop the space handle
   H5Sclose( id_space_tmp );
   ierr = 0;
   for(n=0;n<nds;++n) {
      if( idims[n] != idimd[n] ) ierr = ierr + 1;
   }
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_MAX, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  Dataspace dimension different \n");
      free( idims );
      ierr = 3;
      goto return_point;
   }

   // create a transfer property list with parallel I/O
   // create a transfer property list with or without parallel I/O
   if( iop == 0 ) {
      id_plist = H5Pcreate( H5P_DATASET_XFER );
      H5Pset_dxpl_mpio( id_plist, H5FD_MPIO_COLLECTIVE );
   } else if( iop == 1 ) {
      id_plist = H5Pcreate( H5P_DATASET_XFER );
      H5Pset_dxpl_mpio( id_plist, H5FD_MPIO_INDEPENDENT );
   } else {
      id_plist = H5P_DEFAULT;
   }

   // create a memory dataspace
   id_mem = H5Screate_simple( nd, icnt, icnt );

   // select the area where the data will be written
   for(n=0;n<nds;++n) {
      istd[n] = 1;
      iblk[n] = 1;
   }
   H5Sselect_hyperslab( id_space, H5S_SELECT_SET, ioff, istd, icnt, iblk );

   // perform the parallel write operations
   ierr = H5Dread( id_data, id_type, id_mem, id_space, id_plist, data );
   if( ierr != 0 ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) fprintf( stdout, " e  There was a problem writing \n");
      H5Pclose( id_plist );
      ierr = 100;
      goto return_point;
   }

   // drop handles no longer needed
   H5Sclose( id_mem );
   if( iop == 0 || iop == 1 ) H5Pclose( id_plist );

   //drop the sizes array
   free( idims );

return_point:
   // restore previous error handler
   H5Eset_auto( id_estack, old_err_func, old_err_client_data );
   // drop error stack no longer needed
   H5Eclear( id_estack );
   H5Eclose_stack( id_estack );

   return ierr;
}

// the high level call
int inUtils_HDF_DatasetRead( MPI_Comm comm, hid_t id_data,
                             hid_t id_type, hid_t id_space,
                             int nd,
                             hsize_t* ioff,
                             hsize_t* icnt,
                             void* data )
{
   return inUtils_HDF_DatasetReadCall( 0, comm, id_data, id_type, id_space,
                                       nd, ioff, icnt, data );
}

