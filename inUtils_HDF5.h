/******************************************************************************
  *** Function prototypes for an example of how to use the INHDF5 API
  *** 
  *** Copyright (c) 2020 Ioannis Nompelis
  *** Ioannis Nompelis <nompelis@nobelware.com>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <hdf5.h>


int inUtils_HDF_CreateFile( char* filename, MPI_Comm comm,
                            hid_t* id_alist, 
                            hid_t* id_file );

int inUtils_HDF_OpenFile( char* filename, MPI_Comm comm,
                          hid_t* id_alist, 
                          hid_t* id_file );

int inUtils_HDF_CloseFile( MPI_Comm comm,
                           hid_t* id_alist,
                           hid_t* id_file );


int inUtils_HDF_CreateGroup( char* grpname, MPI_Comm comm,
                             hid_t id_tree,
                             hid_t* id_grp );

int inUtils_HDF_CloseGroup( MPI_Comm comm, hid_t* id_grp );


int inUtils_HDF_OpenGroup( char* grpname, MPI_Comm comm,
                           hid_t id_tree,
                           hid_t* id_grp );


int inUtils_HDF_CreateDataset( char* dataname, MPI_Comm comm,
                               hid_t id_tree,
                               hid_t id_type,
                               hid_t id_space,
                               hid_t* id_data );

int inUtils_HDF_CloseDataset( MPI_Comm comm, hid_t* id_data );

int inUtils_HDF_OpenDataset( char* dataname, MPI_Comm comm, hid_t id_tree,
                             hid_t* id_data,
                             hid_t* id_type,
                             hid_t* id_space );

int inUtils_HDF_DatasetWrite( MPI_Comm comm, hid_t id_data,
                              hid_t id_type, hid_t id_space,
                              int nd,
                              hsize_t* ioff,
                              hsize_t* icnt,
                              void* data );

int inUtils_HDF_DatasetRead( MPI_Comm comm, hid_t id_data,
                             hid_t id_type, hid_t id_space,
                             int nd,
                             hsize_t* ioff,
                             hsize_t* icnt,
                             void* data );

