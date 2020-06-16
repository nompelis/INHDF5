ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *** An example of how to use the INHDF5 API
c *** 
c *** Copyright (c) 2020 Ioannis Nompelis
c *** Ioannis Nompelis <nompelis@nobelware.com>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Program inHDF_example
       Use MPI
       Use HDF5
       Implicit none
       Character*100 :: filename, groupname, datasetname
       Integer*4 :: ierr, irank, nrank
       Integer*4 :: icom
       Integer(HID_T) :: id_alist, id_file, id_grp, id_data
       Integer(HID_T) :: id_type, id_space
       Integer(HSIZE_T) :: idims(1),ioff(1),icnt(1)
       Integer*4 :: hdf_err
       Integer*4,dimension(10), TARGET :: idum    !--- needs to be imovable!
       Type(C_PTR) :: ptr


       !--- start the MPI API
       call MPI_INIT( ierr )
       icom = MPI_COMM_WORLD
       call MPI_COMM_RANK( icom, irank, ierr )
       call MPI_COMM_SIZE( icom, nrank, ierr )

       if( nrank .gt. 2 ) then
          PRINT*,'This demo is meant to run with 2 processes'
          call MPI_FINALIZE( ierr )
          STOP
       endif

       !--- start the HDF5 API
       call h5open_f( hdf_err )

       !--- create a file, group, dataset
       filename = 'test.h5'
       call inUtils_HDF_CreateFile( filename, icom,
     &                        id_alist, id_file, ierr )

       groupname = 'my_group'
       call inUtils_HDF_CreateGroup( groupname, icom,
     &                        id_alist, id_file, id_grp, ierr )

       datasetname = 'my_dataset'
       id_type = H5T_NATIVE_INTEGER
       idims(1) = 20
       call h5screate_simple_f( 1, idims, id_space, hdf_err ) 

       call inUtils_HDF_CreateDataset( datasetname, icom,
     &                        id_grp, id_type, id_space, id_data, ierr )

       !--- write to dataset
       icnt(1) = 10    !--- all proceses write 10 integers
       if( irank .eq. 0 ) then
          ioff(1) = 0    !---- offsets for process in file
          icnt(1) = 10   !---- count for process
          idum(1:10) = 1
       endif
       if( irank .eq. 1 ) then
          ioff(1) = 10   !---- offsets for process in file
          icnt(1) = 10   !---- count for process
          idum(1:10) = 2
       endif
       ptr = C_LOC( idum(1) )         !--- get the C pointer
       call inUtils_HDF_DatasetWrite( icom, id_data, id_type, id_space,
     &                                      1, ioff, icnt, ptr, ierr )

       !-- close the space no longer needed
       call h5sclose_f( id_space, hdf_err ) 

       !-- close dataset, group, file
       call inUtils_HDF_CloseDataset( icom, id_data, ierr )
       call inUtils_HDF_CloseGroup( icom,
     &                        id_alist, id_file, id_grp, ierr )
       call inUtils_HDF_CloseFile( icom,
     &                        id_alist, id_file, ierr )

       !--- close the HDF5 and MPI APIs
       call h5close_f( hdf_err )
       call MPI_FINALIZE( ierr )

       End program

