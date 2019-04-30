cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to create a file in parallel
c **** The routine accepts a communicator over which to make the parallel
c **** file creation and return error codes. It accepts a filename string
c **** and returns valid file/prop-list handles on success.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CreateFile( filename, icom,
     &                        id_alist, id_file, ierr )

       Use MPI
       Use HDF5
       Implicit None

       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom
       Character*(*),intent(INOUT) :: filename
       Integer(hid_t),intent(INOUT) :: id_file, id_alist
       Integer(hid_t) :: id_clist
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Creating file in parallel: ',trim(filename)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       !--- create creation list and file access property list with parallel I/O
       call h5pcreate_f(H5P_FILE_CREATE_F, id_clist, iher)
       call h5pcreate_f(H5P_FILE_ACCESS_F, id_alist, iher)
       call h5pset_fapl_mpio_f(id_alist, icom, MPI_INFO_NULL, iher)

       call h5fcreate_f( filename, H5F_ACC_TRUNC_F, id_file, iher,
     &                   id_clist, id_alist )
       call MPI_ALLREDUCE(iher,ierr,1,MPI_INTEGER,MPI_SUM,icom,ier)
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully created file in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to create file in parallel'
          ierr = 1
          if(iher.eq.0) then
             call h5fclose_f( id_file, iher )
          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
          call h5pclose_f(id_alist, iher)
       endif

       call h5pclose_f(id_clist, iher)

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to close an open parallel file and its prop. lists
c **** The routine accepts a communicator over which the file is opened
c **** and the associated property lists.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CloseFile( icom,
     &                        id_alist, id_file, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_file, id_alist
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Closing file'
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5fclose_f( id_file, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully closed file in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to close file in parallel'
          ierr = 1
          if(iher.eq.0) then

          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       call h5pclose_f(id_alist, iher)

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to open an existing file in parallel
c **** The routine accepts a communicator over which to make the parallel
c **** file creation and return error codes. It accepts a filename string
c **** and returns valid file/prop-list handles on success.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_OpenFile( filename, icom,
     &                        id_alist, id_file, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom
       Character*(*),intent(INOUT) :: filename
       Integer(hid_t),intent(INOUT) :: id_file, id_alist
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Opening file in parallel: ',trim(filename)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       !--- create creation list and file access property list with parallel I/O
       call h5pcreate_f(H5P_FILE_ACCESS_F, id_alist, iher)
       call h5pset_fapl_mpio_f(id_alist, icom, MPI_INFO_NULL, iher)

       call h5fopen_f( filename, H5F_ACC_RDWR_F, id_file, iher, id_alist )
       call MPI_ALLREDUCE(iher,ierr,1,MPI_INTEGER,MPI_SUM,icom,ier)
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully opened RDWR file in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to open RDWR file in parallel'
          ierr = 1
          if(iher.eq.0) then
             call h5fclose_f( id_file, iher )
          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
          call h5pclose_f(id_alist, iher)
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to create a group inside an HDF5 handle
c **** The routine accepts a communicator over which the handle is open
c **** and the associated property lists.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CreateGroup( grpname, icom,
     &                        id_alist, id_tree, id_grp, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Character*(*),intent(INOUT) :: grpname
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_tree, id_alist, id_grp
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Creating group in parallel: ',trim(grpname)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5gcreate_f( id_tree, trim(grpname), id_grp, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully created group in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to create group in parallel'
          ierr = 1
          if(iher.eq.0) then

          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to close an open group handle in a parallel handle.
c **** The routine accepts a communicator over which the file is opened
c **** and the associated property lists.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CloseGroup( icom,
     &                        id_alist, id_tree, id_grp, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_tree, id_alist, id_grp
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher

       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Closing group'
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5gclose_f( id_grp, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully closed group in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to close group in parallel'
          ierr = 1
          if(iher.eq.0) then

          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to open an existing group handle in a parallel handle.
c **** The routine accepts a communicator over which the file is opened
c **** and the associated property lists.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190423
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190423
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_OpenGroup( grpname, icom,
     &                        id_alist, id_tree, id_grp, ierr )

       Use MPI
       Use HDF5
       Use mgds_mpivars
       Use mgds_hdfvars
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Character*(*),intent(INOUT) :: grpname
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_tree, id_alist, id_grp
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Opening existing group: ',trim(grpname)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5gopen_f( id_tree, trim(grpname), id_grp, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully opened group in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to open group in parallel'
          ierr = 1
          if(iher.eq.0) then
             call h5gclose_f( id_grp, iher )
          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to create a dataset inside an HDF5 handle
c **** The routine accepts a communicator over which the handle is open
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190424
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190424
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CreateDataset( dataname, icom,
     &                        id_tree, id_type, id_space, id_data, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Character*(*),intent(INOUT) :: dataname
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_tree, id_type, id_space, id_data
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Creating dataset in parallel: ',trim(dataname)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5dcreate_f( id_tree, trim(dataname),
     &                   id_type, id_space, id_data, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully created dataset in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to create dataset in parallel'
          ierr = 1
          if(iher.eq.0) then
             call h5dclose_f( id_data, iher )
          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to close an open dataset
c **** The routine accepts a communicator over which the handle is open
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190424
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190424
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_CloseDataset( icom, id_data, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_data
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Closing dataset in parallel'
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5dclose_f( id_data, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully closed dataset in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to close dataset in parallel'
          ierr = 1
          if(iher.eq.0) then

          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to open an existing dataset inside an HDF5 handle
c **** The routine accepts a communicator over which the handle is open
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190424
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190424
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_OpenDataset( dataname, icom,
     &                        id_tree, id_data, id_type, id_space, ierr )

       Use MPI
       Use HDF5
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Character*(*),intent(INOUT) :: dataname
       Integer*4,intent(IN) :: icom
       Integer(hid_t),intent(INOUT) :: id_tree, id_type, id_space, id_data
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Opening existing dataset in parallel: ',
     &                        trim(dataname)
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       call h5dopen_f( id_tree, trim(dataname), id_data, iher )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if(ierr.eq.0) then
#ifdef _DEBUG_
          if(irank.eq.0) PRINT*,'Successfully opened dataset in parallel'
#endif
       else
          if(irank.eq.0) PRINT*,'Failed to open dataset in parallel'
          ierr = 1
          if(iher.eq.0) then
             call h5dclose_f( id_data, iher )
          else
             !-- when errors are not trapped, we clean the HDF5 error stack
             call h5eclear_f( iher )
          endif
       endif

       !-- retrieve the space of the dataset
       call h5dget_space_f( id_data, id_space, iher )

       !-- retrieve the datatype in the dataset
       call h5dget_type_f( id_data, id_type, iher )

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to write data to an open dataset in parallel.
c **** The routine requires the pointer to the data, counts and offsets,
c **** and two data spaces describing the data on disc and in memory that
c **** have been pre-defined.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190424
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190424
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_DatasetWrite( icom, id_data, id_type, id_space,
     &                                      nd, ioff, icnt, data, ierr )

       Use MPI
       Use HDF5
       Use ISO_C_BINDING
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom, nd
       Integer(hid_t),intent(IN) :: id_data, id_type, id_space
       Integer(hsize_t),intent(IN) :: icnt(nds), ioff(nds)
       Type(C_PTR),intent(IN) :: data
       Integer(hid_t) :: id_dtype, id_dspace, id_mem, id_plist
       Integer*4 :: nds,ndd                                         ! rank check
       Integer(hsize_t) :: idims(nd),idimd(nd), idim1(nd),idim2(nd) ! dims check
       Logical :: flag
       Integer*4 :: n
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Writing to existing dataset in parallel'
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       !-- retrieve the datatype in the dataset
       call h5dget_type_f( id_data, id_dtype, iher )

       !-- check for concistency of data types
       call h5tequal_f( id_type, id_dtype, flag, iher )
       !-- immediatelly release the datatype
       call h5tclose_f( id_dtype, iher )
       if( flag .eqv. .FALSE. ) then
          if(irank.eq.0) write(*,*) 'Error: datatype is inconcistent'
          ierr = 1
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype check passed'
#endif

       !-- retrieve the space of the dataset
       call h5dget_space_f( id_data, id_dspace, iher )

       !-- check for concistency of data space rank
       call h5sget_simple_extent_ndims_f( id_space, nds, iher )
       call h5sget_simple_extent_ndims_f( id_dspace, ndd, iher )
       if( nds .ne. ndd .OR. nds .ne. nd .OR. ndd .ne. nd ) then
          if(irank.eq.0) write(*,*) 'Error: dataspace ranks are different'
          if(irank.eq.0) write(*,*) '- sizes: nd/nds/ndd =',nd,nds,ndd
          call h5sclose_f( id_dspace, iher )
          ierr = 2
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype rank passed'
#endif

       !-- check for concistency of data space dimensions
       call h5sget_simple_extent_dims_f( id_space, idims, idim1, iher )
       call h5sget_simple_extent_dims_f( id_dspace, idimd, idim2, iher )
       !-- immediatelly drop the space handle
       call h5sclose_f( id_dspace, iher )
       ierr = 0
       do n = 1,nds
          if( idims(n) .ne. idimd(n) ) ierr = ierr + 1
       enddo
       if( ierr.ne.0 ) then
          if(irank.eq.0) write(*,*) 'Error: dataspace dimension are different'
          if(irank.eq.0) write(*,*) '       Number of diffs: ',ierr
          ierr = 3
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype dimenions passed'
#endif

       !--- create a transfer property list with parallel I/O
       call h5pcreate_f( H5P_DATASET_XFER_F, id_plist, iher )
       call h5pset_dxpl_mpio_f( id_plist, H5FD_MPIO_COLLECTIVE_F, iher )

       !--- create a memory dataspace
       call h5screate_simple_f( nd, icnt, id_mem, iher )

       !--- select the area where the data will be written
       call h5sselect_hyperslab_f( id_space, H5S_SELECT_SET_F, ioff, icnt, iher)

       !--- perform the parallel write operations
       call h5dwrite_f( id_data, id_type, data, iher,
     &                  id_mem, id_space, id_plist )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if( ierr.ne.0 ) then
          if(irank.eq.0) write(*,*) 'Error: there was a problem writing: ',ierr
          !-- when errors are not trapped, we clean the HDF5 error stack
          call h5eclear_f( iher )
          ierr = 100
#ifdef _DEBUG_
       else
          if(irank.eq.0) write(*,*) '- Completed writing operation'
#endif
       endif

       !-- drop the memory dataspace
       call h5sclose_f( id_mem, iher )

       !-- drop the property list
       call h5pclose_f( id_plist, iher )

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c **** Subroutine to read data from an open dataset in parallel.
c **** The routine requires the pointer to a buffer, counts and offsets,
c **** and two data spaces describing the data on disc and in memory that
c **** have been pre-defined.
c Ioannis Nompelis <nompelis@nobelware.com>     Created:       20190429
c Ioannis Nompelis <nompelis@nobelware.com>     Last modified: 20190429
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       Subroutine inUtils_HDF_DatasetRead( icom, id_data, id_type, id_space,
     &                                     nd, ioff, icnt, data, ierr )

       Use MPI
       Use HDF5
       Use ISO_C_BINDING
       Implicit None
       Integer*4,intent(OUT) :: ierr
       Integer*4,intent(IN) :: icom, nd
       Integer(hid_t),intent(IN) :: id_data, id_type, id_space
       Integer(hsize_t),intent(IN) :: icnt(nds), ioff(nds)
       Type(C_PTR),intent(INOUT) :: data
       Integer(hid_t) :: id_dtype, id_dspace, id_mem, id_plist
       Integer*4 :: nds,ndd                                         ! rank check
       Integer(hsize_t) :: idims(nd),idimd(nd), idim1(nd),idim2(nd) ! dims check
       Logical :: flag
       Integer*4 :: n
       Integer*4 :: irank,nrank
       Integer*4 :: ier,iher


       call MPI_COMM_RANK( icom, irank, ier )
       call MPI_COMM_SIZE( icom, nrank, ier )
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) 'Reading from dataset in parallel'
#endif

       !-- turn off HDF5 error printing
       call  h5eset_auto_f( 0, iher )

       !-- retrieve the datatype in the dataset
       call h5dget_type_f( id_data, id_dtype, iher )

       !-- check for concistency of data types
       call h5tequal_f( id_type, id_dtype, flag, iher )
       !-- immediatelly release the datatype
       call h5tclose_f( id_dtype, iher )
       if( flag .eqv. .FALSE. ) then
          if(irank.eq.0) write(*,*) 'Error: datatype is inconcistent'
          ierr = 1
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype check passed'
#endif

       !-- retrieve the space of the dataset
       call h5dget_space_f( id_data, id_dspace, iher )

       !-- check for concistency of data space rank
       call h5sget_simple_extent_ndims_f( id_space, nds, iher )
       call h5sget_simple_extent_ndims_f( id_dspace, ndd, iher )
       if( nds .ne. ndd .OR. nds .ne. nd .OR. ndd .ne. nd ) then
          if(irank.eq.0) write(*,*) 'Error: dataspace ranks are different'
          if(irank.eq.0) write(*,*) '- sizes: nd/nds/ndd =',nd,nds,ndd
          call h5sclose_f( id_dspace, iher )
          ierr = 2
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype rank passed'
#endif

       !-- check for concistency of data space dimensions
       call h5sget_simple_extent_dims_f( id_space, idims, idim1, iher )
       call h5sget_simple_extent_dims_f( id_dspace, idimd, idim2, iher )
       !-- immediatelly drop the space handle
       call h5sclose_f( id_dspace, iher )
       ierr = 0
       do n = 1,nds
          if( idims(n) .ne. idimd(n) ) ierr = ierr + 1
       enddo
       if( ierr.ne.0 ) then
          if(irank.eq.0) write(*,*) 'Error: dataspace dimension are different'
          if(irank.eq.0) write(*,*) '       Number of diffs: ',ierr
          ierr = 3
          !-- restore default HDF5 library functionality before returning
          call  h5eset_auto_f( 1, iher )
          call MPI_BARRIER( icom, ier )
          return
       endif
#ifdef _DEBUG_
       if(irank.eq.0) write(*,*) '- Datatype dimenions passed'
#endif

       !--- create a transfer property list with parallel I/O
       call h5pcreate_f( H5P_DATASET_XFER_F, id_plist, iher )
       call h5pset_dxpl_mpio_f( id_plist, H5FD_MPIO_COLLECTIVE_F, iher )

       !--- create a memory dataspace
       call h5screate_simple_f( nd, icnt, id_mem, iher )

       !--- select the area where the data will be written
       call h5sselect_hyperslab_f( id_space, H5S_SELECT_SET_F, ioff, icnt, iher)

       !--- perform the parallel write operations
       call h5dread_f( id_data, id_type, data, iher,
     &                 id_mem, id_space, id_plist )
       call MPI_ALLREDUCE( iher, ierr, 1, MPI_INTEGER, MPI_SUM, icom, ier )
       if( ierr.ne.0 ) then
          if(irank.eq.0) write(*,*) 'Error: there was a problem reading: ',ierr
          !-- when errors are not trapped, we clean the HDF5 error stack
          call h5eclear_f( iher )
          ierr = 100
#ifdef _DEBUG_
       else
          if(irank.eq.0) write(*,*) '- Completed reading operation'
#endif
       endif

       !-- drop the memory dataspace
       call h5sclose_f( id_mem, iher )

       !-- drop the property list
       call h5pclose_f( id_plist, iher )

       !-- restore default HDF5 library functionality before returning
       call  h5eset_auto_f( 1, iher )

       return
       end

