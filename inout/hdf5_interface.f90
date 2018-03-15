!$
!===================================================================================================
!
!   interface for HDF5 file input & output
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    hdf5_write_data
!                               hdf5_read_data
!                               hdf5_set_attr
!                               hdf5_get_attr
!                               hdf5_open_group
!                               hdf5_close_group
!                               hdf5_file_create
!                               hdf5_file_open
!                               hdf5_file_close
!
!   Public type lists:          No
!
!===================================================================================================
module hdf5_interface

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    use, intrinsic  :: ISO_C_BINDING

    use HDF5
    use H5LT

    implicit none
    private 
    public  :: hdf5_write_data, hdf5_read_data
    public  :: hdf5_set_attr, hdf5_get_attr

    public  :: hdf5_open_group, hdf5_close_group
    public  :: hdf5_file_create, hdf5_file_open, hdf5_file_close
    
    ! --------------------------------------------------------------------------
    ! local parameter
    integer           :: hdf5_err                                               ! HDF5 error code
    integer           :: hdf5_stat                                              ! existence state
    integer           :: hdf5_rank                                              ! rank of data
    integer(SIZE_T)   :: hdf5_size                                              ! size of the 1D array
    integer(HID_T)    :: dset                                                   ! data set handle
    integer(HID_T)    :: dspace                                                 ! data or file space handle
    integer(HID_T)    :: memspace                                               ! data space handle for individual procs
    integer(HID_T)    :: plist                                                  ! property list handle
    integer(HSIZE_T)  :: dims1(1)                                               ! dims type for 1-D array
    integer(HSIZE_T)  :: dims2(2)                                               ! dims type for 2-D array
    integer(HSIZE_T)  :: dims3(3)                                               ! dims type for 3-D array
    integer(HSIZE_T)  :: dims4(4)                                               ! dims type for 4-D array

    ! Generic HDF5 write procedure interface
    interface hdf5_write_data
        module procedure hdf5_write_real
        module procedure hdf5_write_real_1Darray
        module procedure hdf5_write_real_2Darray
        module procedure hdf5_write_real_3Darray
        module procedure hdf5_write_real_4Darray
        module procedure hdf5_write_integer
        module procedure hdf5_write_integer_1Darray
        module procedure hdf5_write_integer_2Darray
        module procedure hdf5_write_integer_3Darray
        module procedure hdf5_write_integer_4Darray
    end interface hdf5_write_data

    ! Generic HDF5 read procedure interface
    interface hdf5_read_data
        module procedure hdf5_read_real
        module procedure hdf5_read_real_1Darray
        module procedure hdf5_read_real_2Darray
        module procedure hdf5_read_real_3Darray
        module procedure hdf5_read_real_4Darray
        module procedure hdf5_read_integer
        module procedure hdf5_read_integer_1Darray
        module procedure hdf5_read_integer_2Darray
        module procedure hdf5_read_integer_3Darray
        module procedure hdf5_read_integer_4Darray
    end interface hdf5_read_data

    ! Generic HDF5 set attr procedure interface
    interface hdf5_set_attr
        module procedure hdf5_set_attr_string
        module procedure hdf5_set_attr_int
        module procedure hdf5_set_attr_int_1Darray
        module procedure hdf5_set_attr_real
        module procedure hdf5_set_attr_real_1Darray
    end interface hdf5_set_attr
    
    ! Generic HDF5 get attr procedure interface
    interface hdf5_get_attr
        module procedure hdf5_get_attr_string
        module procedure hdf5_get_attr_int
        module procedure hdf5_get_attr_int_1Darray
        module procedure hdf5_get_attr_real
        module procedure hdf5_get_attr_real_1Darray
    end interface hdf5_get_attr
    
contains
    !$
    !===============================================================================================
    ! HDF5_FILE_CREATE creates HDF5 file
    !===============================================================================================
    subroutine hdf5_file_create(filename, file_id)

        character(len=*),   intent(in)    :: filename                           ! name of file
        integer(HID_T), intent(inout) :: file_id                                ! file handle
    
        ! Create the file
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)

    end subroutine hdf5_file_create

    !$
    !===============================================================================================
    ! HDF5_FILE_OPEN opens HDF5 file
    !===============================================================================================
    subroutine hdf5_file_open(filename, file_id, mode)

        character(len=*),  intent(in)      :: filename                          ! name of file
        integer(HID_T), intent(inout)  :: file_id                               ! file handle
        character(len=*),  intent(in)      :: mode                              ! access mode to file
    
        integer :: open_mode                                                    ! HDF5 open mode
    
        ! Determine access type
        open_mode = H5F_ACC_RDONLY_F
        if (trim(mode) == 'w') then
            open_mode = H5F_ACC_RDWR_F
        end if
    
        ! Open file
        call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)

    end subroutine hdf5_file_open

    !$
    !===============================================================================================
    ! HDF5_FILE_CLOSE closes HDF5 file
    !===============================================================================================
    subroutine hdf5_file_close(file_id)

        integer(HID_T), intent(inout) :: file_id                                ! file handle
    
        ! Close the file
        call h5fclose_f(file_id, hdf5_err)

    end subroutine hdf5_file_close
    
    !$
    !===============================================================================================
    ! HDF5_OPEN_GROUP creates/opens HDF5 group to temp_group
    !===============================================================================================
    subroutine hdf5_open_group(hdf5_fh, group, hdf5_grp)

        integer(HID_T), intent(in)    :: hdf5_fh                                ! file handle of main output file
        character(len=*),   intent(in)    :: group                              ! name of group
        integer(HID_T), intent(inout) :: hdf5_grp                               ! handle for group
    
        logical :: status                                                       ! does the group exist
    
        ! Check if group exists
        call h5ltpath_valid_f(hdf5_fh, trim(group), .true., status, hdf5_err) 
    
        ! Either create or open group
        if (status) then
            call h5gopen_f(hdf5_fh, trim(group), hdf5_grp, hdf5_err)
        else
            call h5gcreate_f(hdf5_fh, trim(group), hdf5_grp, hdf5_err)
        end if

    end subroutine hdf5_open_group

    !$
    !===============================================================================================
    ! HDF5_CLOSE_GROUP closes HDF5 temp_group
    !===============================================================================================
    subroutine hdf5_close_group(hdf5_grp)

        integer(HID_T), intent(inout) :: hdf5_grp
    
        ! Close the group
        call h5gclose_f(hdf5_grp, hdf5_err)
    
    end subroutine hdf5_close_group
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! HDF5_WRITE_INTEGER writes integer scalar data
    !===============================================================================================
    subroutine hdf5_write_integer(group, name, buffer)

        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        integer,        intent(in) :: buffer                                    ! data to write
    
        ! Set rank and dimensions
        hdf5_rank = 1
        dims1(1) = 1
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, (/ buffer /), hdf5_err)

    end subroutine hdf5_write_integer

    !$
    !===============================================================================================
    ! HDF5_READ_INTEGER reads integer scalar data
    !===============================================================================================
    subroutine hdf5_read_integer(group, name, buffer)

        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        integer,        intent(inout) :: buffer                                 ! read data to here 
    
        integer :: buffer_copy(1)                                               ! need an array for read
    
        ! Set up dimensions
        dims1(1) = 1
    
        ! Read data
        call h5ltread_dataset_int_f(group, name, buffer_copy, dims1, hdf5_err)
        buffer = buffer_copy(1)

    end subroutine hdf5_read_integer

    !$
    !===============================================================================================
    ! HDF5_WRITE_INTEGER_1DARRAY writes integer 1-D array
    !===============================================================================================
    subroutine hdf5_write_integer_1Darray(group, name, buffer, length)

        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        integer,        intent(in) :: buffer(:)                                 ! data to write
        integer,        intent(in) :: length(1)                                 ! length of array to write
    
        ! Set rank and dimensions of data
        hdf5_rank = 1
        dims1(1) = length(1)
    
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, buffer, hdf5_err)

    end subroutine hdf5_write_integer_1Darray

    !$
    !===============================================================================================
    ! HDF5_READ_INTEGER_1DARRAY reads integer 1-D array
    !===============================================================================================
    subroutine hdf5_read_integer_1Darray(group, name, buffer, length)

        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        integer,        intent(inout) :: buffer(:)                              ! read data to here
        integer,        intent(in)    :: length(1)                              ! length of array

        ! Set dimensions
        dims1(1) = length(1)
    
        ! Read data
        call h5ltread_dataset_int_f(group, name, buffer, dims1, hdf5_err)
    
    end subroutine hdf5_read_integer_1Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_INTEGER_2DARRAY writes integer 2-D array
    !===============================================================================================
    subroutine hdf5_write_integer_2Darray(group, name, buffer, length)

        integer,        intent(in) :: length(2)                                 ! length of array dimensions
        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        integer,        intent(in) :: buffer(length(1),length(2))               ! data to write
    
        ! Set rank and dimensions
        hdf5_rank = 2
        dims2 = length
    
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims2, buffer, hdf5_err)

    end subroutine hdf5_write_integer_2Darray

    !$
    !===============================================================================================
    ! HDF5_READ_INTEGER_2DARRAY reads integer 2-D array
    !===============================================================================================
    subroutine hdf5_read_integer_2Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(2)                              ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        integer,        intent(inout) :: buffer(length(1),length(2))            ! data to read
    
        ! Set rank and dimensions
        dims2 = length
    
        ! Write data
        call h5ltread_dataset_int_f(group, name, buffer, dims2, hdf5_err)

    end subroutine hdf5_read_integer_2Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_INTEGER_3DARRAY writes integer 3-D array
    !===============================================================================================
    subroutine hdf5_write_integer_3Darray(group, name, buffer, length)

        integer,        intent(in) :: length(3)                                 ! length of array dimensions
        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        integer,        intent(in) :: buffer(length(1),length(2),length(3))     ! data to write
    
        ! Set rank and dimensions
        hdf5_rank = 3
        dims3 = length
    
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims3, buffer, hdf5_err)

    end subroutine hdf5_write_integer_3Darray

    !$
    !===============================================================================================
    ! HDF5_READ_INTEGER_3DARRAY reads integer 3-D array
    !===============================================================================================
    subroutine hdf5_read_integer_3Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(3)                              ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        integer,        intent(inout) :: buffer(length(1),length(2), length(3)) ! data to read
        
        ! Set rank and dimensions
        dims3 = length
        
        ! Write data
        call h5ltread_dataset_int_f(group, name, buffer, dims3, hdf5_err)

    end subroutine hdf5_read_integer_3Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_INTEGER_4DARRAY writes integer 4-D array
    !===============================================================================================
    subroutine hdf5_write_integer_4Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(4)                                                  ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                                      ! name of group
        character(len=*),   intent(in)    :: name                                                   ! name of data
        integer,        intent(in)    :: buffer(length(1),length(2),length(3),length(4))            ! data to write
        
        ! Set rank and dimensions
        hdf5_rank = 4
        dims4 = length
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims4, buffer, hdf5_err)
        
    end subroutine hdf5_write_integer_4Darray

    !$
    !===============================================================================================
    ! HDF5_READ_INTEGER_4DARRAY reads integer 4-D array
    !===============================================================================================
    subroutine hdf5_read_integer_4Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(4)                                                  ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                                      ! name of group
        character(len=*),   intent(in)    :: name                                                   ! name of data
        integer,        intent(inout) :: buffer(length(1),length(2), length(3),length(4))           ! data to read
        
        ! Set rank and dimensions
        dims4 = length
        
        ! Write data
        call h5ltread_dataset_int_f(group, name, buffer, dims4, hdf5_err)

    end subroutine hdf5_read_integer_4Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_REAL writes integer scalar data
    !===============================================================================================
    subroutine hdf5_write_real(group, name, buffer)

        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        real(KREAL),intent(in) :: buffer                                        ! data to write
        real(8)  :: local_bufffer
        
        ! Set rank and dimensions
        hdf5_rank = 1
        dims1(1) = 1
        local_bufffer = buffer
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, (/ local_bufffer /), hdf5_err)

    end subroutine hdf5_write_real

    !$
    !===============================================================================================
    ! HDF5_READ_REAL reads real scalar data
    !===============================================================================================
    subroutine hdf5_read_real(group, name, buffer)

        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        real(KREAL),intent(inout) :: buffer                                     ! read data to here 
        
        real(8) :: buffer_copy(1)                                               ! need an array for read
        
        ! Set up dimensions
        dims1(1) = 1
        
        ! Read data
        call h5ltread_dataset_double_f(group, name, buffer_copy, dims1, hdf5_err)
        buffer = buffer_copy(1)

    end subroutine hdf5_read_real

    !$
    !===============================================================================================
    ! HDF5_WRITE_REAL_1DARRAY writes real 1-D array
    !===============================================================================================
    subroutine hdf5_write_real_1Darray(group, name, buffer, length)

        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        real(KREAL),intent(in) :: buffer(:)                                     ! data to write
        integer,        intent(in) :: length(1)                                 ! length of array to write
        real(8)  :: local_buffer(SIZE(buffer))
        
        ! Set rank and dimensions of data
        hdf5_rank = 1
        dims1(1) = length(1)
        local_buffer = buffer
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, local_buffer, hdf5_err)

    end subroutine hdf5_write_real_1Darray

    !$
    !===============================================================================================
    ! HDF5_READ_REAL_1DARRAY reads real 1-D array
    !===============================================================================================
    subroutine hdf5_read_real_1Darray(group, name, buffer, length)

        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        real(KREAL),intent(inout) :: buffer(:)                                  ! read data to here
        integer,        intent(in)    :: length(1)                              ! length of array
        real(8)  :: local_buffer(SIZE(buffer))
        
        ! Set dimensions
        dims1(1) = length(1)
        
        ! Read data
        call h5ltread_dataset_double_f(group, name, local_buffer, dims1, hdf5_err)
        buffer = local_buffer

    end subroutine hdf5_read_real_1Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_REAL_2DARRAY writes real 2-D array
    !===============================================================================================
    subroutine hdf5_write_real_2Darray(group, name, buffer, length)

        integer,        intent(in) :: length(2)                                 ! length of array dimensions
        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        real(KREAL),intent(in) :: buffer(length(1),length(2))                   ! data to write
        real(8)  :: local_buffer(length(1),length(2))
        
        ! Set rank and dimensions
        hdf5_rank = 2
        dims2 = length
        local_buffer = buffer
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims2, local_buffer, hdf5_err)

    end subroutine hdf5_write_real_2Darray

    !$
    !===============================================================================================
    ! HDF5_READ_REAL_2DARRAY reads real 2-D array
    !===============================================================================================
    subroutine hdf5_read_real_2Darray(group, name, buffer, length)
    
        integer,        intent(in)    :: length(2)                              ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        real(KREAL),intent(inout) :: buffer(length(1),length(2))                ! data to read
        real(8)  :: local_buffer(length(1),length(2))
        
        ! Set rank and dimensions
        dims2 = length
        
        ! Write data
        call h5ltread_dataset_double_f(group, name, local_buffer, dims2, hdf5_err)
        buffer = local_buffer

    end subroutine hdf5_read_real_2Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_REAL_3DARRAY writes real 3-D array
    !===============================================================================================
    subroutine hdf5_write_real_3Darray(group, name, buffer, length)

        integer,        intent(in) :: length(3)                                 ! length of array dimensions
        integer(HID_T), intent(in) :: group                                     ! name of group
        character(len=*),   intent(in) :: name                                  ! name of data
        real(KREAL),intent(in) :: buffer(length(1),length(2),length(3))         ! data to write
        real(8)  :: local_buffer(length(1),length(2),length(3))
        
        ! Set rank and dimensions
        hdf5_rank = 3
        dims3 = length
        local_buffer = buffer
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims3, local_buffer, hdf5_err)

    end subroutine hdf5_write_real_3Darray

    !$
    !===============================================================================================
    ! HDF5_READ_REAL_3DARRAY reads real 3-D array
    !===============================================================================================
    subroutine hdf5_read_real_3Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(3)                              ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                  ! name of group
        character(len=*),   intent(in)    :: name                               ! name of data
        real(KREAL),intent(inout) :: buffer(length(1),length(2), length(3))     ! data to read
        real(8)  :: local_buffer(length(1),length(2),length(3))
        
        ! Set rank and dimensions
        dims3 = length
        
        ! Write data
        call h5ltread_dataset_double_f(group, name, local_buffer, dims3, hdf5_err)
        buffer = local_buffer

    end subroutine hdf5_read_real_3Darray

    !$
    !===============================================================================================
    ! HDF5_WRITE_REAL_4DARRAY writes real 4-D array
    !===============================================================================================
    subroutine hdf5_write_real_4Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(4)                                                  ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                                      ! name of group
        character(len=*),   intent(in)    :: name                                                   ! name of data
        real(KREAL),intent(in)    :: buffer(length(1),length(2), length(3),length(4))               ! data to write
        real(8)  :: local_buffer(length(1),length(2),length(3),length(4))
        
        ! Set rank and dimensions
        hdf5_rank = 4
        dims4 = length
        local_buffer = buffer
        
        ! if already exist, delete
        hdf5_stat = h5ltfind_dataset_f(group, name)
        if(hdf5_stat==1) then
            call h5Ldelete_f(group, name, hdf5_err)
        end if
        
        ! write dataset
        call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims4, local_buffer, hdf5_err)

    end subroutine hdf5_write_real_4Darray

    !$
    !===============================================================================================
    ! HDF5_READ_REAL_4DARRAY reads real 4-D array
    !===============================================================================================
    subroutine hdf5_read_real_4Darray(group, name, buffer, length)

        integer,        intent(in)    :: length(4)                                                  ! length of array dimensions
        integer(HID_T), intent(in)    :: group                                                      ! name of group
        character(len=*),   intent(in)    :: name                                                   ! name of data
        real(KREAL),intent(inout) :: buffer(length(1),length(2), length(3),length(4))               ! data to read
        real(8)  :: local_buffer(length(1),length(2),length(3),length(4))
        
        ! Set rank and dimensions
        dims4 = length
    
        ! Write data
        call h5ltread_dataset_double_f(group, name, local_buffer, dims4, hdf5_err)
        buffer = local_buffer

    end subroutine hdf5_read_real_4Darray
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! HDF5_SET_ATTR_STRING set string attribute
    !===============================================================================================
    subroutine hdf5_set_attr_string(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        character(len=*),   intent(in)    :: buffer                                                 ! string
        
        ! Set attr
        call h5ltset_attribute_string_f(grp_id, obj_name, attr_name, buffer, hdf5_err)

    end subroutine hdf5_set_attr_string
    
    !$
    !===============================================================================================
    ! HDF5_SET_ATTR_INT set INT attribute
    !===============================================================================================
    subroutine hdf5_set_attr_int(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group or file id
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        integer,        intent(in)    :: buffer                                                     ! integer
        
        ! Set attr
        hdf5_size = 1
        call h5ltset_attribute_int_f(grp_id, obj_name, attr_name, (/ buffer /) , hdf5_size, hdf5_err)

    end subroutine hdf5_set_attr_int
     
    !$
    !===============================================================================================
    ! HDF5_SET_ATTR_INT set INT 1D array attribute
    !===============================================================================================
    subroutine hdf5_set_attr_int_1Darray(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group or file id
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        integer,        intent(in)    :: buffer(:)                                                  ! integer array
        
        ! Set attr
        hdf5_size = SIZE(buffer)
        call h5ltset_attribute_int_f(grp_id, obj_name, attr_name, buffer, hdf5_size, hdf5_err)

    end subroutine hdf5_set_attr_int_1Darray
     
    !$
    !===============================================================================================
    ! HDF5_SET_ATTR_REAL set REAL attribute
    !===============================================================================================
    subroutine hdf5_set_attr_real(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group or file id
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        real(KREAL),intent(in)    :: buffer                                                         ! integer array
        real(8)  :: local_buffer
        
        local_buffer = buffer
        ! Set attr
        hdf5_size = 1
        call h5ltset_attribute_double_f(grp_id, obj_name, attr_name, (/ local_buffer /), hdf5_size, hdf5_err)

    end subroutine hdf5_set_attr_real
     
    !$
    !===============================================================================================
    ! HDF5_SET_ATTR_REAL_1Darray set REAL 1D array attribute
    !===============================================================================================
    subroutine hdf5_set_attr_real_1Darray(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group or file id
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        real(KREAL),intent(in)    :: buffer(:)                                                      ! integer array
        real(8),allocatable  :: local_buffer(:)
        
        ! Set attr
        hdf5_size = size(buffer)
        allocate(local_buffer(hdf5_size))
        local_buffer = buffer
        call h5ltset_attribute_double_f(grp_id, obj_name, attr_name, local_buffer, hdf5_size, hdf5_err)

    end subroutine hdf5_set_attr_real_1Darray
    
    !$
    !===============================================================================================
    ! HDF5_GET_ATTR_STRING get string attribute
    !===============================================================================================
    subroutine hdf5_get_attr_string(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        character(len=*),   intent(inout) :: buffer                                                 ! string
        
        ! Set attr
        call h5ltget_attribute_string_f(grp_id, obj_name, attr_name, buffer, hdf5_err)

    end subroutine hdf5_get_attr_string
    
    !$
    !===============================================================================================
    ! HDF5_GET_ATTR_INT get INT attribute
    !===============================================================================================
    subroutine hdf5_get_attr_int(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        integer,        intent(inout) :: buffer                                                     ! integer
        
        integer :: buffer_copy(1)
        
        ! Set attr
        call h5ltget_attribute_int_f(grp_id, obj_name, attr_name, buffer_copy , hdf5_err)
        buffer = buffer_copy(1)

    end subroutine hdf5_get_attr_int
    
    !$
    !===============================================================================================
    ! HDF5_GET_ATTR_INT_1Darray get INT 1D array attribute
    !===============================================================================================
    subroutine hdf5_get_attr_int_1Darray(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        integer,        intent(inout) :: buffer(:)                                                  ! integer
        
        ! Set attr
        call h5ltget_attribute_int_f(grp_id, obj_name, attr_name, buffer , hdf5_err)

    end subroutine hdf5_get_attr_int_1Darray
       
    !$
    !===============================================================================================
    ! HDF5_GET_ATTR_REAL get REAL attribute
    !===============================================================================================
    subroutine hdf5_get_attr_real(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        real(KREAL),intent(inout) :: buffer                                                         ! real
        
        real(8) :: buffer_copy(1)                                               ! need an array for read
        
        ! Set attr
        call h5ltget_attribute_double_f(grp_id, obj_name, attr_name, buffer_copy , hdf5_err)

        buffer = buffer_copy(1)

    end subroutine hdf5_get_attr_real
    
    !$
    !===============================================================================================
    ! HDF5_GET_ATTR_REAL_1Darray get REAL 1D array attribute
    !===============================================================================================
    subroutine hdf5_get_attr_real_1Darray(grp_id, obj_name, attr_name, buffer)

        integer(HID_T), intent(in)    :: grp_id                                                     ! id of upper group
        character(len=*),   intent(in)    :: obj_name                                               ! name of object
        character(len=*),   intent(in)    :: attr_name                                              ! name of attribute
        real(KREAL),intent(inout) :: buffer(:)                                                      ! real 1D array
        
        real(8),allocatable :: buffer_copy(:)                                   ! need an array for read
        
        hdf5_size = size(buffer)
        allocate(buffer_copy(hdf5_size))
        
        ! Set attr
        call h5ltget_attribute_double_f(grp_id, obj_name, attr_name, buffer_copy , hdf5_err)

        buffer = buffer_copy

    end subroutine hdf5_get_attr_real_1Darray
    
end module hdf5_interface
