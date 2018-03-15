!$
!===================================================================================================
!
!   read & write array from/to binary file  
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    snap_array_read
!                               snap_array_write
!                               
!   Public type lists:          No
!
!===================================================================================================
module snap_interface 

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    implicit none
    private
    public  :: snap_array_read, snap_array_write
    
    interface  snap_array_read
        module procedure  read_array1D_real4
        module procedure  read_array1D_real8
        module procedure  read_array1D_int
        module procedure  read_array1D_log
        module procedure  read_array1D_chr
        module procedure  read_array2D_real4
        module procedure  read_array2D_real8
        module procedure  read_array2D_int
        module procedure  read_array2D_log
        module procedure  read_array2D_chr
        module procedure  read_array3D_real4
        module procedure  read_array3D_real8
        module procedure  read_array3D_int
        module procedure  read_array3D_log
        module procedure  read_array3D_chr
        module procedure  read_array4D_real4
        module procedure  read_array4D_real8
        module procedure  read_array4D_int
        module procedure  read_array4D_log
        module procedure  read_array4D_chr
        module procedure  read_array5D_real8
        module procedure  read_array6D_real8
    end interface snap_array_read
    
    interface snap_array_write
        module procedure  write_array1D_real4
        module procedure  write_array1D_real8
        module procedure  write_array1D_int
        module procedure  write_array1D_log
        module procedure  write_array1D_chr
        module procedure  write_array2D_real4
        module procedure  write_array2D_real8
        module procedure  write_array2D_int
        module procedure  write_array2D_log
        module procedure  write_array2D_chr
        module procedure  write_array3D_real4
        module procedure  write_array3D_real8
        module procedure  write_array3D_int
        module procedure  write_array3D_log
        module procedure  write_array3D_chr
        module procedure  write_array4D_real4
        module procedure  write_array4D_real8
        module procedure  write_array4D_int
        module procedure  write_array4D_log
        module procedure  write_array4D_chr
        module procedure  write_array5D_real8
        module procedure  write_array6D_real8
    end interface snap_array_write
    
    integer  :: NUM(10)                                                         ! size of array per dimension
    integer  :: LOWER(10)                                                       ! lower bould of array per dimension
    integer  :: UPPER(10)                                                       ! upper bould of array per dimension
    integer  :: i_allocate 
    
contains
    !$
    !===============================================================================================
    ! read array-1D
    !===============================================================================================
    subroutine read_array1D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in out), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array1D_real4
    
    subroutine read_array1D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array1D_real8
    
    subroutine read_array1D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in out), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array1D_int
    
    subroutine read_array1D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in out), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array1D_log
    
    subroutine read_array1D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in out), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array1D_chr
    
    !$
    !===============================================================================================
    ! read array-2D
    !===============================================================================================
    subroutine read_array2D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in out), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array2D_real4
    
    subroutine read_array2D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array2D_real8
    
    subroutine read_array2D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in out), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array2D_int
    
    subroutine read_array2D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in out), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array2D_log
    
    subroutine read_array2D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in out), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array2D_chr
    
    !$
    !===============================================================================================
    ! read array-3D
    !===============================================================================================
    subroutine read_array3D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in out), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array3D_real4

    subroutine read_array3D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array3D_real8

    subroutine read_array3D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in out), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array3D_int

    subroutine read_array3D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in out), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array3D_log

    subroutine read_array3D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in out), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array3D_chr

    !$
    !===============================================================================================
    ! read array-4D
    !===============================================================================================
    subroutine read_array4D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in out), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array4D_real4

    subroutine read_array4D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array4D_real8

    subroutine read_array4D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in out), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array4D_int

    subroutine read_array4D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in out), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array4D_log

    subroutine read_array4D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in out), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4)), stat=i_allocate)
            read(unit=unit_)  array
        end if 
    end subroutine read_array4D_chr

    !$
    !===============================================================================================
    ! read array-5D & 6D
    !===============================================================================================
    subroutine read_array5D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:, :, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        NUM(5) = SIZE(array, dim=5); LOWER(5) = LBOUND(array, dim=5); UPPER(5) = UBOUND(array, dim=5);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        read(unit=unit_)  NUM(5), LOWER(5), UPPER(5)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4), LOWER(5):UPPER(5)), stat=i_allocate)
            read(unit=unit_)  array
        end if
    end subroutine read_array5D_real8
    
    subroutine read_array6D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in out), allocatable  :: array(:, :, :, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        NUM(5) = SIZE(array, dim=5); LOWER(5) = LBOUND(array, dim=5); UPPER(5) = UBOUND(array, dim=5);
        NUM(6) = SIZE(array, dim=6); LOWER(6) = LBOUND(array, dim=6); UPPER(6) = UBOUND(array, dim=6);
        read(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        read(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        read(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        read(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        read(unit=unit_)  NUM(5), LOWER(5), UPPER(5)
        read(unit=unit_)  NUM(6), LOWER(6), UPPER(6)
        if (NUM(1) > 0)  then
            if (allocated(array))   deallocate(array)
            allocate(array(LOWER(1):UPPER(1), LOWER(2):UPPER(2), LOWER(3):UPPER(3), LOWER(4):UPPER(4), LOWER(5):UPPER(5), LOWER(6):UPPER(6)), stat=i_allocate)
            read(unit=unit_)  array
        end if
    end subroutine read_array6D_real8

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! write array-1D
    !===============================================================================================
    subroutine write_array1D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in), allocatable  :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array1D_real4
    
    subroutine write_array1D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable      :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array1D_real8
    
    subroutine write_array1D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in), allocatable      :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array1D_int
    
    subroutine write_array1D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in), allocatable      :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array1D_log
    
    subroutine write_array1D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in), allocatable      :: array(:)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array1D_chr
    
    !$
    !===============================================================================================
    ! write array-2D
    !===============================================================================================
    subroutine write_array2D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array2D_real4
    
    subroutine write_array2D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array2D_real8
    
    subroutine write_array2D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array2D_int
    
    subroutine write_array2D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array2D_log
    
    subroutine write_array2D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in), allocatable  :: array(:, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array2D_chr
    
    !$
    !===============================================================================================
    ! write array-3D
    !===============================================================================================
    subroutine write_array3D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array3D_real4

    subroutine write_array3D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array3D_real8

    subroutine write_array3D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array3D_int

    subroutine write_array3D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array3D_log

    subroutine write_array3D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in), allocatable  :: array(:, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array3D_chr

    !$
    !===============================================================================================
    ! write array-4D
    !===============================================================================================
    subroutine write_array4D_real4 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(4), intent(in), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array4D_real4

    subroutine write_array4D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array4D_real8

    subroutine write_array4D_int (unit_, array)
        integer, intent(in)                   :: unit_ 
        integer, intent(in), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array4D_int

    subroutine write_array4D_log (unit_, array)
        integer, intent(in)                   :: unit_ 
        logical, intent(in), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array4D_log

    subroutine write_array4D_chr (unit_, array)
        integer, intent(in)                   :: unit_ 
        character(len=*), intent(in), allocatable  :: array(:, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array4D_chr

    !$
    !===============================================================================================
    ! write array-5D & 6D
    !===============================================================================================
    subroutine write_array5D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable  :: array(:, :, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        NUM(5) = SIZE(array, dim=5); LOWER(5) = LBOUND(array, dim=5); UPPER(5) = UBOUND(array, dim=5);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        write(unit=unit_)  NUM(5), LOWER(5), UPPER(5)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array5D_real8
    
    subroutine write_array6D_real8 (unit_, array)
        integer, intent(in)                   :: unit_ 
        real(8), intent(in), allocatable  :: array(:, :, :, :, :, :)
        NUM(1) = SIZE(array, dim=1); LOWER(1) = LBOUND(array, dim=1); UPPER(1) = UBOUND(array, dim=1);
        NUM(2) = SIZE(array, dim=2); LOWER(2) = LBOUND(array, dim=2); UPPER(2) = UBOUND(array, dim=2);
        NUM(3) = SIZE(array, dim=3); LOWER(3) = LBOUND(array, dim=3); UPPER(3) = UBOUND(array, dim=3);
        NUM(4) = SIZE(array, dim=4); LOWER(4) = LBOUND(array, dim=4); UPPER(4) = UBOUND(array, dim=4);
        NUM(5) = SIZE(array, dim=5); LOWER(5) = LBOUND(array, dim=5); UPPER(5) = UBOUND(array, dim=5);
        NUM(6) = SIZE(array, dim=6); LOWER(6) = LBOUND(array, dim=6); UPPER(6) = UBOUND(array, dim=6);
        write(unit=unit_)  NUM(1), LOWER(1), UPPER(1)
        write(unit=unit_)  NUM(2), LOWER(2), UPPER(2)
        write(unit=unit_)  NUM(3), LOWER(3), UPPER(3)
        write(unit=unit_)  NUM(4), LOWER(4), UPPER(4)
        write(unit=unit_)  NUM(5), LOWER(5), UPPER(5)
        write(unit=unit_)  NUM(6), LOWER(6), UPPER(6)
        if (NUM(1) > 0)  write(unit=unit_)  array
    end subroutine write_array6D_real8
    
end module snap_interface 
