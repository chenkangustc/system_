!$
!===================================================================================================
!
!   class for code compile-time selection
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Check_system
!                               Check_compiler
!                               Set_openMP
!
!   Public type lists:          No
!
!===================================================================================================
module environment
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    use omp_lib
    
    implicit none
    private
    public  :: Check_system, Check_compiler, Set_openMP
    
    ! --------------------------------------------------------------------------
    ! logical for compile option, currently is only for openMP
    type, private  ::  openmp_switch_tp
        logical, public          :: is_openmp        = .FALSE.                  ! active openMP ?
        integer, public          :: threads_limit    = 1 
        integer, public          :: direction_limit  = 24                       ! lower limit of direction number for parallel with openMP
    end type openmp_switch_tp
    
    type(openmp_switch_tp)  :: a_openmp
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Check_system ()
    
        character(len=MAX_WORD_LEN)  :: system_with

        ! check system
    
    end subroutine Check_system
 
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Check_compiler ()
        
        character(len=MAX_WORD_LEN)  :: compiled_by 
        character(len=MAX_WORD_LEN)  :: compiled_with
        
        ! check for version 
    
    end subroutine Check_compiler

    !$
    !===============================================================================================
    ! set executive parameters for openMP environment
    !===============================================================================================
    subroutine Set_openMP (nthread)
        
        integer, intent(in)  :: nthread
        
        integer  :: reserve_thread = 0                                          ! number of reserve thread for not parallel execute
        logical  :: is_nested      = .TRUE.                                     ! is neseted parallel ?
        integer  :: max_level      = 2                                          ! level of nested parallel
        logical  :: is_dynamic     = .TRUE.                                     ! is dynamic adjust ?
        
        if (nthread <= 0)  then
            a_openmp%is_openmp = .FALSE.
        else
            a_openmp%is_openmp = .TRUE.
        end if 
        
        if (a_openmp%is_openmp)  then
            a_openmp%threads_limit = MIN(omp_get_max_threads () - reserve_thread, nthread)
            call omp_set_num_threads (a_openmp%threads_limit)                   ! set threads number for parallel region
            call omp_set_nested (is_nested)                                     ! set nested parallel logocal
            call omp_set_dynamic (is_dynamic)                                   ! set dynamic adjust logocal
            call omp_set_max_active_levels (max_level)                          ! set max active level
        else
            a_openmp%threads_limit = 1 
            call omp_set_num_threads (a_openmp%threads_limit)                   ! set threads number for parallel region
            call omp_set_nested (is_nested)                                     ! set nested parallel logocal
            call omp_set_dynamic (is_dynamic)                                   ! set dynamic adjust logocal
            call omp_set_max_active_levels (max_level)                          ! set max active level
        end if 
        
        write(*, '(1x, A, I4)') 'threads_limit = ', a_openmp%threads_limit 
    
    end subroutine Set_openMP
    
end module environment
