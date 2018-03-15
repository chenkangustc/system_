!$
!===================================================================================================
!
!   control subroutine for adjoint calculation at initial condition, with various type define
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_adjoint_gpt
!
!   Public type lists:          No
!
!===================================================================================================
module driver_adjoint_gpt
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    
    use transit_to_solver,          only : Transit_xsec_gpt
    use iteration_initialize,       only : Init_iteration_variable
    use iteration_control,          only : Driving_iteration
    use reactivity,                 only : Get_gpt_source
    
    use output_hdf5,                only : Print_binary_hdf5
    use output_visit,               only : Print_vtk_files
    
    implicit none
    private
    public  :: Run_adjoint_gpt
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_adjoint_gpt ()
        
        logical  :: storage_eigen                                               ! keep eigenvlaue type before adjoint iteration
    
        ! ----------------------------------------------------------------------
        ! perform some prepare before adjoint solve
        call quad%negate ()
        call sweep%set (mesh, geom, bound, quad)
        storage_eigen = ns%flag%is_eigen
        
        ns%flag%is_eigen = .FALSE.
        call Get_gpt_source ()
        
        ! ----------------------------------------------------------------------
        ! prepare xsec used in iteration
        call Transit_xsec_gpt (iter_count_unpert%eigenvalue)
        
        ! initialize iteration parameter
        call Init_iteration_variable (is_adjoint=.TRUE.)
                
        ! perform interation
        call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.TRUE.)
        
        ! ----------------------------------------------------------------------
        ! reset changing after adjoint solve
        call quad%negate ()
        call sweep%set (mesh, geom, bound, quad)
        ns%flag%is_eigen = storage_eigen
        
        call Print_binary_hdf5 (is_adjoint=.TRUE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        call Print_vtk_files (is_adjoint=.TRUE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        
    end subroutine Run_adjoint_gpt
    
end module driver_adjoint_gpt
