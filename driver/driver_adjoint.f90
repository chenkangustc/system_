!$
!===================================================================================================
!
!   control subroutine for adjoint calculation at initial condition, with various type define
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_adjoint
!
!   Public type lists:          No
!
!===================================================================================================
module driver_adjoint
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    
    use transit_to_solver,          only : Transit_xsec_adjoint
    use iteration_initialize,       only : Init_iteration_variable
    use iteration_control,          only : Driving_iteration
    use output_hdf5,                only : Print_binary_hdf5
    use output_visit,               only : Print_vtk_files
    
    implicit none
    private
    public  :: Run_adjoint
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_adjoint ()
        
        logical  :: storage_eigen                                               ! keep eigenvlaue type before adjoint iteration
        integer  :: ig 
        
        ! ----------------------------------------------------------------------
        ! perform some prepare before adjoint solve
        call quad%negate ()
        call sweep%set (mesh, geom, bound, quad)
        storage_eigen = ns%flag%is_eigen
        
        select case(nt%flag%adjoint_type)
        case('HOMOGENEOUS', 'homogeneous')
            ns%flag%is_eigen = .TRUE.
            call iter_adjoint%homo ()
            
        case('ONE', 'one')
            ns%flag%is_eigen = .FALSE.
            call iter_adjoint%one ()
            
        case('SIGMA_F_NU', 'sigma_f_nu')
            ns%flag%is_eigen = .FALSE.
            call iter_adjoint%neutron (xsec)
            
        case('SIGMA_F_KAPPA', 'sigma_f_kappa')
            ns%flag%is_eigen = .FALSE.
            call iter_adjoint%energy (xsec)
            
        case('ACTUAL', 'actual')
            ns%flag%is_eigen = .FALSE.
            call iter_adjoint%actual (Q_ext)
        end select
        
        ! ----------------------------------------------------------------------
        ! prepare xsec used in iteration
        call Transit_xsec_adjoint (is_eigen=ns%flag%is_eigen)
        call Init_iteration_variable (is_adjoint=.TRUE.)
        call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.TRUE.)
        
        ! ----------------------------------------------------------------------
        ! reset changing after adjoint solve
        call quad%negate ()
        call sweep%set (mesh, geom, bound, quad)
        ns%flag%is_eigen = storage_eigen
        
        ! change direction
        call flux_adjoint%fix ()
!        do ig = 1, SIZE(flux_adjoint%ngs, dim=1)
!            flux_adjoint%ngs(ig)%scalar = 1.0 
!            flux_adjoint%ngs(ig)%angular = 1.0 
!        end do 
        
        call Print_binary_hdf5 (is_adjoint=.TRUE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        call Print_vtk_files (is_adjoint=.TRUE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        
    end subroutine Run_adjoint
    
end module driver_adjoint
