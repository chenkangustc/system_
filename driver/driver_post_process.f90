!$
!===================================================================================================
!
!   Post process after the calculation completed: print tail, close file and so on.
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_post_process
!
!   Public type lists:          No
!
!===================================================================================================
module driver_post_process
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use coefficient_iteration
    
    use files_dirs
    use output_main,            only : Print_tail_main
	
	use imp_driving_post_process
    
    implicit none
    private
    public   :: Run_post_process

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_post_process ()
        
        ! --------------------------------------------------------------------------
        call time_event%print (FILES%MAIN)
        
        call Print_tail_main(FILES%MAIN, TRIM(ns%flag%case_title))
        call Print_tail_main(OUTPUT_UNIT, TRIM(ns%flag%case_title))
        !--------------------------------------------------------------------------
		call Run_output()
		call Free_imp_thermal()
        ! --------------------------------------------------------------------------
        call Free_steady_memory ()
        
        if (nt%flag%is_perturb )  then
            call Free_perturb_memory ()
        end if
        
        if (nt%flag%is_transient )  then
            call Free_transient_memory ()
        end if
        
        ! delete the tmp file and close the input and output files
        call Close_file (FILES%CASENAME, is_keep=.FALSE.)
        
        ! output file
        call Close_file (FILES%MAIN, is_keep=.TRUE.)
        call Close_file (FILES%MEMORY, is_keep=.TRUE.)
        call Close_file (FILES%TIMELIST, is_keep=.TRUE.)
        call Close_file (FILES%DET, is_keep=.TRUE.)
        
        call Close_file (FILES%PT, is_keep=nt%flag%is_perturb)
        call Close_file (FILES%REACTIVITY, is_keep=nt%flag%is_transient)
        
    end subroutine Run_post_process
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! deallocate memory used in steady calculation
    !===============================================================================================
    subroutine Free_steady_memory ()
    
        integer  :: i
        
        ! for geometry
        call mesh%clean ()
        call geom%clean ()
        call bound%clean ()
        call mesh_vtk%clean ()
        
        ! for material
        call mat_info%clean ()
        
        call xsec_inp%clean ()
        call xsec_init%clean ()
        call xsec_iter%clean ()
        call xsec%clean ()
            
        call param_inp%clean ()
        call param%clean ()
        
        call Q_ext%clean ()
        
        ! for quadrature set & sweeping
        call quad%clean ()
        call sweep%clean ()
        
        ! for iteration parameter
        call iter_q%clean ()
        call iter_flux%clean ()
        call iter_flux%dist%clean ()
        call iter_count%clean ()
        
        call flux_scat%clean ()
        
        call iter_adjoint%clean ()
        
        ! for result container
        call timelist%clean ()
        
        call flux_forward%clean ()
        call flux_adjoint%clean ()
        
        call dist_power%clean ()
        call dist_flux%clean ()
        call dist_fission_rate%clean ()
        
        if (allocated(dist_dnps))  then
            do i = 1, SIZE(dist_dnps)
                call dist_dnps(i)%clean ()
            end do
        end if
        
        ! for control rod information
        call cr_bank%clean ()
        
        ! for perturbation
        call pert_xsec%clean ()
        call pert_q%clean ()
        call pert_cr%clean ()
        call pert_mat%clean ()
        call pert_case%clean ()
        call pert_th%clean ()
        
        ! for iteration
        call coeff_surface%clean ()
        call coeff_nodal%clean ()
        call coeff_source%clean ()
        
        ! for SOR
        call accele_SOR%clean ()
        
        ! for feedback
        call self_fdbk%clean ()
        call self_lra%clean ()
        call det%clean ()
        
    end subroutine Free_steady_memory
    
    !$
    !===============================================================================================
    ! deallocate memory used in perturbation calculation
    !===============================================================================================
    subroutine Free_perturb_memory ()
        
        call xsec_unpert%clean ()
        
        call flux_forward_unpert%clean ()
        
        call flux_adjoint_pt%clean ()
        call flux_adjoint_gpt%clean ()
    
    end subroutine Free_perturb_memory
    
    !$
    !===============================================================================================
    ! deallocate memory used in transient calculation
    !===============================================================================================
    subroutine Free_transient_memory ()

        ! for time step
        call time_step%clean ()
        
        ! for transient process
        call shape_last%clean ()
        call shape_current%clean ()
        call shape_predict%clean ()

        call amplitude%clean ()
        call pk_parameter%clean ()
        
        call dnp_solver%clean ()
        call pk_param%clean ()
        call pk_solver%clean ()
        
    end subroutine Free_transient_memory
    
end module driver_post_process
