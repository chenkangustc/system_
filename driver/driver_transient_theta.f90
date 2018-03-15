!$
!===================================================================================================
!
!   control subroutine for transient execution (theta method), advance for every time step
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_transient_theta
!
!   Public type lists:          No
!
!===================================================================================================
module driver_transient_theta

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use transit_to_solver,          only : Transit_xsec_theta
    use transit_from_solver,        only : Transit_precursor
    use iteration_initialize,       only : Init_iteration_variable
    use iteration_control,          only : Driving_iteration
    use feedback,                   only : Check_feedback_transient, Check_xsec_transient

    use time_advancing,             only : Get_initial_precursor, Update_precursor
    use perturbation,               only : Perform_perturbation
    use process_pcqs
    
    use output_timelist,            only : Print_timelist
    use output_hdf5,                only : Print_binary_hdf5
    use output_visit,               only : Print_vtk_files
    
    use timestep_header,            only : TimeStepInfo
    
    implicit none
    private
    public  :: Run_transient_theta
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_transient_theta ()
    
        type(TimeStepInfo)  :: step_Macro, step_Medial
        integer             :: i_Macro, i_Medial
        integer             :: id_Macro, id_Medial
        integer             :: id_Total                                         ! time step counter for the whole transient precess
        logical             :: is_pass
        real(KREAL)         :: ctime
        
        ! ----------------------------------------------------------------------
        ! obtain initial flux and precursor concentration
        step_Macro = time_step%info (0, is_step=.FALSE.)
        step_Medial = time_step%info (0, is_step=.TRUE.)
        
        call dnp_solver%init (mesh, geom, xsec, param, flux_forward, dist_dnps, timelist)
!        call Get_initial_precursor ()
!        call Transit_precursor ()

        call Normal_adjoint_flux ()
        
        call Generate_pk_parameter (0)
        call Generate_reactivity (0)
        
        ! ----------------------------------------------------------------------
        ! cycle for macro time step
        id_Total = 0
        id_Macro = 0
        macro: do i_Macro = 1, time_step%get_section_count ()
            id_Macro = id_Macro + 1
            step_Macro = time_step%info (id_Macro, is_step=.FALSE.)
        
            ! ------------------------------------------------------------------
            ! cycle for micro time step
            id_Medial = 0
            medial: do i_Medial = 1, time_step%step_per_section(i_Macro)
                id_Medial = id_Medial + 1
                id_Total = id_Total + 1
                step_Medial = time_step%info (id_Total, is_step=.TRUE.)     
        
!                tidx = tidx + 1
                ctime = step_Medial%right
                
                call Perform_perturbation (id_Total, ctime)
                call Check_xsec_transient ()
                
                nested: do 
                    call Transit_xsec_theta (step_Medial%pace)
                    
!                    if (id_Total == 1)  then
!                        call xsec%print (301)
!                        call xsec_iter%print (302)
!                        call Q_ext%print (303)
!                        stop(0)
!                    end if 
                    
                    call Init_iteration_variable (is_adjoint=.FALSE.)
                    call Driving_iteration (is_eigen=.FALSE., is_adjoint=.FALSE. , is_transient=.TRUE., tidx=id_Total, ctime=ctime)
                    
                    if (.NOT. ns%feedback%is_feedback)  then
                        exit nested
                    end if
                    
                    ! treatment feedback
                    call Check_feedback_transient (is_pass, step_Medial%left, step_Medial%right)
                    
                    if (is_pass)  then
                        exit nested
                    end if
                    
                end do nested
                
                ! ------------------------------------------------------------------
                ! post transient treatment
                ! calculate the map of precursor concentration per delayed groups
!                call Update_precursor (step_Medial)
!                call Transit_precursor () 
                call dnp_solver%advance (mesh, xsec, param, flux_forward, dist_dnps, step_Medial) 
                
                call Print_transient_theta_timelist (id_Total, ctime)
            end do medial
            
            call Print_transient_theta_dist (id_Total, ctime)
        end do macro
            
    end subroutine Run_transient_theta
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_transient_theta_timelist (tidx, ctime)
        
        integer, intent(in)  :: tidx
        real(KREAL), intent(in)  :: ctime
        
        call timelist%set (geom, mesh, dist_power, self_fdbk, pk_parameter%rho, pk_parameter%beta)
        call Print_timelist (timelist, tidx, ctime, FILES%TIMELIST)
        call det%get (ns, geom, mesh, flux_forward)
        call det%pvalue (FILES%DET, tidx, ctime)
        
    end subroutine Print_transient_theta_timelist
        
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_transient_theta_dist (tidx, ctime)
        
        integer, intent(in)  :: tidx
        real(KREAL), intent(in)  :: ctime
        
        call Print_binary_hdf5 (is_adjoint=.FALSE., is_transient=.TRUE., tidx=tidx, ctime=ctime)
        call Print_vtk_files (is_adjoint=.FALSE., is_transient=.TRUE., tidx=tidx, ctime=ctime)
    
    end subroutine Print_transient_theta_dist
        
end module driver_transient_theta
