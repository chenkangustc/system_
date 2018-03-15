!$
!===================================================================================================
!
!   control subroutine for transient execution (pcqs method), advance for every time step
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_transient_pcqs
!
!   Public type lists:          No
!
!===================================================================================================
module driver_transient_pcqs
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use driver_adjoint,             only : Run_adjoint
    use transit_to_solver,          only : Transit_xsec_pcqs
    use transit_from_solver,        only : Transit_precursor
    use iteration_initialize,       only : Init_iteration_variable
    use iteration_control,          only : Driving_iteration   
    use feedback,                   only : Check_feedback_transient, Check_xsec_transient

    use time_advancing,             only : Get_initial_precursor, Update_precursor
    use perturbation,               only : Perform_perturbation
    use coefficient_pk_header
    use process_pcqs
    use process_pk,                 only : Pk_solver_RBDF, Pk_solver_VODE, Print_step_info
    
    use output_timelist,            only : Print_timelist
    use output_visit,               only : Print_vtk_files
    use output_hdf5,                only : Print_binary_hdf5
    
    use timestep_header,            only : TimeStepInfo
    
    implicit none
    private
    public  :: Run_transient_pcqs
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_transient_pcqs ()
        
        type(TimeStepInfo)  :: step_Macro, step_Medial
        integer             :: i_Macro, i_Medial
        integer             :: id_Macro, id_Medial
        real(KREAL)         :: ctime, ltime
        integer             :: id_Total
        logical             :: is_pass
        
        ! ----------------------------------------------------------------------
        ! obtain initial flux and precursor concentration      
        step_Macro = time_step%info (0, is_step=.FALSE.)
        step_Medial = time_step%info (0, is_step=.TRUE.)
        
        call dnp_solver%init (mesh, geom, xsec, param, flux_forward, dist_dnps, timelist)
!        call Get_initial_precursor ()
!        call Transit_precursor ()
        
        ! obtain initial shape function
        call Normal_adjoint_flux ()
        call Get_initial_flux_shape ()
            
        ! generate point kinetics parameters
        call Generate_pk_parameter (0)
        call Generate_reactivity (0)
        
        ! initial condition for pk equation
        call Get_initial_amplitude ()
        call Get_initial_precursor_shape ()
        
        pk_param%neutron = 1.0D0
        pk_param%precursor = pk_param%partial_beta / pk_param%partial_lambda / pk_param%generation_time 
        
        ! ----------------------------------------------------------------------
        ! cycle for macro time step
        id_Total = 0
        id_Macro = 0
        macro: do i_Macro = 1, time_step%get_section_count ()
            id_Macro = id_Macro + 1
            step_Macro = time_step%info (id_Macro, is_step=.FALSE.)    
            ctime = step_Macro%right
            
            ! get cross section for current state, mainly for external source
            call Perform_perturbation (id_Total, step_Macro%right)
            call Check_xsec_transient ()
            
            ! predict flux per macro time step
            call Init_iteration_variable (is_adjoint=.FALSE.)
            call Transit_xsec_pcqs (step_Macro%pace)
            call Driving_iteration (is_eigen=.FALSE., is_adjoint=.FALSE. , is_transient=.TRUE., tidx=id_Macro, ctime=step_Macro%right)
            call Get_flux_shape (id_Macro)
            
            ! cycle for micro time step
            ! the point of index=1 need not shape interpolation
            id_Medial = 0
            medial: do i_Medial = 1, time_step%step_per_section(i_Macro)
                id_Medial = id_Medial + 1
                id_Total = id_Total + 1
                step_Medial = time_step%info (id_Total, is_step=.TRUE.)  
                
                ctime = step_Medial%right
                ltime = ctime - step_Medial%pace
                
                ! get cross section for current state
                call Perform_perturbation (id_Total, ctime)
                call Check_xsec_transient ()
                
                ! generate reactivity per time step                
                call Generate_pk_parameter (id_Total)
                call Generate_reactivity (id_Total)
                call pk_solver%advance (pk_param, left=ltime, right=ctime)
                
                amplitude%flux = pk_param%neutron
                amplitude%precursor = pk_param%precursor
                
                ! the last point is the macro time point, so the information output should be in the macro time cycle
                if (i_Medial == time_step%step_per_section(i_Macro))  then
                    call Transit_shape_function (id_Macro)
                else 
                    call Shape_interpolation (step_Macro, step_Medial)
                end if
                
                ! update precusor shape
!                call Get_precursor_shape (step_Medial)
                
                ! generation power level
!                call Get_correct_flux (id_Total)
!                call Get_correct_precursor (id_Total)
                call Regenerate_power (id_Total)
                call dnp_solver%advance (mesh, xsec, param, flux_forward, dist_dnps, step_Medial) 
                
                ! perform feedback
                call Check_feedback_transient (is_pass, step_Medial%left, step_Medial%right)
                call Print_transient_pcqs_timelist (id_Total, ctime)
            end do medial
                
            call Print_transient_pcqs_dist (id_Total, ctime)
        end do macro
        
        call pk_param%clean ()
        call pk_solver%clean ()
        
    end subroutine Run_transient_pcqs
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_transient_pcqs_timelist (id_Total, ctime)
        
        integer, intent(in)  :: id_Total
        real(KREAL), intent(in)  :: ctime
        
        call timelist%set (geom, mesh, dist_power, self_fdbk, pk_parameter%rho, pk_parameter%beta)
        call Print_timelist (timelist, id_Total, ctime, FILES%TIMELIST)
        call det%get (ns, geom, mesh, flux_forward)
        call det%pvalue (FILES%DET, id_Total, ctime)
        
    end subroutine Print_transient_pcqs_timelist
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_transient_pcqs_dist (id_Total, ctime)
        
        integer, intent(in)  :: id_Total
        real(KREAL), intent(in)  :: ctime
        
        call Print_binary_hdf5 (is_adjoint=.FALSE., is_transient=.TRUE., tidx=id_Total, ctime=ctime)
        call Print_vtk_files (is_adjoint=.FALSE., is_transient=.TRUE., tidx=id_Total, ctime=ctime)
    
    end subroutine Print_transient_pcqs_dist
    
end module driver_transient_pcqs
