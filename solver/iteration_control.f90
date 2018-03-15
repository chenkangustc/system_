!$
!===================================================================================================
!
!   whole iteration control of this transport solver
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_iteration
!
!   Public type lists:          No
!
!===================================================================================================
module iteration_control

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use iteration_header,           only : IterationCriterion
    
    use global
    use coefficient_iteration
    
    use transit_from_solver,        only : Transit_flux
    use iteration_process
    use iteration_process_adjoint
    use output_main

    implicit none 
    private
    public  :: Driving_iteration
    
    ! 
    type(IterationCriterion), pointer  :: criteria => NULL()

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_iteration (is_eigen, is_adjoint, is_transient, tidx, ctime)
        
        ! intent parameters
        logical, intent(in)  :: is_eigen
        logical, intent(in)  :: is_adjoint
        logical, intent(in), optional      :: is_transient
        integer, intent(in), optional      :: tidx
        real(KREAL), intent(in), optional  :: ctime
        
        if (ns%output%is_log)  then
            call Print_iteration_title(FILES%MAIN, is_eigen, is_adjoint, is_transient, tidx, ctime)
            call Print_iteration_title(OUTPUT_UNIT, is_eigen, is_adjoint, is_transient, tidx, ctime)
        end if 
        
        ! perform iteration
        call iter_init%import (iter_q, iter_flux, iter_count, flux_scat)
        call Process_iteration (is_eigen, is_adjoint)
        call iter_init%export (iter_q, iter_flux, iter_count, flux_scat)
        
        ! multiply normal factor after iteration only for FSP
        if (.NOT. is_eigen .and. ns%method%is_Ks  .and. xsec_iter%is_fission () )  then
                iter_flux%info%moment(0, :, :) = iter_flux%info%moment(0, :, :) * Q_ext%iter_normal * (1.0D0/iter_count%coeff_FSP)
                iter_flux%dist%nodal = iter_flux%dist%nodal * Q_ext%iter_normal * (1.0D0/iter_count%coeff_FSP)
        end if
        
        ! transit from solver
        call Transit_flux (is_eigen, is_adjoint, is_transient, tidx, ctime)
        
        ! print iteration informtaion
        if (ns%output%is_log)  then
            call time_program%update (FILES%MAIN)
            call Print_iteration_summary (FILES%MAIN, is_eigen=is_eigen, is_adjoint=is_adjoint)
        end if 
        
        ! reset back for extrapolation
        if (.NOT. is_eigen .and. ns%method%is_Ks  .and. xsec_iter%is_fission () )  then
                iter_flux%info%moment(0, :, :) = iter_flux%info%moment(0, :, :) / (Q_ext%iter_normal * (1.0D0/iter_count%coeff_FSP))
                iter_flux%dist%nodal = iter_flux%dist%nodal / (Q_ext%iter_normal * (1.0D0/iter_count%coeff_FSP))
        end if
        
    end subroutine Driving_iteration
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Process_iteration (is_eigen, is_adjoint)
        
        logical, intent(in)  :: is_eigen
        logical, intent(in)  :: is_adjoint
        
        integer  :: ig, iig, ia, ir, il, id, ml, iz, i_cycle
        
        logical, parameter  :: is_rubost   = .FALSE.                            ! ks update, more rubost or faster
        logical, parameter  :: is_scaling  = .FALSE.
        logical             :: is_passing
        
        real(KREAL)  :: rho                                                     ! extrapolation factor for LW
        real(KREAL)  :: error_inner                                             ! error for iteration parameter
        real(KREAL)  :: error_flux
        real(KREAL)  :: error_fission_rate
        real(KREAL)  :: error_eigen
        real(KREAL)  :: eigenvalue_new
        
        real(KREAL), save  :: fq_old, fq_new
        integer, save      :: n_iter = 0                                        ! number of transport solver evoked, Note: save
        
        if (.NOT. is_eigen .and. ns%method%is_Ks  .and. xsec_iter%is_fission ())  then
            n_iter = n_iter + 1
        end if
        call Select_criteria (is_eigen, i_cycle=1)
        
        if (ns%state%scat_order /= 0) then 
            call coeff_source%set (quad)
        end if
        
        if (.FALSE.)  then
            do il = 1, SIZE(quad%is_symmetry)
                write(102, fmt='(1x, I3, A, I3)') il, '-->', quad%is_symmetry(il)
            end do 
            do il = 1, SIZE(quad%is_order)
                write(103, fmt='(1x, I3, A, I3)') il, '-->', quad%is_order(il)
            end do 
            
            do il = 1, SIZE(quad%directions)
                write(104, fmt='(1x, I3, *(ES12.5, TR3))') il, quad%directions(il)%wmu, quad%directions(il)%xmu
            end do 
        end if 
        
        ! @adjoint
        if (is_adjoint)  then
            call Update_adjoint_fission ()
        end if
        
        ! initial guess for eigenvalue
        if (n_iter == 1)  then
            fq_old = 0.0D0
            fq_new = 0.0D0
        end if
        
        ! fq_new as initial condition 
        if (.NOT. is_eigen .and. ns%method%is_Ks .and. xsec_iter%is_fission ())  then
            if (.NOT. is_adjoint) then 
                call Get_forward_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            ! @adjoint
            else 
                call Get_adjoint_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            end if 
        end if
        
        ! @adjoint
        if (is_adjoint)  then
            call iter_count%set_downscatter (xsec_iter)
        else
            call iter_count%set_upscatter (xsec_iter)
        end if
        
        ! begin outer (fission source) iterations
        iter_count%out = 0
    10  iter_count%out = iter_count%out + 1
        
        if (.NOT. is_eigen .and. ns%method%is_Ks  .and. xsec_iter%is_fission () )  then
            iter_count%coeff_FSP = fq_new * (1.0D0-iter_count%ks) / iter_count%ks                ! NOTE:
        else
            iter_count%coeff_FSP = 1.0D0
        end if
        
        ! ----------------------------------------------------------------------
        ! energy group iterations
        upscatter: do i_cycle = 1, SIZE(iter_count%ng_start)
        call Select_criteria (is_eigen, i_cycle)
        energy_group: do ig = iter_count%ng_start(i_cycle), iter_count%ng_end(i_cycle), iter_count%ng_step(i_cycle)
            
            call time_event%begin(2)
            call coeff_nodal%set (ig, geom, mesh, xsec_iter, quad)
            call coeff_surface%set (ig, geom, mesh, xsec_iter, quad)
            call time_event%end(2)

            ! generate source moments, except in-group source
            call Generate_source_moments (ig, is_adjoint)
            
            ! ------------------------------------------------------------------
            ! begin inner iteration
            iter_count%in = 0
    20      iter_count%in = iter_count%in + 1
        
            if (ns%method%is_LW) then
                call accele_LW%LW (iter_flux, iter_count%in, ig)
            end if
        
            ! add in-group source
            call Add_in_group_source (ig)
            
            ! transfer to old iteration flux
            call iter_flux%transit (ig)
            
            ! sweep space-angle mesh 
            call time_event%begin(3)
            call Get_average_flux (ig, is_adjoint=is_adjoint)
            call time_event%end(3)
            
            ! maxial rms for flux for all nodals
            call iter_flux%cycle_in (criteria, ig, error_inner)
                        
            if (ns%method%is_LW ) then
                call accele_LW%set_LW (iter_flux, iter_count%in, ig)
            end if
        
            ! end of inner iteration when eigenvalue problem
            call criteria%check_inner (iter_count%in, error_inner, is_passing)
            if (is_eigen .AND. .NOT. is_passing)  then
                goto 20
            end if
            
            ! perform high-order flux moments calculation
            ! For conventional FSP, generate moments during the inner iteration
            ! For FSP with fission material, generate moments out the innera iteration
            if (iter_count%out == 1) then
                call time_event%begin(4)
                call Re_generate_flux_moments (ig)
                call time_event%end(4)
            end if
            
            ! cycle inner iteration when FSP problem, add in-group-source to total source
            if (.NOT. is_eigen .AND. .NOT. is_passing)  then
                call Add_in_group_source_moments (ig)
                goto 20
            end if
            
            ! perform high-order flux moments calculation
            if (iter_count%out > 1) then
                call time_event%begin(4)
                call Re_generate_flux_moments (ig)
                call time_event%end(4)
            end if
            
            ! calculate source moment by flux moments
            if (.NOT. is_adjoint)  then
                call Update_source_moment ()
            else 
                call Update_adjoint_moment ()
            end if
            
        end do energy_group
        end do upscatter
        
        ! ----------------------------------------------------------------------
        ! cycle outer iteration when FSP problem
        iter_exit: if (.NOT. is_eigen) then

            if (.NOT. xsec_iter%is_fission ()) then
                return
            end if
            
            if (.NOT. is_adjoint) then 
                call Get_forward_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            ! @adjoint
            else 
                call Get_adjoint_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            end if 
            
            ! scaling method for acceleration
            if (.NOT. ns%method%is_Ks .AND.  is_scaling )  then
                call Q_ext%strength (geom, accele_Scaling%source)
                accele_Scaling%factor = accele_Scaling%source / (accele_Scaling%source + fq_old - fq_new)
                iter_flux%info%moment(0,:,:) = iter_flux%info%moment(0,:,:) * accele_Scaling%factor
                iter_q%fission%moment(0,:) = iter_q%fission%moment(0,:) * accele_Scaling%factor
                iter_q%fission%old(:) = iter_q%fission%old(:) * accele_Scaling%factor
            end if
            
            if (ns%method%is_Ks )  then
                if (.NOT. is_rubost)  then
                    eigenvalue_new = iter_count%ks * fq_new/fq_old
                else 
                    eigenvalue_new = iter_count%ks * (fq_new + fq_old) / (2.0*fq_old)
                end if
                
                error_eigen = ABS(eigenvalue_new-iter_count%ks)/eigenvalue_new
                iter_count%ks = eigenvalue_new
            else 
                error_eigen = 0.0
                iter_count%ks = 1.0
            end if
            
            if (ns%output%is_log)  then
                call Print_iteration_step_FSP (OUTPUT_UNIT, iter_count, error_eigen, error_flux, error_inner)
                call Print_iteration_step_FSP (FILES%MAIN, iter_count, error_eigen, error_flux, error_inner)
            end if 
            
            ! end cycle when fission material with FSP
            call criteria%check_outer (iter_count%out, error_eigen, error_flux, error_fission_rate, is_passing)
            if (.NOT. is_passing)  then
                goto 10
            end if
            
        ! cycle outer iteration when eigenvalue problem
        else iter_exit
            if (.NOT. is_adjoint) then 
                call Get_forward_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            ! @adjoint cycle critieration    
            else
                call Get_adjoint_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
            end if 
            
            ! for eigenvalue update
            eigenvalue_new = iter_count%eigenvalue * fq_new/fq_old
            error_eigen = ABS(eigenvalue_new-iter_count%eigenvalue)/eigenvalue_new
            
            iter_count%eigenvalue = eigenvalue_new
            
            if (ns%output%is_log)  then
                call Print_iteration_step_eigen (OUTPUT_UNIT, iter_count, error_eigen, error_flux, error_inner)
                call Print_iteration_step_eigen (FILES%MAIN, iter_count, error_eigen, error_flux, error_inner)
            end if 
            
            ! end of outer iteration
            call criteria%check_outer (iter_count%out, error_eigen, error_flux, error_fission_rate, is_passing)
            if (.NOT. is_passing)  then
                goto 10
            end if

        end if iter_exit
        
    end subroutine Process_iteration
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Select_criteria (is_eigen, i_cycle)
        
        logical, intent(in)  :: is_eigen
        integer, intent(in)  :: i_cycle
        
        type(IterationCriterion), target  :: criteria_tmp
            
        if (i_cycle == 1)  then
            if (is_eigen )  then
                criteria => criteria_eigen
            else
                criteria => criteria_fsp
                if (.NOT. xsec_iter%is_fission ())  then
                    criteria_tmp = criteria_fsp
                    criteria_tmp%max_inner = criteria_tmp%max_inner * 100
                    criteria => criteria_tmp
                end if
            end if
        else
            if (is_eigen )  then
                criteria_tmp = criteria_eigen
                criteria_tmp%max_inner = MAX(4, criteria_tmp%max_inner/4)
                criteria => criteria_tmp
            else
                criteria_tmp = criteria_fsp
                criteria_tmp%max_inner = MAX(4, criteria_tmp%max_inner/4)
                criteria => criteria_tmp
            end if
        end if
    
    end subroutine Select_criteria

end module iteration_control
