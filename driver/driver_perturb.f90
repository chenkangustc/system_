!$
!===================================================================================================
!
!   control subroutine for perturbation execution
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_perturb
!
!   Public type lists:          No
!
!===================================================================================================
module driver_perturb

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use driver_adjoint_gpt,         only : Run_adjoint_gpt
    use driver_steady,              only : Run_steady
    use transit_to_solver,          only : Transit_xsec_steady
    use reactivity,                 only : Static_reactivity
    use perturbation,               only : Move_CRWorth
    
    implicit none
    private
    public  :: Run_perturb
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_perturb ()
        
        integer  :: ipert, istep
        
        ! save inital xsec & flux
        call Transit_perturbation ()
        
        if ((worth_cr%pos_type == worth_cr%BY_INP) .OR. (worth_cr%pos_type == worth_cr%BY_STEP))  then
            do ipert = 1, worth_cr%getNstep()
                call Move_CRWorth (ipert, istep)
                call Run_steady ()
                
                write(OUTPUT_UNIT, fmt='(1x, A, I4, A, I4, A, F9.6)')  &
                    &   'Bank=', worth_cr%bank, ';  Step=', istep, ';  Keff=', iter_count%eigenvalue
                write(110, fmt='(1x, A, I4, A, I4, A, F9.6)')  &
                    &   'Bank=', worth_cr%bank, ';  Step=', istep, ';  Keff=', iter_count%eigenvalue
            end do 
            
        else
            call pert_case%fix_conf (mat_info)
            do ipert = 1, pert_case%n_pert
                call pert_case%get_xsec (ipert, xsec)
                ! perform gpt once only
                if (ipert == 1)  then
                    call Run_adjoint_gpt ()
                    flux_adjoint_gpt = flux_adjoint
                end if
                
                call Run_steady ()
                call Static_reactivity (FILES%PT, ipert)
            end do
        end if 
        
    end subroutine Run_perturb
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_perturbation ()
        
        integer  :: i, j
        
        ! xsec
        do i = 1, SIZE(xsec%matrixs, dim=1)
            do j = 1, SIZE(xsec%matrixs, dim=2)
                xsec_unpert%matrixs(i, j) = xsec%matrixs(i, j)
            end do
        end do
        
        ! counter
        iter_count_unpert = iter_count
        
        ! flux
        flux_forward_unpert = flux_forward
        flux_adjoint_pt = flux_adjoint
    
    end subroutine Transit_perturbation
        
end module driver_perturb
