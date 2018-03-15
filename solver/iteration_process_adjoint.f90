!$
!===================================================================================================
!
!   module for iteration calculation kernel for adjoint flux
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Update_adjoint_moment
!                               Update_adjoint_fission
!                               Get_adjoint_error
!
!   Public type lists:          No
!
!===================================================================================================
module iteration_process_adjoint
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use coefficient_iteration
    
    use iteration_header,               only : IterationCriterion
    use vector_operation
    
    implicit none
    private
    public  :: Update_adjoint_fission, Update_adjoint_moment, Get_adjoint_error

contains    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_adjoint_fission ()
    
        integer  :: ia, ir, iz, il, ig
    
        iter_adjoint%sigma_f = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = (ia-1)*ns%state%nodal + ir
                do ig = 1, ns%state%ng
                    iter_adjoint%sigma_f(il) = iter_adjoint%sigma_f(il) + xsec_iter%matrixs(iz, ia)%sigma_f_nu(ig)
                end do 
            end do
        end do
        
    end subroutine Update_adjoint_fission

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_adjoint_moment ()
        
        integer  :: ia, ml, ir, il, iz, id, iig
        
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                il = ml + ir
                iz = mesh%zone(ir)
                do id = 0, 8
                    iter_adjoint%q_moments(id,il) = 0.0
                    do iig = 1, ns%state%ng
                        iter_adjoint%q_moments(id,il) = iter_adjoint%q_moments(id,il) + xsec_iter%matrixs(iz,ia)%chi_steady(iig) * iter_flux%info%moment(id,il,iig)
                    end do
                end do
            end do
        end do 
        
    end subroutine Update_adjoint_moment
    
    !$
    !===============================================================================================
    !
    !===============================================================================================    
    subroutine Get_adjoint_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
        
        type(IterationCriterion), intent(in)  :: criteria
        real(KREAL), intent(out)  :: error_fission_rate
        real(KREAL), intent(out)  :: error_flux
        real(KREAL), intent(out)  :: fq_new
        real(KREAL), intent(out)  :: fq_old
        
        real(KREAL)  :: old_rate(ns%state%layer * ns%state%nodal)
        real(KREAL)  :: new_rate(ns%state%layer * ns%state%nodal)
        real(KREAL)  :: division
        integer  :: ia, ml, ir, il
        
        fq_new = 0.0
        fq_old = 0.0
        new_rate = 0.0
        old_rate = 0.0
        
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                il = ml + ir
                new_rate(il) = iter_adjoint%q_moments(0,il)*iter_adjoint%sigma_f(il)*geom%area(ir)*geom%height(ia)
                old_rate(il) = iter_adjoint%q_old(il)*iter_adjoint%sigma_f(il)*geom%area(ir)*geom%height(ia)
                fq_new = fq_new + new_rate(il)
                fq_old = fq_old + old_rate(il)
            end do
        end do
        
        error_flux = get_vector_error (iter_adjoint%q_old, iter_adjoint%q_moments(0,:), criteria%error_type)
        error_fission_rate = get_vector_error (old_rate, new_rate, criteria%error_type)
        
        do il = 1, ns%deduce%nodal_total
            iter_adjoint%q_old(il) = iter_adjoint%q_moments(0,il)
        end do
    
    end subroutine Get_adjoint_error

end module iteration_process_adjoint
