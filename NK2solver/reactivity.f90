!$
!===================================================================================================
!
!   module for reactivity calculation by perturbation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Static_reactivity
!
!   Public type lists:          No
!
!===================================================================================================
module reactivity
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use contain_header,             only : GroupsFlux
    
    implicit none 
    private
    public  :: Static_reactivity, Get_gpt_source
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Static_reactivity (file_out, i_pert)
        
        integer, intent(in)  :: file_out
        integer, intent(in)  :: i_pert
        
        ! local variables
        real(KREAL)  :: rho_ept                                             ! reactivity eastimated by exact perturbation
        real(KREAL)  :: rho_dp                                              !                       by delta keff
        real(KREAL)  :: rho_pt                                              !                       by first-order perturbation
        real(KREAL)  :: rho_gpt                                             !                       by gpt based second order perturbation
        
        ! by exact perturbation
        call Reactivity_by_EPT (rho_ept)
        
        ! by delta keff
        rho_dp = (iter_count%eigenvalue - iter_count_unpert%eigenvalue) / (iter_count%eigenvalue * iter_count_unpert%eigenvalue)
        
        ! by first-order
        call Reactivity_by_PT (rho_pt)
        
        ! by second-order
        call Reactivity_by_GPT (rho_pt, rho_gpt)
        
        ! ----------------------------------------------------------------------
        ! result output
        write(file_out, "(1x, A)")          '_______________________________________________________________________________'
        write(file_out, "(1x, A, I3)")      'Perturbation index  is     :', i_pert
        write(file_out, "(1x, A)")          'reactivity estimation      :'
        
        write(file_out, "(1x, A, ES13.6)")  'current keff        : ', iter_count%eigenvalue
        write(file_out, "(1x, A, ES13.6)")  '    by exact pt     : ', rho_ept
        write(file_out, "(1x, A, ES13.6)")  '    by keff         : ', rho_dp
        write(file_out, "(1x, A, ES13.6)")  '    by first-order  : ', rho_pt
        write(file_out, "(1x, A, ES13.6)")  '    by second-order : ', rho_gpt
        write(file_out, "(1x, A)")          ' '
        
    end subroutine Static_reactivity
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_gpt_source ()
    
        integer  :: ig, iig, ir, ia, im, iz
        
        ! local variables
        real(KREAL)  :: denominator
        real(KREAL)  :: rho_t, rho_s, rho_f_nu, rho
        real(KREAL)  :: tmp_1, tmp_2
        
        ! get integration
        call Get_integration (flux_forward_unpert, flux_adjoint_pt, denominator, rho_t, rho_s, rho_f_nu)
        rho = rho_t - rho_s - rho_f_nu/iter_count_unpert%eigenvalue
        
        ! get source
        iter_adjoint%source = 0.0
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    rho_t = 0.0
                    rho_s = 0.0
                    rho_f_nu = 0.0
                    
                    rho_t = flux_adjoint_pt%ngs(ig)%scalar(ir,ia) * (xsec%matrixs(iz,ia)%sigma_t(ig) - xsec_unpert%matrixs(iz,ia)%sigma_t(ig))
                    do iig = 1, ns%state%ng
                        rho_s = rho_s + flux_adjoint_pt%ngs(ig)%scalar(ir,ia) * (xsec%matrixs(iz,ia)%sigma_s(iig,ig,1) - xsec_unpert%matrixs(iz,ia)%sigma_s(iig,ig,1))
                        rho_f_nu = rho_f_nu + flux_adjoint_pt%ngs(ig)%scalar(ir,ia) * ((xsec%matrixs(iz,ia)%sigma_f_nu(iig) - xsec_unpert%matrixs(iz,ia)%sigma_f_nu(iig)))
                    end do
                    
                    tmp_1 = rho_t - rho_s - rho_f_nu/iter_count_unpert%eigenvalue
                    tmp_2 = xsec%matrixs(iz,ia)%sigma_f_nu(ig) * flux_forward_unpert%ngs(ig)%scalar(ir,ia)
                    
                    iter_adjoint%source(iz, ia, ig) = iter_adjoint%source(iz, ia, ig) + tmp_1/rho - tmp_2/denominator
                end do
            end do
        end do
    
    end subroutine Get_gpt_source
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Reactivity_by_EPT (rho)
        
        real(KREAL), intent(in out)  :: rho
        
        ! local variables
        real(KREAL)  :: denominator
        real(KREAL)  :: rho_t, rho_s, rho_f_nu
        
        call Get_integration (flux_forward, flux_adjoint_pt, denominator, rho_t, rho_s, rho_f_nu)
        
        rho = - (rho_t - rho_s - rho_f_nu/iter_count_unpert%eigenvalue) / denominator
        
        if (.FALSE.)  then
            write(110, *)  ''
            write(110, *)  '______________________'
            write(110, *)  'by EPT:'
            write(110, *)  'rho_t       = ', rho_t
            write(110, *)  'rho_s       = ', rho_s
            write(110, *)  'rho_f_nu    = ', rho_f_nu/iter_count_unpert%eigenvalue
            write(110, *)  'denominator = ', denominator
        end if 
        
    end subroutine Reactivity_by_EPT
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Reactivity_by_PT (rho)
        
        real(KREAL), intent(in out)  :: rho
        
        ! local variables
        real(KREAL)  :: denominator
        real(KREAL)  :: rho_t, rho_s, rho_f_nu
        
        call Get_integration (flux_forward_unpert, flux_adjoint_pt, denominator, rho_t, rho_s, rho_f_nu)
        
        rho = - (rho_t - rho_s - rho_f_nu/iter_count_unpert%eigenvalue) / denominator
        
        write(110, *)  ''
        write(110, *)  '______________________'
        write(110, *)  'by PT:'
        write(110, *)  'rho_t       = ', rho_t
        write(110, *)  'rho_s       = ', rho_s
        write(110, *)  'rho_f_nu    = ', rho_f_nu/iter_count_unpert%eigenvalue
        write(110, *)  'denominator = ', denominator
    
    end subroutine Reactivity_by_PT
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Reactivity_by_GPT (rho_in, rho)
        
        real(KREAL), intent(in)      :: rho_in
        real(KREAL), intent(in out)  :: rho
        
        ! local variables
        real(KREAL)  :: denominator
        real(KREAL)  :: rho_t, rho_s, rho_f_nu
        real(KREAL)  :: factor
        
        call Get_integration (flux_forward_unpert, flux_adjoint_gpt, denominator, rho_t, rho_s, rho_f_nu)
        
        ! correct rho
        factor = 0.0
        factor = (rho_t - rho_s - rho_f_nu/iter_count_unpert%eigenvalue)
        factor = factor - rho_f_nu*((1.0/iter_count%eigenvalue) - (1.0/iter_count_unpert%eigenvalue))
        
        if (factor < EPS_ZERO)  then
            rho = rho_in * (1.0 - factor)
        else    
            rho = rho_in * (1.0 - factor/(1.0 + factor))
        end if
        
        write(110, *)  ''
        write(110, *)  '______________________'
        write(110, *)  'by GPT:'
        write(110, *)  'factor     = ', factor
        write(110, *)  'correction = ', rho_in / rho
    
    end subroutine Reactivity_by_GPT
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_integration (forward, adjoint, denominator, rho_t, rho_s, rho_f_nu)
        
        type(GroupsFlux), intent(in)     :: forward
        type(GroupsFlux), intent(in)     :: adjoint
        real(KREAL), intent(in out)  :: denominator
        real(KREAL), intent(in out)  :: rho_t
        real(KREAL), intent(in out)  :: rho_s
        real(KREAL), intent(in out)  :: rho_f_nu
        
        ! local variables
        integer  :: iig, ig, ir, ia, im, iz, is
        
        ! obtains the denominator
        denominator = 0.0
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do iig = 1, ns%state%ng
                        denominator = denominator + xsec%matrixs(iz,ia)%chi_steady(ig) * adjoint%ngs(ig)%scalar(ir,ia)   &
                            &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * forward%ngs(iig)%scalar(ir,ia)                   &
                            &   * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! obtains reactivity
        rho_t = 0.0        
        ! by sigma_t -- by scalar
        do ig = 1, ns%state%ng  
            do ia = 1, ns%state%layer 
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    rho_t = rho_t + adjoint%ngs(ig)%scalar(ir,ia) * forward%ngs(ig)%scalar(ir,ia)                   &
                        &   * (xsec%matrixs(iz,ia)%sigma_t(ig) - xsec_unpert%matrixs(iz,ia)%sigma_t(ig))            &
                        &   * geom%area(ir) * geom%height(ia)
                end do
            end do
        end do
        
        if (.FALSE.)  then
            write(111, *) ''
            write(111, *) '____________________________________'
            write(111, *) 'rho total by scalar  :',  rho_t
        end if 
        
        rho_t = 0.0
        do ig = 1, ns%state%ng  
            do is = 1, ns%deduce%direction
                do ia = 1, ns%state%layer 
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        rho_t = rho_t + adjoint%ngs(ig)%angular(ir,ia,is) * forward%ngs(ig)%angular(ir,ia,is)           &
                            &   * (xsec%matrixs(iz,ia)%sigma_t(ig) - xsec_unpert%matrixs(iz,ia)%sigma_t(ig))            &
                            &   * geom%area(ir) * geom%height(ia) * quad%directions(is)%wmu
                    end do
                end do
            end do
        end do
        write(111, *) 'rho total by angular :',  rho_t
        
        ! by sigma_s
        rho_s = 0.0
        do ig = 1, ns%state%ng 
            do ia = 1, ns%state%layer 
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do iig = 1, ns%state%ng
                        rho_s = rho_s + adjoint%ngs(ig)%scalar(ir,ia) * forward%ngs(iig)%scalar(ir,ia)                        &
                            &   * (xsec%matrixs(iz,ia)%sigma_s(iig,ig,1) - xsec_unpert%matrixs(iz,ia)%sigma_s(iig,ig,1))      &
                            &   * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! @TODO--harmonic expansion
        if (ns%state%scat_order /= 0)  then
        end if
        
        ! by sigma_f_nu
        rho_f_nu = 0.0
        do ig = 1, ns%state%ng 
            do ia = 1, ns%state%layer 
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do iig = 1, ns%state%ng
                        rho_f_nu = rho_f_nu + adjoint%ngs(ig)%scalar(ir,ia) * forward%ngs(iig)%scalar(ir,ia)              &
                            &   * (xsec%matrixs(iz,ia)%sigma_f_nu(iig) - xsec_unpert%matrixs(iz,ia)%sigma_f_nu(iig))      &
                            &   * geom%area(ir) * geom%height(ia) * xsec_unpert%matrixs(iz,ia)%chi_steady(ig)
                    end do
                end do
            end do
        end do
        
    end subroutine Get_integration
    
end module reactivity
