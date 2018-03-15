!$
!===================================================================================================
!
!   initialize the parameter used in iteration process
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Init_iteration_variable
!
!   Public type lists:          No
!
!===================================================================================================
module iteration_initialize
 
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    implicit none
    private
    public  :: Init_iteration_variable

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Init_iteration_variable (is_adjoint)
    
        logical, intent(in)  :: is_adjoint
        
        ! loop index
        integer  :: ig, ia, ir, il, iz, iy
        integer  :: iy_count
        
        ! local variables
        real(KREAL)   :: tmp_sum, tmp_flux                                  ! temporary for initializing
        integer           :: ml
        integer, save     :: n_iter = 0                                         ! number of transport solver evoked, Note: save
        
        if (.NOT. is_adjoint)  then
            n_iter = n_iter + 1
        end if
        
        ! for adjoint or the first one none-adjoint iteration
        iter: if (n_iter==1 .or. is_adjoint)  then
            iter_count%eigenvalue = 1.0D0
            iter_count%ks = 1.0D0
        
            ! initialize the interation source
            iter_q%fission%moment  = 0.0
            iter_q%fission%old     = 0.0
            
            iter_q%info%total_moment     = 0.0
            iter_q%info%out_group_moment = 0.0
            
            ! initialize the iteration flux
            iter_flux%info%moment      = 0.0
            iter_flux%info%moment_omp  = 0.0
            iter_flux%info%old         = 0.0
            
            iter_flux%dist%surface    = 1.0                                     ! here not zero
            iter_flux%dist%nodal      = 1.0
            iter_flux%dist%point      = 1.0
            
            iter_flux%dist%axi_surf = 1.0
            iter_flux%dist%rad_surf = 1.0
            
            ! initialize anisotropic scatter flux variables
            do ig = 1, ns%state%ng
                flux_scat%ngs(ig)%aniso_zero  = 0.0
                flux_scat%ngs(ig)%aniso_cos   = 0.0
                flux_scat%ngs(ig)%aniso_sin   = 0.0
            end do
            
            ! initialize fission-source and flux by sigma_f_nu 
            tmp_flux = 1.0
            do ia = 1, ns%state%layer
                ml = (ia-1) * ns%state%nodal
                do ir = 1, ns%state%nodal
                    il = ml + ir
                    iz = mesh%zone(ir)
                    tmp_sum = 0.0
                    do ig = 1, ns%state%ng
                        iter_flux%info%moment(0,il,ig) = tmp_flux
                        tmp_sum = tmp_sum + xsec%matrixs(iz,ia)%sigma_f_nu(ig) * tmp_flux
                    end do
                    iter_q%fission%moment(0,il) = tmp_sum
                    iter_q%fission%old(il) = tmp_sum
                end do
            end do
            
            ! initialize surface flux per nodal by flux (all surface and all direction)
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    il = (ia-1)*ns%state%nodal + ir
                    iter_flux%dist%surface(:,ir,ia,:) = iter_flux%info%moment(0, il, 1)
                end do
            end do
           
            ! initialize surface flux per nodal by flux (all surface and all direction)
            do ia = 1, ns%state%layer
                iy_count = 0
                do ir = 1, ns%state%nodal
                    il = (ia-1)*ns%state%nodal + ir
                    
                    ! for radial
                    do ig = 1, ns%state%ng
                        iter_flux%dist%rad_surf(:, ia, :, ig) = iter_flux%info%moment(0, il, ig)
                    end do
                    
                    ! for axial
                    do ig = 1, ns%state%ng
                        if (ia == 1)  then
                            iter_flux%dist%axi_surf(ir, 1, :, ig) = iter_flux%info%moment(0, il, ig)
                        end if
                        if (ia == ns%state%layer)  then
                            iter_flux%dist%axi_surf(ir, 2, :, ig) = iter_flux%info%moment(0, il, ig)
                        end if
                    end do
                end do
            end do
        
        ! for any other iteration, use value from last time
        else iter
        
        end if iter
        
        ! initialize fission source term for adjoint calculation
        if (is_adjoint) then 
            iter_adjoint%q_moments(0, :) = iter_q%fission%moment(0, :)
            iter_adjoint%q_old(:) = iter_q%fission%old(:)
        end if 
                
    end subroutine Init_iteration_variable
    
end module iteration_initialize
