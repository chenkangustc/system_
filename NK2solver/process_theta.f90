!$
!===================================================================================================
!
!   module for performing perturbation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Generate_TFSP
!
!   Public type lists:          No
!
!===================================================================================================
module process_theta
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    implicit none 
    private
    public  :: Generate_TFSP
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Generate_TFSP (time_pace)
    
        real(KREAL), intent(in)  :: time_pace
        
        ! local variables
        real(KREAL)  :: q_first(ns%state%nodal, ns%state%layer, ns%state%ng)                  ! first part of modified exteral source
        real(KREAL)  :: q_second(ns%state%nodal, ns%state%layer, ns%state%ng)                 ! second part of modified exteral source
!        real(KREAL)  :: time_pace                                                             ! time step length per time step
        real(KREAL)  :: tmp
        real(KREAL)  :: factor
        real(KREAL)  :: strength
        integer  :: im, ig, ip, ia, ir, iz, is
        
        ! select time step length
!        if (present(is_macro) .AND. is_macro )  then
!            time_pace = time_step%get_current(tidx, is_step=.FALSE.)
!        else 
!            time_pace = time_step%get_current(tidx, is_step=.TRUE.)
!        end if
        
        ! ----------------------------------------------------------------------
        ! construct TFSP by three step
        
        ! obtain modified total cross section
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    associate(velocity => param%matrixs(iz, ia)%velocity(ig))
                    if (velocity > EPS_ZERO)  then
                        xsec_iter%matrixs(iz, ia)%sigma_t(ig) = xsec_iter%matrixs(iz, ia)%sigma_t(ig) + 1.0/(time_pace*velocity)
                    else 
                        xsec_iter%matrixs(iz, ia)%sigma_t(ig) = xsec_iter%matrixs(iz, ia)%sigma_t(ig) + 0.0D0
                    end if
                    end associate
                end do
            end do
        end do
        
        ! obtain modified fission spectrum
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    tmp = 0.0
                    do ip = 1, nt%state%dg
                        ! NOTE: steady spectrum is chi_steady, this is default
                        factor = 1.0 / (1.0+time_pace*param%matrixs(iz, ia)%lambda(ip))
                        tmp = tmp - param%matrixs(iz, ia)%beta(ip)*param%matrixs(iz, ia)%chi_delay(ip, ig)*factor

                        ! steady spectrum is chi_prompt
!                        factor = (param%matrixs(iz, ia)%beta(ip) * time_pace) / (1.0+time_pace*param%matrixs(iz, ia)%lambda(ip)) 
!                        xsec_iter%matrixs(iz, ia)%chi_steady(ig) = xsec_iter%matrixs(iz, ia)%chi_steady(ig) + param%matrixs(iz, ia)%lambda(ip)*param%matrixs(iz, ia)%chi_delay(ip, ig)*factor  &
!                            &   - xsec_iter%matrixs(iz, ia)%chi_steady(ig)*param%matrixs(iz, ia)%beta(ip)
                    end do
                    
                    xsec_iter%matrixs(iz, ia)%chi_steady(ig) = xsec_iter%matrixs(iz, ia)%chi_steady(ig) + tmp
                end do
            end do
        end do
        
        ! obtain modified external source (scalar value, consist of three parts)
        ! the first part--precursor
        q_first = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                do ig = 1, ns%state%ng
                    tmp = 0.0
                    do ip = 1, nt%state%dg
                        factor = 1.0 / (1.0+time_pace*param%matrixs(iz, ia)%lambda(ip))
                        tmp = tmp + param%matrixs(iz, ia)%chi_delay(ip, ig)*dnp_solver%Lprecursor(ip, ir, ia)*param%matrixs(iz, ia)%lambda(ip)*factor
                    end do
                    q_first(ir, ia, ig) = tmp
                end do
            end do 
        end do
        
        ! the second part--velocity
        q_second = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ig = 1, ns%state%ng
                    iz = mesh%zone(ir)
                    associate(velocity => param%matrixs(iz, ia)%velocity(ig))
                    if (velocity > EPS_ZERO)  then
                        q_second(ir, ia, ig) = q_second(ir, ia, ig) + flux_forward%ngs(ig)%scalar(ir, ia) * (1.0/(time_pace*velocity))
                    else
                        q_second(ir, ia, ig) = q_second(ir, ia, ig) + 0.0D0
                    end if
                    end associate
                end do
            end do
        end do
        
        call Q_ext%strength (geom, strength)
!        write(850, *)  "tidx", tidx, "strength", strength
        
        ! the total result
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ig = 1, ns%state%ng
                    Q_ext%iter_scalar(ir,ia)%intensity(ig) = q_first(ir,ia,ig) + q_second(ir, ia, ig) + Q_ext%iter_scalar(ir,ia)%intensity(ig)
                end do
            end do
        end do
        
        call Q_ext%strength (geom, strength)
!        write(851, *)  "tidx", tidx, "strength", strength
    
    end subroutine Generate_TFSP
    
end module process_theta
