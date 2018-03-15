!$
!===================================================================================================
!
!   module for theta differential method
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Get_initial_flux
!                               Get_initial_precursor
!                               Update_precursor
!
!   Public type lists:          No
!
!===================================================================================================
module time_advancing
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use transit_from_solver,        only : Transit_precursor
    use timestep_header,            only : TimeStepInfo
    
    implicit none 
    private
    public  :: Get_initial_precursor, Update_precursor, Get_physics_kinetics_parameter
    
    integer, parameter  :: FISSION_RATE_CONSTANT = 0
    integer, parameter  :: FISSION_RATE_LINEAR = 1
    integer, parameter  :: METHOD_ = 0
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_initial_precursor ()
    
        integer  :: ig, ia, ir, ip, iz
        real(KREAL)  :: reaction_rate
    
        ! obtain fission reaction rate to generate precursor concentration
        
        dnp_solver%Rprecursor = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                
                reaction_rate = 0.0
                do ig = 1, ns%state%ng
                    reaction_rate = reaction_rate + flux_forward%ngs(ig)%scalar(ir, ia)*xsec%matrixs(iz, ia)%sigma_f_nu(ig)
                end do
                
                do ip = 1, nt%state%dg
                    associate(beta => param%matrixs(iz, ia)%beta(ip), lambda => param%matrixs(iz, ia)%lambda(ip))
                    if (lambda > EPS_ZERO)  then
                        dnp_solver%Rprecursor(ip, ir, ia) = reaction_rate * (beta/lambda)
                    else
                        dnp_solver%Rprecursor(ip, ir, ia) = 0.0D0
                    end if
                    end associate
                end do
            end do
        end do
        
        ! update precursor concentration for the previous time point
        dnp_solver%Lprecursor = dnp_solver%Rprecursor
    
    end subroutine Get_initial_precursor

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_precursor (a_step)
    
        type(TimeStepInfo), intent(in)  :: a_step
!        integer, intent(in)  :: tidx
!        logical, intent(in), optional  :: is_macro
        
        ! loop index
        integer  :: ig, ia, ir, ip, iz
        
        real(KREAL)  :: reaction_rate
        real(KREAL)  :: factor
        real(KREAL)  :: time_pace
        real(KREAL)  :: ctime
        
        ! distinguish section and step
        time_pace = a_step%pace
        ctime = a_step%right
        
!        if (present(is_macro) .AND. is_macro)  then
!            time_pace = time_step%get_current (tidx, is_step=.FALSE.)
!            ctime = time_step%get_point (tidx, is_step=.FALSE.)
!        else
!            time_pace = time_step%get_current (tidx, is_step=.TRUE.)
!            ctime = time_step%get_point (tidx, is_step=.TRUE.)
!        end if
        
        ! obtain fission reaction rate to generate precursor concentration
        select case(METHOD_)
        case (FISSION_RATE_CONSTANT)
            dnp_solver%Rprecursor = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    
                    reaction_rate = 0.0
                    do ig = 1, ns%state%ng
                        reaction_rate = reaction_rate + flux_forward%ngs(ig)%scalar(ir, ia)*xsec%matrixs(iz,ia)%sigma_f_nu(ig)
                    end do
                    
                    do ip = 1, nt%state%dg
                        associate(beta => param%matrixs(iz,ia)%beta(ip), lambda => param%matrixs(iz,ia)%lambda(ip))
                        if (lambda > EPS_ZERO)  then
                            factor = 1.0 / (1.0 + time_pace*lambda)
                            dnp_solver%Rprecursor(ip, ir, ia) = dnp_solver%Lprecursor(ip, ir, ia)*factor + reaction_rate*time_pace*beta*factor
                        else 
                            dnp_solver%Rprecursor(ip, ir, ia) = 0.0D0
                        end if 
                        end associate
                    end do
                end do
            end do
        
        case (FISSION_RATE_LINEAR)
            dnp_solver%Rprecursor = 0.0
            
        end select 
        
        ! transit precursor concentration to last time point
        dnp_solver%Lprecursor = dnp_solver%Rprecursor
    
    end subroutine Update_precursor

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_physics_kinetics_parameter ()
        
        real(KREAL)  :: beta(nt%state%dg)
        real(KREAL)  :: lambda(nt%state%dg)
        real(KREAL)  :: denominator
        real(KREAL)  :: numerator
        
        integer  :: ig, ia, ir, ip, iz, im
        
        beta = 0.0
        lambda = 0.0
        
        do ip = 1, nt%state%dg
            ! for beta
            numerator = 0.0
            denominator = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do ig = 1, ns%state%ng
                        numerator = numerator + param%matrixs(iz, ia)%sigma_bvf(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)
                        denominator = denominator + xsec%matrixs(iz, ia)%sigma_f_nu(ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
            
            beta(ip) = numerator / denominator
            
            ! for lambda
            numerator = 0.0
            denominator = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do ig = 1, ns%state%ng
                        numerator = numerator + param%matrixs(iz, ia)%sigma_bvf(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)
                        denominator = denominator + param%matrixs(iz, ia)%sigma_bvl(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
            
            lambda(ip) = numerator / denominator
        end do
        
        ! set this set to all the material and nodal
        do im = 1, ns%state%mat
            param_inp%mats(im)%beta = beta
            param_inp%mats(im)%lambda = lambda
        end do
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                param%matrixs(iz, ia)%beta = beta
                param%matrixs(iz, ia)%lambda = lambda
            end do
        end do
        
        write(FILES%REACTIVITY, *) ' '
        write(FILES%REACTIVITY, "(1x, A)")  'physics kinetics parameter:'
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))")  'partial beta   = ', beta
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))")  'partial lambda = ', lambda
    
    end subroutine Get_physics_kinetics_parameter
    
end module time_advancing
