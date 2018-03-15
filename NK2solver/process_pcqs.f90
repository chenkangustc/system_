!$
!===================================================================================================
!
!   module for macro-step of predict-correct quasi-static method
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Get_initial_flux_shape                          
!                               Get_initial_precursor_shape
!                               Get_initial_amplitude
!                               Get_flux_shape
!                               Get_precursor_shape
!                               Transit_shape_function
!                               Normal_adjoint_flux
!
!                               Shape_interpolation
!                               Generate_pk_parameter
!                               Generate_reactivity
!                               Regenerate_power
!                               Get_correct_flux
!                               Get_correct_precursor
!
!   Public type lists:          No
!
!===================================================================================================
module process_pcqs
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use timestep_header,    only : TimeStepInfo
    
    implicit none 
    private
    
    ! for macro step
    public  :: Get_initial_flux_shape, Get_initial_precursor_shape, Get_initial_amplitude
    public  :: Get_flux_shape, Get_precursor_shape
    public  :: Transit_shape_function, Normal_adjoint_flux
    
    ! for media step
    public  :: Shape_interpolation, Generate_pk_parameter, Generate_reactivity
    public  :: Regenerate_power, Get_correct_flux, Get_correct_precursor
    
contains
    !$
    !===============================================================================================
    ! at time point zero, equal flux shape to flux
    !===============================================================================================
    subroutine Get_initial_flux_shape ()
        
        integer  :: ig
        
        do ig = 1, ns%state%ng
            shape_current%flux_angular(:, :, :, ig) = flux_forward%ngs(ig)%angular(:, :, :)
            shape_current%flux_scalar(:, :, ig) = flux_forward%ngs(ig)%scalar(:, :)
        end do
    
        shape_last%flux_angular = shape_current%flux_angular
        shape_last%flux_scalar = shape_current%flux_scalar
        
        shape_predict%flux_angular = shape_current%flux_angular
        shape_predict%flux_scalar = shape_current%flux_scalar
        
    end subroutine Get_initial_flux_shape
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_initial_precursor_shape ()
    
        integer  :: ip, ia, ir, im, ig, is, iz
        real(KREAL)  :: normal_factor(nt%state%dg)
        
        do ip = 1, nt%state%dg
            shape_current%precursor(ip, :, :) = dnp_solver%Lprecursor(ip, :, :) / (pk_parameter%partial_beta(ip)/(pk_parameter%partial_lambda(ip)*pk_parameter%generation_time))
        end do
        
        ! normalize 
        do ip = 1, nt%state%dg
            normal_factor(ip) = 0.0
            
            do is = 1, ns%deduce%direction
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            
                            normal_factor(ip) = normal_factor(ip) + shape_current%precursor(ip,ir,ia) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu    &
                                &   * param%matrixs(iz,ia)%chi_delay(ip,ig)
                        end do
                    end do
                end do
            end do
            
            shape_current%precursor(ip,:,:) = shape_current%precursor(ip,:,:) / normal_factor(ip)
        end do
        
        shape_last%precursor = shape_current%precursor
        shape_predict%precursor = shape_current%precursor
    
    end subroutine Get_initial_precursor_shape
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_initial_amplitude ()
    
        integer  :: ip
        
        amplitude%flux = 1.0
        
        do ip = 1, nt%state%dg
            amplitude%precursor(ip) = pk_parameter%partial_beta(ip) / (pk_parameter%partial_lambda(ip)*pk_parameter%generation_time)
        end do
    
    end subroutine Get_initial_amplitude
    
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_flux_shape (tidx)
    
        integer, intent(in)  :: tidx
    
        ! loop index
        integer  :: ia, ir, im, ig, is, iz
        
        ! local variables
        real(KREAL)  :: normal_factor
        
        ! ----------------------------------------------------------------------
        ! get normal factor
        normal_factor = 0.0
        
        do is = 1, ns%deduce%direction
            do ig = 1, ns%state%ng
                do ia = 1, ns%state%layer
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        associate(velocity => param%matrixs(iz,ia)%velocity(ig))
                        if (velocity > EPS_ZERO)  then
                            normal_factor = normal_factor + flux_forward%ngs(ig)%angular(ir,ia,is) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu   &
                                &   * geom%area(ir)*geom%height(ia) / velocity
                        else
                            normal_factor = normal_factor + 0.0D0
                        end if 
                        end associate
                    end do
                end do
            end do
        end do
        
        do ig = 1, ns%state%ng
            shape_predict%flux_angular(:, :, :, ig) = flux_forward%ngs(ig)%angular(:, :, :) / normal_factor
            shape_predict%flux_scalar(:, :, ig) = flux_forward%ngs(ig)%scalar(:, :) / normal_factor
        end do
        
    end subroutine Get_flux_shape
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_precursor_shape (step_Medial)
        
        type(TimeStepInfo)  :: step_Medial
    
        ! loop index
        integer  :: ia, ir, im, ig, is, ip, iz
        
        ! local variables
        real(KREAL)  :: normal_factor(nt%state%dg)
        real(KREAL)  :: factor
        real(KREAL)  :: denominator
        real(KREAL)  :: time_pace
        
        ! ----------------------------------------------------------------------
        time_pace = step_Medial%right
!        time_pace = time_step%get_current (tidx, is_step=.TRUE.)
    
        ! obtains the factor first
        factor = 0.0
        
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    factor = factor + xsec%matrixs(iz,ia)%sigma_f_nu(ig) * shape_current%flux_scalar(ir,ia,ig)     &
                        &   * geom%area(ir) * geom%height(ia)
                end do
            end do
        end do
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                do ip = 1, nt%state%dg
                    denominator = 1 + time_pace * (pk_parameter%partial_beta(ip)/pk_parameter%generation_time) * (amplitude%flux/amplitude%precursor(ip))
                    
                    shape_current%precursor(ip, ir, ia) = shape_last%precursor(ip, ir, ia) + factor*param%matrixs(iz,ia)%beta(ip)*(time_pace*amplitude%flux/amplitude%precursor(ip))
                    shape_current%precursor(ip, ir, ia) = shape_current%precursor(ip, ir, ia) / denominator
                end do
            end do
        end do
        
        ! normalize 
        do ip = 1, nt%state%dg
            normal_factor(ip) = 0.0
            
            do is = 1, ns%deduce%direction
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            
                            normal_factor(ip) = normal_factor(ip) + shape_current%precursor(ip,ir,ia) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu    &
                                &   * param%matrixs(iz,ia)%chi_delay(ip,ig)
                        end do
                    end do
                end do
            end do
            
            shape_current%precursor(ip,:,:) = shape_current%precursor(ip,:,:) / normal_factor(ip)
        end do
        
        ! transfer to current to last
        shape_last%precursor = shape_current%precursor
    
    end subroutine Get_precursor_shape
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_shape_function (tidx)
        
        integer, intent(in)  :: tidx
        
        shape_current = shape_predict
        shape_last = shape_predict
        
    end subroutine Transit_shape_function
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Normal_adjoint_flux ()
    
        ! loop index
        integer  :: ia, ir, im, ig, is, iz
        
        ! local variables
        real(KREAL)  :: normal_factor
        
        ! ----------------------------------------------------------------------
        ! get normal factor
        normal_factor = 0.0
        do is = 1, ns%deduce%direction
            do ig = 1, ns%state%ng
                do ia = 1, ns%state%layer
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        associate(velocity => param%matrixs(iz,ia)%velocity(ig))
                        if (velocity > EPS_ZERO)  then
                            normal_factor = normal_factor + flux_forward%ngs(ig)%angular(ir,ia,is) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu   &
                                &   * geom%area(ir)*geom%height(ia) / velocity
                        else 
                            normal_factor = normal_factor + 0.0D0
                        end if
                        end associate
                    end do
                end do
            end do
        end do
        
        do ig = 1, ns%state%ng
            flux_adjoint%ngs(ig)%angular(:, :, :) = flux_adjoint%ngs(ig)%angular(:, :, :) / normal_factor
            flux_adjoint%ngs(ig)%scalar(:, :) = flux_adjoint%ngs(ig)%scalar(:, :) / normal_factor
        end do
    
    end subroutine Normal_adjoint_flux 
    
    ! --------------------------------------------------------------------------
    ! following is for media step
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! re-generate absolute power and precursor, integreate and distribution
    ! same function as <Get_correct_flux & Get_correct_precursor>
    !===============================================================================================
    subroutine Regenerate_power (tidx)
        
        integer, intent(in)  :: tidx
        
        integer  :: ig, ia, ir, iz, ip, il
        real(KREAL)  :: matrix(ns%state%nodal, ns%state%layer)              ! transfer into distribution
        
        ! ----------------------------------------------------------------------
        ! for nuetron
        do ig = 1, ns%state%ng
            flux_forward%ngs(ig)%angular(:, :, :) = 0.0
            flux_forward%ngs(ig)%angular(:, :, :) = shape_current%flux_angular(:, :, :, ig) * amplitude%flux
        end do
        call flux_forward%set_scalar (quad)
        
        matrix = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ig = 1, ns%state%ng
                    matrix(ir, ia) = matrix(ir, ia) + flux_forward%ngs(ig)%scalar(ir, ia)
                end do
            end do
        end do
        call dist_flux%set (matrix)

        matrix = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                do ig = 1, ns%state%ng
                    matrix(ir,ia) = matrix(ir,ia) + flux_forward%ngs(ig)%scalar(ir, ia) * xsec_iter%matrixs(iz,ia)%sigma_f_kappa(ig)
                end do
            end do
        end do
        call dist_power%set (matrix)
        
        matrix = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                do ig = 1, ns%state%ng
                    matrix(ir,ia) = matrix(ir,ia) + flux_forward%ngs(ig)%scalar(ir, ia) * xsec_iter%matrixs(iz,ia)%sigma_f_nu(ig)
                end do
            end do
        end do
        call dist_fission_rate%set (matrix)
          
        ! power level
        timelist%power = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = (ia-1) * ns%state%nodal + ir
                do ig = 1, ns%state%ng
                    timelist%power = timelist%power + flux_forward%ngs(ig)%scalar(ir,ia)*xsec%matrixs(iz,ia)%sigma_f_kappa(ig)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
        
        ! ----------------------------------------------------------------------
        ! for precursor
        do ip = 1, nt%state%dg
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    dist_dnps(ip)%matrix(ir, ia) = shape_current%precursor(ip, ir, ia) * amplitude%precursor(ip)
                end do
            end do
        end do
        
        timelist%precursor = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ip = 1, nt%state%dg
                    timelist%precursor(ip) = timelist%precursor(ip) + dist_dnps(ip)%matrix(ir,ia)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
    
    end subroutine Regenerate_power

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_correct_flux (tidx)
        
        integer, intent(in)  :: tidx
        
        integer  :: ig, ia, ir, iz, ip, il
        
        ! ----------------------------------------------------------------------
        ! for nuetron
        do ig = 1, ns%state%ng
            flux_forward%ngs(ig)%scalar(:, :) = 0.0
            flux_forward%ngs(ig)%scalar(:, :) = shape_current%flux_scalar(:, :, ig) * amplitude%flux
            
            flux_forward%ngs(ig)%angular(:, :, :) = 0.0
            flux_forward%ngs(ig)%angular(:, :, :) = shape_current%flux_angular(:, :, :, ig) * amplitude%flux
        end do
        
        ! power level
        timelist%power = 0.0
            
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = (ia-1) * ns%state%nodal + ir
                do ig = 1, ns%state%ng
                    timelist%power = timelist%power + flux_forward%ngs(ig)%scalar(ir,ia)*xsec%matrixs(iz,ia)%sigma_f_kappa(ig)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
    
    end subroutine Get_correct_flux
    
    !$
    !===============================================================================================
    ! equal to subroutine 'Update_precursor'
    !===============================================================================================    
    subroutine Get_correct_precursor (tidx)
        
        integer, intent(in)  :: tidx
        
        integer  :: ig, ia, ir, iz, ip, il
    
        ! ----------------------------------------------------------------------
        ! for precursor
        do ip = 1, nt%state%dg
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    dist_dnps(ip)%matrix(ir, ia) = shape_current%precursor(ip, ir, ia) * amplitude%precursor(ip)
                end do
            end do
        end do
        
        timelist%precursor = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ip = 1, nt%state%dg
                    timelist%precursor(ip) = timelist%precursor(ip) + dist_dnps(ip)%matrix(ir,ia)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
        
    end subroutine Get_correct_precursor
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Shape_interpolation (step_Macro, step_Medial)
    
        type(TimeStepInfo), intent(in)  :: step_Macro                                     ! current macro step index
        type(TimeStepInfo), intent(in)  :: step_Medial                                     ! current mdedia step index
    
        integer  :: ig, ip
        real(KREAL)  :: time_begin
        real(KREAL)  :: time_end
        real(KREAL)  :: time_current
        real(KREAL)  :: frac_last, frac_next
        
        ! ----------------------------------------------------------------------
        ! get time interval
        time_end = step_Macro%right
        time_begin = step_Macro%left 
        time_current = step_Medial%right
         
!        time_end = time_step%get_point(index_macro, is_step=.FALSE.)
!        time_begin = time_end - time_step%get_current(index_macro, is_step=.FALSE.)
!        time_current = time_step%get_point(tidx, is_step=.TRUE.)
        
        frac_last = (time_end - time_current) / (time_end - time_begin)
        frac_next = (time_current - time_begin) / (time_end - time_begin)
        
        ! ----------------------------------------------------------------------
        ! interpolation shape function
        
        ! for flux
        do ig = 1, ns%state%ng
            shape_current%flux_angular(:, :, :, ig) = frac_last*shape_last%flux_angular(:, :, :, ig) + frac_next*shape_predict%flux_angular(:, :, :, ig)
            shape_current%flux_scalar(:, :, ig) = frac_last*shape_last%flux_scalar(:, :, ig) + frac_next*shape_predict%flux_scalar(:, :, ig)
        end do
        
        ! for precursor
        do ip = 1, nt%state%dg
            shape_current%precursor(ip, :, :) = frac_last*shape_last%precursor(ip, :, :) + frac_next*shape_predict%precursor(ip, :, :)
        end do
    
    end subroutine Shape_interpolation
    
    !$
    !===============================================================================================
    ! generate parameter for point kinetics approximate
    ! ----------------------------------------------------------------------------------------------
    !
    ! NOTE: total weight is normal to 1
    !===============================================================================================
    subroutine Generate_pk_parameter (tidx)
        
        integer, intent(in)  :: tidx
        
        ! loop index
        integer  :: is, ig, iig, ir, ia, il, im, ip, iz
        
        ! local variables
        real(KREAL)  :: denominator, numerator
        real(KREAL)  :: lambda_denominator, lambda_numerator
        
        real(KREAL)  :: beta(nt%state%dg)
        real(KREAL)  :: lambda(nt%state%dg)
       
        if (nt%method%scheme == 'THETA')  then
            ! obtains the denominator
            denominator = 0.0
            do ig = 1, ns%state%ng
                do iig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            denominator = denominator + xsec%matrixs(iz,ia)%chi_steady(ig) * flux_adjoint%ngs(ig)%scalar(ir,ia)    &
                                &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * flux_forward%ngs(iig)%scalar(ir,ia)                    &
                                &   * geom%area(ir) * geom%height(ia)
                        end do
                    end do
                end do
            end do
            
            ! obtains partial beta fraction
            pk_parameter%beta = 0.0
            pk_parameter%partial_beta = 0.0
            
            do ip = 1, nt%state%dg
                do ig = 1, ns%state%ng
                    do iig = 1, ns%state%ng
                        do ia = 1, ns%state%layer
                            do ir = 1, ns%state%nodal
                                iz = mesh%zone(ir)
                                pk_parameter%partial_beta (ip) = pk_parameter%partial_beta(ip)                                                                  &
                                    &   + param%matrixs(iz,ia)%chi_delay(ip,ig) * param%matrixs(iz,ia)%beta(ip) * flux_adjoint%ngs(ig)%scalar(ir,ia)    &
                                    &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * flux_forward%ngs(iig)%scalar(ir,ia)                                            &
                                    &   * geom%area(ir) * geom%height(ia)
                            end do
                        end do
                    end do
                end do
            end do
            
            pk_parameter%partial_beta = pk_parameter%partial_beta / denominator
            pk_parameter%beta = SUM(pk_parameter%partial_beta)
            
            ! obtains generation time
            pk_parameter%generation_time = 0.0
            
            do is = 1, ns%deduce%direction
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            associate(velocity => param%matrixs(iz,ia)%velocity(ig))
                            if (velocity > EPS_ZERO)  then
                                pk_parameter%generation_time = pk_parameter%generation_time + flux_forward%ngs(ig)%angular(ir,ia,is) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu * velocity
                            else
                                pk_parameter%generation_time = pk_parameter%generation_time + 0.0D0
                            end if 
                            end associate
                        end do
                    end do
                end do
            end do
            
            pk_parameter%generation_time = pk_parameter%generation_time / denominator
            
            ! obtains delay constant
            pk_parameter%partial_lambda = 0.0
            
            do ip = 1, nt%state%dg
                lambda_numerator = 0.0
                lambda_denominator = 0.0
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            
                            ! use concentration as weighting function instead of precursor shape function
                            ! Note--change back to shape in the future
                            lambda_numerator = lambda_numerator + flux_adjoint%ngs(ig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,ig)    &
                                &   * dnp_solver%Rprecursor(ip,ir,ia) * param%matrixs(iz,ia)%lambda(ip) * geom%area(ir) * geom%height(ia)
                                
                            lambda_denominator = lambda_denominator + flux_adjoint%ngs(ig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,ig)  &
                                &   * dnp_solver%Rprecursor(ip,ir,ia) * geom%area(ir) * geom%height(ia)
                            
                        end do
                    end do
                end do
                pk_parameter%partial_lambda(ip) = lambda_numerator / lambda_denominator
            end do

        else 
            !-----------------------------------------------------------------------
            ! obtains the denominator
            denominator = 0.0
            do ig = 1, ns%state%ng
                do iig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            denominator = denominator + xsec%matrixs(iz,ia)%chi_steady(ig) * flux_adjoint%ngs(ig)%scalar(ir,ia)    &
                                &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * shape_current%flux_scalar(ir,ia,iig)                   &
                                &   * geom%area(ir) * geom%height(ia)
                        end do
                    end do
                end do
            end do
            
            ! ----------------------------------------------------------------------
            ! obtains partial beta fraction
            pk_parameter%beta = 0.0
            pk_parameter%partial_beta = 0.0
            
            do ip = 1, nt%state%dg
                do ig = 1, ns%state%ng
                    do iig = 1, ns%state%ng
                        do ia = 1, ns%state%layer
                            do ir = 1, ns%state%nodal
                                iz = mesh%zone(ir)
                                pk_parameter%partial_beta (ip) = pk_parameter%partial_beta(ip)                                                                  &
                                    &   + param%matrixs(iz,ia)%chi_delay(ip,ig) * param%matrixs(iz,ia)%beta(ip) * flux_adjoint%ngs(ig)%scalar(ir,ia)    &
                                    &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * shape_current%flux_scalar(ir,ia,iig)                                            &
                                    &   * geom%area(ir) * geom%height(ia)
                            end do
                        end do
                    end do
                end do
            end do
            
            pk_parameter%partial_beta = pk_parameter%partial_beta / denominator
            pk_parameter%beta = SUM(pk_parameter%partial_beta)
            
            ! ----------------------------------------------------------------------
            ! obtains generation time
            pk_parameter%generation_time = 0.0
            
            do is = 1, ns%deduce%direction
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            associate(velocity => param%matrixs(iz,ia)%velocity(ig))
                            if (velocity > EPS_ZERO)  then
                                pk_parameter%generation_time = pk_parameter%generation_time + shape_current%flux_angular(ir,ia,is,ig) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu   &
                                    &   * geom%area(ir)*geom%height(ia) / velocity
                            else
                                pk_parameter%generation_time = pk_parameter%generation_time + 0.0D0
                            end if 
                            end associate
                        end do
                    end do
                end do
            end do
            
            pk_parameter%generation_time = pk_parameter%generation_time / denominator
            
            ! ----------------------------------------------------------------------
            ! obtains delay constant
            pk_parameter%partial_lambda = 0.0
            
            do ip = 1, nt%state%dg
                lambda_numerator = 0.0
                lambda_denominator = 0.0
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            
                            ! use concentration as weighting function instead of precursor shape function
                            ! Note--change back to shape in the future
                            lambda_numerator = lambda_numerator + flux_adjoint%ngs(ig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,ig)    &
                                &   * dnp_solver%Rprecursor(ip,ir,ia) * param%matrixs(iz,ia)%lambda(ip) * geom%area(ir) * geom%height(ia)
                                
                            lambda_denominator = lambda_denominator + flux_adjoint%ngs(ig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,ig)  &
                                &   * dnp_solver%Rprecursor(ip,ir,ia) * geom%area(ir) * geom%height(ia)
                        end do
                    end do
                end do
                pk_parameter%partial_lambda(ip) = lambda_numerator / lambda_denominator
            end do
        end if
        
        write(FILES%REACTIVITY, *) ' '
        write(FILES%REACTIVITY, *) 'effective kinetics parameter by thoery method:'
        write(FILES%REACTIVITY, "(1x, A, I6)") 'index = ', tidx
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))") 'denominator      =', denominator
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))") 'beta             =', pk_parameter%beta
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))") 'generation time  =', pk_parameter%generation_time
        write(FILES%REACTIVITY, *) ' '
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))") 'partial beta   = ', pk_parameter%partial_beta
        write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))") 'partial lambda = ', pk_parameter%partial_lambda
        
        ! by citation method
        if (nt%flag%kinetics_type == 'SRAC')  then
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
                            do iig = 1, ns%state%ng
                                numerator = numerator + param%matrixs(iz, ia)%sigma_bvf(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)     &
                                    &   * param%matrixs(iz,ia)%chi_delay(ip,iig) * flux_adjoint%ngs(iig)%scalar(ir,ia)
                                denominator = denominator + xsec%matrixs(iz, ia)%sigma_f_nu(ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)     &
                                    &   * xsec%matrixs(iz,ia)%chi_steady(iig) * flux_adjoint%ngs(iig)%scalar(ir,ia)
                            end do
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
                            do iig = 1, ns%state%ng
                                numerator = numerator + param%matrixs(iz, ia)%sigma_bvf(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia)     &
                                    &   * flux_adjoint%ngs(iig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,iig)
                                denominator = denominator + param%matrixs(iz, ia)%sigma_bvl(ip, ig) * flux_forward%ngs(ig)%scalar(ir, ia) * geom%area(ir) * geom%height(ia) &
                                    &   * flux_adjoint%ngs(iig)%scalar(ir,ia) * param%matrixs(iz,ia)%chi_delay(ip,iig)
                            end do
                        end do
                    end do
                end do
                
                lambda(ip) = numerator / denominator
            end do
            
            write(FILES%REACTIVITY, *) ' '
            write(FILES%REACTIVITY, "(1x, A)")  'effective kinetics parameter by citation method:'
            write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))")  'partial beta   = ', beta
            write(FILES%REACTIVITY, "(1x, A, *(TR3, ES11.4))")  'partial lambda = ', lambda
        end if
        
        pk_param%beta = pk_parameter%beta
        pk_param%partial_beta = pk_parameter%partial_beta
        pk_param%partial_lambda = pk_parameter%partial_lambda
        pk_param%generation_time = pk_parameter%generation_time
        
    end subroutine Generate_pk_parameter
    
    !$
    !===============================================================================================
    ! generate reactivity and source term
    !===============================================================================================
    subroutine Generate_reactivity (tidx)
        
        integer, intent(in)  :: tidx
    
        ! loop index
        integer  :: is, ig, iig, ir, ia, im, ip, iz
        
        ! local variables
        real(KREAL)  :: denominator
        real(KREAL)  :: rho_t, rho_s, rho_f_nu                              ! reactivity by total, scatter, fission xsec
        real(KREAL), save  :: normal_source
        
        ! ----------------------------------------------------------------------
        ! obtains the denominator
        denominator = 0.0
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    do iig = 1, ns%state%ng
                        denominator = denominator + xsec%matrixs(iz,ia)%chi_steady(ig) * flux_adjoint%ngs(ig)%scalar(ir,ia)    &
                            &   * xsec%matrixs(iz,ia)%sigma_f_nu(iig) * shape_current%flux_scalar(ir,ia,iig)                   &
                            &   * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! ----------------------------------------------------------------------
        ! obtains reactivity
        rho_t = 0.0
        rho_s = 0.0
        rho_f_nu = 0.0
        
        ! by sigma_t
        do ig = 1, ns%state%ng 
            do is = 1, ns%deduce%direction 
                do ia = 1, ns%state%layer 
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        im = mat_info%loading(iz, ia)
                        rho_t = rho_t + flux_adjoint%ngs(ig)%angular(ir,ia,is) * shape_current%flux_angular(ir,ia,is,ig)          &
                            &   * (xsec%matrixs(iz,ia)%sigma_t(ig) - xsec_init%matrixs(iz,ia)%sigma_t(ig)) * quad%directions(is)%wmu     &
                            &   * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! by sigma_s
        do ig = 1, ns%state%ng 
            do ia = 1, ns%state%layer 
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    im = mat_info%loading(iz, ia)
                    do iig = 1, ns%state%ng
                        rho_s = rho_s + flux_adjoint%ngs(ig)%scalar(ir,ia) * shape_current%flux_scalar(ir,ia,iig)               &
                            &   * (xsec%matrixs(iz,ia)%sigma_s(iig,ig,1) - xsec_init%matrixs(iz,ia)%sigma_s(iig,ig,1))                 &
                            &   * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! @TODO--harmonic expansion
        if (ns%state%scat_order /= 0)  then
        end if
        
        ! by sigma_f_nu
        do ig = 1, ns%state%ng 
            do ia = 1, ns%state%layer 
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    im = mat_info%loading(iz, ia)
                    do iig = 1, ns%state%ng
                        rho_f_nu = rho_f_nu + flux_adjoint%ngs(ig)%scalar(ir,ia) * shape_current%flux_scalar(ir,ia,iig)     &
                            &   * (xsec%matrixs(iz,ia)%sigma_f_nu(iig) - xsec_init%matrixs(iz,ia)%sigma_f_nu(iig))                 &
                            &   * geom%area(ir) * geom%height(ia) * xsec_init%matrixs(iz,ia)%chi_steady(ig)
                    end do
                end do
            end do
        end do
        
        pk_parameter%rho = - (rho_t - rho_s - rho_f_nu) / denominator
        
        ! ----------------------------------------------------------------------
        ! obtains source
        pk_parameter%source = 0.0
        
        do ig = 1, ns%state%ng 
            do is = 1, ns%deduce%direction 
                do ia = 1, ns%state%layer 
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        pk_parameter%source = pk_parameter%source + flux_adjoint%ngs(ig)%angular(ir,ia,is) * Q_ext%matrixs(iz,ia)%intensity(ig)  &
                            &   * quad%directions(is)%wmu * geom%area(ir) * geom%height(ia)
                    end do
                end do
            end do
        end do
        
        ! NOTE--this line add for improved
        if (.FALSE.)  then
            do ig = 1, ns%state%ng 
                do is = 1, ns%deduce%direction 
                    do ia = 1, ns%state%layer 
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            pk_parameter%source = pk_parameter%source - shape_current%flux_angular(ir,ia,is,ig) * iter_adjoint%source(iz,ia,ig)  &
                                &   * quad%directions(is)%wmu * geom%area(ir) * geom%height(ia)
                        end do
                    end do
                end do
            end do
        end if
        
        ! for initial rho
        if (tidx == 0)  then
            normal_source = 0.0
            do is = 1, ns%deduce%direction
                do ig = 1, ns%state%ng
                    do ia = 1, ns%state%layer
                        do ir = 1, ns%state%nodal
                            iz = mesh%zone(ir)
                            associate(velocity => param%matrixs(iz,ia)%velocity(ig))
                            if (velocity > EPS_ZERO)  then
                                normal_source = normal_source + flux_forward%ngs(ig)%angular(ir,ia,is) * flux_adjoint%ngs(ig)%angular(ir,ia,is) * quad%directions(is)%wmu   &
                                    &   * geom%area(ir)*geom%height(ia) / velocity
                            else 
                                normal_source = normal_source + 0.0D0
                            end if 
                            end associate
                        end do
                    end do
                end do
            end do
            pk_parameter%rho_init = - (pk_parameter%source / normal_source) * pk_parameter%generation_time
        end if
        
        pk_parameter%rho = pk_parameter%rho + pk_parameter%rho_init
        pk_parameter%source = pk_parameter%source / normal_source

        pk_param%rho = pk_parameter%rho
        pk_param%q = pk_parameter%source
        
    end subroutine Generate_reactivity
    
end module process_pcqs
