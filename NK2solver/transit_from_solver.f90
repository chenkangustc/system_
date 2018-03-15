!$
!===================================================================================================
!
!   dumping calculation result into distribution matrix, generate timelist parameter sepecific
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Transit_flux
!                               Transit_precursor
!
!   Public type lists:          No
!
!===================================================================================================
module transit_from_solver
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global

    implicit none 
    private
    public  :: Transit_flux, Transit_precursor
    
contains
    !$
    !===============================================================================================
    ! dumming flux, reaction rate, and power
    !===============================================================================================
    subroutine Transit_flux (is_eigen, is_adjoint, is_transient, tidx, ctime)
    
        ! intent parameters
        logical, intent(in)            :: is_eigen
        logical, intent(in)            :: is_adjoint
        logical, intent(in), optional  :: is_transient
        integer, intent(in), optional  :: tidx
        real(KREAL), intent(in), optional  :: ctime
        
        ! loop index
        integer  :: ia, ir, im, il, ig, is, iz
        
        ! local variables
        real(KREAL)  :: power                                               ! power level for FSP
        real(KREAL)  :: matrix(ns%state%nodal, ns%state%layer)              ! transfer into distribution
        real(KREAL)  :: zone_(ns%state%zone, ns%state%layer)
        
        ! ----------------------------------------------------------------------
        ! obtains absolute power level
        power = 0.0
        do is = 1, ns%deduce%direction
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    il = (ia-1) * ns%state%nodal + ir
                    do ig = 1, ns%state%ng
                        power = power + iter_flux%dist%nodal(ir, ia, is, ig) * quad%directions(is)%wmu     &
                            &   * xsec_iter%matrixs(iz,ia)%sigma_f_kappa(ig)*geom%area(ir)*geom%height(ia)
                    end do
                end do
            end do
        end do

        ! ----------------------------------------------------------------------
        ! get normal factor
        timelist%initial_power = ns%flag%power_level
        
        if (.NOT. nt%flag%is_transient .and. .NOT. nt%flag%is_perturb)  then
            if (is_adjoint )  then
!                timelist%power = ns%flag%power_level
!                timelist%normal_factor = ns%flag%power_level / power
            else
                if (xsec_iter%is_fission ())  then
                    timelist%power = ns%flag%power_level
                    timelist%normal_factor = ns%flag%power_level / power
                else
                    timelist%power = power
                    timelist%normal_factor = 1.0
                end if
            end if
        
        ! for perturbation
        else if (nt%flag%is_perturb)  then
            if (is_adjoint )  then
!                timelist%power = power
!                timelist%normal_factor = 1.0
            else
                timelist%power = power
                timelist%normal_factor = 1.0
            end if
        
        ! for transient
        else 
            if (present(is_transient) .and. is_transient)  then
                if (tidx == 0)  then
                    if (is_eigen )  then
                        timelist%power = ns%flag%power_level
                        timelist%normal_factor = ns%flag%power_level / power
                    else
                        timelist%power = power
                        timelist%normal_factor = 1.0
                    end if
                else 
                    timelist%power = power
                    timelist%normal_factor = 1.0
                end if
            else if (is_adjoint)  then
!                timelist%power = ns%flag%power_level
!                timelist%normal_factor = ns%flag%power_level / power
            else
                timelist%power = ns%flag%power_level
                timelist%normal_factor = ns%flag%power_level / power
            end if
        end if
        
        ! ----------------------------------------------------------------------
        ! transfer flux into GroupsFlux type
        if (.NOT. is_adjoint)  then
            do ig = 1, ns%state%ng
                flux_forward%ngs(ig)%angular(:, :, :) = iter_flux%dist%nodal(:, :, :, ig) * timelist%normal_factor
            end do
            call flux_forward%set_scalar (quad)
            if (.NOT. ns%misc%is_angwise)  then
                call flux_forward%sca2ang ()
            end if 
        else
            do ig = 1, ns%state%ng
                flux_adjoint%ngs(ig)%angular(:, :, :) = iter_flux%dist%nodal(:, :, :, ig) * timelist%normal_factor
            end do
            call flux_adjoint%set_scalar (quad)
            if (.NOT. ns%misc%is_angwise)  then
                call flux_adjoint%sca2ang ()
            end if 
        end if
        
        ! ----------------------------------------------------------------------
        ! transfer flux into DistributionParameter type
        ! for flux
        matrix = 0.0
        if (.NOT. is_adjoint)  then
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    do ig = 1, ns%state%ng
                        matrix(ir, ia) = matrix(ir, ia) + flux_forward%ngs(ig)%scalar(ir, ia)
                    end do
                end do
            end do
        else
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    do ig = 1, ns%state%ng
                        matrix(ir, ia) = matrix(ir, ia) + flux_adjoint%ngs(ig)%scalar(ir, ia)
                    end do
                end do
            end do
        end if
        call dist_flux%set (matrix)
        
        ! for power
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
        
        ! write xxx
        if (.False.)  then
            do ig = 1, ns%state%ng
                call dist_fission_rate%set (flux_forward%ngs(ig)%scalar)
                call dist_fission_rate%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=zone_)
                do iz = 1, ns%state%zone
                    write(unit=600+ig, fmt="(1x, *(TR3, ES13.6))")  zone_(iz, :)
                end do
            end do
        end if
        
        ! for fission rate
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
          
    end subroutine Transit_flux
 
    !$
    !===============================================================================================
    ! dumming precursor into matrix distribution and generate absolute integral value
    !===============================================================================================
    subroutine Transit_precursor ()

        integer  :: ir, ia, ip
        
        ! precursor per nodal per layer
        do ip = 1, nt%state%dg
            dist_dnps(ip)%matrix(:, :) = dnp_solver%Rprecursor(ip, :, :)
        end do
        
        ! get integral information
        timelist%precursor = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ip = 1, nt%state%dg
                    timelist%precursor(ip) = timelist%precursor(ip) + dist_dnps(ip)%matrix(ir,ia)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
    
    end subroutine Transit_precursor
    
end module transit_from_solver
