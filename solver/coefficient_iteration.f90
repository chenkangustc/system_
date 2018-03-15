!$
!===================================================================================================
!
!   module for global parameters of transport iteration
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module coefficient_iteration
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use coefficient_header_transverse,      only : SurfaceCoefficient, NodalCoefficient
    use coefficient_header_anisotropic,     only : AnisotropicSourceCoefficient
    use coefficient_header_accelarate,      only : LWExtrapolation, SORIterationInner, FSPScalingFactor
    
    implicit none 
    public
    
    ! --------------------------------------------------------------------------
    ! transverse coefficient 
    type(SurfaceCoefficient)                :: coeff_surface                    ! transverse coefficient of surface
    type(NodalCoefficient)                  :: coeff_nodal                      ! transverse coefficient of nodal
    
    ! anisotropic expansion coefficent
    type(AnisotropicSourceCoefficient)      :: coeff_source
    
    ! --------------------------------------------------------------------------
    ! non-coefficient parameter
    type(LWExtrapolation)                   :: accele_LW
    type(SORIterationInner)                 :: accele_SOR
    type(FSPScalingFactor)                  :: accele_Scaling
    
end module coefficient_iteration
