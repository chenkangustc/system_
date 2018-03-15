!$
!===================================================================================================
!
!    this module is for decay heat parameter definition;
!    use fitting method by PARCS & NESTLE
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          DecayHeat
!
!===================================================================================================
module decay_heat_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use geometry_header,        only : Meshing
    use material_header,        only : CrossSection
    
    implicit none 
    private
    public  :: DecayHeat
    
    ! --------------------------------------------------------------------------
    ! type for decay heat
    type  DecayHeat
        real(KREAL), public  :: alpha                                       ! sum of fractional powers
        real(KREAL), public, allocatable  :: partial_alpha(:)               ! decay heat precursor group fractional power
        real(KREAL), public, allocatable  :: partial_lambda(:)              ! decay heat precursor decay constants [sec^-1]
        
        real(KREAL), public, allocatable  :: precursor(:, :, :)             ! decay heat precursor concentrations at each node [J/cm^3]
        real(KREAL), public, allocatable  :: power(:, :)                    ! real power density 
    contains
        procedure, public  :: alloc => Alloc_DecayHeat
        procedure, public  :: clean => Free_DecayHeat
        procedure, public  :: steady => Steady_DecayHeat
        procedure, public  :: transient => Transient_DecayHeat
        procedure, public  :: get_power => Get_DecayHeat_power
    end type DecayHeat
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_DecayHeat (this)
        
        class(DecayHeat), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocatd status first
        call this%clean ()
        
        allocate(this%partial_alpha(nt%state%hg), stat=i_allocate)
        allocate(this%partial_lambda(nt%state%hg), stat=i_allocate)
        
        allocate(this%precursor(nt%state%hg, ns%state%nodal, ns%state%layer), stat=i_allocate)
        allocate(this%power(ns%state%nodal, ns%state%layer), stat=i_allocate)
        
        this%partial_alpha  = REAL_ZERO
        this%partial_lambda = REAL_ZERO
        
        this%precursor = REAL_ZERO
        this%power = REAL_ZERO
        
        ! set default
        this%partial_alpha =  [2.35402E-2, 1.89077E-2, 1.39236E-2, 6.90315E-3, 3.56888E-3, 3.31633E-3]
        this%partial_lambda = [1.05345E-1, 8.37149E-3, 5.20337E-4, 4.73479E-5, 3.28153E-6, 1.17537E-11]
        this%alpha = SUM(this%partial_alpha)
        
    end subroutine Alloc_DecayHeat
    
    !$
    !===============================================================================================
    ! finalizer for class of DecayHeat
    !===============================================================================================
    subroutine Free_DecayHeat (this)
    
        class(DecayHeat), intent(in out)  :: this
        
        if (allocated(this%partial_alpha))          deallocate(this%partial_alpha)
        if (allocated(this%partial_lambda))         deallocate(this%partial_lambda)
        
        if (allocated(this%precursor))              deallocate(this%precursor)
        if (allocated(this%power))                  deallocate(this%power)
           
    end subroutine Free_DecayHeat
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Steady_DecayHeat (this, mesh, xsec, flux)
        
        class(DecayHeat), intent(in out)  :: this
        type(Meshing), intent(in)         :: mesh
        type(CrossSection), intent(in)    :: xsec
        real(KREAL), intent(in)       :: flux(:, :, :)
        
        real(KREAL)  :: rate
        integer          :: ih, ia, ir, ig, iz
        
        do ih = 1, nt%state%hg
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    
                    rate = 0.0
                    do ig = 1, ns%state%ng
                        rate = rate + xsec%matrixs(iz,ia)%sigma_f_kappa(ig) * flux(ir,ia,ig)
                    end do
                    
                    this%precursor(ih, ir, ia) = rate * this%partial_alpha(ih) / this%partial_lambda(ih)
                end do
            end do
        end do
    
    end subroutine Steady_DecayHeat
    
    !$
    !===============================================================================================
    ! update decay heat precursor concentrations based on new flux
    !===============================================================================================
    subroutine Transient_DecayHeat (this, mesh, xsec, flux, step_length)
        
        class(DecayHeat), intent(in out)  :: this
        type(Meshing), intent(in)         :: mesh
        type(CrossSection), intent(in)    :: xsec
        real(KREAL), intent(in)       :: flux(:, :, :)
        real(KREAL), intent(in)       :: step_length
        
        real(KREAL)  :: base, factor, rate
        integer          :: ia, ir, ih, iz, ig
        
        do ih = 1, nt%state%hg
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    
                    rate = 0.0
                    do ig = 1, ns%state%ng
                        rate = rate + xsec%matrixs(iz,ia)%sigma_f_kappa(ig) * flux(ir,ia,ig)
                    end do
                    
                    base = this%precursor(ih,ir,ia) * EXP(-this%partial_lambda(ih)*step_length)
                    factor = (1.0-EXP(-this%partial_lambda(ih)*step_length)) * (this%partial_alpha(ih)/this%partial_lambda(ih))
                    
                    this%precursor(ih,ir,ia) = base + factor * rate
                end do
            end do
        end do
        
    end subroutine Transient_DecayHeat
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Get_DecayHeat_power (this, power)
        
        class(DecayHeat), intent(in out)  :: this
        real(KREAL), intent(in out)   :: power(:, :)
        
        real(KREAL)  :: fract
        integer          :: ih, ir, ia
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                fract = 0.0
                do ih = 1, nt%state%hg
                    fract = fract + this%partial_lambda(ih) * this%precursor(ih,ir,ia)
                end do
                
                power(ir, ia) = (1.0-this%alpha) * power(ir, ia) + fract
            end do
        end do
    
    end subroutine Get_DecayHeat_power

end module decay_heat_header
