!$
!===================================================================================================
!
!   class for transient process parameter
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          DelayedNeutron
!                               ShapeFunction
!                               AmplitudeFunction
!                               PointKineticsParameter
!
!===================================================================================================
module transient_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
        
    implicit none 
    private
    public  :: DelayedNeutron
    public  :: ShapeFunction, AmplitudeFunction, PointKineticsParameter
    
    ! --------------------------------------------------------------------------
    ! type for delayed neutron precursor concentration
    type  DelayedNeutron
        real(KREAL), allocatable  :: concentration(:, :, :)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_PrecursorConcentration
        procedure, public  :: clean => Free_PrecursorConcentration
    end type DelayedNeutron
    
    ! --------------------------------------------------------------------------
    ! type for shape function
    type  ShapeFunction
        real(KREAL), allocatable  :: flux_angular(:, :, :, :)               ! angular flux
        real(KREAL), allocatable  :: flux_scalar(:, :, :)                   ! scalar flux
        real(KREAL), allocatable  :: precursor(:, :, :)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_ShapeFunction
        procedure, public  :: clean => Free_ShapeFunction
        procedure, private :: Equal_ShapeFunction
        generic, public    :: assignment(=) => Equal_ShapeFunction
    end type ShapeFunction
    
    ! type for amplitude function
    type  AmplitudeFunction
        real(KREAL)               :: flux
        real(KREAL)               :: flux_init                              ! initial value for flux amplitude 
        real(KREAL)               :: flux_condition                         ! flux normalization condition

        real(KREAL), allocatable  :: precursor(:)
        real(KREAL), allocatable  :: precursor_init(:)
        real(KREAL), allocatable  :: precursor_condition(:)
    contains
        procedure, public  :: alloc => Allocate_AmplitudeFunction
        procedure, public  :: clean => Free_AmplitudeFunction
        procedure, public  :: update => Update_AmplitudeFunction
    end type AmplitudeFunction
    
    ! type for point kinetics parameter
    type  PointKineticsParameter
        real(KREAL)               :: rho_init                               ! initial reactivity
        real(KREAL)               :: rho                                    ! current reactivity
        real(KREAL)               :: source                                 ! external source intensity
        real(KREAL)               :: generation_time                        ! neutron generation time
        real(KREAL)               :: beta                                   ! total delayed neutron fraction
        real(KREAL), allocatable  :: partial_beta(:)                        ! group delayed neutron fraction
        real(KREAL), allocatable  :: partial_lambda(:)                      ! group decay constant
    contains 
        procedure, public  :: alloc => Allocate_PointKineticsParameter
        procedure, public  :: clean => Free_PointKineticsParameter
        procedure, public  :: print => Print_PointKineticsParameter
    end type PointKineticsParameter
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_PrecursorConcentration (this)
        
        class(DelayedNeutron), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%concentration(nt%state%dg, ns%state%nodal, ns%state%layer), stat=i_allocate)
        
        this%concentration = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%concentration)
        
    end subroutine Allocate_PrecursorConcentration
    
    !$
    !===============================================================================================
    ! finalizer for class of DelayedNeutron
    !===============================================================================================
    subroutine Free_PrecursorConcentration (this)
    
        class(DelayedNeutron), intent(in out)  :: this
        
        if (allocated(this%concentration))          deallocate(this%concentration)
        
        this%memory = REAL_ZERO
        
    end subroutine Free_PrecursorConcentration 
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_ShapeFunction (this)
        
        class(ShapeFunction), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%flux_angular(ns%state%nodal, ns%state%layer, ns%deduce%direction, ns%state%ng), stat=i_allocate)
        allocate(this%flux_scalar(ns%state%nodal, ns%state%layer, ns%state%ng), stat=i_allocate)
        allocate(this%precursor(nt%state%dg, ns%state%nodal, ns%state%layer), stat=i_allocate)
        
        this%flux_angular  = REAL_ZERO
        this%flux_scalar   = REAL_ZERO
        this%precursor     = REAL_ZERO

        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%flux_angular)
        this%memory = this%memory + REAL_BYTE * SIZE(this%flux_scalar)
        this%memory = this%memory + REAL_BYTE * SIZE(this%precursor)
        
    end subroutine Allocate_ShapeFunction
    
    !$
    !===============================================================================================
    ! finalizer for class of ShapeFunction
    !===============================================================================================
    subroutine Free_ShapeFunction (this)
        
        class(ShapeFunction), intent(in out)  :: this
        
        if (allocated(this%flux_angular))       deallocate(this%flux_angular)
        if (allocated(this%flux_scalar))        deallocate(this%flux_scalar)
        if (allocated(this%precursor))          deallocate(this%precursor)
        
        this%memory = REAL_ZERO
        
    end subroutine Free_ShapeFunction
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Equal_ShapeFunction (left, right)
        
        class(ShapeFunction), intent(in out)  :: left
        type(ShapeFunction), intent(in)       :: right
        
        left%flux_angular = right%flux_angular
        left%flux_scalar = right%flux_scalar
        left%precursor = right%precursor
    
    end subroutine Equal_ShapeFunction
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_AmplitudeFunction (this)
        
        class(AmplitudeFunction), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%precursor(nt%state%dg), stat=i_allocate)
        allocate(this%precursor_init(nt%state%dg), stat=i_allocate)
        allocate(this%precursor_condition(nt%state%dg), stat=i_allocate)
        
        this%flux = REAL_ZERO
        this%flux_init = REAL_ZERO
        this%flux_condition = REAL_ZERO
        
        this%precursor  = REAL_ZERO
        this%precursor_init  = REAL_ZERO
        this%precursor_condition  = REAL_ZERO
        
    end subroutine Allocate_AmplitudeFunction
    
    !$
    !===============================================================================================
    ! finalizer for class of AmplitudeFunction
    !===============================================================================================
    subroutine Free_AmplitudeFunction (this)
        
        class(AmplitudeFunction), intent(in out)  :: this
        
        if (allocated(this%precursor))              deallocate(this%precursor)
        if (allocated(this%precursor_init))         deallocate(this%precursor_init)
        if (allocated(this%precursor_condition))    deallocate(this%precursor_condition)
        
    end subroutine Free_AmplitudeFunction
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_AmplitudeFunction (this, vector)
        
        class(AmplitudeFunction), intent(in out)  :: this
        real(KREAL), intent(in)  :: vector(:)
        
        ! this size not correct
        if (SIZE(vector) /= nt%state%dg+1)  then
        
        else 
            this%flux = vector(1)
            this%precursor = vector(2:nt%state%dg+1)
        end if
    
    end subroutine Update_AmplitudeFunction
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_PointKineticsParameter (this)
    
        class(PointKineticsParameter), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check alloated status first
        call this%clean ()
        
        allocate(this%partial_beta(nt%state%dg), stat=i_allocate)
        allocate(this%partial_lambda(nt%state%dg), stat=i_allocate)
        
        this%partial_beta   = REAL_ZERO
        this%partial_lambda = REAL_ZERO
        
    end subroutine Allocate_PointKineticsParameter
    
    !$
    !===============================================================================================
    ! finalizer for class of PointKineticsParameter
    !===============================================================================================
    subroutine Free_PointKineticsParameter (this)
        
        class(PointKineticsParameter), intent(in out)  :: this
        
        if (allocated(this%partial_beta))       deallocate(this%partial_beta)
        if (allocated(this%partial_lambda))     deallocate(this%partial_lambda)
        
    end subroutine Free_PointKineticsParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_PointKineticsParameter (this, unit_)
        
        class(PointKineticsParameter), intent(in out)  :: this
        integer, intent(in)  :: unit_
        
        integer  :: ip
        
        write(unit=unit_, fmt="(1x, A)")  ' '
        write(unit=unit_, fmt="(1x, A)")  '______________________________________________'
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'rho_init           : ', this%rho_init
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'rho                : ', this%rho
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'source             : ', this%source
        
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'generation_time    : ', this%generation_time
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'beta               : ', this%beta
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'partial_beta       : ', this%partial_beta
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'partial_lambda     : ', this%partial_lambda
    
    end subroutine Print_PointKineticsParameter
        
end module transient_header
