!$
!===================================================================================================
!
!    this module is for abstract class define for thermal property
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          CoolantProperty
!                               CladProperty
!                               GapProperty
!                               FuelProperty
!
!===================================================================================================
module abstract_property_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private
    public  :: CoolantProperty, CladProperty, GapProperty, FuelProperty
    
    ! --------------------------------------------------------------------------
    ! class for coolant property
    type, abstract  :: CoolantProperty
        integer          :: coolant_type                                        ! 
        real(KREAL)      :: density                                             ! density (Kg/M^3)
        real(KREAL)      :: enthalpy                                            ! enthalpy (J/Kg)
        real(KREAL)      :: temperature                                         ! temperature (K)
        real(KREAL)      :: capacity                                            ! heat cappcity (J/Kg-K)
        real(KREAL)      :: conductivity                                        ! thermal conductivity (W/m-K)
        real(KREAL)      :: viscosity                                           ! viscosity (Pa.s)
                         
        real(KREAL)      :: Nusselt_number                                      ! Nusselt_number        
    contains
        procedure(Set_CoolantProperty), deferred, public                      :: set
        procedure(Get_density_by_temperature_coolant), deferred, public       :: get_density
        procedure(Get_enthalpy_by_temperature_coolant), deferred, public      :: get_enthalpy
        procedure(Get_temperature_by_enthalpy_coolant), deferred, public      :: get_temperature
        procedure(Get_capacity_by_temperature_coolant), deferred, public      :: get_capacity
        procedure(Get_conductivity_by_temperature_coolant), deferred, public  :: get_conductivity
        procedure(Get_viscosity_by_temperature_coolant), deferred, public     :: get_viscosity
        procedure(Get_Nusselt_number), deferred, public                       :: get_nusselt
    end type CoolantProperty
    
    ! abstract interface for CoolantProperty class
    abstract interface        
        subroutine Set_CoolantProperty(this, type, option)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            integer, intent(in)  :: type
            integer, intent(in)  :: option
        end subroutine Set_CoolantProperty
        
        function Get_density_by_temperature_coolant(this, t_in) result(density)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: density
        end function Get_density_by_temperature_coolant
        
        function Get_enthalpy_by_temperature_coolant(this, t_in) result(enthalpy)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in) :: t_in
            real(KREAL) :: enthalpy
        end function Get_enthalpy_by_temperature_coolant
        
        function Get_temperature_by_enthalpy_coolant(this, h_in) result(temperature)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: h_in
            real(KREAL)  :: temperature
        end function Get_temperature_by_enthalpy_coolant
        
        function Get_capacity_by_temperature_coolant(this, t_in) result(capacity)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: capacity
        end function Get_capacity_by_temperature_coolant
        
        function Get_conductivity_by_temperature_coolant(this, t_in) result(conductivity)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: conductivity
        end function Get_conductivity_by_temperature_coolant
        
        function Get_viscosity_by_temperature_coolant(this, t_in) result(viscosity)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: viscosity
        end function Get_viscosity_by_temperature_coolant
        
        function Get_Nusselt_number(this, P2D, velocity, Dh, t_in) result(nusselt)
            import  :: CoolantProperty
            import  :: KREAL
            class(CoolantProperty), intent(in out)  :: this
            real(KREAL), intent(in) :: P2D     
            real(KREAL), intent(in) :: velocity
            real(KREAL), intent(in) :: Dh      
            real(KREAL), intent(in) :: t_in    
            real(KREAL) :: nusselt
        end function Get_Nusselt_number
    end interface

    ! --------------------------------------------------------------------------
    ! class for clad property
    type, abstract  :: CladProperty
        integer          :: clad_type                                           !
        real(KREAL)      :: density                                             ! density (Kg/M^3)
        real(KREAL)      :: capacity                                            ! heat cappcity (J/Kg-K)
        real(KREAL)      :: conductivity                                        ! thermal conductivity (W/m-K)
        real(KREAL)      :: expansion                                           ! thermal expansion in radial (1/K)
    contains
        procedure(Set_CladProperty), deferred, public                    :: set
        procedure(Get_density_by_temperature_clad), deferred, public     :: get_density
        procedure(Get_capacity_by_temperature_clad), deferred, public    :: get_capacity
        procedure(Get_conductivity_by_temperature_clad), deferred, public:: get_conductivity
        procedure(Get_expansion_by_temperature_clad), deferred, public   :: get_expansion
    end type CladProperty
    
    ! abstract interface for CladProperty class
    abstract interface
        subroutine Set_CladProperty(this, type)
            import  :: CladProperty
            import  :: KREAL
            class(CladProperty), intent(in out)  :: this
            integer, intent(in)  :: type
        end subroutine Set_CladProperty
        
        function Get_density_by_temperature_clad(this, t_in) result(density)
            import  :: CladProperty
            import  :: KREAL
            class(CladProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: density
        end function Get_density_by_temperature_clad
        
        function Get_capacity_by_temperature_clad(this, t_in) result(capacity)
            import  :: CladProperty
            import  :: KREAL
            class(CladProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: capacity
        end function Get_capacity_by_temperature_clad
        
        function Get_conductivity_by_temperature_clad(this, t_in) result(conductivity)
            import  :: CladProperty
            import  :: KREAL
            class(CladProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: conductivity
        end function Get_conductivity_by_temperature_clad
        
        function Get_expansion_by_temperature_clad(this, t_in)  result(expansion)
            import  :: CladProperty
            import  :: KREAL
            class(CladProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: expansion
        end function Get_expansion_by_temperature_clad
    end interface

    ! --------------------------------------------------------------------------
    ! class for gap property
    type, abstract  :: GapProperty
        integer          :: gas_type                                            !
        real(KREAL)      :: h_transfer                                          ! heat transfer coefficient(W/m^2-K)
    contains
        procedure(Set_GapProperty), deferred, public                    :: set
        procedure(Get_transfer_by_temperature_gap), deferred, public    :: get_transfer
    end type GapProperty
    
    ! abstract interface for GapProperty class
    abstract interface 
        subroutine Set_GapProperty(this, type, x_Xe, X_Kr)
            import  :: GapProperty
            import  :: KREAL        
            class(GapProperty), intent(in out)  :: this
            integer, intent(in)  :: type
            real(KREAL), intent(in), optional  :: x_Xe
            real(KREAL), intent(in), optional  :: x_Kr
        end subroutine Set_GapProperty
        
        function Get_transfer_by_temperature_gap(this, t_in, pellet, gap, is_inner) result(h_transfer)
            import  :: GapProperty
            import  :: KREAL        
            class(GapProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL), intent(in)  :: pellet
            real(KREAL), intent(in)  :: gap
            logical, intent(in)  :: is_inner
            real(KREAL)  :: h_transfer
        end function Get_transfer_by_temperature_gap
    end interface 
    
    ! --------------------------------------------------------------------------
    ! class for fuel property
    type, abstract  :: FuelProperty
        real(KREAL)      :: t_melting                                           ! melting temperature (K)
        real(KREAL)      :: mol_mass                                            ! mol mass (g/mol)
        
        integer          :: fuel_type                                           !
        real(KREAL)      :: density                                             ! density (Kg/M^3)
        real(KREAL)      :: capacity                                            ! heat cappcity (J/Kg-K)
        real(KREAL)      :: conductivity                                        ! thermal conductivity (W/m-K)
        real(KREAL)      :: expansion                                           ! thermal expansion in axial (1/K)
    contains
        procedure(Set_FuelProperty), deferred, public                    :: set
        procedure(Get_density_by_temperature_fuel), deferred, public     :: get_density
        procedure(Get_capacity_by_temperature_fuel), deferred, public    :: get_capacity
        procedure(Get_conductivty_by_temperature_fuel), deferred, public :: get_conductivity
        procedure(Get_expansion_by_temperature_fuel), deferred, public   :: get_expansion
    end type FuelProperty
    
    ! abstract interface for FuelProperty class
    abstract interface
        subroutine Set_FuelProperty(this, type, weight_Zr, x_Pu, x_O2M)
            import  :: FuelProperty
            import  :: KREAL        
            class(FuelProperty), intent(in out)  :: this
            integer, intent(in)  :: type
            real(KREAL), intent(in), optional  :: weight_Zr
            real(KREAL), intent(in), optional  :: x_Pu
            real(KREAL), intent(in), optional  :: x_O2M
        end subroutine Set_FuelProperty
        
        function Get_density_by_temperature_fuel(this, t_in) result(density)
            import  :: FuelProperty
            import  :: KREAL        
            class(FuelProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: density
        end function Get_density_by_temperature_fuel
        
        function Get_conductivty_by_temperature_fuel(this, t_in) result(conductivity)
            import  :: FuelProperty
            import  :: KREAL        
            class(FuelProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: conductivity
        end function Get_conductivty_by_temperature_fuel
        
        function Get_capacity_by_temperature_fuel(this, t_in) result(capacity)
            import  :: FuelProperty
            import  :: KREAL        
            class(FuelProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: capacity
        end function Get_capacity_by_temperature_fuel
        
        function Get_expansion_by_temperature_fuel(this, t_in) result(expansion)
            import  :: FuelProperty
            import  :: KREAL
            class(FuelProperty), intent(in out)  :: this
            real(KREAL), intent(in)  :: t_in
            real(KREAL)  :: expansion
        end function Get_expansion_by_temperature_fuel
    end interface 
    
end module abstract_property_header
