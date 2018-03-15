!$
!===================================================================================================
!
!   module for parameter object use in TH calculation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module th_global

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use coolant_water_parcs_header,     only : WaterProperty_PARCS
    use coolant_water_refprop_header,   only : WaterProperty_REFPROP
    use coolant_HLM_header,             only : CoolantProperty_HLM
    use clad_Zr_header,                 only : CladProperty_Zr
    use clad_steels_header,             only : CladProperty_steels
    use gap_gas_header,                 only : GapProperty_gas
    use fuel_metallic_header,           only : FuelProperty_metallic
    use fuel_caramic_header,            only : FuelProperty_ceramic

    use abstract_property_header,       only : CoolantProperty, CladProperty, GapProperty, FuelProperty
    
    use gth_geometry_header,            only : ThermalScale, ThermalGeometry, ThermalAssemblyGeometry
    use gth_thermal_header,             only : ThermalDesign, ThermalChannel, ThermalHotPoint
    use gth_power_header,               only : LinearPower
    
    implicit none
    public 
    
    ! --------------------------------------------------------------------------
    ! NOTE: type dependence
    !   property <-- geometry <-- thermal <-- power
    ! --------------------------------------------------------------------------
    
    ! coolant & fuel rod material
    type(WaterProperty_PARCS), pointer      :: a_coolant_parcs => NULL()
    type(WaterProperty_REFPROP), pointer    :: a_coolant_refprop => NULL()
    type(CoolantProperty_HLM), pointer      :: a_coolant_HLM => NULL()
    type(CladProperty_Zr), pointer          :: a_clad_Zr => NULL()
    type(CladProperty_steels), pointer      :: a_clad_steels => NULL()
    type(GapProperty_gas), pointer          :: a_gap_gas => NULL()
    type(FuelProperty_metallic), pointer    :: a_fuel_metallic => NULL()
    type(FuelProperty_ceramic), pointer     :: a_fuel_ceramic => NULL()
    
    class(CoolantProperty), pointer         :: a_coolant => NULL()
    class(CladProperty), pointer            :: a_clad => NULL()
    class(GapProperty), pointer             :: a_gap => NULL()
    class(FuelProperty), pointer            :: a_fuel => NULL()

    ! geometry information
    type(ThermalScale)                          :: nth
    type(ThermalGeometry)                       :: geom_th    
    type(ThermalAssemblyGeometry), allocatable  :: geom_assm(:)                 ! assembly by geometry
    
    ! thermal information
    type(ThermalDesign)                     :: design
    type(ThermalHotPoint)                   :: hot_point

    type(ThermalChannel)                    :: avg_channel
    type(ThermalChannel)                    :: hot_channel
    
    ! type for power generation
    type(LinearPower)                       :: th_power
    
end module th_global
