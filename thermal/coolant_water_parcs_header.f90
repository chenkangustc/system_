!$
!===================================================================================================
!
!    this module is for water property class (water properties at 15.5 MPa, from PARCS code);
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          WaterProperty_PARCS
!
!===================================================================================================
module coolant_water_parcs_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : CoolantProperty
    
    implicit none 
    private
    public  :: WaterProperty_PARCS
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! --------------------------------------------------------------------------
    ! type for water proterty
    type, extends(CoolantProperty)  :: WaterProperty_PARCS
    contains
        procedure, public  :: set => Set_WaterProperty
        
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_enthalpy => Get_enthalpy_by_temperature
        procedure, public  :: get_temperature => Get_temperature_by_enthalpy
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivity_by_temperature
        procedure, public  :: get_viscosity => Get_viscosity_by_temperature
        
        procedure, public  :: get_nusselt => Get_Nusselt_number
    end type WaterProperty_PARCS
    
    ! private the real function name
    private :: Set_WaterProperty
    private :: Get_density_by_temperature, Get_enthalpy_by_temperature
    private :: Get_temperature_by_enthalpy, Get_capacity_by_temperature
    private :: Get_conductivity_by_temperature, Get_viscosity_by_temperature, Get_Nusselt_number
    
    real(KREAL), parameter  :: TUPPER = (340.0 + CKELVIN) * 1.05
    real(KREAL), parameter  :: TLOWER = (280.0 + CKELVIN) * 0.95
    
contains
    !$
    !===============================================================================================
    ! set information, do nothing
    !===============================================================================================    
    subroutine Set_WaterProperty (this, type, option)
        
        class(WaterProperty_PARCS), intent(in out)  :: this
        integer, intent(in)  :: type
        integer, intent(in)  :: option
    
    end subroutine Set_WaterProperty
    
    !$
    !===============================================================================================
    ! cubic polynomial for density as a function of temperature, max err=0.0692%
    !  15.5 Mpa,  TLOWER < T < TUPPER, rho in Kg/M^3, T in K
    !===============================================================================================    
    function Get_density_by_temperature(this, t_in)  result(density)
      
        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
        
        character(len=MAX_WORD_LEN)  :: info
        real(KREAL)  :: t
        
        if (t_in > TUPPER)  then
            t  = TUPPER
            write(unit=info, fmt=*) 'coolant water-PARCS, get density, upper exceed  ', t_in
            !  call a_warning%set (INFO_LIST_PROPERTY, info)
            !  call a_warning%print (FILES%TH_WARNING)
        else if (t_in < TLOWER)  then
            t = TLOWER
            write(unit=info, fmt=*) 'coolant water-PARCS, get density, lower exceed  ', t_in
            !  call a_warning%set (INFO_LIST_PROPERTY, info)
            !  call a_warning%print (FILES%TH_WARNING)
        else 
            t = t_in
        end if
        
        t = t - CKELVIN
        
        this%density = 5.9901166D+03+t*(-5.1618182D+01+t*(1.7541848D-01+t*(-2.0613054D-04)))
        density = this%density
      
    end function Get_density_by_temperature
        
    !$
    !===============================================================================================
    ! cubic polynomial for enthalpy as a function of temperature, max err=0.0306%
    !  15.5 Mpa,  TLOWER < T < TUPPER, output h in J/Kg, T in K
    !===============================================================================================
    function Get_enthalpy_by_temperature(this, t_in)  result(enthalpy)

        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: enthalpy
        
        real(KREAL)  :: t
        
        if (t_in > TUPPER)  then
            t  = TUPPER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, upper exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else if (t_in < TLOWER)  then
            t = TLOWER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, lower exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else 
            t = t_in
        end if
        
        t = t - CKELVIN
        
        this%enthalpy = -5.9301427D+03+t*(6.5488800D+01+t*(-2.1237562D-01+t*(2.4941725D-04)))
        
        ! change unit
        this%enthalpy = this%enthalpy * 1000.0D0 
        enthalpy = this%enthalpy
      
    end function Get_enthalpy_by_temperature
    
!!!    !$
!!!    !===============================================================================================
!!!    ! cubic polynomial for temperature as a function of enthalpy, max err=0.0055%
!!!    !  15.5 Mpa,  TLOWER < T < TUPPER, T in K, input h in J/Kg
!!!    !===============================================================================================
!!!    function Get_temperature_by_enthalpy(this, h_in)  result(temperature)
!!!
!!!        class(WaterProperty_PARCS), intent(in out)  :: this
!!!        real(KREAL), intent(in)  :: h_in
!!!        
!!!        real(KREAL)  :: temperature
!!!        real(KREAL)  :: hk
!!!        
!!!        hk = 0.001D0 * h_in
!!!        this%temperature = 1.4851739D+02+hk*(-1.2764991D-01+hk*(3.0781294D-04+hk*(-9.5429959D-08)))
!!!        
!!!        this%temperature = this%temperature + CKELVIN
!!!        temperature = this%temperature
!!!      
!!!    end function Get_temperature_by_enthalpy
    
    !$
    !===============================================================================================
    ! temperaure obtained from enthalpy [k]
    !===============================================================================================
    function Get_temperature_by_enthalpy (this, h_in)  result(temperature)

        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in) :: h_in                                         ! enthalpy, J/kg
        real(KREAL) :: temperature
        
        real(KREAL), parameter     :: delta = 1.0D-6                            ! this value is ok 
        integer, parameter         :: max_iter = 20

        integer     :: niter
        real(KREAL) :: step
        real(KREAL) :: tmelt
        real(KREAL) :: told, tnew
        
        ! if excess, equal max or min
        if (h_in <= this%get_enthalpy(TLOWER))  then
            this%temperature = TLOWER            
            temperature = this%temperature
            return
        else if (h_in >= this%get_enthalpy(TUPPER))  then
            this%temperature = TUPPER 
            temperature = this%temperature
            return
        end if
                
        niter = 0
        tmelt = TLOWER
        tnew  = tmelt + 20.0D0
        
        do
            niter = niter + 1
            told = tnew
            
            step =  (told-tmelt) / (this%get_enthalpy(told)-this%get_enthalpy(tmelt))
            tnew = told + (h_in-this%get_enthalpy (told)) * step
            
            ! normal exit
            if (100.0D0 * ABS(tnew-told)/told < delta) then
                exit
            end if
            
            ! abnormal exit
            if (niter > max_iter) then
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get temperature, iteration exceed')  
                !  call a_warning%print (FILES%TH_WARNING)
                exit
            end if
        end do
        
        this%temperature = tnew
        
        ! set value
        temperature = this%temperature

    end function Get_temperature_by_enthalpy
    
    !$
    !===============================================================================================
    ! quartic polynomial for heat capacity, max error=0.3053%
    !  15.5 Mpa,  TLOWER < T < TUPPER, output h in J/Kg-K, T in K
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: capacity
        
        real(KREAL)  :: t

        if (t_in > TUPPER)  then
            t  = TUPPER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, upper exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else if (t_in < TLOWER)  then
            t = TLOWER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, lower exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else 
            t = t_in
        end if
        
        t = t - CKELVIN
        
        this%capacity = 3.0455749D+03+t*(-4.0684599D+01+t*(2.0411250D-01+t*(-4.5526705D-04+t*(3.8115453D-07))))
        
        ! change unit
        this%capacity = this%capacity * 1000D0
        capacity = this%capacity
      
    end function Get_capacity_by_temperature
    
    !$
    !===============================================================================================
    ! cubic polynomial for thermal conductivity, max err=0.0204%
    !  15.5 Mpa,  TLOWER < T < TUPPER, k in w/m-K, T in K
    !===============================================================================================
    function Get_conductivity_by_temperature(this, t_in)  result(conductivity)
        
        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: conductivity
        
        real(KREAL)  :: t

        if (t_in > TUPPER)  then
            t  = TUPPER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, upper exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else if (t_in < TLOWER)  then
            t = TLOWER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, lower exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else 
            t = t_in
        end if
        
        t = t - CKELVIN
        
        this%conductivity = 8.9182016D-01+t*(-2.1996892D-03+t*(9.9347652D-06+t*(-2.0862471D-08)))
        conductivity = this%conductivity
        
    end function Get_conductivity_by_temperature
    
    !$
    !===============================================================================================
    ! cubic polynomial for dynamics viscosity, max err=0.0641%
    !  15.5 Mpa,  TLOWER < T < TUPPER, mu in Pa-sec, T in K
    !===============================================================================================
    function Get_viscosity_by_temperature(this, t_in)  result(viscosity)
    
        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: viscosity
        
        real(KREAL)  :: t
        
        if (t_in > TUPPER)  then
            t  = TUPPER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, upper exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else if (t_in < TLOWER)  then
            t = TLOWER
            !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant water-PARCS, get enthalpy, lower exceed')
            !  call a_warning%print (FILES%TH_WARNING)
        else 
            t = t_in
        end if
        
        t = t - CKELVIN
        
        this%viscosity = 9.0836878D-04+t*(-7.4542195D-06+t*(2.3658072D-08+t*(-2.6398601D-11)))
        viscosity = this%viscosity
      
    end function Get_viscosity_by_temperature

    !$
    !===============================================================================================
    ! nusselt number
    !===============================================================================================
    function Get_Nusselt_number (this, P2D, velocity, Dh, t_in)  result(nusselt)

        class(WaterProperty_PARCS), intent(in out)  :: this
        real(KREAL), intent(in) :: P2D                                          ! pitch to diameter ratio
        real(KREAL), intent(in) :: velocity                                     ! flow velocity
        real(KREAL), intent(in) :: Dh                                           ! characteristic length
        real(KREAL), intent(in) :: t_in                                         ! coolant temperature
        real(KREAL) :: nusselt                                                  
                                                                                
        ! local varibles                                                        
        real(KREAL)  :: rho                                                     ! density
        real(KREAL)  :: mu                                                      ! viscosity
        real(KREAL)  :: Cp                                                      ! heat capacity
        real(KREAL)  :: Cd                                                      ! heat conductivity
                                                                                
        real(KREAL)  :: Re                                                      ! Reynolds number, Re = (V*L/mu)
        real(KREAL)  :: Pr                                                      ! Prandlt number, Pr = (mu/a)
        real(KREAL)  :: Pe                                                      ! Peclet number, Pe = Pr*Re
        
        ! ----------------------------------------------------------------------
        rho = this%get_density (t_in)
        mu  = this%get_viscosity(t_in) / rho                                    ! transfer dynamics viscosity to kinetics viscosity
        Cp  = this%get_capacity (t_in)
        Cd  = this%get_conductivity (t_in)
        
        Re = velocity*Dh / mu
        Pr = mu*rho*Cp / Cd
        Pe = Re*Pr
        
        nusselt = 0.023D0*(Pr**0.4D0)*(Re**0.8D0)                               ! Dittus-Boelter relationship

    end function Get_Nusselt_number
    
end module coolant_water_parcs_header
