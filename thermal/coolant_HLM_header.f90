!$
!===================================================================================================
!
!   module for thermal property of coolant HLM
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          CoolantProperty_HLM
!
!===================================================================================================
module coolant_HLM_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : CoolantProperty
    
    implicit none
    private 
    public  :: CoolantProperty_HLM
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! --------------------------------------------------------------------------
    ! type for working flow
    type, extends(CoolantProperty)  :: CoolantProperty_HLM
        private
        real(KREAL)  :: t_melting                                               ! melting temperature (K)
        real(KREAL)  :: t_boiling                                               ! boiling temperature (K)
        real(KREAL)  :: mol_mass                                                ! mol mass (g/mol)       
    contains
        procedure, public  :: set => Set_CoolantProperty
        
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_enthalpy => Get_enthalpy_by_temperature
        procedure, public  :: get_temperature => Get_temperature_by_enthalpy
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivity_by_temperature
        procedure, public  :: get_viscosity => Get_viscosity_by_temperature
        
        procedure, public  :: get_nusselt => Get_Nusselt_number
    end type CoolantProperty_HLM
    
    ! --------------------------------------------------------------------------
    ! private real function name
    private :: Set_CoolantProperty
    private :: Get_density_by_temperature, Get_enthalpy_by_temperature
    private :: Get_temperature_by_enthalpy, Get_capacity_by_temperature
    private :: Get_conductivity_by_temperature, Get_viscosity_by_temperature, Get_Nusselt_number

contains
    !$
    !===============================================================================================
    ! set working flow basic information
    !===============================================================================================
    subroutine Set_CoolantProperty (this, type, option)
    
        class(CoolantProperty_HLM), intent(in out)  :: this
        integer, intent(in)  :: type
        integer, intent(in)  :: option
        
        this%coolant_type = type
        
        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            this%t_melting  = 371.0D0
            this%t_boiling  = 1155.0D0
            this%mol_mass   = 23.00D0
            
        case(2)                                                                 ! IAEA-Lead
            this%t_melting  = 600.6D0
            this%t_boiling  = 2021.0D0
            this%mol_mass   = 207.00D0
            
        case(3)                                                                 ! IAEA-LBE
            this%t_melting  = 398.0D0
            this%t_boiling  = 1927.0D0
            this%mol_mass   = 208.00D0
        
        case(4)                                                                 ! OECD/NEA beam trip 
            this%t_melting  = 397.7D0
            this%t_boiling  = 1927.0D0
            this%mol_mass   = 208.00D0
        
        case default
            call a_error%set (INFO_LIST_INPUT, 'HLM type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_CoolantProperty
    
    !$
    !===============================================================================================
    ! obtain density from temperature [kg/m^3]
    !===============================================================================================
    function Get_density_by_temperature (this, t_in)  result(density)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in
        real(KREAL) :: density
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
        
            this%density = 1014.0D0-0.235D0*t
            
        case(2)                                                                 ! IAEA-Lead
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, uppper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%density = 11441.0D0-1.2795D0*t  
            
        case(3)                                                                 ! IAEA-LBE
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%density = 11065.0D0-1.293D0*t
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%density = 11112.0D0-1.375D0*t
            
        end select
            
        ! set value
        density = this%density
        
    end function Get_density_by_temperature

    !$
    !===============================================================================================
    ! specific enthalpy obtained from temperature [J/kg]
    !===============================================================================================
    function Get_enthalpy_by_temperature (this, t_in)  result(enthalpy)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in
        real(KREAL) :: enthalpy
        
        real(KREAL)  :: t
        real(KREAL)  :: tmelt
        
        tmelt   = this%t_melting
        
        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
        
            this%enthalpy = 164.8D0*(t-tmelt)-1.97D-2*(t**2-tmelt**2)+4.167D-4*(t**3-tmelt**3) + 4.56D5*((t**(-1)-tmelt**(-1)))
            
        case(2)                                                                 ! IAEA-Lead
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
        
            this%enthalpy = 176.2D0*(t-tmelt)-2.4615D-2*(t**2-tmelt**2)+5.147D-6*(t**3-tmelt**3) + 1.524D6*((t**(-1)-tmelt**(-1)))
            
        case(3)                                                                 ! IAEA-LBE
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%enthalpy = 164.8D0*(t-tmelt)-1.97D-2*(t**2-tmelt**2)+4.167D-6*(t**3-tmelt**3) + 4.56D5*((t**(-1)-tmelt**(-1)))
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get enthalpy, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
!            this%enthalpy = 164.8D0*(t-tmelt)-1.97D-2*(t**2-tmelt**2)+4.167D-6*(t**3-tmelt**3) + 4.56D5*((t**(-1)-tmelt**(-1)))
            this%enthalpy = 10010.0+33.1025*(t-tmelt)-5.6628E-3*(t-tmelt)**2+1.4823E-6*(t-tmelt)**3
            this%enthalpy = 1000.0 * this%enthalpy / this%mol_mass
            
        end select
        
        ! set value
        enthalpy = this%enthalpy

    end function Get_enthalpy_by_temperature

    !$
    !===============================================================================================
    ! temperaure obtained from enthalpy [k]
    !===============================================================================================
    function Get_temperature_by_enthalpy (this, h_in)  result(temperature)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: h_in                                         ! enthalpy, J/kg
        real(KREAL) :: temperature
        
        real(KREAL), parameter     :: delta = 1.0D-5                            ! this value is ok 
        integer, parameter         :: max_iter = 6

        integer         :: niter
        real(KREAL) :: step
        real(KREAL) :: tmelt
        real(KREAL) :: told, tnew
        
        ! if excess, equal max or min
        if (h_in <= this%get_enthalpy(this%t_melting))  then
            this%temperature = this%t_melting            
            temperature = this%temperature
            return
        else if (h_in >= this%get_enthalpy(this%t_boiling))  then
            this%temperature = this%t_boiling            
            temperature = this%temperature
            return
        end if
                
        niter = 0
        tmelt = this%t_melting
        tnew  = tmelt + 300.0D0
        
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
    ! heat capacity obtained from temperature [J/(kg*K)]
    !===============================================================================================
    function Get_capacity_by_temperature (this, t_in)  result(capacity)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in
        real(KREAL) :: capacity
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%capacity = -3.001D6*t**(-2)+1658.0D0-0.8479D0*t+4.454D-4*t**2
        
        case(2)                                                                 ! IAEA-Lead
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%capacity = 175.1D0-4.961D-2*t+1.985D-5*t**2-2.099D-9*t**3-1.524D6*t**(-2)
            
        case(3)                                                                 ! IAEA-LBE
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%capacity = 164.8D0-3.94D-2*t+1.25D-5*t**2-4.56D5*t**(-2)
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%capacity = 146.5D0
            
        end select
        
        ! set value 
        capacity = this%capacity

    end function Get_capacity_by_temperature

    !$
    !===============================================================================================
    ! heat conductivity coefficient [W/(m*K)]
    !===============================================================================================
    function Get_conductivity_by_temperature (this, t_in)  result(conductivity)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in
        real(KREAL) :: conductivity
        
        real(KREAL)  :: t
         
        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%conductivity = 110.0D0-6.48D-2*t+1.16D-5*t**2
        
        case(2)                                                                 ! IAEA-Lead
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%conductivity = 9.2D0+0.011D0*t
            
        case(3)                                                                 ! IAEA-LBE
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%conductivity = 3.284D0+1.617D-2*t-2.305D-6*t**2
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%conductivity = 3.9021D0+0.0123D0*t
            
        end select
        
        ! set value
        conductivity = this%conductivity
        
    end function Get_conductivity_by_temperature

    !$
    !===============================================================================================
    ! dynamics viscosity obtained from temperature [Pa.s]
    !===============================================================================================
    function Get_viscosity_by_temperature (this, t_in)  result(viscosity)

        class(CoolantProperty_HLM), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in
        real(KREAL) :: viscosity
        
        real(KREAL)  :: t
        
        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case(1)                                                                 ! IAEA-Na
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%viscosity = EXP(556.835D0/t - 6.4406D0) / t**0.3958D0
        
        case(2)                                                                 ! IAEA-Lead
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%viscosity = 4.55D-4 * EXP(1069.0D0/t)
            
        case(3)                                                                 ! IAEA-LBE
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%viscosity = 4.94D-4 * EXP(754.1D0/t)
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in < this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in > this%t_boiling)  then
                t = this%t_boiling
                !  call a_warning%set (INFO_LIST_PROPERTY, 'coolant HLM, get viscosity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else
                t = t_in
            end if
            
            this%viscosity = 4.94D-4 * EXP(754.1D0/t)
            
        end select
        
        ! set value 
        viscosity = this%viscosity

    end function Get_viscosity_by_temperature

    !$
    !===============================================================================================
    ! nusselt number
    !===============================================================================================
    function Get_Nusselt_number (this, P2D, velocity, Dh, t_in)  result(nusselt)

        class(CoolantProperty_HLM), intent(in out)  :: this
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
        
        integer          :: correlation = 2                                     ! select a correlation 
        
        ! ----------------------------------------------------------------------
        rho = this%get_density (t_in)
        mu  = this%get_viscosity(t_in) / rho                                    ! transfer dynamics viscosity to kinetics viscosity
        Cp  = this%get_capacity (t_in)
        Cd  = this%get_conductivity (t_in)
        
        Re = velocity*Dh / mu
        Pr = mu*rho*Cp / Cd
        Pe = Re*Pr
        
        ! ----------------------------------------------------------------------
        select case (this%coolant_type)
        case (1, 2, 3)                                                          ! IAEA-Na, Lead & LBE
            ! 1--IAEA book; 2--PSI FAST project
            if (correlation == 1)  then
                if (1.2D0<=P2D .AND. P2D<=2.0D0 .AND. 1.0D0<=Pe .AND. Pe<=4000.0D0) then
                    nusselt = 7.55D0*P2D-20.0D0*P2D**(-13)+3.67D0*Pe**(0.56D0+0.19D0*P2D)/(90.0D0*P2D**2)
                else
                    call a_error%set (INFO_LIST_INPUT, 'nusselt number not obtained, P/D or Pe is not in validity range')  
                    call a_error%print (FILES%MAIN)
                end if
            else
                if (1.1D0<=P2D .AND. P2D<=1.95D0 .AND. 30.0D0<=Pe .AND. Pe<=5000.0D0) then
                    nusselt = 0.047D0 * (1.0D0 - EXP(-3.8D0*(P2D-1.0D0))) * (250.0D0 + Pe**0.77D0)
                else
                    call a_error%set (INFO_LIST_INPUT, 'nusselt number not obtained, P/D or Pe is not in validity range')  
                    call a_error%print (FILES%MAIN)
                end if
            end if
            
        case (4)
            Pe = velocity*rho * (Dh*Cp/Cd)
            nusselt = 4.0+0.16*P2D**5+0.33*(P2D**3.8)*((Pe/100.0)**0.86)
            
        end select
        
        ! set value
        this%Nusselt_number = nusselt

    end function Get_Nusselt_number

end module coolant_HLM_header
