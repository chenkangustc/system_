!$
!===================================================================================================
!
!    this module is for property class ceramic with 95% theory density;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          FuelProperty_ceramic
!
!===================================================================================================
module fuel_caramic_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : FuelProperty
    
    implicit none 
    private
    public  :: FuelProperty_ceramic
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! type for fuel property
    type, extends(FuelProperty)  :: FuelProperty_ceramic
        private
        real(KREAL)      :: TD                                                  ! theory density
        real(KREAL)      :: x_Pu                                                ! weight fraction Pu/(Pu+U)
        real(KREAL)      :: x_O2M                                               ! oxygen excess (U.Pu)O2+x
    contains
        procedure, public  :: set => Set_FuelProperty
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivty_by_temperature
        procedure, public  :: get_expansion => Get_expansion_by_temperature
    end type FuelProperty_ceramic
    
    ! private the real function name   
    private :: Set_FuelProperty
    private :: Get_density_by_temperature, Get_capacity_by_temperature, Get_conductivty_by_temperature, Get_expansion_by_temperature
    
contains
    !$
    !===============================================================================================
    ! for UO2, recommand 7
    ! for MOX, recommand 8
    !===============================================================================================
    subroutine Set_FuelProperty (this, type, weight_Zr, x_Pu, x_O2M)
    
        class(FuelProperty_ceramic), intent(in out)  :: this
        integer, intent(in)  :: type
        real(KREAL), intent(in), optional      :: weight_Zr                     ! NOTE: not use in this type                           
        real(KREAL), intent(in), optional      :: x_Pu                           
        real(KREAL), intent(in), optional      :: x_O2M                         ! atomic weight input
        
        this%fuel_type = type
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! IAEA-UO2
            this%t_melting = 3120.0D0
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case(2)                                                                 ! IAEA-PuO2
            this%t_melting = 2663.0D0
            this%mol_mass  = 276.045D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 1.0D0
            this%x_O2M = 2.0D0
        
        case(3)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            this%t_melting = 3023.0D0
            this%mol_mass  = 271.2D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.2D0*276.045D0 / (0.2D0*276.045D0 + 0.8D0*270.3D0)
            this%x_O2M = 2.0D0
            
            if (PRESENT(x_Pu))  then
                this%x_Pu  = x_Pu*276.045D0 / (x_Pu*276.045D0 + (1.0D0-x_Pu)*270.3D0)
            end if
            if (PRESENT(x_O2M))  then
                this%x_O2M  = x_O2M
            end if
        
        case(4)                                                                 ! PARCS-UO2, same as 1
            this%t_melting = 3120.0D0
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case(5)                                                                 ! PARCS-MOX, same as 3
            this%t_melting = 3023.0D0
            this%mol_mass  = 271.2D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.2D0*276.045D0 / (0.2D0*276.045D0 + 0.8D0*270.3D0)
            this%x_O2M = 2.0D0
            
        case(6)                                                                 ! COBRA-MATPRO-UO2, same as 1
            this%t_melting = 3120.0D0                                             
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case(7)                                                                 ! COBRA-NEA-UO2, same as 1
            this%t_melting = 3120.0D0
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case(8)                                                                 ! RELAP5-UO2, same as 1
            this%t_melting = 3120.0D0
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case(9)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x, same as 3
            this%t_melting = 3023.0D0
            this%mol_mass  = 271.2D0
            this%TD        = 0.95D0
                        
            this%x_Pu  = 0.2D0
            this%x_O2M = 2.0D0
            
            if (PRESENT(x_Pu))  then
                this%x_Pu  = x_Pu*276.045D0 / (x_Pu*276.045D0 + (1.0D0-x_Pu)*270.3D0)
            end if
            if (PRESENT(x_O2M))  then
                this%x_O2M  = x_O2M
            end if

        case(10)                                                                ! OECD/NEA beam trip, same as 1
            this%t_melting = 3120.0D0
            this%mol_mass  = 270.3D0
            this%TD        = 0.95D0
            
            this%x_Pu  = 0.0D0
            this%x_O2M = 2.0D0
            
        case default
            call a_error%set (INFO_LIST_INPUT, 'fuel ceramic type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_FuelProperty
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_density_by_temperature (this, t_in)  result(density)
    
        class(FuelProperty_ceramic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
        
        real(KREAL)  :: t
        
        real(KREAL)  :: k_1, k_2, k_3, ED, k
        real(KREAL)  :: expansion_U, expansion_Pu, expansion
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! IAEA-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (273.0D0<=t .AND. t<=923.0D0)  then
                expansion_U = 0.99734D0+9.802D-6*t-2.705D-10*t**2+4.291D-13*t**3
            else if (t<=this%t_melting)  then
                expansion_U = 0.99672D0+1.179D-5*t-2.429D-9*t**2+1.219D-12*t**3
            end if

            this%density = 10960.0D0
            this%density = this%density * expansion_U**(-3)
            this%density = this%density * this%TD
            
        case(2)                                                                 ! IAEA-PuO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (273.0D0<=t .AND. t<=923.0D0)  then
                expansion_Pu = 0.99734D0+9.802D-6*t-2.705D-10*t**2+4.291D-13*t**3
            else if (t<=this%t_melting)  then
                expansion_Pu = 0.99672D0+1.179D-5*t-2.429D-9*t**2+1.219D-12*t**3
            end if

            this%density = 11460.0D0
            this%density = this%density * expansion_Pu**(-3)
            this%density = this%density * this%TD
        
        case(3)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 11080.0D0 / (1.0D0+2.04D-5*(t-273.15D0)+8.7D-9*(t-273.15D0)**2)
            this%density = this%density * this%TD
        
        case(4)                                                                 ! PARCS-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 10282.0D0
            this%density = 10412.0D0
            
        case(5)                                                                 ! PARCS-MOX
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 11080.0D0
            
        case(6)                                                                 ! COBRA-MATPRO-UO2, same as PARCS
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 10282.0D0
            
        case(7)                                                                 ! COBRA-NEA-UO2, same as PARCS
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 10282.0D0
            
        case(8)                                                                 ! RELAP5-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 10980.0D0
            this%density = this%density * this%TD
            
        case(9)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            k_1 = 1.0D-5
            k_2 = 3.0D-3
            k_3 = 4.0D-2
            ED  = 6.9D-20
            k   = 1.38D-23
            expansion_U = k_1*t - K_2 + K_3*EXP(-ED/(k*t))
            
            k_1 = 9.0D-6
            k_2 = 2.7D-3
            k_3 = 7.0D-2
            ED  = 7.0D-20
            k   = 1.38D-23
            expansion_Pu = k_1*t - K_2 + K_3*EXP(-ED/(k*t))
            
            expansion = this%x_Pu*expansion_Pu + (1.0D0-this%x_Pu)*expansion_U
            
            this%density = (1.0D0-this%x_Pu)*10980.0D0 + this%x_Pu*11460.0D0
            this%density = this%density * (1.0D0-3.0D0*expansion)
            this%density = this%density * this%TD
            
        case(10)                                                                ! OECD/NEA beam trip 
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 10354.0D0
            
        end select
        
        ! set value
        density = this%density
    
    end function Get_density_by_temperature

    !$
    !===============================================================================================
    ! functions for thermal properties of fuel regions thermal conductivity in w/m-K, t in K
    !===============================================================================================
    function Get_conductivty_by_temperature(this, t_in)  result(conductivity)
        
        class(FuelProperty_ceramic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: conductivity
        
        real(KREAL)  :: Cv_U, Cv_Pu
        real(KREAL)  :: degree
        real(KREAL)  :: t
        
        real(KREAL)  :: A, B, e_th, Cv, D, t_1, t_2    
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! IAEA-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t / 1000.0D0
            
            this%conductivity = 100.0D0/(7.5408D0+17.692D0*t+3.6142D0*t**2) + (6400.0D0/t**2.5D0)*EXP(-16.35D0/t)
            
        case(2)                                                                 ! IAEA-PuO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 8.441D0 - 7.445D-3*t + 2.236D-6*t**2
            this%conductivity = this%conductivity * (1.0D0 - 2.5D0*((1.0D0-this%TD)/this%TD))
        
        case(3)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            degree = 2.0D0 - this%x_O2M
            this%conductivity = 1.0D0/(1.528D0*SQRT(degree+0.00931D0)-0.1055D0+2.885D-4*t) + 76.38D-12*t**3
        
        case(4)                                                                 ! PARCS-UO2  
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 1.05D0+t*(0.0D0+t*(0.0D0+t*0.0D0))+2150.0D0/(t-73.15D0)
            
        case(5)                                                                 ! PARCS-MOX
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 1.05D0+t*(0.0D0+t*(0.0D0+t*0.0D0))+2150.0D0/(t-73.15D0)
            this%conductivity = this%conductivity * 0.90D0
            
        case(6)                                                                 ! COBRA-MATPRO-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%conductivity = COBRA_fuel_conductivity (t, real(1.0, KREAL))
            this%conductivity = this%conductivity * 6.23D3
            
        case(7)                                                                 ! COBRA-NEA-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%conductivity = COBRA_fuel_conductivity (t, real(-1.0, KREAL))
            this%conductivity = this%conductivity * 6.23D3
            
        case(8)                                                                 ! RELAP5-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            A    = 0.339D0
            B    = 0.06867D0
            e_th = 1.0D-5*t - 3.0D-3 + 4.0D-2*EXP(-6.9D-20/(1.38D-23*t))
            Cv   = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2)
            D    = 0.95D0
            
            if (t < 1364.0D0)  then
                t_1 = 6.5D0 - 0.00469D0*t
            else if (t > 1834.0D0)  then
                t_1 = -1.0D0
            else 
                t_1 = 0.10284D0 + (-0.1D0 - 0.10284D0) * (t-1364.0D0) / (1834.0D0-1364.0D0)
            end if
            
            if (t < 1800.0D0)  then
                t_2 = t
            else if (t > 2300.0D0)  then
                t_2 = 2050.0D0
            else 
                t_2 = 1800.0D0 + (2050.0D0-1800.0D0) * (t-1800.0D0) / (2300.0D0-1800.0D0)
            end if
            
            this%conductivity = (D/(1.0D0+t_1*(1.0D0-D))) * (Cv/(A+B*t_2)/(1.0D0+3*e_th)) + 5.2997D-3*t*EXP(-13358.0D0/t)*(1.0D0+0.169D0*((13358.0D0/t)+2.0D0)**2)
            
        case(9)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            A    = 0.339D0 + 12.6D0*ABS(2.0D0-this%x_O2M)
            B    = 0.06867D0 * (1.0D0 + 0.6238D0*this%x_Pu)
            e_th = 1.0D-5*t - 3.0D-3 + 4.0D-2*EXP(-6.9D-20/(1.38D-23*t))
            D    = this%TD

            Cv_U   = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2)
            Cv_Pu  = (347.4D0*571.000D0**2*EXP(571.000D0/t)) / (t**2*(EXP(571.000D0/t)-1)**2)
            Cv = Cv_Pu*this%x_Pu + Cv_U*(1-this%x_Pu)
            
            if (t < 1364.0D0)  then
                t_1 = 6.5D0 - 0.00469D0*t
            else if (t > 1834.0D0)  then
                t_1 = -1.0D0
            else 
                t_1 = 0.10284D0 + (-0.1D0 - 0.10284D0) * (t-1364.0D0) / (1834.0D0-1364.0D0)
            end if
            
            if (t < 1800.0D0)  then
                t_2 = t
            else if (t > 2300.0D0)  then
                t_2 = 2050.0D0
            else 
                t_2 = 1800.0D0 + (2050.0D0-1800.0D0) * (t-1800.0D0) / (2300.0D0-1800.0D0)
            end if
            
            this%conductivity = (D/(1.0D0+t_1*(1.0D0-D))) * (Cv/(A+B*t_2)/(1.0D0+3*e_th)) + 5.2997D-3*t*EXP(-13358.0D0/t)*(1.0D0+0.169D0*((13358.0D0/t)+2.0D0)**2)
            
        case(10)                                                                ! OECD/NEA beam trip
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 1.0/(0.042+2.71D-4*t) + 6.9D-11*t**3
           
        end select
        
        ! set value
        conductivity = this%conductivity
     
    end function Get_conductivty_by_temperature
    
    !$
    !===============================================================================================
    ! heat capacity in J/Kg-K, t in K
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(FuelProperty_ceramic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: capacity
        
        real(KREAL)  :: t
        real(KREAL)  :: rho
        real(KREAL)  :: Cp_U, Cp_Pu
        real(KREAL)  :: Cp_UO2, Cp_PuO2
        real(KREAL)  :: k_1, k_2, k_3, theta, ED, R
         
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! IAEA-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t / 1000.0D0
            
            this%capacity = 52.1743D0+87.951D0*t-82.2411D0*t**2+31.542D0*t**3-2.6334D0*t**4-0.71391D0*t**(-2)
            this%capacity = this%capacity * 1000.0D0/ this%mol_mass
            
        case(2)                                                                 ! IAEA-PuO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%capacity = -4.243D-6*t**2 + 2.366D-3*t + 293.1D0
        
        case(3)                                                                 ! IAEA-MOX (U0.8Pu0.2)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t / 1000.0D0
            
            Cp_U = 52.1743D0+87.951D0*t-82.2411D0*t**2+31.542D0*t**3-2.6334D0*t**4-0.71391D0*t**(-2)
            Cp_U = Cp_U * 1000.0D0/ this%mol_mass
            
            t = t_in
            Cp_Pu = -4.243D-6*t**2 + 2.366D-3*t + 293.1D0
        
            this%capacity = 0.2D0*Cp_Pu + 0.8D0*Cp_U
        
        case(4)                                                                 ! PARCS-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if

            rho = this%get_density (t)
            
            this%capacity = 162.3D0*rho+t*(0.3038D0*rho+t*(-2.391D-4*rho+t*6.404D-8*rho))
            this%capacity = this%capacity / rho
            
        case(5)                                                                 ! PARCS-MOX
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if

            rho = this%get_density (t)
            
            this%capacity = 162.3D0*rho+t*(0.3038D0*rho+t*(-2.391D-4*rho+t*6.404D-8*rho))
            this%capacity = this%capacity / rho
            
        case(6)                                                                 ! COBRA-MATPRO-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%capacity = COBRA_fuel_capacity (t, real(1.0, KREAL))
            this%capacity = this%capacity * 4186.3D0
            
        case(7)                                                                 ! COBRA-NEA-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%capacity = COBRA_fuel_capacity (t, real(-1.0, KREAL))
            this%capacity = this%capacity * 4186.3D0
            
        case(8)                                                                 ! RELAP5-UO2
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%capacity = (296.7D0*535.285D0**2*EXP(535.285D0/t)) / (t**2*(EXP(535.285D0/t)-1)**2) + 2.43D-2*t    &
                &   + ((2*8.745D7*1.577D5) / (2*8.3143D0*t**2)) * EXP(-1.577D5/(8.3143D0*t))
            
        case(9)                                                                 ! RELAP5-MOX (U1-c.Puc)O2+x
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            k_1   = 296.7D0
            k_2   = 2.43D-2
            k_3   = 8.745D7
            theta = 535.285D0
            ED    = 1.577D5
            R     = 8.3143D0
            Cp_U = (k_1*theta**2*EXP(theta/t)) / (t**2*(EXP(theta/t)-1)**2) + k_2*t    &
                &   + ((this%x_O2M*k_3*ED) / (2*R*t**2)) * EXP(-ED/(R*t))
                
            k_1   = 347.4D0
            k_2   = 3.95D-4
            k_3   = 3.860D7
            theta = 571.000D0
            ED    = 1.967D5
            R     = 8.3143D0
            Cp_Pu = (k_1*theta**2*EXP(theta/t)) / (t**2*(EXP(theta/t)-1)**2) + k_2*t    &
                &   + ((this%x_O2M*k_3*ED) / (2*R*t**2)) * EXP(-ED/(R*t))
                
            this%capacity = (this%x_Pu*276.045D0*Cp_Pu + (1-this%x_Pu)*270.3D0*Cp_U) / (this%x_Pu*276.045D0 + (1-this%x_Pu)*270.3D0)
            
        case(10)                                                                ! OECD/NEA beam trip
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel ceramic, get capacity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if

            Cp_UO2 = 81.825 + 0.78695*t - 1.1552E-3*t**2 + 9.9037E-7*t**3 - 5.1982E-10*t**4 + 1.5241E-13*t**5 - 1.7906E-17*t**6
            Cp_PuO2 = -4.9236E6/(t**2) + 240.89 + 0.32556*t - 3.5398E-4*t**2 + 1.512E-7*t**3 - 1.9707E-11*t**4
            this%capacity = (214.65*Cp_UO2 + 55.56*Cp_PuO2) / 270.21
            
        end select
        
        capacity = this%capacity
        
    end function Get_capacity_by_temperature
    
    !$
    !===============================================================================================
    ! thermal expansion of axial in 1/K, t in K
    !===============================================================================================
    function Get_expansion_by_temperature (this, t_in)  result(expansion)
        
        class(FuelProperty_ceramic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: expansion
    
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case (1)                                                                ! IAEA-UO2
            this%expansion = 1.0D-2
        
        case (2)                                                                ! IAEA-PuO2
            this%expansion = 1.0D-2
        
        case (3)                                                                ! IAEA-MOX (U0.8Pu0.2)O2+x
            this%expansion = 1.0D-2
        
        case (4)                                                                ! PARCS-UO2
            this%expansion = 1.0D-2
        
        case (5)                                                                ! PARCS-MOX
            this%expansion = 1.0D-2
        
        case (6)                                                                ! COBRA-MATPRO-UO2, same as PARCS
            this%expansion = 1.0D-2
        
        case (7)                                                                ! COBRA-NEA-UO2, same as PARCS
            this%expansion = 1.0D-2
        
        case (8)                                                                ! RELAP5-UO2
            this%expansion = 1.0D-2
        
        case (9)                                                                ! RELAP5-MOX (U1-c.Puc)O2+x
            this%expansion = 1.0D-2

        case (10)                                                               ! OECD/NEA beam trip 
            this%expansion = 1.0D-2
        end select
        
        expansion = this%expansion
        
    end function Get_expansion_by_temperature
    
    ! --------------------------------------------------------------------------
    ! following is from COBRA src
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !     input:
    !     t = temperature (f)
    !     flag  = flag to use matpro correlation, if negative >  0.0
    !           = flag to use nea    correlation, if negative <  0.0
    !     output:
    !     conductivity = supplied fuel thermal conductivity (btu/ft/s/f)
    ! ----------------------------------------------------------------------------------------------
    !     1 btu/(ft.s.f) = 6.23E3 W/(m.K)
    !     1 f = (9/5)*c + 32 = (9/5)*(K-273.15) + 32
    !===============================================================================================
    function COBRA_fuel_conductivity (t_in, flag)  result(conductivity)
        
        real(KREAL), intent(in)  :: t_in
        real(KREAL), intent(in)  :: flag
        real(KREAL)  :: conductivity

        real(KREAL)  :: c095
        real(KREAL)  :: c(7)
        real(KREAL)  :: tt
        real(KREAL)  :: w
        
        real(KREAL), parameter  :: tclo = 500.0D0-273.15D0
        real(KREAL), parameter  :: tcup = 3100.0D0-273.15D0
        real(KREAL), parameter  :: tklo = 293.0D0
        real(KREAL), parameter  :: tkup = 3100.0D0

        c095 = 0.01605D0
        c    = [32.0D0, 1.8D0, 273.15D0, 1.05D0, 2150.0D0, 73.15D0, 6226.48D0]
    
        ! ----------------------------------------------------------------------
        if (flag > 0.0D0)  then                                                 ! MATPRO
            tt = (t_in - c(1)) / c(2)
            if (tt < tclo)  then
                tt = tclo
            end if
            if (tt > tcup)  then
                tt = tcup
            end if
            
            w = 40.4D0 / (464.0D0 + tt)
            if (w < 0.0194D0)  then
                w = 0.0194D0
            end if
            w = w + 1.216D-4*EXP(1.867D-3*tt)
            
            conductivity = c095 * w
            
        else                                                                    ! NEA 
            tt = (t_in - c(1))/c(2)+ c(3)
            if (tt < tklo)  then
                tt = tklo
            end if
            if (tt > tkup)  then
                tt = tkup
            end if
            
            w = c(4) + c(5)/(tt - c(6))
            
            conductivity = w / c(7)
        end if
    
    end function COBRA_fuel_conductivity
    
    !$
    !===============================================================================================
    !     input:
    !     t = temperature (f)
    !     flag = flag to use matpro correlation, if negative >  0.0
    !          = flag to use nea    correlation, if negative <  0.0
    !     output:
    !     capacity = supplied fuel specific heat (btu/lb/f)
    ! ----------------------------------------------------------------------------------------------
    !     1 btu/(lb.f) = 4186.3 J/(Kg.K)
    !     1 f = (9/5)*c + 32 = (9/5)*(K-273.15) + 32
    !===============================================================================================
    function COBRA_fuel_capacity (t_in, flag)  result(capacity)
    
        real(KREAL), intent(in)  :: t_in
        real(KREAL), intent(in)  :: flag
        real(KREAL)  :: capacity
    
        real(KDOUBLE)   :: r, teinst, oxymet, ed, rk(3)
        real(KDOUBLE)   :: tt, w1, w2, w3, w
        real(KREAL)     :: c(8)
        
        real(KREAL), parameter  :: tklo = 298.0D0
        real(KREAL), parameter  :: tkup = 3100.0D0
        real(KREAL), parameter  :: uliq = 502.95D0/4184.0D0
        
        r       = 8.3143D0
        teinst  = 535.285D0
        oxymet  = 2.0D0 
        ed      = 1.577D5
        
        rk = [296.7D0, 2.43D-2, 8.745D7]
        c  = [32.0D0, 1.8D0, 273.15D0, 162.3D0, 0.3038D0, 2.391D-4, 6.404D-8, 4184.0D0]
        
        ! ----------------------------------------------------------------------
        if (flag > 0.0D0)  then                                                 ! MATPRO
            tt = (t_in - c(1))/c(2) + c(3)
            if (tt > tkup)  then
                capacity = uliq
                return
            end if
            
            if (tt < tklo)  then
                tt = tklo
            end if
            
            w1 = EXP(teinst/tt)
            w2 = rk(1)*teinst*teinst*w1
            w3 = tt*tt*(w1-1.0D0)*(w1-1.0D0)
            
            w = w2/w3 + rk(2)*tt
            w = w + (oxymet/2.0D0)*rk(3)*ed*EXP(-ed/(r*tt))/(r*tt*tt)
            
            capacity = 2.3889D-4 * w
            
        else                                                                    ! NEA 
            tt = (t_in - c(1))/c(2) + c(3)
            if (tt < tklo)  then
                tt = tklo
            end if
            if (tt > tkup)  then
                tt = tkup
            end if
            
            w = c(4) + ((c(7)*tt-c(6))*tt+c(5))*tt
            
            capacity = w / c(8)
        end if
    
    end function COBRA_fuel_capacity
    
end module fuel_caramic_header
