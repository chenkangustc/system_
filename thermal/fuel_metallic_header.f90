!$
!===================================================================================================
!
!    this module is for property class metallic;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          FuelProperty_metallic
!
!===================================================================================================
module fuel_metallic_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : FuelProperty
    
    implicit none 
    private
    public  :: FuelProperty_metallic
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! type for fuel property
    type, extends(FuelProperty)  :: FuelProperty_metallic
        private       
        real(KREAL)      :: weight_Zr                                           ! weight fraction of Zr
        real(KREAL)      :: weight_Pu                                           ! weight fraction of Pu
        real(KREAL)      :: x_Zr                                                ! atomic fraction of Zr
        real(KREAL)      :: x_Pu                                                ! atomic fraction of Pu
    contains
        procedure, public  :: set => Set_FuelProperty
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivty_by_temperature
        procedure, public  :: get_expansion => Get_expansion_by_temperature
    end type FuelProperty_metallic
    
    ! private the real function name   
    private :: Set_FuelProperty
    private :: Get_density_by_temperature, Get_capacity_by_temperature, Get_conductivty_by_temperature, Get_expansion_by_temperature
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_FuelProperty (this, type, weight_Zr, x_Pu, x_O2M)
    
        class(FuelProperty_metallic), intent(in out)  :: this
        integer, intent(in)  :: type
        real(KREAL), intent(in), optional      :: weight_Zr
        real(KREAL), intent(in), optional      :: x_Pu                          ! NOTE: not use for this type
        real(KREAL), intent(in), optional      :: x_O2M                         ! NOTE: not use for this type
        
        this%fuel_type = type
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! U
            this%t_melting = 1405.0D0
            this%mol_mass  = 238.0D0
            
        case(2)                                                                 ! Pu
            this%t_melting = 913.0D0
            this%mol_mass  = 244.0D0
        
        case(3)                                                                 ! Pu-Zr
            this%weight_Zr = 0.618D0
            if (PRESENT(weight_Zr))  then
                this%weight_Zr = weight_Zr
            end if
            
            this%weight_Pu = 1.0D0 - this%weight_Zr
            
            this%x_Zr = (this%weight_Zr*244.0D0) / (91.22D0 + this%weight_Zr*(244.0D0-91.22D0))
            this%x_Pu = 1.0D0 - this%x_Zr
            
            this%t_melting = this%x_Zr*2128.0D0 + this%x_Pu*913.0D0
            this%mol_mass  = this%x_Zr*91.22D0 + this%x_Pu*244.0D0
            
        case default
            call a_error%set (INFO_LIST_INPUT, 'fuel metallic type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_FuelProperty
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_density_by_temperature (this, t_in)  result(density)
    
        class(FuelProperty_metallic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
        
        real(KREAL)  :: t
        real(KREAL)  :: key_Zr(15), value_Zr(15)
        real(KREAL)  :: key_Pu(22), value_Pu(22)
        real(KREAL)  :: density_Zr, density_Pu
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! U
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (273.0D0<=t .AND. t<=942.0D0)  then
                this%density = 19.36D3 - 1.03347D0*t
            else if (t<=1049.0D0)  then
                this%density = 19.092D3 - 0.9807D0*t
            else
                this%density = 18.447D3 - 0.5166D0*t
            end if
            
        case(2)                                                                 ! Pu
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            ! step-wise approximation
            if (t <= 395.0D0)  then
                this%density = 19816.0D0
            else if (t <= 479.0D0)  then
                this%density = 17770.0D0
            else if (t <= 592.0D0)  then
                this%density = 17140.0D0
            else if (t <= 724.0D0)  then
                this%density = 15920.0D0
            else if (t <= 749.0D0)  then
                this%density = 16010.0D0
            else
                this%density = 16480.0D0
            end if
        
        case(3)                                                                 ! Pu-Zr
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get density, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            density_Zr = 6570.0D0
            density_Pu = 19750.0D0
            
            key_Zr   = [ 293.0D0,  400.0D0,  500.0D0,  600.0D0,  700.0D0,  800.0D0,  900.0D0, 1000.0D0, 1100.0D0, 1135.0D0, 1135.0D0, 1200.0D0, 1400.0D0, 1600.0D0, 1800.0D0]
            value_Zr = [6570.0D0, 6558.0D0, 6546.0D0, 6532.0D0, 6518.0D0, 6503.0D0, 6488.0D0, 6471.0D0, 6456.0D0, 6450.0D0, 6476.0D0, 6465.0D0, 6429.0D0, 6391.0D0, 6351.0D0]
            
            key_Pu   = [  293.0D0,   350.0D0,   397.6D0,   397.6D0,   400.0D0,   450.0D0,   487.9D0,   487.9D0,   500.0D0,   550.0D0,   593.1D0,   593.1D0,   600.0D0,   700.0D0,   736.0D0,   736.0D0,   755.7D0,   755.7D0,   800.0D0,   900.0D0,   913.0D0,  1800.0D0]
            value_Pu = [19750.0D0, 19586.0D0, 19448.0D0, 17724.0D0, 17720.0D0, 17625.0D0, 17553.0D0, 17096.0D0, 17075.0D0, 16991.0D0, 16918.0D0, 15843.0D0, 15846.0D0, 15884.0D0, 15898.0D0, 15921.0D0, 15935.0D0, 16444.0D0, 16369.0D0, 16201.0D0, 16179.0D0, 16179.0D0]
            
            call linear_interpolation (key_Zr, value_Zr, t, density_Zr)
            call linear_interpolation (key_Pu, value_Pu, t, density_Pu)
            
            this%density = density_Zr*this%x_Zr + density_Pu*this%x_Pu
            
        end select
        
        ! set value
        density = this%density
    
    end function Get_density_by_temperature

    !$
    !===============================================================================================
    ! functions for thermal properties of fuel regions thermal conductivity in w/m-K, t in K
    !===============================================================================================
    function Get_conductivty_by_temperature(this, t_in)  result(conductivity)
        
        class(FuelProperty_metallic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: conductivity
        
        real(KREAL)  :: t
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        
        real(KREAL)  :: wdeltap                                                 ! wdeltap, the weight fraction of delta-pu phase
        real(KREAL)  :: vdeltap                                                 ! valphap, the volume fraction of delta-pu phase
        real(KREAL)  :: valphaz                                                 ! vdeltaz, the volume fraction of alpha-zr phase
        real(KREAL)  :: a,b,c,d
        real(KREAL)  :: density_Zr(2,8), density_Pu(2,11)
        real(KREAL)  :: rhozr,rhopu,rhodeltap,rhoalphaz
        real(KREAL)  :: cd,cdalphazr,cddeltapu
        
        integer  :: i, j
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! U
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 22.0D0 + 0.023D0*(t-273.15D0)
            
        case(2)                                                                 ! Pu
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            ! step-wise approximation
            if (t <= 395.0D0)  then
                key   = [273.15D0, 395.0D0]
                value = [5.2D0, 6.6D0]
                call linear_interpolation (key, value, t, this%conductivity)
            else if (t <= 479.0D0)  then
                key   = [395.0D0, 479.0D0]
                value = [7.87D0, 8.67D0]
                call linear_interpolation (key, value, t, this%conductivity)
            else if (t <= 592.0D0)  then
                key   = [479.0D0, 592.0D0]
                value = [8.97D0, 10.5D0]
                call linear_interpolation (key, value, t, this%conductivity)
            else if (t <= 724.0D0)  then
                key   = [592.0D0, 724.0D0]
                value = [10.97D0, 12.1D0]
                call linear_interpolation (key, value, t, this%conductivity)
            else if (t <= 749.0D0)  then
                this%conductivity = 7.72D0
            else
                this%conductivity = 12.35D0
            end if
        
        case(3)                                                                 ! Pu-Zr
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get conductivity, exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            density_Zr(1,1) = 500.0D0;   density_Zr(2,1) = 6.546D0
            density_Zr(1,2) = 600.0D0;   density_Zr(2,2) = 6.532D0
            density_Zr(1,3) = 700.0D0;   density_Zr(2,3) = 6.518D0
            density_Zr(1,4) = 800.0D0;   density_Zr(2,4) = 6.503D0
            density_Zr(1,5) = 900.0D0;   density_Zr(2,5) = 6.488D0
            density_Zr(1,6) = 1000.0D0;  density_Zr(2,6) = 6.471D0
            density_Zr(1,7) = 1100.0D0;  density_Zr(2,7) = 6.456D0
            density_Zr(1,8) = 1135.0D0;  density_Zr(2,8) = 6.450D0
            
            density_Pu(1,1) = 500.0D0;  density_Pu(2,1) = 17.075D0
            density_Pu(1,2) = 550.0D0;  density_Pu(2,2) = 16.991D0
            density_Pu(1,3) = 593.1D0;  density_Pu(2,3) = 15.843D0              ! extend to 540 K;
            density_Pu(1,4) = 600.0D0;  density_Pu(2,4) = 15.846D0
            density_Pu(1,5) = 700.0D0;  density_Pu(2,5) = 15.884D0
            density_Pu(1,6) = 736.0D0;  density_Pu(2,6) = 15.898D0
            density_Pu(1,7) = 736.0D0;  density_Pu(2,7) = 15.921D0
            density_Pu(1,8) = 755.7D0;  density_Pu(2,8) = 15.935D0
            density_Pu(1,9) = 755.7D0;  density_Pu(2,9) = 16.444D0
            density_Pu(1,10) = 800.0D0; density_Pu(2,10) = 16.369D0
            density_Pu(1,11) = 900.0D0; density_Pu(2,11) = 16.201D0
            
            associate (wz => this%weight_Zr, wp => this%weight_Pu)
            temperature: if (540.0D0<=t .AND. t<=870.0D0)  then
                
                ! delta-Pu phase
                if (0.0D0<=wz .AND. wz<=0.45D0)  then
                    a = (1.0D0-SQRT(1.0D0-wp))*3.225D0 + SQRT(1.0D0-wp)*(29.469D0-118.811D0*wp+88.893D0*wp**2)
                    b = (1.0D0-SQRT(1.0D0-wp))*0.0296D0 + SQRT(1.0D0-wp)*(0.0117D0-0.00716D0*wp)
                    c = SQRT(1.0D0-wp)*1.922D-5
                    this%conductivity = a+b*t+c*t**2
                
                ! alpha-Zr phase
                else if (0.8D0<=wz .AND. wz<=1.0D0)  then  
                    a = (1.0D0-SQRT(1.0D0-wz))*8.853D0 + SQRT(1.0D0-wz)*(30.57D0-82.301D0*wz+71.456D0*wz**2)
                    b = (1.0D0-SQRT(1.0D0-wz))*7.082D-3 + SQRT(1.0D0-wz)*(0.01895D0-0.02453D0*wz)
                    c = (1.0D0-SQRT(1.0D0-wz))*2.533D-6 + SQRT(1.0D0-wz)*8.111D-6
                    d = (1.0D0-SQRT(1.0D0-wz))*2.992D3
                    this%conductivity = a+b*t+c*t**2+d/t
                    
                ! delta-pu + alpha-zr phase
                else if (0.45D0<wz .AND. wz<0.8D0)  then                
                    ! bruggeman model -- aaa fuels handbook from ANL
                    ! delta-pu phase contains 45wt% zr and 55wt% pu, alpha-zr phase contains 80wt% zr and 20wt% pu
                    ! according to lever rule
                    wdeltap = (wz-0.8D0)/(0.45D0-0.8D0)
                    rhozr = 6.488D0
                    rhopu = 16.201D0
                    
                    call linear_interpolation (density_Zr(1,:), density_Zr(2,:), t, rhozr)
                    call linear_interpolation (density_Pu(1,:), density_Pu(2,:), t, rhopu)
            
                    rhodeltap = 1.0D0 / (0.45D0/rhozr + 0.55D0/rhopu)
                    rhoalphaz = 1.0D0 / (0.8D0/rhozr + 0.2D0/rhopu)
                    vdeltap = rhoalphaz*wdeltap / (rhodeltap - (rhodeltap-rhoalphaz)*wdeltap)
                    valphaz = 1.0D0 - vdeltap
            
                    ! delta-pu phase, 45wt% zr and 55wt% pu
                    a = (1.0D0-SQRT(1.0D0-0.55D0))*3.225D0 + SQRT(1.0D0-0.55D0)*(29.469D0-118.811D0*0.55D0+88.893D0*0.55D0**2)
                    b = (1.0D0-SQRT(1.0D0-0.55D0))*0.0296D0 + SQRT(1.0D0-0.55D0)*(0.0117D0-0.00716D0*0.55D0)
                    c = SQRT(1.0D0-0.55D0)*1.922D-5
                    cddeltapu = a+b*t+c*t**2
            
                    ! alpha-zr phase, 80wt% zr and 20wt% pu
                    a = (1.0D0-SQRT(1.0D0-0.8D0))*8.853D0 + SQRT(1.0D0-0.8D0)*(30.57D0-82.301D0*0.8D0+71.456D0*0.8D0**2)
                    b = (1.0D0-SQRT(1.0D0-0.8D0))*7.082D-3 + SQRT(1.0D0-0.8D0)*(0.01895D0-0.02453D0*0.8D0)
                    c = (1.0D0-SQRT(1.0D0-0.8D0))*2.533D-6 + SQRT(1.0D0-0.8D0)*8.111D-6
                    d = (1.0D0-SQRT(1.0D0-0.8D0))*2.992D3
                    cdalphazr = a+b*t+c*t**2+d/t
            
                    a = (3.0D0*vdeltap-1.0D0)*cddeltapu+(3.0D0*valphaz-1.0D0)*cdalphazr
            
                    this%conductivity = (a+SQRT(a**2+8.0D0*cdalphazr*cddeltapu))/4.0D0
                end if
                
            else temperature
                a = (1.0D0-SQRT(1.0D0-wz))*8.853D0 + SQRT(1.0D0-wz)*(wz*(-98.806D0+147.895D0*wz-26.883D0*wz**2)+(1.0D0-wz)*9.507D0)
                b = (1.0D0-SQRT(1.0D0-wz))*7.082D-3 + SQRT(1.0D0-wz)*(wz*(0.0512D0-0.0601D0*wz)+(1.0D0-wz)*0.0184D0)
                c = (1.0D0-SQRT(1.0D0-wz))*2.533D-6 + SQRT(1.0D0-wz)*(wz*8.699D-6)
                d = (1.0D0-SQRT(1.0D0-wz))*2.992D3
                
                this%conductivity =  a+b*t+c*t**2+d/t
                
            end if temperature
            end associate
            
        end select
        
        ! set value
        conductivity = this%conductivity
     
    end function Get_conductivty_by_temperature
    
    !$
    !===============================================================================================
    ! heat capacity of in J/m^3-K, t in K
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(FuelProperty_metallic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: capacity
    
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        real(KREAL)  :: rho                                                     ! density of fuel
        real(KREAL)  :: degree
        real(KREAL)  :: t
        real(KREAL)  :: capacity_Zr, capacity_Pu

        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case(1)                                                                 ! U
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (293.0D0<=t .AND. t<=942.0D0)  then
                this%capacity = 104.82D0 + 5.3686D-3*t + 10.1823D-5*t**2
            else if (t<=1049.0D0)  then
                this%capacity = 176.4D0
            else
                this%capacity = 156.8D0
            end if
            
        case(2)                                                                 ! Pu
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            ! step-wise approximation
            if (t <= 395.0D0)  then
                key   = [273.15D0, 395.0D0]
                value = [32.0D0, 34.3D0]
                call linear_interpolation (key, value, t, this%capacity)
            else if (t <= 479.0D0)  then
                key   = [395.0D0, 479.0D0]
                value = [34.3D0, 36.0D0]
                call linear_interpolation (key, value, t, this%capacity)
            else if (t <= 592.0D0)  then
                key   = [479.0D0, 592.0D0]
                value = [34.8D0, 39.8D0]
                call linear_interpolation (key, value, t, this%capacity)
            else if (t <= 724.0D0)  then
                this%capacity = 37.7D0
            else if (t <= 749.0D0)  then
                this%capacity = 37.4D0
            else
                this%capacity = 35.0D0
            end if
            
            this%capacity = this%capacity * 1000.0D0 / this%mol_mass
            
        case(3)                                                                 ! Pu-Zr
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'fuel metallic, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            capacity_Zr = 0.0D0
            if (298.15D0<=t .AND. t<=1135.0D0)  then
                capacity_Zr = 22.839D0+9.091D-3*t-2.132D4*t**(-2)
            else 
                capacity_Zr = 12.885D0+9.976D-3*t+5.518D6*t**(-2)
            end if
            capacity_Zr = capacity_Zr * 1000.0D0 / 91.22D0
           
            capacity_Pu = 0.0D0
            if (298.15D0<=t .AND. t<=397.6D0)  then
                capacity_Pu = 18.126D0 + 4.482D-2*t
            else if (t <= 487.9D0)  then
                capacity_Pu = 27.416D0 + 1.306D-2*t
            else if (t <= 593.1D0)  then
                capacity_Pu = 22.023D0 + 2.296D-2*t
            else if (t <= 736.0D0)  then
                capacity_Pu = 28.478D0 + 1.081D-2*t
            else if (t <= 755.7D0)  then
                capacity_Pu = 35.560D0
            else if (t <= 913.0D0)  then
                capacity_Pu = 33.720D0
            else
                capacity_Pu = 33.720D0
            end if
            capacity_Pu = capacity_Pu * 1000.0D0 / 244.0D0
            
            this%capacity = capacity_Zr*this%x_Zr + capacity_Pu*this%x_Pu
            
        end select
        
        capacity = this%capacity
        
    end function Get_capacity_by_temperature
    
    !$
    !===============================================================================================
    ! thermal expansion of axial in 1/K, t in K
    !===============================================================================================
    function Get_expansion_by_temperature (this, t_in)  result(expansion)
        
        class(FuelProperty_metallic), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: expansion
        
        ! ----------------------------------------------------------------------
        select case(this%fuel_type)
        case (1)                                                                ! U
            this%expansion = 1.0D-2

        case (2)                                                                ! Pu
            this%expansion = 1.0D-2
        
        case (3)                                                                ! Pu-Zr
            this%expansion = 1.0D-2
        
        end select
        
        expansion = this%expansion
    
    end function Get_expansion_by_temperature
    
    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine linear_interpolation (x, y, x_point, y_point)
        
        real(KREAL), intent(in)  :: x(:)
        real(KREAL), intent(in)  :: y(:)
        real(KREAL), intent(in)  :: x_point
        real(KREAL), intent(out) :: y_point
        
        integer  :: i
        
        do i = LBOUND(x, dim=1), UBOUND(x, dim=1)
            if (x_point <= x(i) + EPS_HIGH)  then
                y_point = y(i)
                exit
                
            else if (x_point < x(i) - EPS_HIGH)  then
                y_point = y(i-1) + (x_point-x(i-1)) * ((y(i)-y(i-1))/(x(i)-x(i-1)))
                exit 
            end if
        end do
    
    end subroutine linear_interpolation
    
end module fuel_metallic_header
