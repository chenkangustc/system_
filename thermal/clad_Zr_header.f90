!$
!===================================================================================================
!
!    this module is for clad property class Zr;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          CladProperty_Zr
!
!===================================================================================================
module clad_Zr_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : CladProperty
    
    implicit none 
    private
    public  :: CladProperty_Zr
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! --------------------------------------------------------------------------
    ! type for calding property
    type, extends(CladProperty)  :: CladProperty_Zr
        private
        real(KREAL)  :: t_melting                                               ! melting temperature (K)
        real(KREAL)  :: mol_mass                                                ! mol mass (g/mol)
    contains
        procedure, public  :: set => Set_CladProperty
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivty_by_temperature
        procedure, public  :: get_expansion => Get_expansion_by_temperature
    end type CladProperty_Zr
    
    ! private the real function name   
    private :: Set_CladProperty
    private :: Get_density_by_temperature, Get_capacity_by_temperature, Get_conductivty_by_temperature, Get_expansion_by_temperature
    
contains   
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CladProperty (this, type)
    
        class(CladProperty_Zr), intent(in out)  :: this
        integer, intent(in)  :: type
        
        this%clad_type = type
        
        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA metal-Zr
            this%t_melting = 2128.0D0
            this%mol_mass  = 91.22D0
            
        case(2)                                                                 ! IAEA- Zr + 1% Nb (E-110)
            this%t_melting = 2110.0D0
            this%mol_mass  = 91.24D0
        
        case(3)                                                                 ! IAEA- Zr + 2.5% Nb (E-125)
            this%t_melting = 2100.0D0
            this%mol_mass  = 91.26D0
        
        case(4)                                                                 ! PARCS Zr-alloy
            this%t_melting = 2118.15D0
            this%mol_mass  = 91.38D0
            
        case(5)                                                                 ! MATPRO Zr-alloy, same as 4
            this%t_melting = 2118.15D0
            this%mol_mass  = 91.38D0
            
        case(6)                                                                 ! NEA Zr-alloy, same as 4
            this%t_melting = 2118.15D0
            this%mol_mass  = 91.38D0
            
        case default
            call a_error%set (INFO_LIST_INPUT, 'clad Zr type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_CladProperty
  
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_density_by_temperature (this, t_in)  result(density)
    
        class(CladProperty_Zr), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA metal-Zr
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6550.0D0 - 0.1685D0*t
            
        case(2)                                                                 ! IAEA- Zr + 1% Nb (E-110)
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6636.0D0 - 0.286D0*t
        
        case(3)                                                                 ! IAEA- Zr + 2.5% Nb (E-125)
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6657.0D0 - 0.2861D0*t
        
        case(4)                                                                 ! PARCS Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6600.0D0
            
        case(5)                                                                 ! MATPRO Zr-alloy, same as 4
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6600.0D0
            
        case(6)                                                                 ! NEA Zr-alloy, same as 4
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 6600.0D0
            
        end select
        
        ! set value
        density = this%density
    
    end function Get_density_by_temperature
    
    !$
    !===============================================================================================
    ! volumetric heat capacity in J/m^3-K, t in K
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(CladProperty_Zr), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: capacity
        
        real(KREAL)  :: rho
        real(KREAL)  :: t
        
        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA metal-Zr
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (298.0D0<=t .AND. t<=1100.0D0)  then
                this%capacity = 238.596D0+0.181D0*t-96.1D-6*t**2+36.2D-9*t**3
            else
                this%capacity = 276.462D0+0.0141D0*t-3.08D-6*t**2+10.7D-9*t**3
            end if
            
        case(2)                                                                 ! IAEA- Zr + 1% Nb (E-110)
            if (t_in > 2000.0D0)  then
                t  = 2000.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<=1100.0D0)  then
                this%capacity = 238.0D0 + 0.159D0*t
            else
                this%capacity = 281.0D0 + 0.0663D0*t
            end if
        
        case(3)                                                                 ! IAEA- Zr + 2.5% Nb (E-125)
            if (t_in > 1600.0D0)  then
                t  = 1600.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<=1100.0D0)  then
                this%capacity = 221.0D0 + 0.172D0*t - 5.87D-5*t**2
            else
                this%capacity = 380.0D0
            end if
        
        case(4)                                                                 ! PARCS Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            rho = this%get_density(t)
            this%capacity = 252.54D0*rho+t*(0.11474D0*rho+t*(0.0D0*rho+t*0.0D0*rho))
            this%capacity = this%capacity / rho
        
        case(5)                                                                 ! MATPRO Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%conductivity = COBRA_clad_conductivity (t, real(1.0, KREAL))
            this%conductivity = this%conductivity * 6.23D3
            
        case(6)                                                                 ! NEA Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%conductivity = COBRA_clad_conductivity (t, real(-1.0D0, KREAL))
            this%conductivity = this%conductivity * 6.23D3
        
        end select
        
        ! set value
        capacity = this%capacity
        
    end function Get_capacity_by_temperature

    !$
    !===============================================================================================
    ! thermal conductivity  in w/m-K, t in K
    !===============================================================================================
    function Get_conductivty_by_temperature(this, t_in)  result(conductivity)

        class(CladProperty_Zr), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: conductivity
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA metal-Zr
            if (t_in > 2000.0D0)  then
                t  = 2000.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 8.8527D0+7.0820D-3*t+2.5329D-6*t**2+2.9918D3*t**(-1)
            
        case(2)                                                                 ! IAEA- Zr + 1% Nb (E-110)
            if (t_in > 1600.0D0)  then
                t  = 1600.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<=1100.0D0)  then
                this%conductivity = 23.5D0 - 0.0192D0*t + 1.68D-5*t**2
            else
                this%conductivity = 1.5D0 + 0.02D0*t
            end if
        
        case(3)                                                                 ! IAEA- Zr + 2.5% Nb (E-125)
            if (t_in > 1100.0D0)  then
                t  = 1100.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 14.0D0 + 0.0115D0*t
        
        case(4)                                                                 ! PARCS Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 7.51D0+t*(2.09D-2+t*(-1.45D-5+t*7.67D-9))
        
        case(5)                                                                 ! MATPRO Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%capacity = COBRA_clad_capacity (t, real(1.0, KREAL))
            this%capacity = this%capacity * 4186.3D0
            
        case(6)                                                                 ! NEA Zr-alloy
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad Zr, get conductivity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            t = t - CKELVIN
            t = 32.0D0 + (9.0D0/5.0D0)*t
            
            this%capacity = COBRA_clad_capacity (t, real(-1.0, KREAL))
            this%capacity = this%capacity * 4186.3D0
        
        end select
        
        ! set value
        conductivity = this%conductivity
        
    end function Get_conductivty_by_temperature
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_expansion_by_temperature (this, t_in)  result(expansion)
        
        class(CladProperty_Zr), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: expansion
    
        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case (1)                                                                ! IAEA metal-Zr
            this%expansion = 1.0D-2
            
        case (2)                                                                ! IAEA- Zr + 1% Nb (E-110)
            this%expansion = 1.0D-2
            
        case (3)                                                                ! IAEA- Zr + 2.5% Nb (E-125)
            this%expansion = 1.0D-2
            
        case (4)                                                                ! PARCS Zr-alloy
            this%expansion = 1.0D-2
            
        case (5)                                                                ! MATPRO Zr-alloy
            this%expansion = 1.0D-2
            
        case (6)                                                                ! NEA Zr-alloy
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
    !     flag = flag to use matpro correlation, if negative >  0.0
    !          = flag to use nea    correlation, if negative <  0.0
    !     output:
    !     conductivity=supplied clad thermal conductivity (btu/ft/s/f)
    ! ----------------------------------------------------------------------------------------------
    !     1 btu/(ft.s.f) = 6.23E3 W/(m.K)
    !     1 f = (9/5)*c + 32 = (9/5)*(K-273.15) + 32
    !===============================================================================================
    function COBRA_clad_conductivity (t_in, flag)  result(conductivity)
        
        real(KREAL), intent(in)  :: t_in
        real(KREAL), intent(in)  :: flag
        real(KREAL)  :: conductivity
        
        real(KDOUBLE)  :: coef(5), tt, w
        real(KREAL)    :: c(8)
        
        real(KREAL), parameter  :: tklo = 298.0D0
        real(KREAL), parameter  :: tkup = 3100.0D0
        
        coef = [7.51D0, 0.0209D0, 1.45D-5, 7.67D-9, 1.6056D-4]
        c    = [32.0D0, 1.8D0, 273.15D0, 7.51D0, 2.09D-2, 1.45D-5, 7.67D-9, 6226.48D0]
    
        ! ----------------------------------------------------------------------
        if (flag > 0.0D0)  then                                                 ! MATPRO
            tt = (t_in - c(1))/c(2) + c(3)
            if (tt > tkup)  then
                tt = tkup
            end if
            if (tt < tklo)  then
                tt = tklo
            end if
            
            w = ((coef(4)*tt-coef(3))*tt+coef(2))*tt + coef(1)
            
            conductivity = coef(5) * w
            
        else                                                                    ! NEA
            tt = (t_in - c(1))/c(2) + c(3)
            if (tt < tklo)  then
                tt = tklo
            end if
            if (tt > tkup)  then
                tt = tkup
            end if
            
            w = c(4) + ((c(7)*tt-c(6))*tt+c(5))*tt
            
            conductivity = w / c(8)
        end if
    
    end function COBRA_clad_conductivity
    
    !$
    !===============================================================================================
    !     input:
    !     t = temperature (f)
    !     flag = flag to use matpro correlation, if negative > 0.0
    !          = flag to use nea    correlation, if negative < 0.0
    !     output:
    !     capacity=supplied clad specific heat (btu/lb/f)
    ! ----------------------------------------------------------------------------------------------
    !     1 btu/(lb.f) = 4186.3 J/(Kg.K)
    !     1 f = (9/5)*c + 32 = (9/5)*(K-273.15) + 32
    !===============================================================================================
    function COBRA_clad_capacity (t_in, flag)  result(capacity)
        
        real(KREAL), intent(in)  :: t_in
        real(KREAL), intent(in)  :: flag
        real(KREAL)  :: capacity
    
        real(KREAL)  :: tt,w
        real(KREAL)  :: td(13)
        real(KREAL)  :: cd(13)
        real(KREAL)  :: c(6)

        real(KREAL), parameter  :: tklo = 293.0D0
        real(KREAL), parameter  :: tkup = 3100.0D0
        
        td = [300.0D0, 400.0D0, 640.0D0, 1090.0D0, 1093.0D0, 1113.0D0, 1133.0D0, 1153.0D0, 1173.0D0, 1193.0D0, 1213.0D0, 1233.0D0, 1248.0D0]
        cd = [0.0671D0, 0.07212D0, 0.07904D0, 0.08955D0, 0.11988D0, 0.14089D0, 0.14686D0, 0.1717D0, 0.1949D0, 0.18388D0, 0.1478D0, 0.1120D0, 0.0850D0]
        c  = [32.0D0, 1.8D0, 273.15D0, 252.54D0, 0.11474D0, 4184.0D0]
        
        ! ----------------------------------------------------------------------
        if (flag > 0.0D0)  then                                                 ! MATPRO
            tt = (t_in - c(1))/c(2) + c(3)
            
            if (tt <= td(1))  then
                capacity = cd(1)
            else if (tt >= td(13))  then
                capacity = cd(13)
            else
                call linear_interpolation (td, cd, tt, w)
                capacity = w
            end if
            
        else                                                                    ! NEA
            tt = (t_in - c(1))/c(2) + c(3)
            if (tt < tklo)  then
                tt = tklo
            end if
            if (tt > tkup)  then
                tt = tkup
            end if
            
            w = c(4) + c(5)*tt
            capacity = w / c(6)
        end if
    
    end function COBRA_clad_capacity
    
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
    
end module clad_Zr_header
