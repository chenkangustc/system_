!$
!===================================================================================================
!
!    this module is for clad property class steels;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          CladProperty_steels
!
!===================================================================================================
module clad_steels_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : CladProperty
    
    implicit none 
    private
    public  :: CladProperty_steels
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! --------------------------------------------------------------------------
    ! type for calding property
    type, extends(CladProperty)  :: CladProperty_steels
        private
        real(KREAL)  :: t_melting                                               ! melting temperature (K)
        real(KREAL)  :: mol_mass                                                ! mol mass (g/mol)
    contains
        procedure, public  :: set => Set_CladProperty
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivty_by_temperature
        procedure, public  :: get_expansion => Get_expansion_by_temperature
    end type CladProperty_steels
    
    ! private the real function name   
    private :: Set_CladProperty
    private :: Get_density_by_temperature, Get_capacity_by_temperature, Get_conductivty_by_temperature, Get_expansion_by_temperature
    
contains 
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CladProperty (this, type)
    
        class(CladProperty_steels), intent(in out)  :: this
        integer, intent(in)  :: type
        
        this%clad_type = type
        
        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA, 316
            this%t_melting = 1703.0D0
            this%mol_mass  = 55.9354D0
            
        case(2)                                                                 ! RELAP5, 304, as 1
            this%t_melting = 1703.0D0
            this%mol_mass  = 55.9354D0
            
        case(3)                                                                 ! HT9, as 1
            this%t_melting = 1703.0D0
            this%mol_mass  = 55.9354D0
            
        case(4)                                                                 ! OECD/NEA beam trip , as 1
            this%t_melting = 1703.0D0
            this%mol_mass  = 55.9354D0
            
        case default
            call a_error%set (INFO_LIST_INPUT, 'clad steels type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_CladProperty
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_density_by_temperature (this, t_in)  result(density)
    
        class(CladProperty_steels), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
        
        real(KREAL)  :: t
        real(KREAL)  :: expansion

        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA, 316
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%density = 8084.0D0 - 0.4209D0*t - 3.894D-5*t**2
            
        case(2)                                                                 ! RELAP5, 304
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<=1671.0D0)  then
                expansion = 1.57D-5*t + 1.69D-9*t**2
            else 
                expansion = -2.986634D-1 + 1.972573D-4*t
            end if
        
            this%density = 7800.0D0
            this%density = this%density * (1.0D0-3.0D0*expansion)
            
        case(3)                                                                 ! HT9, same as 430 series
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
        
            this%density = 7700.0D0
            
        case(4)                                                                 ! OECD/NEA beam trip, same as 430 series
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get density, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
        
            this%density = 7924.0D0
            
        end select
        
        ! set value
        density = this%density
    
    end function Get_density_by_temperature
    
    !$
    !===============================================================================================
    ! volumetric heat capacity in J/m^3-K, t in K
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(CladProperty_steels), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: capacity
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA, 316
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%capacity = 462.0D0 + 0.134D0*t
            
        case(2)                                                                 ! RELAP5, 304
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<= 1671.0D0)  then
                this%capacity = 326.0D0 - 0.242D0*t + 3.71D0*(t**0.719D0)
            else 
                this%capacity = 691.98D0
            end if
            
        case(3)                                                                 ! HT9, same as 1
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
        
            this%capacity = 460.0D0 + 0.134D0*t
            
        case(4)                                                                 ! OECD/NEA beam trip, same as 1
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get capacity, upper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
        
            this%capacity = 620.0
            
        end select
        
        ! set value
        capacity = this%capacity
        
    end function Get_capacity_by_temperature
    
    !$
    !===============================================================================================
    ! cladding conducivity [(J/Kg-K)]
    !===============================================================================================
    function Get_conductivty_by_temperature(this, t_in)  result(conductivity)

        class(CladProperty_steels), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in                                     
        real(KREAL) :: conductivity
        
        real(KREAL)  :: t

        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case(1)                                                                 ! IAEA, 316
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
        
            this%conductivity = 9.248D0 + 0.01571D0*t
            
        case(2)                                                                 ! RELAP5, 304
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (300.0D0<=t .AND. t<= 1671.0D0)  then
                this%conductivity = 7.58D0 + 0.0189D0*t
            else if (t <= 1727.0D0)  then
                this%conductivity = 610.9393D0 - 0.342176D0*t
            else 
                this%conductivity = 20.0D0
            end if
            
        case(3)                                                                 ! HT9
            if (t_in > 1200.0D0)  then
                t  = 1200.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else if (t_in < 500.0D0)  then
                t  = 500.0D0
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, lower exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            if (500.0D0<=t .AND. t<=1027.0D0) then
                this%conductivity = 17.622D0+2.428D-2*t-1.696D-5*t**2
            else
                this%conductivity = 12.027D0+1.218D-2*t
            end if
            
        case(4)                                                                 ! OECD/NEA beam trip 
            if (t_in > this%t_melting)  then
                t  = this%t_melting
                !  call a_warning%set (INFO_LIST_PROPERTY, 'clad steels, get conductivity, uppper exceed')
                !  call a_warning%print (FILES%TH_WARNING)
            else 
                t = t_in
            end if
            
            this%conductivity = 15.4767D0 + t*3.448D-3
            
        end select
        
        ! set value
        conductivity = this%conductivity
        
    end function Get_conductivty_by_temperature
    
    !$
    !===============================================================================================
    ! thermal expansion of radial in 1/K, t in K
    !===============================================================================================
    function Get_expansion_by_temperature (this, t_in)  result(expansion)
        
        class(CladProperty_steels), intent(in out)  :: this
        real(KREAL), intent(in) :: t_in                                     
        real(KREAL) :: expansion
        
        ! ----------------------------------------------------------------------
        select case(this%clad_type)
        case (1)                                                                ! IAEA, 316
            this%expansion = 1.0D-2
            
        case (2)                                                                ! RELAP5, 304
            this%expansion = 1.0D-2
        
        case (3)                                                                ! HT9
            this%expansion = 1.0D-2
        
        case (4)                                                                ! OECD/NEA beam trip 
            this%expansion = 1.0D-2
        
        end select
        
        expansion = this%expansion
        
    end function Get_expansion_by_temperature
    
end module clad_steels_header
