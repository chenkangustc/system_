!$
!===================================================================================================
!
!    this module is for gap property class;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          GapProperty_gas
!
!===================================================================================================
module gap_gas_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector, WarningCollector
    use abstract_property_header,   only : GapProperty
    
    implicit none
    private
    public  :: GapProperty_gas
    
    type(ErrorCollector)    :: a_error                                          ! to print error information
    type(WarningCollector)  :: a_warning                                        ! to print warning information when excess define field
    
    ! --------------------------------------------------------------------------
    ! type for calding property
    type, extends(GapProperty)  :: GapProperty_gas
        private
        real(KREAL)      :: mol_mass                                            ! mol mass (g/mol)
        real(KREAL)      :: conductivity                                        ! thermal conductivity (W/m-K)
                         
        real(KREAL)      :: x_He                                                ! atomic fraction of He
        real(KREAL)      :: x_Kr                                                ! atomic fraction of Kr
        real(KREAL)      :: x_Xe                                                ! atomic fraction of Xe
    contains
        procedure, public  :: set => Set_GapProperty
        procedure, public  :: get_transfer => Get_transfer_by_temperature
    end type GapProperty_gas
    
    ! private the real function name   
    private :: Set_GapProperty, Get_transfer_by_temperature
    
contains 
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_GapProperty (this, type, x_Xe, x_Kr)
    
        class(GapProperty_gas), intent(in out)  :: this
        integer, intent(in)  :: type
        real(KREAL), intent(in), optional  :: x_Xe
        real(KREAL), intent(in), optional  :: x_Kr
        
        this%gas_type = type
        
        ! ----------------------------------------------------------------------
        select case(this%gas_type)
        case(1)                                                                 ! RELAP5, He+Kr+Xe
            this%x_Xe = 0.0D0                                                   ! 0.7594, for PWR
            this%x_Kr = 0.0D0                                                   ! 0.1340, for PWR
            if (PRESENT(x_Xe))  then
                this%x_Xe = x_Xe
            end if
            if (PRESENT(x_Kr))  then
                this%x_Kr = x_Kr
            end if
            
            this%x_He = 1.0D0 - this%x_Xe - this%x_Kr
            this%mol_mass = this%x_He*4.003D0 + this%x_Kr*83.8D0 + this%x_Xe*131.3D0
             
        case(2)                                                                 ! IAEA-Na
            this%mol_mass = 23.00D0

        case(3)                                                                 ! IAEA-LBE
            this%mol_mass = 208.0D0
            
        case(4)                                                                 ! beam interruption benchmark, He
            this%mol_mass = 4.003D0

        case(5)                                                                 ! PWR MOX benchmark
            this%mol_mass = 16.0D0
           
        case(6)                                                                 ! BN-600
            this%mol_mass = 16.0D0
           
        case default
            call a_error%set (INFO_LIST_INPUT, 'gas type is not pre-defined')  
            call a_error%print (FILES%MAIN)
        end select
    
    end subroutine Set_GapProperty
    
    !$
    !===============================================================================================
    ! heat transfer coefficient [(W/K-m^2)]
    !===============================================================================================
    function Get_transfer_by_temperature (this, t_in, pellet, gap, is_inner)  result(h_transfer)
    
        class(GapProperty_gas), intent(in out)  :: this
        real(KREAL), intent(in)      :: t_in                                    ! average temperature
        real(KREAL), intent(in)      :: pellet                                  ! fuel pellet radius
        real(KREAL), intent(in)      :: gap                                     ! gap thickness
        logical, intent(in)          :: is_inner                                ! inner face no outer face of gap
        real(KREAL)  :: h_transfer
        
        real(KREAL)  :: t
        real(KREAL)  :: multi
        
        real(KREAL)  :: y(3), x(3)
        real(KREAL)  :: k(3), m(3)
        real(KREAL)  :: phi(3,3), psi(3,3)
        integer  :: i, j
    
        ! ----------------------------------------------------------------------
        select case(this%gas_type)
        case(1)                                                                 ! RELAP5, He+Kr+Xe
            if (is_inner )  then
                multi = (pellet+0.5D0*gap) / (gap*pellet)
            else 
                multi = (pellet+0.5D0*gap) / (gap*(pellet+gap))
            end if
            
            t = t_in
            
            x(1) = this%x_He;  m(1) = 4.003D0;  k(1) = 2.639D-3*t**0.7085D0
            x(2) = this%x_Kr;  m(2) = 83.3D0 ;  k(2) = 8.247D-5*t**0.8363D0
            x(3) = this%x_Xe;  m(3) = 131.3D0;  k(3) = 4.351D-5*t**0.8616D0
            
            phi = 0.0D0
            do i = 1, 3
                do j = 1, 3
                    phi(i, j) = (1+(k(i)/k(j))**0.50D0*(m(i)/m(j))**0.25D0)**2 / (2**1.5D0*(1+m(i)/m(j))**0.5D0)
                end do
            end do
            
            psi = 0.0D0
            do i = 1, 3
                do j = 1, 3
                    psi(i, j) = phi(i, j) * (1+2.41D0*(m(i)-m(j))*(m(i)-0.142D0*m(j))/(m(i)+m(j))**2)
                end do
            end do
            
            this%conductivity = 0.0D0
            do i = 1, 3
                y(i) = 0.0D0
                do j = 1, 3
                    if (j == i)  then
                        y(i) = y(i) + (1-1)*psi(i,j)*x(j)
                    else 
                        y(i) = y(i) + (1-0)*psi(i,j)*x(j)
                    end if
                end do
                
                this%conductivity = this%conductivity + k(i)*x(i)/(x(i)+y(i))
            end do
            
            this%h_transfer = this%conductivity * multi
            
        case(2)                                                                 ! IAEA-Na
            if (is_inner )  then
                multi = (pellet+0.5D0*gap) / (gap*pellet)
            else 
                multi = (pellet+0.5D0*gap) / (gap*(pellet+gap))
            end if
            
            t = t_in
        
            this%conductivity = 110.0D0-6.48D-2*t+1.16D-5*t**2
            this%h_transfer = this%conductivity * multi
            
        case(3)                                                                 ! IAEA-LBE
            if (is_inner )  then
                multi = (pellet+0.5D0*gap) / (gap*pellet)
            else 
                multi = (pellet+0.5D0*gap) / (gap*(pellet+gap))
            end if
            
            t = t_in
        
            this%conductivity = 3.284D0+1.617D-2*t-2.305D-6*t**2
            this%h_transfer = this%conductivity * multi
            
        case(4)                                                                 ! beam interruption benchmark, He
            t = t_in
        
            this%h_transfer = 3.623D-3 * t**0.66D0 / gap
            
        case(5)                                                                 ! PWR MOX benchmark
            t = t_in
            this%h_transfer = 1.0D4
            
        case(6)                                                                 ! BN-600
            t = t_in
            this%h_transfer = 1.0D4
            
        end select
        
        ! set value
        h_transfer = this%h_transfer
    
    end function Get_transfer_by_temperature
    
end module gap_gas_header
