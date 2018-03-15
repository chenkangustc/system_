!$
!===================================================================================================
!
!    this module is for water property class (water properties at 15.5 MPa, from NIST REFPROP);
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          WaterProperty_REFPROP
!
!===================================================================================================
module coolant_water_refprop_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use abstract_property_header,   only : CoolantProperty
    
    implicit none 
    private
    public  :: WaterProperty_REFPROP
    
    ! --------------------------------------------------------------------------
    ! for state parameter, the name is the same as NIST REFPROP code
    type, private  :: stat_parameter_tp
        double precision  :: wmm                                                ! molecular weight [g/mol]
        double precision  :: t                                                  ! temperature [K]
        double precision  :: p                                                  ! pressure [kPa]
        double precision  :: rho                                                ! molar density [mol/L]
        double precision  :: eta                                                ! viscosity (uPa.s)
        double precision  :: tcx                                                ! thermal conductivity (W/m.K)
        double precision  :: e                                                  ! internal energy [J/mol]
        double precision  :: h                                                  ! enthalpy [J/mol]
        double precision  :: s                                                  ! entropy [J/mol-K]
        double precision  :: Cv                                                 ! isochoric heat capacity [J/mol-K]
        double precision  :: Cp                                                 ! isobaric heat capacity [J/mol-K]
        double precision  :: w                                                  ! speed of sound [m/s]
        double precision  :: hjt                                                ! isenthalpic Joule-Thompson coefficient [K/kPa]
        double precision  :: D                                                  ! overall (bulk) molar density [mol/L]
        double precision  :: Dl                                                 ! molar density [mol/L] of the liquid phase
        double precision  :: Dv                                                 ! molar density [mol/L] of the vapor phase
        integer           :: q                                                  ! vapor quality [basis specified by kq]
    end type stat_parameter_tp
    
    ! type for initial parameter
    type, private  :: initial_parameter_tp
        real(KREAL)                      :: pressure                            ! pressure (Mpa)

        integer                          :: nc_max
        integer                          :: nc
        integer                          :: ierr
        
        double precision, allocatable    :: x(:) 
        double precision, allocatable    :: x_liquid(:)
        double precision, allocatable    :: x_vapor(:)
                                         
        character(len=3)                 :: hrf
        character(len=255)               :: hfmix
        character(len=255)               :: herr
        character(len=255), allocatable  :: hfiles(:)
    end type initial_parameter_tp
    
    ! --------------------------------------------------------------------------
    ! type for water proterty
    type, extends(CoolantProperty)  :: WaterProperty_REFPROP
        private        
        type(stat_parameter_tp)          :: point
        type(initial_parameter_tp)       :: init
    contains
        procedure, public  :: set => Set_WaterProperty
        
        procedure, public  :: get_density => Get_density_by_temperature
        procedure, public  :: get_enthalpy => Get_enthalpy_by_temperature
        procedure, public  :: get_temperature => Get_temperature_by_enthalpy
        procedure, public  :: get_capacity => Get_capacity_by_temperature
        procedure, public  :: get_conductivity => Get_conductivity_by_temperature
        procedure, public  :: get_viscosity => Get_viscosity_by_temperature
        
        procedure, public  :: get_nusselt => Get_Nusselt_number
    end type WaterProperty_REFPROP
    
    ! private the real function name
    private :: Set_WaterProperty
    private :: Get_density_by_temperature, Get_enthalpy_by_temperature
    private :: Get_temperature_by_enthalpy, Get_capacity_by_temperature
    private :: Get_conductivity_by_temperature, Get_viscosity_by_temperature, Get_Nusselt_number
    
contains
    !$
    !===============================================================================================
    ! set information, blank temperary (option=1, 15.5Mpa-water; option=2, 7.0MPa-steam)
    !===============================================================================================
    subroutine Set_WaterProperty (this, type, option)
        
        class(WaterProperty_REFPROP), intent(in out)  :: this
        integer, intent(in)  :: type
        integer, intent(in)  :: option
        
        double precision  :: ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas
        integer  :: i_allocate
        
        ! prepare
        this%init%nc_max    = 20                                                ! max number of components
        this%init%nc        = 1                                                 ! real number of components
        
        allocate(this%init%x(this%init%nc_max), stat=i_allocate)
        allocate(this%init%x_liquid(this%init%nc_max), stat=i_allocate)
        allocate(this%init%x_vapor(this%init%nc_max), stat=i_allocate)
        allocate(this%init%hfiles(this%init%nc_max), stat=i_allocate)
        
        this%init%hfiles(1) = 'WATER.FLD'                                       ! Fluid name
        this%init%hfmix     = 'HMX.BNC'                                         ! Mixture file name
        this%init%hrf       = 'DEF'                                             ! Reference state (DEF means default)
        
        if (option == 1)  then
            this%init%pressure    = 15.50                                       ! pressure pre-define, [MPa]
            this%init%x(1)        = 1.0
            this%init%x_liquid(1) = 1.0
            this%init%x_vapor(1)  = 0.0
        else if (option == 2)  then
            this%init%pressure    = 7.00                                        ! pressure pre-define, [MPa]
            this%init%x(1)        = 1.0
            this%init%x_liquid(1) = 0.0
            this%init%x_vapor(1)  = 1.0
        end if 
        
        call SETUP (this%init%nc, this%init%hfiles, this%init%hfmix, this%init%hrf, this%init%ierr, this%init%herr)
        if (this%init%ierr /= 0)   then
            write (OUTPUT_UNIT, *) this%init%herr
        end if
          
        call INFO (1, this%point%wmm, ttrp, tnbpt, tc, pc, Dc, Zc, acf, dip, Rgas)
    
    end subroutine Set_WaterProperty
    
    !$
    !===============================================================================================
    !
    !===============================================================================================    
    function Get_density_by_temperature(this, t_in)  result(density)
      
        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)  :: density
            
        this%point%t = t_in
        this%point%p = this%init%pressure * 1000.0D0
        
        call TPFLSH (this%point%t, this%point%p, this%init%x, this%point%D, this%point%Dl, this%point%Dv,       &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%h,                &
            &   this%point%s, this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        this%density = this%point%D * this%point%wmm
        density = this%density
      
    end function Get_density_by_temperature
        
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_enthalpy_by_temperature(this, t_in)  result(enthalpy)

        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)    :: enthalpy

        this%point%t = t_in
        this%point%p = this%init%pressure * 1000.0D0
        
        call TPFLSH (this%point%t, this%point%p, this%init%x, this%point%D, this%point%Dl, this%point%Dv,       &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%h,                &
            &   this%point%s, this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        
        this%enthalpy = this%point%h*1000.0D0 / this%point%wmm
        enthalpy = this%enthalpy
      
    end function Get_enthalpy_by_temperature
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_temperature_by_enthalpy(this, h_in)  result(temperature)

        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: h_in
        real(KREAL)  :: temperature
        
        this%point%p = this%init%pressure * 1000.0D0
        this%point%h = h_in * this%point%wmm / 1000.0D0
        
        call PHFLSH (this%point%p, this%point%h, this%init%x, this%point%t, this%point%D, this%point%Dl, this%point%Dv,     &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%s,                            &
            &   this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        this%temperature = this%point%t
        temperature = this%temperature
      
    end function Get_temperature_by_enthalpy
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_capacity_by_temperature(this, t_in)  result(capacity)

        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)    :: capacity

        this%point%t = t_in
        this%point%p = this%init%pressure * 1000.0D0
        
        call TPFLSH (this%point%t, this%point%p, this%init%x, this%point%D, this%point%Dl, this%point%Dv,       &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%h,                &
            &   this%point%s, this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        ! change unit
        this%capacity = this%point%cp*1000.0D0 / this%point%wmm                                
        capacity = this%capacity
      
    end function Get_capacity_by_temperature
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_conductivity_by_temperature(this, t_in)  result(conductivity)
        
        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)    :: conductivity

        this%point%t = t_in
        this%point%p = this%init%pressure * 1000.0D0
        
        call TPFLSH (this%point%t, this%point%p, this%init%x, this%point%D, this%point%Dl, this%point%Dv,       &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%h,                &
            &   this%point%s, this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        this%point%t = t_in
        call TRNPRP (this%point%t, this%point%D, this%init%x, this%point%eta, this%point%tcx, this%init%ierr, this%init%herr)
            
        this%conductivity = this%point%tcx
        conductivity = this%conductivity
        
    end function Get_conductivity_by_temperature
    
    !$
    !===============================================================================================
    ! dynamics viscosity [Pa.s]
    !===============================================================================================
    function Get_viscosity_by_temperature(this, t_in)  result(viscosity)
    
        class(WaterProperty_REFPROP), intent(in out)  :: this
        real(KREAL), intent(in)  :: t_in
        real(KREAL)    :: viscosity
        
        this%point%t = t_in
        this%point%p = this%init%pressure * 1000.0D0
        
        call TPFLSH (this%point%t, this%point%p, this%init%x, this%point%D, this%point%Dl, this%point%Dv,       &
            &   this%init%x_liquid, this%init%x_vapor, this%point%q, this%point%e, this%point%h,                &
            &   this%point%s, this%point%cv, this%point%cp, this%point%w, this%init%ierr, this%init%herr)
        
        this%point%t = t_in
        call TRNPRP (this%point%t, this%point%D, this%init%x, this%point%eta, this%point%tcx, this%init%ierr, this%init%herr)
            
        this%viscosity = this%point%eta * 1.0D-6
        viscosity = this%viscosity
      
    end function Get_viscosity_by_temperature

    !$
    !===============================================================================================
    ! nusselt number
    !===============================================================================================
    function Get_Nusselt_number (this, P2D, velocity, Dh, t_in)  result(nusselt)

        class(WaterProperty_REFPROP), intent(in out)  :: this
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
    
end module coolant_water_refprop_header
