!$
!===================================================================================================
!
!   class for reativity format
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          PKReactivity                               
!
!===================================================================================================
module pkreactivity_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: PKReactivity
    
    ! --------------------------------------------------------------------------
    ! type for reactivity function parameter
    type  PKReactivity
        character(len=MAX_WORD_LEN)     :: type                                 ! reactivity type
        real(KREAL)                 :: first                                    ! first parameter for type
        real(KREAL)                 :: second                                   
        real(KREAL)                 :: third                                    
        real(KREAL)                 :: initial                                  ! initial reactivity before the transient start
        real(KREAL)                 :: rho_start                                ! start time of this reactivity
        real(KREAL)                 :: rho_end                                  ! end time of this reactivity
                                                                                
        logical                     :: is_feedback  = .FALSE.                   ! for feedback
        real(KREAL)                 :: T0_fuel      = 800.0                     
        real(KREAL)                 :: T_fuel       = 800.0                     
        real(KREAL)                 :: T0_coolant   = 600.0                     
        real(KREAL)                 :: T_coolant    = 600.0                     
        real(KREAL)                 :: Rho0_coolant   = 10.0                    
        real(KREAL)                 :: Rho_coolant    = 10.0                    
        real(KREAL)                 :: AD       = -4.8684E-4                    ! feedback cefficient
        real(KREAL)                 :: Bv       = -3.9083E-4
        real(KREAL)                 :: AXIAL    = -0.3826E-5
        real(KREAL)                 :: RADIAL   = -0.6537E-5
    contains
        procedure, public  :: get => Get_reactivity
        procedure, public  :: set_init_fdbk => Set_initial_feedback
        procedure, public  :: set_fdbk => Set_feedback_parameter
        procedure, public  :: set_coeff => Set_feedback_coefficient
        
        procedure, private  :: Reactivity_step
        procedure, private  :: Reactivity_ramp
        procedure, private  :: Reactivity_oscillatory
        procedure, private  :: Reactivity_special
        procedure, private  :: Reactivity_feedback
    end type PKReactivity
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Get_reactivity, Set_initial_feedback, Set_feedback_parameter, Set_feedback_coefficient
    
contains
    !$
    !===============================================================================================
    ! input a time point, return a reactivity rho according to the reactivity function type.
    !===============================================================================================
    subroutine Get_reactivity(this, tin, rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in)      :: tin
        real(KREAL), intent(in out)  :: rho
        
        real(KREAL)  :: rho_tmp
        
        select case(TRIM(this%type))
        case ('STEP')
            rho = this%Reactivity_step (tin)
            
        case ('RAMP')
            rho = this%Reactivity_ramp (tin)
            
        case ('OSCILLATORY')
            rho = this%Reactivity_oscillatory (tin)
            
        case ('SPECIAL')
            rho = this%Reactivity_special (tin)
            
        end select
        
        if (this%is_feedback)  then
            rho_tmp = this%Reactivity_feedback (tin)
            rho = rho + rho_tmp
        end if
    
    end subroutine Get_reactivity

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !===============================================================================================
    ! step reactivity function, use only one parameter
    ! rho = constant
    !===============================================================================================
    function Reactivity_step(this, tin) result(rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in)  :: tin
        real(KREAL)              :: rho
        
        if (tin < this%rho_start)  then
            rho = this%initial 
        else if (tin <= this%rho_end) then
            rho = this%initial + this%first
        else
            rho = this%initial
        end if
        
    end function Reactivity_step
    
    !$
    !===============================================================================================
    ! ramp reactivity function, use two parameter
    ! rho = a*t + b
    !===============================================================================================
    function Reactivity_ramp(this, tin) result(rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in) :: tin
        real(KREAL)             :: rho
        
        real(KREAL)  :: t
        
        if (tin < this%rho_start)  then
            rho = this%initial 
        else if (tin <= this%rho_end) then
            t = tin - this%rho_start
            rho = this%initial + (this%first * t) + this%second
        else
            t = this%rho_end - this%rho_start
            rho = this%initial + (this%first * t) + this%second
        end if
    
    end function Reactivity_ramp
    
    !$
    !===============================================================================================
    ! oscillatory reactivity function, use three parameter
    ! rho = a*sin(b*t) + c
    !===============================================================================================
    function Reactivity_oscillatory(this, tin) result(rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in) :: tin
        real(KREAL)             :: rho
        
        real(KREAL)  :: t
        
        if (tin < this%rho_start)  then
            rho = this%initial 
        else if (tin <= this%rho_end) then
            t = tin - this%rho_start
            rho = this%initial + this%first * SIN(this%second*t) + this%third
        else
            rho = this%initial
        end if

    end function Reactivity_oscillatory
    
    !$
    !===============================================================================================
    ! self define, used for extendsion
    ! rho = ----
    !===============================================================================================
    function Reactivity_special(this, tin) result(rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in) :: tin
        real(KREAL)             :: rho
        
        real(KREAL)  :: t
        
        if (tin < this%rho_start)  then
            rho = this%initial
        else if (tin <= this%rho_end) then
            t = tin - this%rho_start
            rho = this%initial
        else
            rho = this%initial
        end if

    end function Reactivity_special
    
    !$
    !===============================================================================================
    ! reactivity define by feedback
    ! rho = ----
    !===============================================================================================
    function Reactivity_feedback(this, tin) result(rho)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in) :: tin
        real(KREAL)             :: rho
        
        real(KREAL)  :: rho_fuel
        real(KREAL)  :: rho_coolant
        real(KREAL)  :: rho_axial
        real(KREAL)  :: rho_radial
        
        real(KREAL)  :: AD
        real(KREAL)  :: Bv
        real(KREAL)  :: AXIAL
        real(KREAL)  :: RADIAL
        
        AD = this%AD
        Bv = this%Bv
        AXIAL  = this%AXIAL
        RADIAL = this%RADIAL
        
        rho_fuel = (AD / this%T_fuel) * (this%T_fuel - this%T0_fuel)
        rho_coolant = (Bv * 100.0 / this%Rho0_coolant) * (this%Rho_coolant - this%Rho0_coolant)
        rho_axial = AXIAL * (this%T_fuel - this%T0_fuel)
        rho_radial = RADIAL * (this%T_coolant - this%T0_coolant)
        
        if (tin < this%rho_start)  then
            rho = 0.0
        else if (tin <= this%rho_end) then
            rho = rho_fuel + rho_coolant + rho_axial + rho_radial
        else
            rho = rho_fuel + rho_coolant + rho_axial + rho_radial
        end if
        
!        write(152, fmt="(1x, F13.6, *(TR3, ES13.6))")  tin, rho_fuel, rho_coolant, rho_axial, rho_radial, rho

    end function Reactivity_feedback
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_initial_feedback (this, T0_fuel, T0_coolant, Rho0_coolant)
    
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in)  :: T0_fuel
        real(KREAL), intent(in)  :: T0_coolant
        real(KREAL), intent(in)  :: Rho0_coolant
        
        this%T0_fuel = T0_fuel
        this%T0_coolant = T0_coolant
        this%Rho0_coolant = Rho0_coolant
        
        this%T_fuel = T0_fuel
        this%T_coolant = T0_coolant
        this%Rho_coolant = Rho0_coolant
    
    end subroutine Set_initial_feedback
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_feedback_parameter (this, T_fuel, T_coolant, Rho_coolant)
    
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in)  :: T_fuel
        real(KREAL), intent(in)  :: T_coolant
        real(KREAL), intent(in)  :: Rho_coolant
        
        this%T_fuel = T_fuel
        this%T_coolant = T_coolant
        this%Rho_coolant = Rho_coolant
    
    end subroutine Set_feedback_parameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_feedback_coefficient (this, rho_dim)
        
        class(PKReactivity), intent(in out)  :: this
        real(KREAL), intent(in)   :: rho_dim(:)
        
        this%AD = rho_dim(1)
        this%Bv = rho_dim(2)
        this%AXIAL = rho_dim(3)
        this%RADIAL = rho_dim(4)
    
    end subroutine Set_feedback_coefficient
    
end module pkreactivity_header
