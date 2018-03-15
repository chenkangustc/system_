!$
!===================================================================================================
!
!   class for link parameter define
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          LinkParameter
!
!===================================================================================================
module link_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use Link
    
    implicit none 
    private
    public  :: LinkParameter
    
    ! --------------------------------------------------------------------------
    ! type for link information 
    type, private  :: link_parameter_tp
        character(len=MAX_WORD_LEN)          :: name          = ''              ! name of this parameter
        real(KREAL)                          :: reference     = 0.0D0           ! reference value
        real(KREAL)                          :: max           = 0.0D0           ! max value of case matrix for this parameter
        real(KREAL)                          :: min           = 0.0D0           ! min value
    end type link_parameter_tp
        
    ! --------------------------------------------------------------------------
    ! type for parameter used in linking:
    ! 1-BU, burnup (MW.D/t)
    ! 2-CB, boron concentration (ppm)
    ! 3-Tf, fuel temperature (K)
    ! 4-Tm, coolant temperature (K)
    ! 5-rho_m, coolant density (kg/m^3)
    ! --------------------------------------------------------------------------
    type  LinkParameter
        integer, public                               :: n_parameter            ! number of parameter above used
        type(link_parameter_tp), public, allocatable  :: info(:)                ! info
        integer, public                               :: link_method = 1
        integer, public                               :: npoint = 1             ! number of link section for burnup interpolation
        real(KREAL), allocatable                      :: burnup_point(:)        ! upper value of burnup section for interpolation
    contains
        procedure, public  :: set_burnup => Set_LinkParameter_point
        procedure, public  :: get_point => Get_LinkParameter_point
        procedure, public  :: set_type => Set_SelfLinking_type
        procedure, public  :: set_value => Set_SelfLinking_value
        procedure, public  :: set_max => Set_max_limitation
        procedure, public  :: set_min => Set_min_limitation
        procedure, public  :: is_active => Is_LinkParameter_used
        procedure, public  :: clean => Free_LinkParameter
    end type LinkParameter
        
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Set_LinkParameter_point, Get_LinkParameter_point, Set_SelfLinking_type, Set_SelfLinking_value
    private  :: Set_max_limitation, Set_min_limitation, Is_LinkParameter_used, Free_LinkParameter
        
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_LinkParameter_point (this, burnup_point)
        
        class(LinkParameter), intent(in out)  :: this
        real(KREAL), intent(in)  :: burnup_point(:)
        
        integer  :: i_allocate
        
        this%npoint = SIZE(burnup_point)
        if (allocated(this%burnup_point))   deallocate(this%burnup_point)
        allocate(this%burnup_point(this%npoint), stat=i_allocate)
        
        this%burnup_point = burnup_point
    
    end subroutine Set_LinkParameter_point

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function Get_LinkParameter_point (this, burnup)  result(ipoint)
        
        class(LinkParameter), intent(in)  :: this
        real(KREAL), intent(in)  :: burnup
        integer  :: ipoint
        integer  :: i 
        
        ! if NO burnup section input, ipoint = 1
        if (.NOT. allocated(this%burnup_point))  then 
            ipoint = 1
            return
        end if 
        
        if ((burnup <= 0.0) .OR. (burnup > this%burnup_point(this%npoint)))  then 
        end if 
        
        do i = 1, this%npoint
            if (this%burnup_point(i) >= burnup)  then 
                ipoint = i 
                exit 
            end if 
        end do 
    
    end function Get_LinkParameter_point
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_SelfLinking_type (this, types)
        
        class(LinkParameter), intent(in out)       :: this
        character(len=MAX_WORD_LEN), intent(in)    :: types(:)
        
        integer  :: i
        integer  :: i_allocate
        
        ! check allocated status
        call this%clean ()
        this%n_parameter = SIZE(types)
        allocate(this%info(this%n_parameter), stat=i_allocate)
        
        do i = 1, SIZE(types)
            this%info(i)%name = TRIM(ADJUSTL(types(i)))
        end do
    
    end subroutine Set_SelfLinking_type
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_SelfLinking_value (this, values)
        
        class(LinkParameter), intent(in out)   :: this
        real(KREAL), intent(in)                :: values(:)
        
        integer  :: i
               
        if (this%n_parameter /= SIZE(values))  then
        end if
        
        do i = 1, SIZE(values)
            this%info(i)%reference = values(i)
        end do
    
    end subroutine Set_SelfLinking_value
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_max_limitation (this, values)
        
        class(LinkParameter), intent(in out)   :: this
        real(KREAL), intent(in)                :: values(:)
        
        integer  :: i
               
        if (this%n_parameter /= SIZE(values))  then
        end if
        
        do i = 1, SIZE(values)
            this%info(i)%max = values(i)
        end do
    
    end subroutine Set_max_limitation
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_min_limitation (this, values)
        
        class(LinkParameter), intent(in out)   :: this
        real(KREAL), intent(in)                :: values(:)
        
        integer  :: i
               
        if (this%n_parameter /= SIZE(values))  then
        end if
        
        do i = 1, SIZE(values)
            this%info(i)%min = values(i)
        end do
    
    end subroutine Set_min_limitation
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function Is_LinkParameter_used (this, symbol)  result(is_true)
        
        class(LinkParameter), intent(in out)   :: this
        character(len=*) , intent(in)          :: symbol
        logical  :: is_true
        
        integer  :: i
        
        is_true = .FALSE.
        do i = 1, this%n_parameter
            if (TRIM(ADJUSTL(this%info(i)%name)) == symbol)  then
                is_true = .TRUE.
                exit
            end if
        end do
        
    end function Is_LinkParameter_used
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_LinkParameter (this)
        
        class(LinkParameter), intent(in out)   :: this
        
        if (allocated(this%info))           deallocate(this%info)
    
    end subroutine Free_LinkParameter
        
end module link_header
