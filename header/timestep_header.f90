!$
!===================================================================================================
!
!   class of time step for performing transient
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          TimeSteper
!                               EqualSteper
!                               TimeStepInfo
!
!===================================================================================================
module timestep_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : WarningCollector
    
    implicit none
    private 
    public :: TimeSteper, EqualSteper, TimeStepInfo
    
    type(WarningCollector)  :: a_warning
    
    ! --------------------------------------------------------------------------
    ! type for time interval
    type, private  :: section_tp
        real(KREAL), public    :: point                                         ! ending point of a time interval
        real(KREAL), public    :: length                                        ! step length of this interval
    end type section_tp
    
    type, private  :: interval_tp
        real(KREAL), public    :: tend                                          ! ending point of a time interval 
        real(KREAL), public    :: dmacro                                        ! macro step length of this interval 
        real(KREAL), public    :: dmiddle                                       ! middle step length of this interval 
    end type interval_tp 
    
    ! --------------------------------------------------------------------------
    ! type for time step 
    type  TimeSteper
        real(KREAL), public                   :: start_time     = REAL_ZERO     ! start time of the whole transient process
        real(KREAL), public                   :: end_time       = REAL_ZERO     ! final time of the whole transient process
        integer, public, allocatable          :: step_per_section(:)            ! how many time step per section
        integer, public                       :: n_section      = INT_ZERO      ! how many section for the whole transient process
        integer, public                       :: n_interval     = INT_ZERO      ! how many interval for the whole transient process
        integer, public                       :: n_step         = INT_ZERO      ! how many step for the whole transient process
        type(section_tp), public, allocatable :: sections(:)                    ! information per section
        integer, private                      :: i_section      = INT_ZERO
        integer, private                      :: i_step         = INT_ZERO
        
        type(interval_tp), public, allocatable :: intervals(:)
    contains
        procedure, public  :: alloc => Allocate_TimeSteper
        procedure, public  :: clean => Free_TimeSteper
        procedure, public  :: set_section => Set_TimeSteper_by_section
        procedure, public  :: set_interval => Set_TimeSteper_by_interval 
        
        procedure, public  :: get_start => Get_start_time
        procedure, public  :: get_end => Get_end_time
        procedure, public  :: get_section_count => Get_section_number
        procedure, public  :: get_step_count => Get_step_number
        procedure, public  :: info => Get_TimeSteper_info
        procedure, public  :: check => Check_step_index
    end type TimeSteper
    
    ! type for equal interval
    type  EqualSteper
        integer, public                       :: nstep     = 1                  ! step number of the whole transient
        real(KREAL), public                   :: stepsize  = REAL_ZERO          ! step size of the whole transient
        real(KREAL), public                   :: total     = REAL_ZERO          ! current ccumulative time
        real(KREAL), public                   :: left      = REAL_ZERO          ! 
        real(KREAL), public                   :: right     = REAL_ZERO          ! time point to next step
        logical, public                       :: is_init   = .FALSE.            ! is initialized ?
    contains
        procedure, public  :: add => Add_EqualSteper_step
    end type EqualSteper
    
    ! type for a step info
    type  TimeStepInfo
        integer, public          :: index                                       ! step total index
        real(KREAL), public      :: left                                        ! left point of a step
        real(KREAL), public      :: right                                       ! right point of a step
        real(KREAL), public      :: pace                                        ! pace length of a step
    contains
        procedure, public  :: print => Print_TimeStepInfo
    end type TimeStepInfo
    
    ! --------------------------------------------------------------------------
    ! private the real function name
contains
    !$
    !===============================================================================================
    ! allocate the section information of TimeSteper
    !===============================================================================================
    subroutine Allocate_TimeSteper (this)
        
        class(TimeSteper), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%step_per_section(this%n_section), stat=i_allocate)
        allocate(this%sections(this%n_section), stat=i_allocate)
        allocate(this%intervals(this%n_interval), stat=i_allocate)
        
        this%step_per_section = INT_ZERO
        this%sections = section_tp(REAL_ZERO, REAL_ZERO)
        this%intervals = interval_tp(REAL_ZERO, REAL_ZERO, REAL_ZERO)
        
    end subroutine Allocate_TimeSteper
    
    !$
    !===============================================================================================
    ! finalizer for class of TimeSteper
    !===============================================================================================
    subroutine Free_TimeSteper (this)
        
        class(TimeSteper), intent(in out)  :: this
        
        if (allocated(this%step_per_section))   deallocate(this%step_per_section)
        if (allocated(this%sections))           deallocate(this%sections)
        if (allocated(this%intervals))          deallocate(this%intervals)
    
    end subroutine Free_TimeSteper
    
    !$
    !===============================================================================================
    ! set total step
    !===============================================================================================
    subroutine Set_TimeSteper_by_section (this)
        
        class(TimeSteper), intent(in out)  :: this
        real(KREAL), parameter  :: EPS_ = 5.0D-6                                ! fix the last digit space
        integer                 :: it
        
        ! obtain step number per section and the total step number
        this%n_step = 0
        do it = 1, this%n_section
            if (it == 1) then
                this%step_per_section(it) = FLOOR((this%sections(it)%point-this%start_time)*(1.0D0+EPS_)/this%sections(it)%length)
                this%n_step = this%n_step + this%step_per_section(it)
            else
                this%step_per_section(it) = FLOOR((this%sections(it)%point-this%sections(it-1)%point)*(1.0D0+EPS_)/this%sections(it)%length)
                this%n_step = this%n_step + this%step_per_section(it)
            end if
        end do
    
    end subroutine Set_TimeSteper_by_section
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Set_TimeSteper_by_interval (this)
    
        class(TimeSteper), intent(in out)  :: this
        
        real(KREAL)  :: point, length 
        integer  :: it, iit, idx 
        integer  :: per_interval 
        integer  :: i_allocate 
        
        this%n_step = 0; this%n_section = 0;
        do it = 1, this%n_interval 
            if (it == 1)  then
                this%n_section = this%n_section + NINT((this%intervals(it)%tend-this%start_time)/this%intervals(it)%dmacro)
                this%n_step = this%n_step + NINT((this%intervals(it)%tend-this%start_time)/this%intervals(it)%dmiddle)
            else
                this%n_section = this%n_section + NINT((this%intervals(it)%tend-this%intervals(it-1)%tend)/this%intervals(it)%dmacro)
                this%n_step = this%n_step + NINT((this%intervals(it)%tend-this%intervals(it-1)%tend)/this%intervals(it)%dmiddle)
            end if 
        end do 
        
        if (allocated(this%step_per_section))   deallocate(this%step_per_section)
        if (allocated(this%sections))           deallocate(this%sections)
        allocate(this%step_per_section(this%n_section), stat=i_allocate)
        allocate(this%sections(this%n_section), stat=i_allocate)
        
        idx = 0 
        point = 0.0D0; length = 0.0D0;
        do it = 1, this%n_interval
            if (it == 1)  then
                per_interval = NINT((this%intervals(it)%tend-this%start_time)/this%intervals(it)%dmacro)
            else 
                per_interval = NINT((this%intervals(it)%tend-this%intervals(it-1)%tend)/this%intervals(it)%dmacro)
            end if 
            
            do iit = 1, per_interval
                idx = idx + 1 
                point = point + this%intervals(it)%dmacro
                length = this%intervals(it)%dmiddle
                
                this%sections(idx)%point = point
                this%sections(idx)%length = length 
            end do 
        end do 
        
        call this%set_section ()
        
    end subroutine Set_TimeSteper_by_interval
    
    
    !$
    !===============================================================================================
    ! get start time
    !===============================================================================================
    function Get_start_time (this)  result(point)
        
        class(TimeSteper), intent(in)  :: this
        real(KREAL)  :: point
        
        point = this%start_time
    
    end function Get_start_time
    
    !$
    !===============================================================================================
    ! get end time
    !===============================================================================================
    function Get_end_time (this)  result(point)
        
        class(TimeSteper), intent(in)  :: this
        real(KREAL)  :: point
        
        point = this%end_time
    
    end function Get_end_time
    
    !$
    !===============================================================================================
    ! get time section number
    !===============================================================================================
    function Get_section_number(this)  result(cnt)
        
        class(TimeSteper), intent(in)  :: this
        integer  :: cnt
    
        cnt = SIZE(this%sections)
        
    end function Get_section_number
    
    !$
    !===============================================================================================
    ! get time step number
    !===============================================================================================
    function Get_step_number(this)  result(cnt)
        
        class(TimeSteper), intent(in)  :: this
        integer  :: cnt
    
        cnt = this%n_step
        
    end function Get_step_number
    
    !$
    !===============================================================================================
    ! default 'is_step' is .TRUE.
    !===============================================================================================
    function Get_TimeSteper_info (this, tid, is_step)  result(a_step)
        
        class(TimeSteper), intent(in out)  :: this
        integer, intent(in)                :: tid
        logical, intent(in), optional      :: is_step
        type(TimeStepInfo)  :: a_step
        
        integer  :: lower, upper, excess
        integer  :: it
        logical  :: is_ending
        
        if (tid == 0)  then               
            a_step%index = tid
            a_step%right = this%start_time
            a_step%left = this%start_time
            a_step%pace = a_step%right - a_step%left 
            return
        end if
        
        if (PRESENT(is_step) .and. (.NOT. is_step))  then
            call this%check (tid, is_current=.TRUE., is_step=.FALSE.)
            this%i_section = tid
            a_step%right = this%sections(tid)%point
            if (tid == 1)  then
                a_step%pace = this%sections(tid)%point - this%start_time
            else
                a_step%pace = this%sections(tid)%point - this%sections(tid-1)%point
            end if
            
        else
            call this%check (tid, is_current=.TRUE., is_step=.TRUE.)
            this%i_step = tid
            
            upper = 0; lower = 0
            do it = 1, this%n_section
                upper = upper + this%step_per_section(it)
                lower = upper - this%step_per_section(it)
                if (tid < upper)  then
                    is_ending = .FALSE.
                    excess = tid - lower
                    exit
                else if (tid == upper) then
                    is_ending = .TRUE.
                    excess = 0
                    exit
                end if
            end do
            
            if (it == 1)  then
                if (is_ending)  then
                    a_step%right = this%sections(it)%point
                    if (this%step_per_section(it) == 1)  then
                        a_step%pace = this%sections(it)%point - this%start_time
                    else
                        a_step%pace = (this%sections(it)%point - this%start_time) - (this%step_per_section(it)-1)*this%sections(it)%length
                    end if
                else
                    a_step%right = this%start_time + this%sections(it)%length * excess
                    a_step%pace = this%sections(it)%length
                end if
            else
                if (is_ending)  then
                    a_step%right = this%sections(it)%point
                    if (this%step_per_section(it) == 1)  then
                        a_step%pace = this%sections(it)%point - this%sections(it-1)%point
                    else
                        a_step%pace = (this%sections(it)%point - this%sections(it-1)%point) - (this%step_per_section(it)-1)*this%sections(it)%length
                    end if
                else
                    a_step%right = this%sections(it-1)%point + this%sections(it)%length * excess
                    a_step%pace = this%sections(it)%length  
                end if
            end if
        end if
        
        a_step%index = tid
        a_step%left = a_step%right - a_step%pace
        
    end function Get_TimeSteper_info
    
    !$
    !===============================================================================================
    ! check for time step index
    !===============================================================================================
    subroutine Check_step_index (this, tid, is_current, is_step)
        
        class(TimeSteper), intent(in)  :: this
        integer, intent(in)  :: tid                                             ! index for current time step
        logical, intent(in)  :: is_current                                      ! return current or next step information
        logical, intent(in), optional  :: is_step                               ! return step or section information
        
        ! for current interval
        if (is_current)  then
            if (PRESENT(is_step) .AND. .NOT. is_step)  then
                if (tid < 0  .or. tid > this%n_section)  then
                    call a_warning%set (INFO_LIST_TIMESTEP, 'time step index out of range')
                    call a_warning%print (OUTPUT_UNIT)
                end if 
            else
                if (tid < 0  .or. tid > this%n_step)  then
                    call a_warning%set (INFO_LIST_TIMESTEP, 'time section index out of range')
                    call a_warning%print (OUTPUT_UNIT)
                end if 
            end if
            
        ! for next interval
        else 
            if (PRESENT(is_step) .AND. .NOT. is_step)  then
                if (tid < 0  .or. tid >= this%n_section)  then
                    call a_warning%set (INFO_LIST_TIMESTEP, 'time step index out of range')
                    call a_warning%print (OUTPUT_UNIT)
                end if 
            else
                if (tid < 0  .or. tid >= this%n_step)  then
                    call a_warning%set (INFO_LIST_TIMESTEP, 'time section index out of range')
                    call a_warning%print (OUTPUT_UNIT)
                end if 
            end if
        end if 
    
    end subroutine Check_step_index
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Add_EqualSteper_step (this, in_step, out_step, is_advance)
        
        class(EqualSteper), intent(in out) :: this
        type(TimeStepInfo), intent(in)     :: in_step
        type(TimeStepInfo), intent(out)    :: out_step
        logical, intent(out)               :: is_advance
        
        if (.NOT. this%is_init)  then
            this%total = REAL_ZERO
            this%left = REAL_ZERO
            this%right = this%stepsize
            this%is_init = .TRUE.
        end if
        
        this%total = this%total + (in_step%right - in_step%left)
        
        is_advance = .FALSE.
        out_step%left = this%left
        out_step%right = this%right
        
        if (ABS(this%total - this%right) <= EPS_EQUAL)  then
            is_advance = .TRUE.
            this%left = this%right
            this%right = this%right + this%stepsize
           
        ! warning
        else if (this%total > this%right)  then
            
        end if
    
    end subroutine Add_EqualSteper_step
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_TimeStepInfo (this, unit_)
        
        class(TimeStepInfo), intent(in)  :: this
        integer, intent(in)              :: unit_
        
        write(unit=unit_, fmt="(1x, I4, *(TR3, F11.6))")  this%index, this%left, this%right, this%pace
    
    end subroutine Print_TimeStepInfo
    
end module timestep_header
