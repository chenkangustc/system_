!$
!===================================================================================================
!
!   define a class for time counting
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          EventTimeCounter
!                               TotalTimeCounter
!                               CPUTimeCounter
!                               DateHolder
!
!===================================================================================================
module timer_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private
    public  :: EventTimeCounter, TotalTimeCounter, CPUTimeCounter, DateHolder
    
    ! --------------------------------------------------------------------------
    type  EventTimeCounter
        integer, public                         :: idx(100)        = -1
        character(MAX_WORD_LEN), public         :: name(100)       = ' '
        integer, public                         :: Ncount(100)     = 0 
        integer(KINT), public                   :: ibeg(100)       = 1
        integer(KINT), public                   :: iend(100)       = 1 
        real(KDOUBLE), public                   :: time(100)       = 0.0D0
        logical, public                         :: is_running(100) = .FALSE.
    contains
        procedure, public  :: set => Set_EventTimeCounter
        procedure, public  :: begin => Begin_EventTimeCounter
        procedure, public  :: end => End_EventTimeCounter
        procedure, public  :: print => Print_EventTimeCounter
    end type EventTimeCounter
    
    ! type for time counting by total time
    type  TotalTimeCounter
        logical                     :: is_running    = .FALSE.                  ! whether the timer is run ?
        real(KDOUBLE), public       :: total         = REAL_ZERO                ! total time since start
        real(KDOUBLE)               :: step          = REAL_ZERO                ! time for current step
        integer(KINT)               :: start_count   = INT_ZERO                 ! count for the start point of time step
        integer(KINT)               :: end_count     = INT_ZERO                 ! count for the end point of time step
    contains 
        procedure, public  :: start => Start_TimeCounter
        procedure, public  :: update => Update_TimeCounter
        procedure, public  :: stop => Stop_TimeCounter
        procedure, public  :: reset => Reset_TimeCounter
    end type TotalTimeCounter
    
    ! type for time counting by CPU time
    type  CPUTimeCounter
        logical                     :: is_running    = .FALSE.                  ! whether the timer is run ?
        real(KDOUBLE), public       :: total         = REAL_ZERO                ! total time since start
        real(KDOUBLE)               :: step          = REAL_ZERO                ! time for current step
        real(KDOUBLE)               :: start_point   = REAL_ZERO                ! start point of time step
        real(KDOUBLE)               :: end_point     = REAL_ZERO                ! end point of time step
    contains
        procedure, public  :: start => Start_CPUTimeCounter
        procedure, public  :: update => Update_CPUTimeCounter
        procedure, public  :: stop => Stop_CPUTimeCounter
        procedure, public  :: reset => Reset_CPUTimeCounter
    end type CPUTimeCounter
    
    ! type for keeping current date and time
    type DateHolder
        character(len=MAX_WORD_LEN)  :: date                                    ! current date
        character(len=MAX_WORD_LEN)  :: time                                    ! current time
    contains
        procedure, public  :: print => Print_DateHolder
    end type DateHolder
    
contains
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Set_EventTimeCounter (this, idx, name)
        
        class(EventTimeCounter), intent(in out)  :: this
        integer, intent(in)                      :: idx 
        character(len=*), intent(in)             :: name 
        
        if ((idx >= SIZE(this%idx)) .OR. (idx <= 0))  then
            return
        end if 
        
        this%idx(idx) = idx 
        this%name(idx) = name 
        this%Ncount(idx) = 0 
        this%ibeg(idx) = 0
        this%iend(idx) = 0
        this%time(idx) = 0.0D0
        this%is_running(idx) = .FALSE.
        
    end subroutine Set_EventTimeCounter

    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Begin_EventTimeCounter (this, idx)
        
        class(EventTimeCounter), intent(in out)  :: this
        integer, intent(in)                      :: idx 
        integer(KINT)  :: count_rate
        
        call SYSTEM_CLOCK (this%ibeg(idx), count_rate)
        this%Ncount(idx) = this%Ncount(idx) + 1 
        this%is_running(idx) = .TRUE.
        
    end subroutine Begin_EventTimeCounter

    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine End_EventTimeCounter (this, idx)
    
        class(EventTimeCounter), intent(in out)  :: this
        integer, intent(in)                      :: idx 
        integer(KINT)  :: count_rate
        
        call SYSTEM_CLOCK (this%iend(idx), count_rate)
        
        if (this%is_running(idx))  then
            this%time(idx) = this%time(idx) + (real(this%iend(idx), KDOUBLE) - real(this%ibeg(idx), KDOUBLE)) / real(count_rate, KDOUBLE)
            this%is_running(idx) = .FALSE.
        end if 
        
    end subroutine End_EventTimeCounter

    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Print_EventTimeCounter (this, unit_)
        
        class(EventTimeCounter), intent(in)  :: this
        integer, intent(in)  :: unit_ 
        
        integer  :: i
        
        write(unit=unit_, fmt='(1x, A)') ' '
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        write(unit=unit_, fmt='(1x, A)') '  Idx     Name                      Num          Time         Avg-time         '
        write(unit=unit_, fmt='(1x, A)') '-------------------------------------------------------------------------------' 
        do i = 1, SIZE(this%idx)
            if (this%idx(i) > 0)  then
            write(unit=unit_, fmt='(1x, I4, TR3, A26, I6, *(TR5, ES12.5))')  this%idx(i), this%name(i), this%Ncount(i), this%time(i), this%time(i)/this%Ncount(i)
            end if 
        end do 
        write(unit=unit_, fmt='(1x, A)') ' '
        
        if (unit_ /= OUTPUT_UNIT)  then
            write(unit=OUTPUT_UNIT, fmt='(1x, A)') ' '
            write(unit=OUTPUT_UNIT, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
            write(unit=OUTPUT_UNIT, fmt='(1x, A)') '  Idx     Name                      Num          Time         Avg-time         '
            write(unit=OUTPUT_UNIT, fmt='(1x, A)') '-------------------------------------------------------------------------------' 
            do i = 1, SIZE(this%idx)
                if (this%idx(i) > 0)  then
                write(unit=OUTPUT_UNIT, fmt='(1x, I4, TR3, A26, I6, *(TR5, ES12.5))')  this%idx(i), this%name(i), this%Ncount(i), this%time(i), this%time(i)/this%Ncount(i)
                end if 
            end do 
            write(unit=OUTPUT_UNIT, fmt='(1x, A)') ' '
        end if 
        
    end subroutine Print_EventTimeCounter

    !$
    !===============================================================================================
    ! starts a TotalTimeCounter and gets the current time
    !===============================================================================================
    subroutine Start_TimeCounter (this)
    
        class(TotalTimeCounter), intent(in out)  :: this
        
        this%is_running = .TRUE.
        this%total = REAL_ZERO
        this%step  = REAL_ZERO
        
        call SYSTEM_CLOCK(this%end_count)
        this%start_count = this%end_count
        
    end subroutine Start_TimeCounter

    !$
    !===============================================================================================
    ! get the total elapse time and the current step time, print to file & stdout
    !===============================================================================================
    subroutine Update_TimeCounter (this, unit_, silent)
    
        class(TotalTimeCounter), intent(in out)  :: this
        integer, intent(in), optional            :: unit_ 
        logical, intent(in), optional            :: silent 
        integer(KINT)  :: count_rate
        integer        :: to_write
        
        if (this%is_running) then
            call SYSTEM_CLOCK(this%end_count, count_rate)
            this%step = (real(this%end_count, KDOUBLE) - real(this%start_count, KDOUBLE)) / real(count_rate, KDOUBLE)
            this%start_count = this%end_count
        else 
            this%step = REAL_ZERO
        end if 
        this%total = this%total + this%step
        
        if (PRESENT(silent) .and. silent)  then 
            return 
        end if 
        
        ! where to write the information
        if (PRESENT(unit_) .and. unit_ /= OUTPUT_UNIT) then
            to_write = unit_
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
            write(unit=to_write, fmt="(1x, 'Current step total-time elapse is:', F11.3, ' seconds')")  this%step
            write(unit=to_write, fmt="(1x, 'Accumulative total-time elapse is:', F11.3, ' seconds')")  this%total
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        end if
        
        ! print the accumulate time and the nearest step time
        to_write = OUTPUT_UNIT
        write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        write(unit=to_write, fmt="(1x, 'Current step total-time elapse is:', F11.3, ' seconds')")  this%step
        write(unit=to_write, fmt="(1x, 'Accumulative total-time elapse is:', F11.3, ' seconds')")  this%total
        write(unit=to_write, fmt="(1x, '----------------------------------------------------')")        
        
    end subroutine Update_TimeCounter
    
    !$
    !===============================================================================================
    ! update TotalTimeCounter and stop it
    !===============================================================================================
    subroutine Stop_TimeCounter (this, unit_)
    
        class(TotalTimeCounter), intent(in out)  :: this
        integer, intent(in), optional            :: unit_                         ! the optional property is handle in the 'update' subroutine
        
        if (.NOT. this%is_running)  return
        
        call this%update (unit_)
        
        this%is_running = .FALSE.
        
    end subroutine Stop_TimeCounter
    
    !$
    !===============================================================================================
    ! reset TotalTimeCounter to initial status
    !===============================================================================================
    subroutine Reset_TimeCounter (this)
    
        class(TotalTimeCounter), intent(in out)      :: this
        
        this%is_running  = .FALSE.
        this%total = REAL_ZERO
        this%step  = REAL_ZERO
        
        this%start_count = INT_ZERO
        this%end_count   = INT_ZERO
        
    end subroutine Reset_TimeCounter
    
    !$
    !===============================================================================================
    ! starts a CPUTimeCounter and gets the current time
    !===============================================================================================
    subroutine Start_CPUTimeCounter (this)
    
        class(CPUTimeCounter), intent(in out)  :: this
        
        this%is_running = .TRUE.
        this%total = REAL_ZERO
        this%step  = REAL_ZERO
        
        call CPU_TIME(this%end_point)
        
        this%start_point = this%end_point
        
    end subroutine Start_CPUTimeCounter

    !$
    !===============================================================================================
    ! get the total elapse time and the current step time, print to file & stdout
    !===============================================================================================
    subroutine Update_CPUTimeCounter (this, unit_, silent)
    
        class(CPUTimeCounter), intent(in out)  :: this
        integer, intent(in), optional          :: unit_ 
        logical, intent(in), optional          :: silent 
        integer  :: to_write
        
        if (this%is_running) then
            call CPU_TIME(this%end_point)
            this%step = this%end_point - this%start_point
            this%start_point = this%end_point
        else 
            this%step = REAL_ZERO
        end if 
        this%total = this%total + this%step       

        if (PRESENT(silent) .and. silent)  then 
            return 
        end if
        
        ! where to write the information
        if (PRESENT(unit_) .and. unit_ /= OUTPUT_UNIT) then
            to_write = unit_
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
            write(unit=to_write, fmt="(1x, 'Current step cpu-time elapse is: ', F11.3, ' seconds')")  this%step
            write(unit=to_write, fmt="(1x, 'Accumulative cpu-time elapse is: ', F11.3, ' seconds')")  this%total
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        end if
        
        ! print the accumulate time and the nearest step time
        to_write = OUTPUT_UNIT
        write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        write(unit=to_write, fmt="(1x, 'Current step cpu-time elapse is: ', F11.3, ' seconds')")  this%step
        write(unit=to_write, fmt="(1x, 'Accumulative cpu-time elapse is: ', F11.3, ' seconds')")  this%total
        write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        
    end subroutine Update_CPUTimeCounter
    
    !$
    !===============================================================================================
    ! update CPUTimeCounter and stop it
    !===============================================================================================
    subroutine Stop_CPUTimeCounter (this, unit_)
    
        class(CPUTimeCounter), intent(in out)  :: this
        integer, intent(in), optional          :: unit_                           ! the optional property is handle in the 'update' subroutine
        
        if (.NOT. this%is_running)  return
        
        call this%update (unit_)
        
        this%is_running = .FALSE.
        
    end subroutine Stop_CPUTimeCounter
    
    !$
    !===============================================================================================
    ! reset CPUTimeCounter to initial status
    !===============================================================================================
    subroutine Reset_CPUTimeCounter (this)
    
        class(CPUTimeCounter), intent(in out)      :: this
        
        this%is_running  = .FALSE.
        this%total = REAL_ZERO
        this%step  = REAL_ZERO
        
        this%start_point = REAL_ZERO
        this%end_point   = REAL_ZERO
        
    end subroutine Reset_CPUTimeCounter
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! print current date and time, print to file & stdout
    !===============================================================================================
    subroutine Print_DateHolder (this, unit_)
        
        class(DateHolder), intent(in out)  :: this
        integer, intent(in), optional      :: unit_
        integer  :: to_write
        
        call DATE_AND_TIME(this%date, this%time)
        
        ! print the current date and time by format
        if (PRESENT(unit_) .and. unit_ /= OUTPUT_UNIT)  then
            to_write = unit_
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
            write(unit=to_write, fmt="(1x, 6A)") 'Current date is:  ', this%date(1:4), '/', this%date(5:6), '/', this%date(7:8)
            write(unit=to_write, fmt="(1x, 6A)") 'Current time is:  ', this%time(1:2), ':', this%time(3:4), ':', this%time(5:6)
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        else
            to_write = OUTPUT_UNIT
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
            write(unit=to_write, fmt="(1x, 6A)") 'Current date is:  ', this%date(1:4), '/', this%date(5:6), '/', this%date(7:8)
            write(unit=to_write, fmt="(1x, 6A)") 'Current time is:  ', this%time(1:2), ':', this%time(3:4), ':', this%time(5:6)
            write(unit=to_write, fmt="(1x, '----------------------------------------------------')")
        end if
        
    end subroutine Print_DateHolder  
    
end module timer_header
