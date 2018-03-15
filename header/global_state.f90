!$
!===================================================================================================
!
!   module for global state parameters and auxiliary parameters
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module global_state
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use state_header,           only : SteadyState, TransientState
    use timer_header,           only : EventTimeCounter, TotalTimeCounter, CPUTimeCounter, DateHolder
    use exception_header,       only : WarningCollector, ErrorCollector, InformationCollector
    
    use timestep_header,        only : TimeSteper
    
    implicit none
    public
    
    ! --------------------------------------------------------------------------
    ! state parameter
    type(SteadyState)                       :: ns                               ! steady parameter and flag
    type(TransientState)                    :: nt                               ! transient parameter and flag
    
    ! --------------------------------------------------------------------------
    ! count for program time
    type(TotalTimeCounter)                  :: time_program                     ! real time for whole problem
    type(CPUTimeCounter)                    :: cputime_program                  ! cpu time for whole problem
    type(DateHolder)                        :: current_date                     ! hold for current data

    type(EventTimeCounter)                  :: time_event
    
    ! --------------------------------------------------------------------------
    ! warning, error and information
    type(WarningCollector)                  :: this_warning
    type(ErrorCollector)                    :: this_error
    type(InformationCollector)              :: this_information
    
    ! --------------------------------------------------------------------------
    type(TimeSteper)                        :: time_step
    
end module global_state
