!$
!===================================================================================================
!
!   class for exception (error, warning & information) treatment
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          ErrorCollector
!                               WarningCollector
!                               InformationCollector
!
!===================================================================================================
module exception_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private
    public  :: ErrorCollector, WarningCollector, InformationCollector
    
    ! --------------------------------------------------------------------------
    ! type for collection error
    type  ErrorCollector
        character(len=MAX_WORD_LEN)    :: info_list                             ! identification for warning type
        character(len=MAX_WORD_LEN)    :: message                               ! message information for warning type
        integer                        :: counting   = 0                        ! count for message seting
        integer                        :: cnt_limit  = 1                        ! upper limit of message setting
    contains
        procedure, public  :: set => Set_ErrorCollector
        procedure, public  :: print => Print_ErrorCollector
    end type ErrorCollector
    
    ! type for collectiong warning
    type  WarningCollector
        character(len=MAX_WORD_LEN)    :: info_list                             
        character(len=MAX_WORD_LEN)    :: message                               
        integer                        :: counting   = 0                      
        integer                        :: cnt_limit  = 50000                    
    contains
        procedure, public  :: set  => Set_WarningCollector
        procedure, public  :: print => Print_WarningCollector 
    end type WarningCollector
    
    ! type for collection information
    type  InformationCollector
        character(len=MAX_WORD_LEN)    :: info_list                             
        character(len=MAX_WORD_LEN)    :: message                               
        integer                        :: counting   = 0                      
        integer                        :: cnt_limit  = 50000                    
    contains
        procedure, public  :: set => Set_InformationCollector
        procedure, public  :: print => Print_InformationCollector
    end type InformationCollector
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Set_ErrorCollector, Print_ErrorCollector
    private  :: Set_WarningCollector, Print_WarningCollector
    private  :: Set_InformationCollector, Print_InformationCollector
    
contains
    !$
    !===============================================================================================
    ! set error collector information (info_list and message)
    !===============================================================================================
    subroutine Set_ErrorCollector (this, info_list, message)
        
        class(ErrorCollector), intent(in out)   :: this
        character(len=*), intent(in)  :: info_list
        character(len=*), intent(in)  :: message
        
        this%counting = this%counting + 1
        this%info_list = TRIM(ADJUSTL(info_list))
        this%message = TRIM(ADJUSTL(message))
    
    end subroutine Set_ErrorCollector
    
    !$
    !===============================================================================================
    ! set warning collector information (info_list and message)
    !===============================================================================================
    subroutine Set_WarningCollector (this, info_list, message)
        
        class(WarningCollector), intent(in out)  :: this
        character(len=*), intent(in)  :: info_list
        character(len=*), intent(in)  :: message
        
        this%counting = this%counting + 1
        this%info_list = TRIM(ADJUSTL(info_list))
        this%message = TRIM(ADJUSTL(message))
    
    end subroutine Set_WarningCollector

    !$
    !===============================================================================================
    ! set information collector information (info_list and message)
    !===============================================================================================
    subroutine Set_InformationCollector (this, info_list, message)
        
        class(InformationCollector), intent(in out)  :: this
        character(len=*), intent(in)  :: info_list
        character(len=*), intent(in)  :: message
        
        this%counting = this%counting + 1
        this%info_list = TRIM(ADJUSTL(info_list))
        this%message = TRIM(ADJUSTL(message))
    
    end subroutine Set_InformationCollector
    
    !$
    !===============================================================================================
    ! Print error message (file & stdout) and abort the program
    !===============================================================================================
    subroutine Print_ErrorCollector (this, unit_)
    
        ! intent parmeters
        class(ErrorCollector), intent(in)  :: this
        integer, intent(in)  :: unit_
        
        integer  :: to_write
        to_write = unit_
        
        ! output error signal
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        write(unit=to_write, fmt="(1X, '@-Error:   ', A)")       TRIM(this%info_list)
        write(unit=to_write, fmt="(1X, '           ', A, ' !')") TRIM(this%message) 
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        
        if (unit_ /= OUTPUT_UNIT)  then
            to_write = OUTPUT_UNIT
            write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
            write(unit=to_write, fmt="(1X, '@-Error:   ', A)")       TRIM(this%info_list)
            write(unit=to_write, fmt="(1X, '           ', A, ' !')") TRIM(this%message) 
            write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        end if
        
        if (this%counting >= this%cnt_limit)  then
            write(unit=to_write, fmt="(1x, A)") ' '
            write(unit=to_write, fmt="(1x, A)") 'Too many error-s, please check carefully'
            write(unit=to_write, fmt="(1x, A)") ' '
            stop
        end if
        
    end subroutine Print_ErrorCollector
        
    !$
    !===============================================================================================
    ! Print warning message (file only) and go on the program if within limit
    !===============================================================================================
    subroutine Print_WarningCollector (this, unit_)
        
        ! intent parmeters
        class(WarningCollector), intent(in)  :: this
        integer, intent(in)  :: unit_
        
        integer  :: to_write
        to_write = unit_
        
        ! output warning signal
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        write(unit=to_write, fmt="(1X, '@-Warning: ', A)")       TRIM(this%info_list)
        write(unit=to_write, fmt="(1X, '           ', A, ' !')") TRIM(this%message)
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        
        if (this%counting >= this%cnt_limit)  then
            write(unit=to_write, fmt="(1x, A)") ' '
            write(unit=to_write, fmt="(1x, A)") 'To many warning-s, please check carefully'
            write(unit=to_write, fmt="(1x, A)") ' '
            stop
        end if
        
    end subroutine Print_WarningCollector
    
    !$
    !===============================================================================================
    ! Print information message (file only) and go on the program if within limit
    !===============================================================================================
    subroutine Print_InformationCollector (this, unit_)
    
        ! intent parmeters
        class(InformationCollector), intent(in)  :: this
        integer, intent(in)  :: unit_
        
        integer  :: to_write
        to_write = unit_
        
        ! output error signal
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        write(unit=to_write, fmt="(1X, '@-Info:    ', A)")       TRIM(this%info_list)
        write(unit=to_write, fmt="(1X, '           ', A, ' !')") TRIM(this%message) 
        write(unit=to_write, fmt="(1X, '-----------------------------------------------')")
        
        if (this%counting >= this%cnt_limit)  then
            write(unit=to_write, fmt="(1x, A)") ' '
            write(unit=to_write, fmt="(1x, A)") 'To many info-s, please check carefully'
            write(unit=to_write, fmt="(1x, A)") ' '
            stop
        end if
        
    end subroutine Print_InformationCollector
    
end module exception_header
    