!$
!===================================================================================================
!
!   parsing command line parameters
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Parse_command_line
!
!   Public type lists:          No
!
!===================================================================================================
module command_line
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
        
    implicit none 
    private
    public  :: Parse_command_line

contains
    !$
    !===============================================================================================
    ! parsing command line parameters, if command line error, then stop
    !===============================================================================================
    subroutine Parse_command_line (is_getting, case_file)
        
        logical, intent(out)  :: is_getting                                     ! is getting input file name from command line
        character(len=MAX_WORD_LEN), intent(in out)  :: case_file

        character(len=MAX_WORD_LEN), allocatable  :: argv(:)                    ! holder for arguments
        integer  :: argc                                                        ! number of argument
        integer  :: i
        
        ! initialize to not getting
        is_getting = .FALSE.
        
        argc = COMMAND_ARGUMENT_COUNT ()
        if (argc == 0)  then
            return
        end if
        
        ! get command line parameter
        allocate(argv(argc))
        do i = 1, argc
            call GET_COMMAND_ARGUMENT (i, argv(i))
        end do
        
        ! parse argument
        i = 1
        argument: do 
            if (INDEX(TRIM(ADJUSTL(argv(i))), '-') == 1)  then
                select case (TRIM(ADJUSTL(argv(i))))
                case ('-v', '-version', '--version')
                    call Print_version (OUTPUT_UNIT)
                    stop
                case ('-?', '-help', '--help')
                    call Print_usage (OUTPUT_UNIT)
                    stop
                
                ! unknown arguments, stop
                case default
                    call Print_unknow (OUTPUT_UNIT, argv(i))
                    stop
                end select
                
            else if (INDEX(TRIM(ADJUSTL(argv(i))), 'input.') == 1)  then
                is_getting = .TRUE.
                case_file = TRIM(ADJUSTL(argv(i)))
                
            ! error arguments, stop
            else 
                call Print_unknow (OUTPUT_UNIT, argv(i))
                stop
            end if

            i = i + 1
            if (i > argc) then
                exit argument
            end if
        end do argument
        
        deallocate(argv)
                
    end subroutine Parse_command_line
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! print current version and copyright information
    !===============================================================================================
    subroutine Print_version (unit_)
        
        integer, intent(in)  :: unit_
        
        integer  :: to_write
        
        to_write = unit_
        
        write(unit=to_write, fmt="(1x, /)")
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, A, ' version: ', I1, '.', I1, I1)") TRIM(CODENAME), VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
        write(unit=to_write, fmt="(1x, 2x, A)") TRIM(COPYRIGHT)
        write(unit=to_write, fmt="(1x, 2x, A)") TRIM(LABORATORY)
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, /)")
    
    end subroutine Print_version
    
    !$
    !===============================================================================================
    ! print the usage of command line parameters
    !===============================================================================================
    subroutine Print_usage (unit_)
        
        integer, intent(in)  :: unit_
        
        integer  :: to_write
        
        to_write = unit_
        
        write(unit=to_write, fmt="(1x, /)")
        write(unit=to_write, fmt="(1x, 'Usage: ', A, ' [option]')")  TRIM(CODENAME)
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, A)")  'Options:'
        write(unit=to_write, fmt="(1x, A)")  'input.xxx     :       use input.xxx as input files'
        write(unit=to_write, fmt="(1x, A)")  '-v, --version :       version information'
        write(unit=to_write, fmt="(1x, A)")  '-?, --help    :       help message'
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, /)")
    
    end subroutine Print_usage
    
    !$
    !===============================================================================================
    ! treat unknow command line parameters
    !===============================================================================================
    subroutine Print_unknow (unit_, argu)
        
        integer, intent(in)           :: unit_
        character(len=*), intent(in)  :: argu
        
        integer  :: to_write
        
        to_write = unit_
        
        write(unit=to_write, fmt="(1x, /)")
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, 'Unknown command line option: ', A)")  TRIM(argu)
        write(unit=to_write, fmt="(1x, 'type [', A, ' -help] for usage!')")  TRIM(CODENAME)
        write(unit=to_write, fmt="(1x, '------------------------------------------------------')")
        write(unit=to_write, fmt="(1x, /)")
    
    end subroutine Print_unknow
    
end module command_line
