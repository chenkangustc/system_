!$
!===================================================================================================
!
!   Pre process before the calculation begin: print tile, open file and so on.
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Run_pre_process
!
!   Public type lists:          No
!
!===================================================================================================
module driver_pre_process
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use command_line
    use files_dirs
    use environment,                only : Set_openMP
    use input_driver,               only : Driving_input
    
    implicit none
    private
    public  :: Run_pre_process

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_pre_process ()
    
        character(len=MAX_WORD_LEN)  :: file_main
        character(len=MAX_WORD_LEN)  :: case_file
        logical  :: is_getting
        
        ! initialize time counter
        call time_program%start ()
        call cputime_program%start ()
        
        call time_event%set (1, 'get Build Info')
        call time_event%set (2, 'get Nodal Coeff')
        call time_event%set (3, 'get Avg Flux')
        call time_event%set (4, 'get Out Moment')
        
        ! parsing command line argument, redict main input file
        call Parse_command_line (is_getting, case_file)
        if (is_getting )  then
            file_main = case_file
        else
            file_main = TRIM(PREFIX_INPUT) // TRIM(EXTEND_CASENAME)
        end if
        
        ! prepare files
        call Open_IO_files (file_main)
        
        ! read input file and allocate memory for variables
        call time_event%begin(1)
        call Driving_input (file_main)
        call time_event%end(1)
        
        ! set openMP environment
        call Set_openMP (ns%misc%nthread)
        
    end subroutine Run_pre_process
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Open_IO_files (file_main)
        
        character(len=*), intent(in)  :: file_main
        
        ! local variables
        character(len=MAX_WORD_LEN)  :: file_in
        character(len=MAX_WORD_LEN)  :: file_out
        character(len=MAX_WORD_LEN)  :: extension
        integer  :: file_unit
        
        logical  :: is_input
        logical  :: is_exist                                                    ! logical for file existence
        integer  :: io_error                                                    ! for file I/O error
        
        ! input file open
        call Is_file_existence (file_main, is_input=.TRUE.)
        call Get_file_extension (file_main, extension)
        
        if (TRIM(ADJUSTL(extension)) /= EXTEND_XML)  then
            call Preprocess_input (file_main, PREFIX_INPUT, PREFIX_TMP, file_out)
            call Open_file (file_out, is_input=.TRUE., file_unit=file_unit)
            FILES%CASENAME = file_unit
        end if
        
        ! output file open
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_MAIN)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%MAIN = file_unit
        
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_TIMELIST)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%TIMELIST = file_unit
        
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_DET)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%DET = file_unit
        
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_MEMORY)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%MEMORY = file_unit
        
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_REACTIVITY)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%REACTIVITY = file_unit
        
!        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_POINTKINETICS)
!        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
!        FILES%POINTKINETICS = file_unit
    
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(EXTEND_PT)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%PT = file_unit
    
    end subroutine Open_IO_files
    
end module driver_pre_process