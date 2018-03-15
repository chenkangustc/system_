!$
!===================================================================================================
!
!   subroutine related to file operation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Preprocess_input
!                               Preprocess_ansys
!                               Get_line_count
!                               Get_file_extension
!                               Get_file_basename
!                               Get_file_dir
!                               Remove_dir_in_path
!                               Open_file
!                               Close_file
!                               Is_file_existence
!                               Create_blank_file
!                               Delete_file
!                               Skip_lines
!                               Get_lines
!
!   Public type lists:          No
!
!===================================================================================================
module files_dirs

    use constants
    use, intrinsic :: ISO_FORTRAN_ENV
    
    use string
    use exception_header,       only : WarningCollector, ErrorCollector
    
    implicit none    
    private
    public  :: Preprocess_input, Preprocess_ansys, Get_line_count
    public  :: Get_file_extension, Get_file_basename, Get_file_dir, Remove_dir_in_path
    public  :: Open_file, Close_file, Is_file_existence, Create_blank_file, Delete_file
    public  :: Skip_lines, Get_lines
    
    interface  Get_lines
        module procedure  Get_one_lines
        module procedure  Get_multi_lines
    end interface Get_lines
    
    type(WarningCollector)  :: a_warning
    type(ErrorCollector)    :: a_error
    
contains
    !$
    !===============================================================================================
    ! preprocess the main input file of the problem
    !   file_out = file_in - char_sub + char_add (sub & add is prefix)
    !   file_out & file_in should not be the same
    !===============================================================================================
    subroutine Preprocess_input (file_in, char_sub, char_add, file_out)
        
        character(len=*), intent(in)  :: file_in
        character(len=*), intent(in)  :: char_sub
        character(len=*), intent(in)  :: char_add
        character(len=MAX_WORD_LEN), intent(out)  :: file_out
        
        ! local variables
        character(len=MAX_LINE_LEN)   :: a_line                                 ! a line of input file
        character(len=MAX_LINE_LEN)   :: line_write                             ! a line to be write after preprocess
        character(len=MAX_LINE_LEN)   :: line_tmp
        
        character(len=MAX_WORD_LEN)  :: words(MAX_WORDS)                        ! hold for words after a line splited
        integer  :: n_word                                                      ! how many words of a line
        integer  :: i_position                                                  ! position of comment character
        integer  :: i
        integer  :: char_len
        integer  :: io_error
        
        integer  :: TMP_IN                                                      ! input unit for temparary usage
        integer  :: TMP_OUT
        
        character(len=MAX_WORD_LEN)  :: char_mark = ''
        character(len=MAX_WORD_LEN)  :: char_keep_start = 'MCNP_GEOM'
        character(len=MAX_WORD_LEN)  :: char_keep_end = 'MCNP_GEOM_END'
        logical                      :: is_keep_start = .FALSE.
        logical                      :: is_keep_end = .FALSE.
        
        if (TRIM(ADJUSTL(char_sub)) == TRIM(ADJUSTL(char_add)))  then
            call a_error%set (INFO_LIST_FRAMEWORK, 'the input and output file is the same')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        file_out = TRIM(char_add) // TRIM(file_in(LEN_TRIM(char_sub)+1: ))
        open (newunit=TMP_IN, file=file_in, access='sequential', form='formatted', status='old', action='read', iostat=io_error)
        open (newunit=TMP_OUT, file=file_out, access='sequential', form='formatted', status='replace', action='write', iostat=io_error)
        
        do
            read (unit=TMP_IN, fmt="(A)", iostat=io_error)  a_line
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read (unit=a_line, fmt="(A)") char_mark
            call Upper_case(char_mark)
            
            if (TRIM(ADJUSTL(char_mark)) == TRIM(ADJUSTL(char_keep_start)))  then
                is_keep_start = .TRUE.
            end if 
            if (TRIM(ADJUSTL(char_mark)) == TRIM(ADJUSTL(char_keep_end)))  then
                is_keep_end = .TRUE.
            end if 
            
            if (is_keep_start .and. (.NOT. is_keep_end))  then
                if (TRIM(a_line) == '')  then
                    write(unit=TMP_OUT, fmt="(A)")  '        '
                else
                    do i = MAX_LINE_LEN, 1, -1
                        if (a_line(i:i) /= ' ') exit
                    end do 
                    write(unit=TMP_OUT, fmt="(A)")  a_line(:i)
                end if
                cycle
            end if 
            
            ! replace the "TAB" in line content with "SPACE"
            call Replace_tab (a_line)
            ! transfer to upper case
            call Upper_case (a_line)
            ! skip blank line
            if (TRIM(a_line) == '')  cycle
            
            ! treat the comment
            i_position = INDEX(ADJUSTL(a_line), '!', back=.FALSE.)
            if (i_position == 1)  cycle
            if (i_position /= 0)  then
                char_len = i_position - 1
            else
                char_len = LEN(TRIM(a_line))
            end if 
            
            line_tmp = a_line(1: char_len)
            
            ! change logical symbol: F--> .false.  T--> .true.
            call Split_string(line_tmp, words, n_word)
            do i = 1, n_word
                if (TRIM(words(i)) == 'F')  then
                    words(i) = '.FALSE.'
                end if
                if (TRIM(words(i)) == 'T')  then
                    words(i) = '.TRUE.'
                end if
            end do
            call Concatenate (line_write, words, n_word)
            
            ! rewrite to file
            write(unit=TMP_OUT, fmt="(A)")  TRIM(line_write)
        end do 
        
        ! post treatment after read the contents
        close(unit=TMP_IN, status='keep', iostat=io_error)
        close(unit=TMP_OUT, status='keep', iostat=io_error)
        
    end subroutine Preprocess_input
    
    !$
    !===============================================================================================
    ! preprocess the ANSYS type input files, remove the ANSYS type note
    !===============================================================================================
    subroutine Preprocess_ansys (file_in, file_out, nline)
        
        character(len=*), intent(in)  :: file_in
        character(len=*), intent(out) :: file_out
        integer, intent(out)          :: nline
        
        ! local variables
        character(len=MAX_LINE_LEN)   :: a_line                                 ! a line of input file
        character(len=MAX_LINE_LEN)   :: line_tmp
        
        character(len=MAX_WORD_LEN)  :: words(MAX_WORDS)                        ! hold for words after a line splited
        integer  :: n_word                                                      ! how many words of a line
        integer  :: i_position                                                  ! position of comment character
        integer  :: i
        integer  :: char_len
        integer  :: io_error
        
        integer  :: TMP_IN                                                      ! input unit for temparary usage
        integer  :: TMP_OUT
        
        file_out = TRIM(file_in) // '_'
        open (newunit=TMP_IN, file=file_in, access='sequential', form='formatted', status='old', action='read', iostat=io_error)
        open (newunit=TMP_OUT, file=file_out, access='sequential', form='formatted', status='replace', action='write', iostat=io_error)
        
        nline = 0 
        do
            read (unit=TMP_IN, fmt="(A)", iostat=io_error)  a_line
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            ! replace the TAB in line content with SPACE
            call Replace_tab (a_line)
            ! transfer to upper case
            call Upper_case (a_line)
            ! skip blank line
            if (TRIM(a_line) == '')  cycle
            
            ! skip useless line of ANSYS file
            if (INDEX(a_line, 'NODE')    /= 0)  cycle
            if (INDEX(a_line, 'ELEM')    /= 0)  cycle
            if (INDEX(a_line, 'LIST')    /= 0)  cycle
            if (INDEX(a_line, 'ANSYS')   /= 0)  cycle
            if (INDEX(a_line, 'VERSION') /= 0)  cycle
            
            nline = nline + 1 
            char_len = LEN(TRIM(a_line))
            line_tmp = a_line(1: char_len)
            
            ! write back to file
            write(unit=TMP_OUT, fmt="(A)") TRIM(line_tmp)
        end do 
        
        ! post treatment after read the contents
        close(unit=TMP_IN, status='keep', iostat=io_error)
        close(unit=TMP_OUT, status='keep', iostat=io_error)
        
    end subroutine Preprocess_ansys 
    
    !$
    !===============================================================================================
    ! count non-blank lines
    !===============================================================================================
    subroutine Get_line_count (file_in, cnt)
    
        character(len=*), intent(in)  :: file_in
        integer, intent(in out)       :: cnt
        
        character(len=MAX_WORD_LEN)  :: a_line
        integer  :: unit_in
        integer  :: io_error
        
        open(newunit=unit_in, file=file_in, action='read', iostat=io_error)
        
        cnt = 0
        do 
            read (unit=unit_in, fmt="(A)", iostat=io_error)  a_line
            if (io_error == IOSTAT_END)  exit
            
            call Replace_tab (a_line)
            if (TRIM(a_line) == '')  cycle
            cnt = cnt + 1
        end do 
        
        close(unit=unit_in, iostat=io_error)
        
    end subroutine Get_line_count
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! with '.'
    !===============================================================================================
    subroutine Get_file_extension (file_in, extension)
        
        character(len=*), intent(in)      :: file_in
        character(len=*), intent(in out)  :: extension
        
        integer  :: ibeg, iend, i, j
        
        extension = ''
        iend = LEN_TRIM(file_in)
        ibeg = INDEX(file_in, '.', back=.TRUE.)
        
        if (ibeg > 0)  then
            do i = ibeg, iend
                j = i + 1 - ibeg
                extension(j:j) = file_in(i:i)
            end do
        end if
    
    end subroutine Get_file_extension
    
    !$
    !===============================================================================================
    ! without '.'
    !===============================================================================================
    subroutine Get_file_basename (file_in, basename)
        
        character(len=*), intent(in)      :: file_in
        character(len=*), intent(in out)  :: basename
        
        integer  :: ibeg, iend, i, j
        
        basename = ''
        iend = INDEX(file_in, '.', back=.TRUE.) - 1
        ibeg = MAX(INDEX(file_in, '\', back=.TRUE.), INDEX(file_in, '/', back=.TRUE.)) + 1
        
        if (iend > -1)  then
            do i = ibeg, iend
                j = i + 1 - ibeg
                basename(j:j) = file_in(i:i)
            end do
        end if
    
    end subroutine Get_file_basename
    
    !$
    !===============================================================================================
    ! with '\' or '/', the directory index is '/' in both windows and linux
    !===============================================================================================
    subroutine Get_file_dir (file_in, dir)
    
        character(len=*), intent(in)      :: file_in
        character(len=*), intent(in out)  :: dir
        
        integer  :: ibeg, iend, i, j
        
        dir = ''
        ibeg = 1
        iend = MAX(INDEX(file_in, '\', back=.TRUE.), INDEX(file_in, '/', back=.TRUE.))
        
        if (iend > 0)  then
            do i = ibeg, iend
                j = i + 1 - ibeg
                dir(j:j) = file_in(i:i)
            end do
        end if
    
    end subroutine Get_file_dir
    
    !$
    !===============================================================================================
    ! get file name without directory
    !===============================================================================================
    subroutine Remove_dir_in_path (file_in, file_out)
        
        character(len=MAX_WORD_LEN), intent(in)      :: file_in
        character(len=MAX_WORD_LEN), intent(in out)  :: file_out
        
        integer  :: i, start
        file_out = ''
        start = MAX(index(file_in, '\', back=.TRUE.), index(file_in, '/', back=.TRUE.))
        
        if (start > 0)  then
            do i = start + 1, len(file_in)
                file_out(i-start: i-start) = file_in(i: i)
            end do
        end if        
        
        file_out = TRIM(ADJUSTL(file_out))
    
    end subroutine Remove_dir_in_path
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! close file
    !===============================================================================================
    subroutine Close_file (file_unit, is_keep)
    
        integer, intent(in)  :: file_unit
        logical, intent(in)  :: is_keep
        
        logical  :: is_open
        integer  :: io_error
        
        inquire(unit=file_unit, opened=is_open)
        if (is_open)  then
            if (is_keep)  then
                close(unit=file_unit, status='keep', iostat=io_error)
            else
                close(unit=file_unit, status='delete', iostat=io_error)
            end if
        end if
    
    end subroutine Close_file
    
    !$
    !===============================================================================================
    ! open file
    !===============================================================================================
    subroutine Open_file (file_name, is_input, file_unit)
        
        character(len=*), intent(in)  :: file_name
        logical, intent(in)           :: is_input
        integer, intent(out)          :: file_unit
        
        integer  :: io_error
        
        ! check existence for input file only
        if (is_input)  then
            call Is_file_existence(file_name, is_input)
        end if
        
        if (is_input)  then
            open(newunit=file_unit, file=file_name, status='old', action='read', iostat=io_error)
        else
            open(newunit=file_unit, file=file_name, status='replace', action='write', iostat=io_error)
        end if
        
        if (io_error /= 0)  then
            call a_error%set (INFO_LIST_FILE, 'open file failed: ' // file_name)  
            call a_error%print (OUTPUT_UNIT)
        end if
    
    end subroutine Open_file
    
    !$
    !===============================================================================================
    ! check a file is existence
    !===============================================================================================
    subroutine Is_file_existence(file_name, is_input)
        
        character(len=*), intent(in) :: file_name
        logical, intent(in)  :: is_input
        logical  :: is_true
        
        inquire(file=file_name, exist=is_true)
        
        if (is_input)  then
            if (.NOT. is_true)  then
                call a_error%set (INFO_LIST_FILE, 'input file does not exist: ' // file_name)  
                call a_error%print (OUTPUT_UNIT)
            end if
        else
            if (is_true)  then
                call a_warning%set (INFO_LIST_FILE, 'output file name has been used: ' // file_name) 
                call a_warning%print (OUTPUT_UNIT)
            end if
        end if
    
    end subroutine Is_file_existence
    
    !$
    !===============================================================================================
    ! create a blank file
    !===============================================================================================
    subroutine Create_blank_file (file_name)
        
        character(len=*), intent(in)  :: file_name
        
        integer  :: file_unit, io_error
        logical  :: is_true
        
        inquire(file=file_name, exist=is_true)
        if (is_true)  then
            call Delete_file (file_name)
        end if
        
        open(newunit=file_unit, file=file_name, status='replace', action='write', iostat=io_error)
        close(unit=file_unit, status='keep', iostat=io_error)
    
    end subroutine Create_blank_file
    
    !$
    !===============================================================================================
    ! delete an exist file
    !===============================================================================================
    subroutine Delete_file (file_name)
        
        character(len=*), intent(in)  :: file_name
        
        integer  :: file_unit, io_error
        logical  :: is_true
        
        inquire(file=file_name, exist=is_true)
        if (is_true)  then
            open(newunit=file_unit, file=file_name, status='old', action='read', iostat=io_error)
            close(unit=file_unit, status='delete', iostat=io_error)
        end if
    
    end subroutine Delete_file
    
    !$
    !===============================================================================================
    ! skip lines of an opened file
    !===============================================================================================
    subroutine Skip_lines (unit_in, n)
        
        integer, intent(in)  :: unit_in
        integer, intent(in)  :: n
        
        character(len=MAX_WORD_LEN)  :: tmp_line
        integer  :: i, io_error
        
        do i = 1, n
            read(unit=unit_in, fmt="(A)", iostat=io_error)  tmp_line
        end do
    
    end subroutine Skip_lines
    
    !$
    !===============================================================================================
    ! get content of one line, with specified column
    !===============================================================================================
    function Get_one_lines (unit_in, y)  result(values)
        
        integer, intent(in)  :: unit_in
        integer, intent(in)  :: y
        real(KREAL), allocatable  :: values(:)
        
        integer  :: i, io_error
        
        allocate(values(y), stat=io_error)
        read(unit=unit_in, fmt=*, iostat=io_error)  values
    
    end function Get_one_lines
    
    !$
    !===============================================================================================
    ! get content of multi lines, with specified column 'y'
    !===============================================================================================
    function Get_multi_lines (unit_in, y, line)  result(values)
        
        integer, intent(in)  :: unit_in
        integer, intent(in)  :: y
        integer, intent(in)  :: line
        real(KREAL), allocatable  :: values(:, :)
        
        integer  :: i, io_error
        allocate(values(line, y), stat=io_error)
        
        do i = 1, line
            read(unit=unit_in, fmt=*, iostat=io_error)  values(i, :)
        end do
    
    end function Get_multi_lines
    
end module files_dirs
    