!===================================================================================================
!
!   module for general string operation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Upper_case
!                               Lower_case
!                               Replace_tab
!                               Split_string
!                               Concatenate
!                               Join_Words
!                               Is_startwith
!                               Is_endwith
!                               Is_number
!                               Int_to_string
!                               String_to_int
!
!   Public type lists:          No
!
!===================================================================================================
module string

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    implicit none
    private
    public  :: Upper_case, Lower_case, Replace_tab, Split_string, Concatenate, Join_Words
    public  :: Is_startwith, Is_endwith, Is_number, Int_to_string, String_to_int
    
    interface  Concatenate
        module procedure  Concatenate_array
        module procedure  Concatenate_words
    end interface Concatenate
    
contains
    !$
    !===============================================================================================
    ! convert string to 'UPPER' case
    !===============================================================================================
    elemental subroutine Upper_case (string)
        
        character(len=*), intent(in out)  :: string
        
        integer, parameter  :: CHAR_DISTANCE = 32
        integer  :: i
        
        do i = 1, LEN(string)
            if (LGE(string(i:i), 'a') .and. (LLE(string(i:i), 'z')))  then 
                string(i:i) = ACHAR(IACHAR(string(i:i)) - CHAR_DISTANCE)
            end if
        end do
        
    end subroutine Upper_case
    
    !$
    !===============================================================================================
    ! convert string to 'lower' case
    !===============================================================================================
    elemental subroutine Lower_case (string)
    
        character(len=*), intent(in out)  :: string
        
        integer, parameter  :: CHAR_DISTANCE = 32
        integer  :: i
        
        do i = 1, LEN(string)
            if (LGE(string(i:i), 'A') .and. (LLE(string(i:i), 'Z'))) then
                string(i:i) = ACHAR(IACHAR(string(i:i)) + CHAR_DISTANCE)
            end if
        end do
    
    end subroutine Lower_case

    !$
    !===============================================================================================
    ! replace tab by one space: IACHR=9 --> IACHR=32, and adjust the left blank
    ! this should be the first routine be invoked for a specific string
    !===============================================================================================
    subroutine Replace_tab (string)
        
        character(len=*), intent(in out)  :: string
        
        integer, parameter  :: CHAR_TAB   = 9
        integer, parameter  :: CHAR_SPACE = 32
        integer  :: i
        
        do i = 1, LEN(string)
            if (IACHAR(string(i:i)) == CHAR_TAB)  then
                string(i:i) = ACHAR(CHAR_SPACE)
            end if
        end do
        
        string = ADJUSTL(string)
        
    end subroutine Replace_tab
    
    !$
    !===============================================================================================
    ! split a string into seperate words
    !===============================================================================================
    subroutine Split_string (string, words, n, sub)
        
        character(len=*), intent(in)     :: string                              ! a string to be splited
        character(len=*), intent(in out) :: words(MAX_WORDS)                    ! hold for words
        integer, intent(out)             :: n                                   ! actual number of words
        character(len=1), intent(in), optional  :: sub                          ! the seperator
        
        character(len=LEN(string))  :: new_string                               ! hold for string
        character(len=1)  :: str                                                ! the actual seperator
        integer           :: ibeg                                             ! start position of words
        integer           :: iend                                               ! end position of words
        integer           :: i  
        
        ! get the separator
        if (PRESENT(sub))  then
            str = sub
        else 
            str = ' '
        end if
        
        new_string = ADJUSTL(string)
        ibeg = 0; iend = 0
        n = 0
        do i = 1, LEN_TRIM(new_string)
            ! get the start position
            if ((ibeg == 0) .and. (new_string(i:i) /= str))  then
                ibeg = i
            end if
                
            if (ibeg > 0) then
                ! get the ending position
                if (new_string(i:i) == str)  then 
                    iend = i - 1
                end if
                ! if the last one and when the seperator is space
                if (i == LEN_TRIM(new_string) .and. str == ' ')  then 
                    iend = i
                end if
                
                ! count for words
                if (iend > 0) then
                    n = n + 1
                    ! word length longer than MAX_WORDS
                    if (iend-ibeg+1 > len(words(n)))  then
                    end if
                    
                    ! get the word and reset the position for next one
                    words(n) = ADJUSTL(new_string(ibeg : iend))
                    ibeg = 0; iend = 0
                end if
            end if
        end do
        
    end subroutine Split_string
    
    !$
    !===============================================================================================
    ! concatenate an array of word with a space
    !===============================================================================================
    subroutine Concatenate_array (string, words, n_word, sub)
        
        character(len=MAX_LINE_LEN), intent(out) :: string
        character(len=*), intent(in)             :: words(:)
        integer, intent(in)                      :: n_word
        character(len=*), intent(in), optional   :: sub
        
        character(len=1)  :: str
        integer  :: cnt 
        integer  :: i
             
        if (PRESENT(sub))  then
            str = sub
        else 
            str = ' '
        end if
             
        string = TRIM(ADJUSTL(words(1)))
        cnt = MIN(n_word, LEN(words))
        if (cnt == 1)  then
            return
        end if
        
        do i = 2, cnt
            string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(words(i)))
        end do
    
    end subroutine Concatenate_array
    
    !$
    !===============================================================================================
    ! concatenate many words with a space, use 'sub' keyword
    !===============================================================================================
    subroutine Concatenate_words (string, a, b, c, d, e, f, g, h, sub)
        
        character(len=MAX_LINE_LEN), intent(in out) :: string
        character(len=*), intent(in)                :: a
        character(len=*), intent(in), optional      :: b
        character(len=*), intent(in), optional      :: c
        character(len=*), intent(in), optional      :: d
        character(len=*), intent(in), optional      :: e
        character(len=*), intent(in), optional      :: f
        character(len=*), intent(in), optional      :: g
        character(len=*), intent(in), optional      :: h
        character(len=*), intent(in), optional      :: sub
        
        character(len=1)  :: str
        
        if (PRESENT(sub))  then
            str = sub
        else 
            str = ' '
        end if
        
        string = TRIM(ADJUSTL(a))
        if (PRESENT(b))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(b))
        if (PRESENT(c))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(c))
        if (PRESENT(d))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(d))
        if (PRESENT(e))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(e))
        if (PRESENT(f))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(f))
        if (PRESENT(g))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(g))
        if (PRESENT(h))   string = TRIM(ADJUSTL(string)) // str // TRIM(ADJUSTL(h))
    
    end subroutine Concatenate_words

    !$
    !===============================================================================================
    ! join many words directely
    !===============================================================================================
    subroutine Join_Words (string, a, b, c, d, e, f, g, h)
        
        character(len=MAX_LINE_LEN), intent(in out) :: string
        character(len=*), intent(in)                :: a
        character(len=*), intent(in), optional      :: b
        character(len=*), intent(in), optional      :: c
        character(len=*), intent(in), optional      :: d
        character(len=*), intent(in), optional      :: e
        character(len=*), intent(in), optional      :: f
        character(len=*), intent(in), optional      :: g
        character(len=*), intent(in), optional      :: h
        
        string = TRIM(ADJUSTL(a))
        if (PRESENT(b))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(b))
        if (PRESENT(c))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(c))
        if (PRESENT(d))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(d))
        if (PRESENT(e))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(e))
        if (PRESENT(f))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(f))
        if (PRESENT(g))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(g))
        if (PRESENT(h))   string = TRIM(ADJUSTL(string)) // TRIM(ADJUSTL(h))
    
    end subroutine Join_Words

    !$
    !===============================================================================================
    ! string start with a sub string without the former bank
    !===============================================================================================
    function Is_startwith (string, sub)  result(is_true)
        
        character(len=*), intent(in)  :: string
        character(len=*), intent(in)  :: sub
        logical  :: is_true

        character(len=LEN_TRIM(ADJUSTL(string)))  :: tmp_string
        character(len=LEN_TRIM(ADJUSTL(sub)))     :: tmp_sub
        
        tmp_string = TRIM(ADJUSTL(string))
        tmp_sub = TRIM(ADJUSTL(sub))
        
        if (INDEX(tmp_string, tmp_sub) == 1) then
            is_true = .TRUE.
        else
            is_true = .FALSE.
        end if
    
    end function Is_startwith
    
    !$
    !===============================================================================================
    ! string end with a sub string exclude the ending bank
    !===============================================================================================
    function Is_endwith (string, sub)  result(is_true)
        
        character(len=*), intent(in)  :: string
        character(len=*), intent(in)  :: sub
        logical  :: is_true

        character(len=LEN_TRIM(ADJUSTL(string)))  :: tmp_string
        character(len=LEN_TRIM(ADJUSTL(sub)))     :: tmp_sub
        integer  :: ibeg
        integer  :: len_string
        integer  :: len_sub
        
        tmp_string = TRIM(ADJUSTL(string))
        tmp_sub = TRIM(ADJUSTL(sub))
        
        len_string = LEN(tmp_string)
        len_sub = LEN(tmp_sub)
        ibeg = len_string - len_sub + 1
        
        if (INDEX(tmp_string, tmp_sub, back=.TRUE.) == ibeg) then
            is_true = .TRUE.
        else 
            is_true = .FALSE.
        end if 
    
    end function Is_endwith
    
    !$
    !===============================================================================================
    ! whether the character of a string in all in 0-9 exclude the blank space
    !===============================================================================================
    function Is_number (string)  result(is_true)
    
        character(len=*), intent(in)  :: string
        logical  :: is_true
        
        character(len=LEN(string))  :: word
        integer  :: i_char
        integer  :: i
        
        word = ADJUSTL(string)
        is_true = .TRUE.
        do i = 1, LEN_TRIM(word)
            i_char = IACHAR(word(i:i))
            if (i_char < IACHAR("0") .or. i_char > IACHAR("9")) then 
                is_true = .FALSE.
                exit
            end if
        end do
      
    end function Is_number

    !$
    !===============================================================================================
    ! convert positive integer to a char string '0000000xx'
    ! 'digit' valid from 1 to 9
    !===============================================================================================
    function Int_to_string (number, digit)  result(str)
        
        integer, intent(in)            :: number
        integer, intent(in), optional  :: digit
        
        character(len=MAX_WORD_LEN)    :: str, line, tmp
        integer, parameter :: PREDEFINE = 9
        
        write(unit=tmp, fmt="(1x, I9)")  number
        if (number <= 9)  then
            line = '00000000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 99)  then
            line = '0000000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 999)  then
            line = '000000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 9999)  then
            line = '00000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 99999)  then
            line = '0000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 999999)  then
            line = '000' // TRIM(ADJUSTL(tmp))        
        else if (number <= 9999999)  then
            line = '00' // TRIM(ADJUSTL(tmp))        
        else if (number <= 99999999)  then
            line = '0' // TRIM(ADJUSTL(tmp))        
        else if (number <= 999999999)  then
            line = '' // TRIM(ADJUSTL(tmp))        
        end if
        
        if (PRESENT(digit) .and. digit <= PREDEFINE)  then
            str = line(PREDEFINE+1-digit: )
        end if
    
    end function Int_to_string

    !$
    !===============================================================================================
    ! convert string to integer
    !===============================================================================================
    function String_to_int (str)  result(number)
        
        character(len=MAX_WORD_LEN), intent(in)  :: str
        integer  :: number
        
        if (Is_number(str) )  then
            read(unit=str, fmt=*)  number
        else
            number = -HUGE(0)
        end if
    
    end function String_to_int
    
end module string 
