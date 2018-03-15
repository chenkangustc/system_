!$
!===================================================================================================
!
!   performs a binary search of an array to find where a specific value lies in the array. 
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    binary_search
!
!   Public type lists:          No
!
!===================================================================================================
module search

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector
    
    implicit none
    private 
    public  :: binary_search
    
    ! --------------------------------------------------------------------------
    interface binary_search
        module procedure  binary_search_real4
        module procedure  binary_search_real8
        module procedure  binary_search_int4
        module procedure  binary_search_int8
    end interface binary_search
    
    type(ErrorCollector)  :: a_error
    integer, parameter    :: MAX_ITERATION = 64

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function binary_search_real4(array, n, value) result(array_index)

        integer, intent(in) :: n
        real(4), intent(in) :: array(n)
        real(4), intent(in) :: value
        integer             :: array_index
        
        integer :: L
        integer :: R
        integer :: n_iteration
        real(4) :: testval
        
        L = 1
        R = n
        if (value < array(L) .or. value > array(R)) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'Value outside of array during binary search')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        n_iteration = 0
        do while (R - L > 1)
            ! Check boundaries
            if (value > array(L) .and. value < array(L+1)) then
                array_index = L
                return
            elseif (value > array(R-1) .and. value < array(R)) then
                array_index = R - 1
                return
            end if
            
            ! Find values at midpoint
            array_index = L + (R - L)/2
            testval = array(array_index)
            if (value >= testval) then
                L = array_index
            elseif (value < testval) then
                R = array_index
            end if
            
            ! check for large number of iterations
            n_iteration = n_iteration + 1
            if (n_iteration == MAX_ITERATION) then
                call a_error%set (INFO_LIST_FRAMEWORK, 'Reached maximum number of iterations on binary search')
                call a_error%print (OUTPUT_UNIT)
                exit
            end if
        end do
        
        array_index = L

    end function binary_search_real4
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function binary_search_real8(array, n, value) result(array_index)

        integer, intent(in) :: n
        real(8), intent(in) :: array(n)
        real(8), intent(in) :: value
        integer             :: array_index
        
        integer :: L
        integer :: R
        integer :: n_iteration
        real(8) :: testval
        
        L = 1
        R = n
        if (value < array(L) .or. value > array(R)) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'Value outside of array during binary search')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        n_iteration = 0
        do while (R - L > 1)
            ! Check boundaries
            if (value > array(L) .and. value < array(L+1)) then
                array_index = L
                return
            elseif (value > array(R-1) .and. value < array(R)) then
                array_index = R - 1
                return
            end if
            
            ! Find values at midpoint
            array_index = L + (R - L)/2
            testval = array(array_index)
            if (value >= testval) then
                L = array_index
            elseif (value < testval) then
                R = array_index
            end if
            
            ! check for large number of iterations
            n_iteration = n_iteration + 1
            if (n_iteration == MAX_ITERATION) then
                call a_error%set (INFO_LIST_FRAMEWORK, 'Reached maximum number of iterations on binary search')
                call a_error%print (OUTPUT_UNIT)
                exit
            end if
        end do
        
        array_index = L

    end function binary_search_real8
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function binary_search_int4(array, n, value) result(array_index)

        integer, intent(in) :: n
        integer, intent(in) :: array(n)
        integer, intent(in) :: value
        integer             :: array_index
        
        integer :: L
        integer :: R
        integer :: n_iteration
        real(8) :: testval
        
        L = 1
        R = n
        if (value < array(L) .or. value > array(R)) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'Value outside of array during binary search')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        n_iteration = 0
        do while (R - L > 1)
            ! Check boundaries
            if (value > array(L) .and. value < array(L+1)) then
                array_index = L
                return
            elseif (value > array(R-1) .and. value < array(R)) then
                array_index = R - 1
                return
            end if
            
            ! Find values at midpoint
            array_index = L + (R - L)/2
            testval = array(array_index)
            if (value >= testval) then
                L = array_index
            elseif (value < testval) then
                R = array_index
            end if
            
            ! check for large number of iterations
            n_iteration = n_iteration + 1
            if (n_iteration == MAX_ITERATION) then
                call a_error%set (INFO_LIST_FRAMEWORK, 'Reached maximum number of iterations on binary search')
                call a_error%print (OUTPUT_UNIT)
                exit 
            end if
        end do
        
        array_index = L

    end function binary_search_int4
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function binary_search_int8(array, n, value) result(array_index)

        integer,    intent(in) :: n
        integer(8), intent(in) :: array(n)
        integer(8), intent(in) :: value
        integer                :: array_index
        
        integer :: L
        integer :: R
        integer :: n_iteration
        real(8) :: testval
        
        L = 1
        R = n
        if (value < array(L) .or. value > array(R)) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'Value outside of array during binary search')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        n_iteration = 0
        do while (R - L > 1)
            ! Check boundaries
            if (value > array(L) .and. value < array(L+1)) then
                array_index = L
                return
            elseif (value > array(R-1) .and. value < array(R)) then
                array_index = R - 1
                return
            end if
            
            ! Find values at midpoint
            array_index = L + (R - L)/2
            testval = array(array_index)
            if (value >= testval) then
                L = array_index
            elseif (value < testval) then
                R = array_index
            end if
            
            ! check for large number of iterations
            n_iteration = n_iteration + 1
            if (n_iteration == MAX_ITERATION) then
                call a_error%set (INFO_LIST_FRAMEWORK, 'Reached maximum number of iterations on binary search')
                call a_error%print (OUTPUT_UNIT)
                exit 
            end if
        end do
        
        array_index = L
        
    end function binary_search_int8

end module search
