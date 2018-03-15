!$
!===================================================================================================
!
!   a linked list class with methods such as append, contains, remove, index, get_item, size, etc.
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          ListInt
!                               ListReal
!                               ListChar
!
!===================================================================================================
module list_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    implicit none
    private
    public  :: ListInt, ListReal, ListChar

    ! --------------------------------------------------------------------------
    ! LISTELEM* types hold one piece of data and a pointer to the next piece of data
    type, private :: ListElemInt
        integer                    :: data
        type(ListElemInt), pointer :: next => null()
        type(ListElemInt), pointer :: prev => null()
    end type ListElemInt

    type, private :: ListElemReal
        real(KREAL)             :: data
        type(ListElemReal), pointer :: next => null()
        type(ListElemReal), pointer :: prev => null()
    end type ListElemReal

    type, private :: ListElemChar
        character(len=MAX_WORD_LEN)     :: data
        type(ListElemChar), pointer :: next => null()
        type(ListElemChar), pointer :: prev => null()
    end type ListElemChar

    ! --------------------------------------------------------------------------
    ! LIST* types contain the linked list with convenience methods. 
    type :: ListInt
        private
        integer                     :: count = 0                                ! Number of elements in list

        ! Used in get_item for fast sequential lookups 
        integer                     :: last_index = HUGE(0)
        type(ListElemInt), pointer  :: last_elem => null()

        ! Pointers to beginning and end of list
        type(ListElemInt), public, pointer :: head => null()
        type(ListElemInt), public, pointer :: tail => null()
    contains
        procedure, public :: append => List_append_int                          ! Add item to end of list
        procedure, public :: clear => List_clear_int                            ! Remove all items
        procedure, public :: contains => List_contains_int                      ! Does list contain?
        procedure, public :: get_item => List_get_item_int                      ! Get i-th item in list
        procedure, public :: index => List_index_int                            ! Determine index of given item
        procedure, public :: insert => List_insert_int                          ! Insert item in i-th position
        procedure, public :: remove => List_remove_int                          ! Remove specified item
        procedure, public :: size => List_size_int                              ! Size of list
    end type ListInt

    type :: ListReal
        private
        integer                     :: count = 0                                ! Number of elements in list

        ! Used in get_item for fast sequential lookups 
        integer                     :: last_index = HUGE(0)
        type(ListElemReal), pointer :: last_elem => null()

        ! Pointers to beginning and end of list
        type(ListElemReal), public, pointer :: head => null()
        type(ListElemReal), public, pointer :: tail => null()
    contains
        procedure, public :: append => List_append_real                         ! Add item to end of list
        procedure, public :: clear => List_clear_real                           ! Remove all items
        procedure, public :: contains => List_contains_real                     ! Does list contain?
        procedure, public :: get_item => List_get_item_real                     ! Get i-th item in list
        procedure, public :: index => List_index_real                           ! Determine index of given item
        procedure, public :: insert => List_insert_real                         ! Insert item in i-th position
        procedure, public :: remove => List_remove_real                         ! Remove specified item
        procedure, public :: size => List_size_real                             ! Size of list
    end type ListReal

    type :: ListChar
        private
        integer                     :: count = 0                                ! Number of elements in list

        ! Used in get_item for fast sequential lookups 
        integer                     :: last_index = HUGE(0)
        type(ListElemChar), pointer :: last_elem => null()

        ! Pointers to beginning and end of list
        type(ListElemChar), public, pointer :: head => null()
        type(ListElemChar), public, pointer :: tail => null()
    contains
        procedure, public :: append => List_append_char                         ! Add item to end of list
        procedure, public :: clear => List_clear_char                           ! Remove all items
        procedure, public :: contains => List_contains_char                     ! Does list contain?
        procedure, public :: get_item => List_get_item_char                     ! Get i-th item in list
        procedure, public :: index => List_index_char                           ! Determine index of given item
        procedure, public :: insert => List_insert_char                         ! Insert item in i-th position
        procedure, public :: remove => List_remove_char                         ! Remove specified item
        procedure, public :: size => List_size_char                             ! Size of list
    end type ListChar

contains
    !$
    !===============================================================================================
    ! LIST_APPEND appends an item to the end of the list. If the list is empty, it becomes the first item
    !===============================================================================================
    subroutine List_append_int(this, data)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: data
        
        type(ListElemInt), pointer :: elem
        
        ! Create element and set dat
        allocate(elem)
        elem % data = data
        
        if (.NOT. ASSOCIATED(this % head)) then
            ! If list is empty, set head and tail to new element
            this % head => elem
            this % tail => elem
        else
            ! Otherwise append element at end of list
            this % tail % next => elem
            elem % prev => this % tail
            this % tail => this % tail % next
        end if
        
        this % count = this % count + 1
    end subroutine List_append_int
    
    !$
    subroutine List_append_real(this, data)
        class(ListReal), intent(in out)  :: this
        real(KREAL), intent(in)      :: data
        
        type(ListElemReal), pointer :: elem
        
        ! Create element and set dat
        allocate(elem)
        elem % data = data
        
        if (.NOT. ASSOCIATED(this % head)) then
            ! If list is empty, set head and tail to new element
            this % head => elem
            this % tail => elem
        else
            ! Otherwise append element at end of list
            this % tail % next => elem
            elem % prev => this % tail
            this % tail => this % tail % next
        end if
        
        this % count = this % count + 1
    end subroutine List_append_real
    
    !$
    subroutine List_append_char(this, data)
        class(ListChar), intent(in out)  :: this
        character(len=*), intent(in)     :: data
        
        type(ListElemChar), pointer :: elem
        
        ! Create element and set dat
        allocate(elem)
        elem % data = data
        
        if (.NOT. ASSOCIATED(this % head)) then
            ! If list is empty, set head and tail to new element
            this % head => elem
            this % tail => elem
        else
            ! Otherwise append element at end of list
            this % tail % next => elem
            elem % prev => this % tail
            this % tail => this % tail % next
        end if
        
        this % count = this % count + 1
    end subroutine List_append_char

    !$
    !===============================================================================================
    ! LIST_CLEAR removes all elements from the list
    !===============================================================================================
    subroutine List_clear_int(this)
        class(ListInt), intent(in out)  :: this
        
        type(ListElemInt), pointer :: current => null()
        type(ListElemInt), pointer :: next => null()
        
        if (this % count > 0) then
            current => this % head
            do while (ASSOCIATED(current))
                ! Set pointer to next element
                next => current % next
        
                ! Deallocate memory for current element
                deallocate(current)
        
                ! Move to next element
                current => next
            end do
        
            nullify(this % head)
            nullify(this % tail)
            this % count = 0
        end if
    end subroutine List_clear_int
    
    !$
    subroutine List_clear_real(this)
        class(ListReal), intent(in out)  :: this
        
        type(ListElemReal), pointer :: current => null()
        type(ListElemReal), pointer :: next => null()
        
        if (this % count > 0) then
            current => this % head
            do while (ASSOCIATED(current))
                ! Set pointer to next element
                next => current % next
        
                ! Deallocate memory for current element
                deallocate(current)
        
                ! Move to next element
                current => next
            end do
        
            nullify(this % head)
            nullify(this % tail)
            this % count = 0
        end if
    end subroutine List_clear_real
    
    !$
    subroutine List_clear_char(this)
        class(ListChar), intent(in out)  :: this
        
        type(ListElemChar), pointer :: current => null()
        type(ListElemChar), pointer :: next => null()
        
        if (this % count > 0) then
            current => this % head
            do while (ASSOCIATED(current))
                ! Set pointer to next element
                next => current % next
        
                ! Deallocate memory for current element
                deallocate(current)
        
                ! Move to next element
                current => next
            end do
        
            nullify(this % head)
            nullify(this % tail)
            this % count = 0
        end if
    end subroutine List_clear_char

    !$
    !===============================================================================================
    ! LIST_CONTAINS determines whether the list contains a specified item. Since it
    ! relies on the index method, it is O(n).
    !===============================================================================================
    function List_contains_int(this, data) result(in_list)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: data
        logical  :: in_list
        
        in_list = (this % index(data) > 0)
    end function List_contains_int
    
    !$
    function List_contains_real(this, data) result(in_list)
        class(ListReal), intent(in out)  :: this
        real(KREAL), intent(in)      :: data
        logical  :: in_list
    
        in_list = (this % index(data) > 0)
    end function List_contains_real
    
    !$
    function List_contains_char(this, data) result(in_list)
        class(ListChar), intent(in out)  :: this
        character(len=*), intent(in)     :: data
        logical  :: in_list
    
        in_list = (this % index(data) > 0)
    end function List_contains_char

    !$
    !===============================================================================================
    ! LIST_GET_ITEM returns the item in the list at position 'i_list'. If the index
    ! is out of bounds, an error code is returned.
    !===============================================================================================
    function List_get_item_int(this, i_list) result(data)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: i_list
        integer  :: data
        
        integer :: last_index
        
        if (i_list < 1 .or. i_list > this % count) then
            ! Check for index out of bounds
            data = INT_INFINITY
        elseif (i_list == 1) then
            data = this % head % data
            this % last_index = 1
            this % last_elem => this % head
        elseif (i_list == this % count) then
            data = this % tail % data
            this % last_index = this % count
            this % last_elem => this % tail
        else
            if (i_list < this % last_index) then
                this % last_index = 1
                this % last_elem => this % head
            end if
        
            do last_index = this % last_index + 1, i_list
                this % last_elem => this % last_elem % next
                this % last_index = last_index
            end do
            data = this % last_elem % data
        end if
    end function List_get_item_int
    
    !$
    function List_get_item_real(this, i_list) result(data)
        class(ListReal), intent(in out)  :: this
        integer, intent(in)              :: i_list
        real(KREAL)  :: data
        
        integer :: last_index
        
        if (i_list < 1 .or. i_list > this % count) then
            ! Check for index out of bounds
            data = REAL_INFINITY
        elseif (i_list == 1) then
            data = this % head % data
            this % last_index = 1
            this % last_elem => this % head
        elseif (i_list == this % count) then
            data = this % tail % data
            this % last_index = this % count
            this % last_elem => this % tail
        else
            if (i_list < this % last_index) then
                this % last_index = 1
                this % last_elem => this % head
            end if
        
            do last_index = this % last_index + 1, i_list
                this % last_elem => this % last_elem % next
                this % last_index = last_index
            end do
            data = this % last_elem % data
        end if
    end function List_get_item_real
    
    !$
    function List_get_item_char(this, i_list) result(data)
        class(ListChar), intent(in out)  :: this
        integer, intent(in)              :: i_list
        character(len=MAX_WORD_LEN) :: data
        
        integer :: last_index
        
        if (i_list < 1 .or. i_list > this % count) then
            ! Check for index out of bounds
            data = ""
        elseif (i_list == 1) then
            data = this % head % data
            this % last_index = 1
            this % last_elem => this % head
        elseif (i_list == this % count) then
            data = this % tail % data
            this % last_index = this % count
            this % last_elem => this % tail
        else
            if (i_list < this % last_index) then
                this % last_index = 1
                this % last_elem => this % head
            end if
        
            do last_index = this % last_index + 1, i_list
                this % last_elem => this % last_elem % next
                this % last_index = last_index
            end do
            data = this % last_elem % data
        end if
    end function List_get_item_char

    !$
    !===============================================================================================
    ! LIST_INDEX determines the first index in the list that contains 'data'. If
    ! 'data' is not present in the list, the return value is -1.
    !===============================================================================================
    function List_index_int(this, data) result(i_list)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: data
        integer  :: i_list
        
        type(ListElemInt), pointer :: elem
        
        i_list = 0
        elem => this % head
        do while (ASSOCIATED(elem))
            i_list = i_list + 1
            if (data == elem % data) exit
            elem => elem % next
        end do
        
        ! Check if we reached the end of the list
        if (.NOT. ASSOCIATED(elem)) i_list = -1
    end function List_index_int
    
    !$
    function List_index_real(this, data) result(i_list)
        class(ListReal), intent(in out)  :: this
        real(KREAL), intent(in)      :: data
        integer  :: i_list
        
        type(ListElemReal), pointer :: elem
        
        i_list = 0
        elem => this % head
        do while (ASSOCIATED(elem))
            i_list = i_list + 1
            if (data == elem % data) exit
            elem => elem % next
        end do
        
        ! Check if we reached the end of the list
        if (.NOT. ASSOCIATED(elem)) i_list = -1
    end function List_index_real
    
    !$
    function List_index_char(this, data) result(i_list)
        class(ListChar), intent(in out)  :: this
        character(len=*), intent(in)     :: data
        integer  :: i_list
        
        type(ListElemChar), pointer :: elem
        
        i_list = 0
        elem => this % head
        do while (ASSOCIATED(elem))
            i_list = i_list + 1
            if (data == elem % data) exit
            elem => elem % next
        end do
        
        ! Check if we reached the end of the list
        if (.NOT. ASSOCIATED(elem)) i_list = -1
    end function List_index_char

    !$
    !===============================================================================================
    ! LIST_INSERT inserts 'data' at index 'i_list' within the list. If 'i_list'
    ! exceeds the size of the list, the data is appends at the end of the list.
    !===============================================================================================
    subroutine List_insert_int(this, i_list, data)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: i_list
        integer, intent(in)             :: data
        
        integer :: i
        type(ListElemInt), pointer :: elem => null()
        type(ListElemInt), pointer :: new_elem => null()
        
        if (i_list > this % count) then
            ! Check whether specified index is greater than number of elements -- if
            ! so, just append it to the end of the list
            call this % append(data)
        
        else if (i_list == 1) then
            ! Check for new head element
            allocate(new_elem)
            new_elem % data = data
            new_elem % next => this % head
            this % head => new_elem
            this % count = this % count + 1
        
        else
            ! Default case with new element somewhere in middle of list
            if (i_list >= this % last_index) then
                i = this % last_index
                elem => this % last_elem
            else
                i = 0
                elem => this % head
            end if
            do while (ASSOCIATED(elem))
                i = i + 1
                if (i == i_list - 1) then
                    ! Allocate new element
                    allocate(new_elem)
                    new_elem % data = data
        
                    ! Put it before the i-th element
                    new_elem % prev => elem % prev
                    new_elem % next => elem
                    new_elem % prev % next => new_elem
                    new_elem % next % prev => new_elem
                    this % count = this % count + 1
                    this % last_index = i_list
                    this % last_elem => new_elem
                    exit
                end if
                i = i + 1
                elem => elem % next
            end do
        end if
    end subroutine List_insert_int
    
    !$
    subroutine List_insert_real(this, i_list, data)
        class(ListReal), intent(in out) :: this
        integer, intent(in)             :: i_list
        real(KREAL), intent(in)     :: data
        
        integer :: i
        type(ListElemReal), pointer :: elem => null()
        type(ListElemReal), pointer :: new_elem => null()
        
        if (i_list > this % count) then
            ! Check whether specified index is greater than number of elements -- if
            ! so, just append it to the end of the list
            call this % append(data)
        
        else if (i_list == 1) then
            ! Check for new head element
            allocate(new_elem)
            new_elem % data = data
            new_elem % next => this % head
            this % head => new_elem
            this % count = this % count + 1
        
        else
            ! Default case with new element somewhere in middle of list
            if (i_list >= this % last_index) then
                i = this % last_index
                elem => this % last_elem
            else
                i = 0
                elem => this % head
            end if
            do while (ASSOCIATED(elem))
                if (i == i_list) then
                    ! Allocate new element
                    allocate(new_elem)
                    new_elem % data = data
        
                    ! Put it before the i-th element
                    new_elem % prev => elem % prev
                    new_elem % next => elem
                    new_elem % prev % next => new_elem
                    new_elem % next % prev => new_elem
                    this % count = this % count + 1
                    this % last_index = i_list
                    this % last_elem => new_elem
                    exit
                end if
                i = i + 1
                elem => elem % next
            end do
        end if
    end subroutine List_insert_real
    
    !$
    subroutine List_insert_char(this, i_list, data)
        class(ListChar), intent(in out) :: this
        integer, intent(in)             :: i_list
        character(len=*), intent(in)    :: data
        
        integer :: i
        type(ListElemChar), pointer :: elem => null()
        type(ListElemChar), pointer :: new_elem => null()
        
        if (i_list > this % count) then
            ! Check whether specified index is greater than number of elements -- if
            ! so, just append it to the end of the list
            call this % append(data)
        
        else if (i_list == 1) then
            ! Check for new head element
            allocate(new_elem)
            new_elem % data = data
            new_elem % next => this % head
            this % head => new_elem
            this % count = this % count + 1
        
        else
            ! Default case with new element somewhere in middle of list
            if (i_list >= this % last_index) then
                i = this % last_index
                elem => this % last_elem
            else
                i = 0
                elem => this % head
            end if
            do while (ASSOCIATED(elem))
                if (i == i_list) then
                    ! Allocate new element
                    allocate(new_elem)
                    new_elem % data = data
        
                    ! Put it before the i-th element
                    new_elem % prev => elem % prev
                    new_elem % next => elem
                    new_elem % prev % next => new_elem
                    new_elem % next % prev => new_elem
                    this % count = this % count + 1
                    this % last_index = i_list
                    this % last_elem => new_elem
                    exit
                end if
                i = i + 1
                elem => elem % next
            end do
        end if
    end subroutine List_insert_char

    !$
    !===============================================================================================
    ! LIST_REMOVE removes the first item in the list that contains 'data'. If 'data'
    ! is not in the list, no action is taken.
    !===============================================================================================
    subroutine List_remove_int(this, data)
        class(ListInt), intent(in out)  :: this
        integer, intent(in)             :: data
        
        type(ListElemInt), pointer :: elem => null()
        
        elem => this % head
        do while (ASSOCIATED(elem))
            ! Check for matching data
            if (elem % data == data) then
        
            ! Determine whether the current element is the head, tail, or a middle
            ! element
            if (ASSOCIATED(elem, this % head)) then
                this % head => elem % next
                if (ASSOCIATED(elem, this % tail)) nullify(this % tail)
                if (ASSOCIATED(this % head)) nullify(this % head % prev)
                deallocate(elem)
            else if (ASSOCIATED(elem, this % tail)) then
                this % tail => elem % prev
                deallocate(this % tail % next)
            else
                elem % prev % next => elem % next
                elem % next % prev => elem % prev
                deallocate(elem)
            end if
        
            ! Decrease count and exit
            this % count = this % count - 1
            exit
            end if
        
            ! Advance pointers
            elem => elem % next
        end do
    end subroutine List_remove_int
    
    !$
    subroutine List_remove_real(this, data)
        class(ListReal), intent(in out)  :: this
        real(KREAL), intent(in)      :: data
        
        type(ListElemReal), pointer :: elem => null()
        
        elem => this % head
        do while (ASSOCIATED(elem))
            ! Check for matching data
            if (elem % data == data) then
        
            ! Determine whether the current element is the head, tail, or a middle
            ! element
            if (ASSOCIATED(elem, this % head)) then
                this % head => elem % next
                if (ASSOCIATED(elem, this % tail)) nullify(this % tail)
                if (ASSOCIATED(this % head)) nullify(this % head % prev)
                deallocate(elem)
            else if (ASSOCIATED(elem, this % tail)) then
                this % tail => elem % prev
                deallocate(this % tail % next)
            else
                elem % prev % next => elem % next
                elem % next % prev => elem % prev
                deallocate(elem)
            end if
        
            ! Decrease count and exit
            this % count = this % count - 1
            exit
            end if
        
            ! Advance pointers
            elem => elem % next
        end do
    end subroutine List_remove_real
    
    !$
    subroutine List_remove_char(this, data)
        class(ListChar), intent(in out)  :: this
        character(len=*), intent(in)     :: data
        
        type(ListElemChar), pointer :: elem => null()
        
        elem => this % head
        do while (ASSOCIATED(elem))
            ! Check for matching data
            if (elem % data == data) then
        
            ! Determine whether the current element is the head, tail, or a middle
            ! element
            if (ASSOCIATED(elem, this % head)) then
                this % head => elem % next
                if (ASSOCIATED(elem, this % tail)) nullify(this % tail)
                if (ASSOCIATED(this % head)) nullify(this % head % prev)
                deallocate(elem)
            else if (ASSOCIATED(elem, this % tail)) then
                this % tail => elem % prev
                deallocate(this % tail % next)
            else
                elem % prev % next => elem % next
                elem % next % prev => elem % prev
                deallocate(elem)
            end if
        
            ! Decrease count and exit
            this % count = this % count - 1
            exit
            end if
        
            ! Advance pointers
            elem => elem % next
        end do
    end subroutine List_remove_char

    !$
    !===============================================================================================
    ! LIST_SIZE returns the number of elements in the list
    !===============================================================================================
    function List_size_int(this) result(size)
        class(ListInt), intent(in out) :: this
        integer        :: size
        
        size = this % count
    end function List_size_int
    
    !$
    function List_size_real(this) result(size)
        class(ListReal), intent(in out) :: this
        integer         :: size
        
        size = this % count
    end function List_size_real
    
    !$
    function List_size_char(this) result(size)
        class(ListChar), intent(in out) :: this
        integer         :: size
        
        size = this % count
    end function List_size_char

end module list_header
