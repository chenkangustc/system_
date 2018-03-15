!$
!===================================================================================================
!
!   a set class based on list class, this results in much worse performance than an implementation 
!   based on hash tables or binary trees
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          SetInt
!                               SetReal 
!                               SetChar
!
!===================================================================================================
module set_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use list_header

    implicit none
    private
    public  :: SetInt, SetReal, SetChar 

    ! --------------------------------------------------------------------------
    ! SET contains a list of elements and methods to add, remove, and perform other basic tasks
    type :: SetInt
        private
        type(ListInt) :: elements
    contains
        procedure, public :: add => Set_add_int
        procedure, public :: clear => Set_clear_int
        procedure, public :: contains => Set_contains_int
        procedure, public :: get_item => Set_get_item_int
        procedure, public :: remove => Set_remove_int
        procedure, public :: size => Set_size_int
    end type SetInt

    type :: SetReal
        private
        type(ListReal) :: elements
    contains
        procedure, public :: add => Set_add_real
        procedure, public :: clear => Set_clear_real
        procedure, public :: contains => Set_contains_real
        procedure, public :: get_item => Set_get_item_real
        procedure, public :: remove => Set_remove_real
        procedure, public :: size => Set_size_real
    end type SetReal

    type :: SetChar
        private
        type(ListChar) :: elements
    contains
        procedure, public :: add => Set_add_char
        procedure, public :: clear => Set_clear_char
        procedure, public :: contains => Set_contains_char
        procedure, public :: get_item => Set_get_item_char
        procedure, public :: remove => Set_remove_char
        procedure, public :: size => Set_size_char
    end type SetChar

contains
    !$
    !===============================================================================================
    ! SET_ADD adds an item to a set if it is not already present in the set
    !===============================================================================================
    subroutine Set_add_int(this, data)
        class(SetInt), intent(in out)  :: this
        integer, intent(in)            :: data
        
        if (.NOT. this % elements % contains(data)) then
            call this % elements % append(data)
        end if
    end subroutine Set_add_int
    
    !$
    subroutine Set_add_real(this, data)
        class(SetReal), intent(in out) :: this
        real(KREAL), intent(in)   :: data
    
        if (.NOT. this % elements % contains(data)) then
            call this % elements % append(data)
        end if
    end subroutine Set_add_real

    !$
    subroutine Set_add_char(this, data)
        class(SetChar), intent(in out) :: this
        character(len=*), intent(in)   :: data
    
        if (.NOT. this % elements % contains(data)) then
            call this % elements % append(data)
        end if
    end subroutine Set_add_char

    !$
    !===============================================================================================
    ! SET_CLEAR removes all items in a set
    !===============================================================================================
    subroutine Set_clear_int(this)
        class(SetInt), intent(in out) :: this

        call this % elements % clear()
    end subroutine Set_clear_int
    
    !$
    subroutine Set_clear_real(this)
        class(SetReal), intent(in out) :: this

        call this % elements % clear()
    end subroutine Set_clear_real

    !$
    subroutine Set_clear_char(this)
        class(SetChar), intent(in out) :: this

        call this % elements % clear()
    end subroutine Set_clear_char

    !$
    !===============================================================================================
    ! SET_CONTAINS determines if a specified item is in a set
    !===============================================================================================
    function Set_contains_int(this, data) result(in_set)
        class(SetInt), intent(in out)  :: this
        integer, intent(in)            :: data
        logical       :: in_set
    
        in_set = this % elements % contains(data)
    end function Set_contains_int
    
    !$
    function Set_contains_real(this, data) result(in_set)
        class(SetReal), intent(in out)  :: this
        real(KREAL), intent(in)    :: data
        logical        :: in_set
    
        in_set = this % elements % contains(data)
    end function Set_contains_real

    !$
    function Set_contains_char(this, data) result(in_set)
        class(SetChar), intent(in out)  :: this
        character(len=*), intent(in)    :: data
        logical        :: in_set
    
        in_set = this % elements % contains(data)
    end function Set_contains_char

    !$
    !===============================================================================================
    ! SET_GET_ITEM returns the i-th item in the set
    !===============================================================================================
    function Set_get_item_int(this, i_list) result(data)
        class(SetInt), intent(in out)  :: this
        integer, intent(in)            :: i_list
        integer :: data
    
        data = this % elements % get_item(i_list)
    end function Set_get_item_int
    
    !$
    function Set_get_item_real(this, i_list) result(data)
        class(SetReal), intent(in out)  :: this
        integer, intent(in)             :: i_list
        real(KREAL) :: data
    
        data = this % elements % get_item(i_list)
    end function Set_get_item_real 

    !$
    function Set_get_item_char(this, i_list) result(data)
        class(SetChar), intent(in out)  :: this
        integer, intent(in)             :: i_list
        character(len=MAX_WORD_LEN) :: data
    
        data = this % elements % get_item(i_list)
    end function Set_get_item_char

    !$
    !===============================================================================================
    ! SET_REMOVE removes the specified item from the set. If it is not in the set, no action is taken.
    !===============================================================================================
    subroutine Set_remove_int(this, data)
        class(SetInt), intent(in out) :: this
        integer, intent(in) :: data
    
        call this % elements % remove(data)
    end subroutine Set_remove_int
    
    !$
    subroutine Set_remove_real(this, data)
        class(SetReal), intent(in out) :: this
        real(KREAL), intent(in) :: data
    
        call this % elements % remove(data)
    end subroutine Set_remove_real

    !$
    subroutine Set_remove_char(this, data)
        class(SetChar), intent(in out) :: this
        character(len=*), intent(in) :: data
    
        call this % elements % remove(data)
    end subroutine Set_remove_char

    !$
    !===============================================================================================
    ! SET_SIZE returns the number of elements in the set
    !===============================================================================================
    function Set_size_int(this) result(size)
        class(SetInt), intent(in out)  :: this
        integer  :: size
    
        size = this % elements % size()
    end function Set_size_int
    
    !$
    function Set_size_real(this) result(size)
        class(SetReal), intent(in out)  :: this
        integer  :: size
   
        size = this % elements % size()
    end function Set_size_real

    !$
    function Set_size_char(this) result(size)
        class(SetChar), intent(in out)  :: this
        integer  :: size
   
        size = this % elements % size()
    end function Set_size_char

end module set_header
