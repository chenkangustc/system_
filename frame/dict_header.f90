!$
!===================================================================================================
!
!   a dictionary class that has (key,value) pairs which used to provide lookup features
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          DictCharInt
!                               DictIntInt
!
!===================================================================================================
module dict_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private
    public  :: DictCharInt, DictIntInt
    
    ! local parameter
    integer, parameter, private :: HASH_SIZE       = 4993
    integer, parameter, private :: HASH_MULTIPLIER = 31
    integer, parameter, private :: DICT_NULL       = -HUGE(0)
    integer, parameter          :: DICT_KEY_LENGTH = 255

    ! --------------------------------------------------------------------------
    ! ELEMKEYVALUE* contains (key,value) pairs and a pointer to the next (key,value) pair
    type, private  :: ElemKeyValueCI
        type(ElemKeyValueCI), pointer   :: next => null()
        character(len=DICT_KEY_LENGTH)  :: key
        integer                         :: value
    end type ElemKeyValueCI

    type, private  :: ElemKeyValueII
        type(ElemKeyValueII), pointer   :: next => null()
        integer                         :: key
        integer                         :: value
    end type ElemKeyValueII

    ! --------------------------------------------------------------------------
    ! HASHLIST* types contain a single pointer to a linked list of (key,value) pairs. 
    ! This type is necesssary so that the Dict types can be dynamically allocated.
    type, private :: HashListCI
        type(ElemKeyValueCI), pointer :: list => null()
    end type HashListCI

    type, private :: HashListII
        type(ElemKeyValueII), pointer :: list => null()
    end type HashListII

    ! --------------------------------------------------------------------------
    ! DICT* is a dictionary of (key,value) pairs with convenience methods as type-bound procedures. 
    ! DictCharInt has character(*) keys and integer values, and DictIntInt has integer keys and values.
    type :: DictCharInt
        private
        type(HashListCI), pointer :: table(:) => null()
    contains
        procedure, public  :: add_key => Dict_add_key_ci
        procedure, public  :: get_key => Dict_get_key_ci
        procedure, public  :: has_key => Dict_has_key_ci
        procedure, public  :: keys => Dict_keys_ci
        procedure, public  :: clear => Dict_clear_ci
        procedure, private :: get_elem => Dict_get_elem_ci
    end type DictCharInt

    type :: DictIntInt
        private
        type(HashListII), pointer :: table(:) => null()
    contains
        procedure, public  :: add_key => Dict_add_key_ii
        procedure, public  :: get_key => Dict_get_key_ii
        procedure, public  :: has_key => Dict_has_key_ii
        procedure, public  :: keys => Dict_keys_ii
        procedure, public  :: clear => Dict_clear_ii
        procedure, private :: get_elem => Dict_get_elem_ii
    end type DictIntInt

contains
    !$
    !===============================================================================================
    ! DICT_ADD_KEY adds a (key,value) entry to a dictionary. If the key is already in the dictionary, 
    !   the value is replaced by the new specified value.
    !===============================================================================================
    subroutine Dict_add_key_ci(this, key, value)
        class(DictCharInt), intent(in out)  :: this
        character(len=*), intent(in)        :: key
        integer, intent(in)                 :: value
        
        integer :: hash
        type(ElemKeyValueCI), pointer :: elem => null()
        type(ElemKeyValueCI), pointer :: new_elem => null()
        
        elem => this % get_elem(key)
        
        if (associated(elem)) then
            elem % value = value
        else
            ! Get hash 
            hash = Dict_hash_key_ci(key)
        
            ! Create new element
            allocate(new_elem)
            new_elem % key = key
            new_elem % value = value
        
            ! Add element to front of list
            new_elem % next => this % table(hash) % list
            this % table(hash) % list => new_elem
        end if
    end subroutine Dict_add_key_ci
    
    !$
    subroutine Dict_add_key_ii(this, key, value)
        class(DictIntInt), intent(in out)  :: this
        integer, intent(in)                :: key
        integer, intent(in)                :: value
        
        integer :: hash
        type(ElemKeyValueII), pointer :: elem => null()
        type(ElemKeyValueII), pointer :: new_elem => null()
        
        elem => this % get_elem(key)
        
        if (associated(elem)) then
            elem % value = value
        else
            ! Get hash 
            hash = Dict_hash_key_ii(key)
        
            ! Create new element
            allocate(new_elem)
            new_elem % key = key
            new_elem % value = value
        
            ! Add element to front of list
            new_elem % next => this % table(hash) % list
            this % table(hash) % list => new_elem
        end if
    end subroutine Dict_add_key_ii

    !$
    !===============================================================================================
    ! DICT_GET_ELEM returns a pointer to the (key,value) pair for a given key. This method is private.
    !===============================================================================================
    function Dict_get_elem_ci(this, key) result(elem)
        class(DictCharInt), intent(in out)  :: this
        character(len=*), intent(in)        :: key
        type(ElemKeyValueCI), pointer :: elem
        
        integer :: hash
        
        ! Check for dictionary not being allocated
        if (.not. associated(this % table)) then
            allocate(this % table(HASH_SIZE))
        end if
        
        hash = Dict_hash_key_ci(key)
        elem => this % table(hash) % list
        do while (associated(elem))
            if (elem % key == key) exit
            elem => elem % next
        end do
    end function Dict_get_elem_ci
    
    !$
    function Dict_get_elem_ii(this, key) result(elem)
        class(DictIntInt), intent(in out)  :: this
        integer, intent(in)                :: key
        type(ElemKeyValueII), pointer :: elem
        
        integer :: hash
        
        ! Check for dictionary not being allocated
        if (.not. associated(this % table)) then
            allocate(this % table(HASH_SIZE))
        end if
        
        hash = Dict_hash_key_ii(key)
        elem => this % table(hash) % list
        do while (associated(elem))
            if (elem % key == key) exit
            elem => elem % next
        end do
    end function Dict_get_elem_ii

    !$
    !===============================================================================================
    ! DICT_GET_KEY returns the value matching a given key. If the dictionary does not contain the key, 
    ! the value DICT_NULL is returned
    !===============================================================================================
    function Dict_get_key_ci(this, key) result(value)
        class(DictCharInt), intent(in out)  :: this
        character(len=*), intent(in)        :: key
        integer  :: value
        
        type(ElemKeyValueCI), pointer :: elem
        
        elem => this % get_elem(key)
        
        if (associated(elem)) then
            value = elem % value
        else
            value = DICT_NULL
        end if
    end function Dict_get_key_ci
    
    !$
    function Dict_get_key_ii(this, key) result(value)
        class(DictIntInt), intent(in out)  :: this
        integer, intent(in)                :: key
        integer  :: value
    
        type(ElemKeyValueII), pointer  :: elem
    
        elem => this % get_elem(key)
    
        if (associated(elem)) then
            value = elem % value
        else
            value = DICT_NULL
        end if
    end function Dict_get_key_ii

    !$
    !===============================================================================================
    ! DICT_HAS_KEY determines whether a dictionary has a (key,value) pair with a given key
    !===============================================================================================
    function Dict_has_key_ci(this, key) result(has)
        class(DictCharInt), intent(in out)  :: this
        character(len=*), intent(in)        :: key
        logical                  :: has
        
        type(ElemKeyValueCI), pointer  :: elem
        
        elem => this % get_elem(key)
        has = associated(elem)
    end function Dict_has_key_ci

    function Dict_has_key_ii(this, key) result(has)
        class(DictIntInt), intent(in out)  :: this
        integer, intent(in)                :: key
        logical             :: has
        
        type(ElemKeyValueII), pointer  :: elem
        
        elem => this % get_elem(key)
        has = associated(elem)
    end function Dict_has_key_ii

    !$
    !===============================================================================================
    ! DICT_HASH_KEY returns the hash value for a given key
    !===============================================================================================
    function Dict_hash_key_ci(key) result(value)
        character(len=*), intent(in) :: key
        integer                      :: value
        
        integer :: i
        
        value = 0
        
        do i = 1, len_trim(key)
            value = HASH_MULTIPLIER * value + ichar(key(i:i))
        end do
        
        ! Added the absolute value on value-1 since the sum in the do loop is
        ! susceptible to integer overflow
        value = 1 + mod(abs(value-1), HASH_SIZE)
    end function Dict_hash_key_ci
    
    !$
    function Dict_hash_key_ii(key) result(value)
        integer, intent(in) :: key
        integer             :: value
        
        value = 0
        
        ! Added the absolute value on value-1 since the sum in the do loop is
        ! susceptible to integer overflow
        value = 1 + mod(abs(key-1), HASH_SIZE)
    end function Dict_hash_key_ii

    !$
    !===============================================================================================
    ! DICT_KEYS returns a pointer to a linked list of all (key,value) pairs
    !===============================================================================================
    function Dict_keys_ci(this) result(keys)
        class(DictCharInt), intent(in out)  :: this
        type(ElemKeyValueCI), pointer :: keys
        
        integer :: i
        type(ElemKeyValueCI), pointer :: current => null()
        type(ElemKeyValueCI), pointer :: elem => null()
        
        keys => null()
        
        do i = 1, size(this % table)
            ! Get pointer to start of bucket i
            elem => this % table(i) % list
        
            do while (associated(elem))
                ! Allocate (key,value) pair
                if (.not. associated(keys)) then
                    allocate(keys)
                    current => keys
                else
                    allocate(current % next)
                    current => current % next
                end if
        
                ! Copy (key,value) pair
                current % key   = elem % key
                current % value = elem % value
        
                ! Move to next element in bucket i
                elem => elem % next
            end do
        end do
    end function Dict_keys_ci
    
    !$
    function Dict_keys_ii(this) result(keys)
        class(DictIntInt), intent(in out)  :: this
        type(ElemKeyValueII), pointer :: keys
        
        integer :: i
        type(ElemKeyValueII), pointer :: current => null()
        type(ElemKeyValueII), pointer :: elem => null()
        
        keys => null()
        
        do i = 1, size(this % table)
            ! Get pointer to start of bucket i
            elem => this % table(i) % list
        
            do while (associated(elem))
                ! Allocate (key,value) pair
                if (.not. associated(keys)) then
                    allocate(keys)
                    current => keys
                else
                    allocate(current % next)
                    current => current % next
                end if
        
                ! Copy (key,value) pair
                current % key   = elem % key
                current % value = elem % value
        
                ! Move to next element in bucket i
                elem => elem % next
            end do
        end do
    end function Dict_keys_ii

    !$
    !===============================================================================================
    ! DICT_CLEAR Deletes and deallocates the dictionary item
    !===============================================================================================
    subroutine Dict_clear_ci(this)
        class(DictCharInt), intent(in out) :: this
        
        integer :: i
        type(ElemKeyValueCI), pointer :: current
        type(ElemKeyValueCI), pointer :: next
        
        if (associated(this % table)) then
            do i = 1, size(this % table)
                current => this % table(i) % list
                do while (associated(current))
                    if (associated(current % next)) then
                        next => current % next
                    else
                        nullify(next)
                    end if
                    deallocate(current)
                    current => next
                end do
                if (associated(this % table(i) % list)) &
                nullify(this % table(i) % list)
            end do
            deallocate(this % table)
        end if
    end subroutine Dict_clear_ci
    
    !$
    subroutine Dict_clear_ii(this)
        class(DictIntInt), intent(in out) :: this
        
        integer :: i
        type(ElemKeyValueII), pointer :: current
        type(ElemKeyValueII), pointer :: next
        
        if (associated(this % table)) then
            do i = 1, size(this % table)
                current => this % table(i) % list
                do while (associated(current))
                    if (associated(current % next)) then
                        next => current % next
                    else
                        nullify(next)
                    end if
                    deallocate(current)
                    current => next
                end do
                if (associated(this % table(i) % list)) &
                nullify(this % table(i) % list)
            end do
            deallocate(this % table)
        end if
    end subroutine Dict_clear_ii

end module dict_header
