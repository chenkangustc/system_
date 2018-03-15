!$
!===================================================================================================
!
!   module of vector operation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    connect_vector
!                               get_vector_error
!                               vector_to_matrix
!                               matrix_to_vector
!                               normalize_vector
!                               sort_vector
!
!   Public type lists:          No
!
!===================================================================================================
module vector_operation
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,       only : WarningCollector, ErrorCollector
    use set_header,             only : SetReal
    use string,                 only : Int_to_string
    
    implicit none 
    private
    public  :: connect_vector, get_vector_error, vector_to_matrix, matrix_to_vector, normalize_vector, sort_vector 
    
    interface  get_vector_error
        module procedure  vector_relative_error
        module procedure  matrix_relative_error
    end interface get_vector_error
    
    interface  vector_to_matrix
        module procedure vector_to_vector
        module procedure vector_to_array_2d
        module procedure vector_to_array_3d
    end interface vector_to_matrix
    
    interface  matrix_to_vector
        module procedure vector_to_vector
        module procedure array2d_to_vector
        module procedure array3d_to_vector
    end interface matrix_to_vector
    
    interface  normalize_vector
        module procedure normalize_vector1d_with_mask
        module procedure normalize_vector2d_with_mask
        module procedure normalize_vector1d_with_zero
        module procedure normalize_vector2d_with_zero
    end interface normalize_vector
    
    ! interpolation scheme
    integer, parameter  :: NORM_1        = 1                                    ! relative error by 1-norm
    integer, parameter  :: NORM_2        = 2                                    ! relative error by 2-norm
    integer, parameter  :: NORM_INF      = 3                                    ! relative error by infinity-norm
    integer, parameter  :: BY_DNTR       = 4                                    ! original method in DNTR
    integer, parameter  :: RMSE          = 5                                    ! relative mean square error
    
    type(WarningCollector) :: a_warning
    type(ErrorCollector)   :: a_error
    
contains
    !$
    !===============================================================================================
    ! join several vectors
    !===============================================================================================
    subroutine connect_vector (summary, a, b, c, d, e, f, g, h)
        
        real(KREAL), intent(in out), allocatable  :: summary(:)
        real(KREAL), intent(in)            :: a(:)
        real(KREAL), intent(in), optional  :: b(:)
        real(KREAL), intent(in), optional  :: c(:)
        real(KREAL), intent(in), optional  :: d(:)
        real(KREAL), intent(in), optional  :: e(:)
        real(KREAL), intent(in), optional  :: f(:)
        real(KREAL), intent(in), optional  :: g(:)
        real(KREAL), intent(in), optional  :: h(:)
        
        integer  :: num
        integer  :: ibeg, iend 
        integer  :: i_allocate
        
        num = SIZE(a)
        if (PRESENT(b))     num = num + SIZE(b)
        if (PRESENT(c))     num = num + SIZE(c)
        if (PRESENT(d))     num = num + SIZE(d)
        if (PRESENT(e))     num = num + SIZE(e)
        if (PRESENT(f))     num = num + SIZE(f)
        if (PRESENT(g))     num = num + SIZE(g)
        if (PRESENT(h))     num = num + SIZE(h)
        
        allocate(summary(num), stat=i_allocate)
        
        iend = 0
        ibeg = iend + 1
        iend = iend + SIZE(a)
        summary(ibeg: iend) = a
        
        if (PRESENT(b))  then
            ibeg = iend + 1
            iend = iend + SIZE(b)
            summary(ibeg: iend) = b
        end if
        if (PRESENT(c))  then
            ibeg = iend + 1
            iend = iend + SIZE(c)
            summary(ibeg: iend) = c
        end if
        if (PRESENT(d))  then
            ibeg = iend + 1
            iend = iend + SIZE(d)
            summary(ibeg: iend) = d
        end if
        if (PRESENT(e))  then
            ibeg = iend + 1
            iend = iend + SIZE(e)
            summary(ibeg: iend) = e
        end if
        if (PRESENT(f))  then
            ibeg = iend + 1
            iend = iend + SIZE(f)
            summary(ibeg: iend) = f
        end if
        if (PRESENT(g))  then
            ibeg = iend + 1
            iend = iend + SIZE(g)
            summary(ibeg: iend) = g
        end if
        if (PRESENT(h))  then
            ibeg = iend + 1
            iend = iend + SIZE(h)
            summary(ibeg: iend) = h
        end if
    
    end subroutine connect_vector
    
    !$
    !===============================================================================================
    ! input is a vector
    !===============================================================================================
    function vector_relative_error (old, new, scheme)  result(output)
        
        real(KREAL), intent(in)  :: old(:)
        real(KREAL), intent(in)  :: new(:)
        integer, intent(in)      :: scheme
        real(KREAL)    :: output
        
        ! local variables
        real(KREAL)  :: delta(SIZE(old))
        real(KREAL)  :: delta_norm
        real(KREAL)  :: old_norm
        real(KREAL)  :: error_max, error_min, ratio
        real(KREAL)  :: tmp
        integer  :: scale, i
        
        scale = SIZE(old)
        delta = new - old
        delta_norm = 0.0D0
        old_norm = 0.0D0
        
        select case(scheme)
        case(NORM_1)
            do i = 1, scale
                delta_norm = delta_norm + ABS(delta(i))
                old_norm = old_norm + ABS(old(i))
            end do
            output = delta_norm / old_norm
            
        case(NORM_2)
            do i = 1, scale
                delta_norm = delta_norm + ABS(delta(i))**2
                old_norm = old_norm + ABS(old(i))**2
            end do
            output = SQRT(delta_norm) / SQRT(old_norm)
        
        case(NORM_INF)
            do i = 1, scale
                if (ABS(delta(i)) >= delta_norm)  then
                    delta_norm = ABS(delta(i))
                end if
                if (ABS(old(i)) >= old_norm)  then
                    old_norm = ABS(old(i))
                end if
            end do
            output = delta_norm / old_norm
            
        case(BY_DNTR)
            error_max = 0.0D0
            error_min = HUGE(1.0_KREAL)
            do i = 1, scale
                if (ABS(old(i)) >= EPS_ZERO)  then
                    ratio = ABS(new(i) / old(i))
                    error_max = MAX(error_max, ratio)
                    error_min = MIN(error_min, ratio)
                end if
            end do
            output = ABS((error_max - error_min) / error_max)
            
        case(RMSE)
            do i = 1, scale
                delta_norm = delta_norm + ABS(delta(i))**2
            end do
            output = SQRT(delta_norm/scale)
    
        case default
            call a_error%set (INFO_LIST_FRAMEWORK, 'normal scheme is out of pre-defined range: ' // Int_to_string(scheme))
            call a_error%print (OUTPUT_UNIT)
        end select
    
    end function vector_relative_error
    
    !$
    !===============================================================================================
    ! input is a matrix
    !===============================================================================================
    function matrix_relative_error (old, new, scheme)  result(output)
        
        real(KREAL), intent(in)  :: old(:, :)
        real(KREAL), intent(in)  :: new(:, :)
        integer, intent(in)      :: scheme
        real(KREAL)    :: output
        
        ! local variables
        real(KREAL)  :: delta(SIZE(old,dim=1), SIZE(old, dim=2))
        integer  :: scale_i, i
        integer  :: scale_j, j
        real(KREAL)  :: delta_norm
        real(KREAL)  :: old_norm
        real(KREAL)  :: error_max, error_min, ratio
        real(KREAL)  :: tmp
        
        scale_i = SIZE(old, dim=1)
        scale_j = SIZE(new, dim=2)
        delta = new - old
        delta_norm = 0.0D0
        old_norm = 0.0D0
        
        select case(scheme)
        case(NORM_1)
            do i = 1, scale_i
                do j = 1, scale_j
                    delta_norm = delta_norm + ABS(delta(i, j))
                    old_norm = old_norm + ABS(old(i, j))
                end do
            end do
            output = delta_norm / old_norm
            
        case(NORM_2)
            do i = 1, scale_i
                do j = 1, scale_j
                    delta_norm = delta_norm + ABS(delta(i, j))**2
                    old_norm = old_norm + ABS(old(i, j))**2
                end do
            end do
            output = SQRT(delta_norm) / SQRT(old_norm)
        
        case(NORM_INF)
            do i = 1, scale_i
                do j = 1, scale_j
                    if (ABS(delta(i, j)) >= delta_norm)  then
                        delta_norm = ABS(delta(i, j))
                    end if
                    if (ABS(old(i, j)) >= old_norm)  then
                        old_norm = ABS(old(i, j))
                    end if
                end do
            end do
            output = delta_norm / old_norm
            
        case(BY_DNTR)
            error_max = 0.0D0
            error_min = HUGE(1.0_KREAL)
            do i = 1, scale_i
                do j = 1, scale_j
                    if (ABS(old(i, j)) >= EPS_ZERO)  then
                        ratio = ABS(new(i, j) / old(i, j))
                        error_max = MAX(error_max, ratio)
                        error_min = MIN(error_min, ratio)
                    end if
                end do
            end do
            output = ABS((error_max - error_min) / error_max)
            
        case(RMSE)
            do i = 1, scale_i
                do j = 1, scale_j
                    delta_norm = delta_norm + ABS(delta(i, j))**2
                end do
            end do
            output = SQRT(delta_norm / (scale_i*scale_j))
    
        case default
            call a_error%set (INFO_LIST_FRAMEWORK, 'normal scheme is out of pre-defined range: ' // Int_to_string(scheme))
            call a_error%print (OUTPUT_UNIT)
        end select
        
    end function matrix_relative_error
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! assign vector to vector
    !===============================================================================================
    subroutine vector_to_vector (left, right)
        
        real(KREAL), intent(in out)  :: left(:)
        real(KREAL), intent(in)      :: right(:)
        integer  :: i, j, idx
        
        if (SIZE(left) /= SIZE(right))  then
            call a_warning%set (INFO_LIST_FRAMEWORK, 'the vector size is not consistent')
            call a_warning%print (OUTPUT_UNIT)
        end if
        
        idx = 0
        do i = 1, SIZE(left)
            idx = idx + 1
            left(i) = right(idx)
        end do
        
    end subroutine vector_to_vector
    
    !$
    !===============================================================================================
    ! assign vector to 2d array (row first)
    !===============================================================================================
    subroutine vector_to_array_2d (left, right)
        
        real(KREAL), intent(in out)  :: left(:, :)
        real(KREAL), intent(in)      :: right(:)
        integer  :: i, j, idx
        
        if (SIZE(left) /= SIZE(right))  then
            call a_warning%set (INFO_LIST_FRAMEWORK, 'the vector size is not consistent')
            call a_warning%print (OUTPUT_UNIT)
        end if
        
        idx = 0
        do i = 1, SIZE(left, dim=1)
            do j = 1, SIZE(left, dim=2)
                idx = idx + 1
                left(i, j) = right(idx)
            end do
        end do
    
    end subroutine vector_to_array_2d
    
    !$
    !===============================================================================================
    ! assign vector to 3d array (row first)
    !===============================================================================================
    subroutine vector_to_array_3d (left, right)
        
        real(KREAL), intent(in out)  :: left(:, :, :)
        real(KREAL), intent(in)      :: right(:)
        integer  :: i, j, k, idx
        
        if (SIZE(left) /= SIZE(right))  then
            call a_warning%set (INFO_LIST_FRAMEWORK, 'the vector size is not consistent')
            call a_warning%print (OUTPUT_UNIT)
        end if
        
        idx = 0
        do i = 1, SIZE(left, dim=1)
            do j = 1, SIZE(left, dim=2)
                do k = 1, SIZE(left, dim=3)
                    idx = idx + 1
                    left(i, j, k) = right(idx)
                end do
            end do
        end do
        
    end subroutine vector_to_array_3d
    
    !$
    !===============================================================================================
    ! assign 2d array to vector (row first)
    !===============================================================================================
    subroutine array2d_to_vector (left, right)
        
        real(KREAL), intent(in out)  :: left(:)
        real(KREAL), intent(in)      :: right(:, :)
        integer  :: i, j, idx
        
        if (SIZE(left) /= SIZE(right))  then
            call a_warning%set (INFO_LIST_FRAMEWORK, 'the vector size is not consistent')
            call a_warning%print (OUTPUT_UNIT)
        end if
        
        idx = 0
        do i = 1, SIZE(right, dim=1)
            do j = 1, SIZE(right, dim=2)
                idx = idx + 1
                left(idx) = right(i, j)
            end do
        end do
    
    end subroutine array2d_to_vector
    
    !$
    !===============================================================================================
    ! assign 3d array to vector (row first)
    !===============================================================================================
    subroutine array3d_to_vector (left, right)
        
        real(KREAL), intent(in out)  :: left(:)
        real(KREAL), intent(in)      :: right(:, :, :)
        integer  :: i, j, k, idx
        
        if (SIZE(left) /= SIZE(right))  then
            call a_warning%set (INFO_LIST_FRAMEWORK, 'the vector size is not consistent')
            call a_warning%print (OUTPUT_UNIT)
        end if
        
        idx = 0
        do i = 1, SIZE(right, dim=1)
            do j = 1, SIZE(right, dim=2)
                do k = 1, SIZE(right, dim=3)
                    idx = idx + 1
                    left(idx) = right(i, j, k)
                end do
            end do
        end do
        
    end subroutine array3d_to_vector
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine normalize_vector1d_with_mask (vector, mask, weight)
        
        real(KREAL), intent(in out)        :: vector(:)
        logical, intent(in), optional      :: mask(:)
        real(KREAL), intent(in), optional  :: weight(:)
        
        real(KREAL)  :: vector_(SIZE(vector))
        logical      :: mask_(SIZE(vector))
        real(KREAL)  :: weight_(SIZE(vector))
        real(KREAL)  :: total, summary
        integer      :: i
        
        mask_ = .TRUE.
        if (PRESENT(mask))  then
            mask_ = mask
        end if
        weight_ = 1.0D0
        if (PRESENT(weight))  then
            weight_ = weight
        end if
        
        vector_ = vector
        total = 0.0D0
        summary = 0.0D0       
        do i = 1, SIZE(vector_, dim=1)
            if (mask_(i))  then
                total = total + vector_(i)*weight_(i)
                summary = summary + weight_(i)
            end if
        end do
        
        vector = 0.0D0
        do i = 1, SIZE(vector_, dim=1)
            vector(i) = vector_(i)*weight_(i) * (summary / total)
        end do
        
    end subroutine normalize_vector1d_with_mask
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine normalize_vector2d_with_mask (matrix, mask, weight)
        
        real(KREAL), intent(in out)        :: matrix(:, :)
        logical, intent(in), optional      :: mask(:, :)
        real(KREAL), intent(in), optional  :: weight(:, :)
        
        real(KREAL)  :: matrix_(SIZE(matrix,dim=1), SIZE(matrix,dim=2))
        logical      :: mask_(SIZE(matrix,dim=1), SIZE(matrix,dim=2))
        real(KREAL)  :: weight_(SIZE(matrix,dim=1), SIZE(matrix,dim=2))
        real(KREAL)  :: total, summary
        integer      :: i, j
        
        mask_ = .TRUE.
        if (PRESENT(mask))  then
            mask_ = mask
        end if
        weight_ = 1.0D0
        if (PRESENT(weight))  then
            weight_ = weight
        end if

        matrix_ = matrix
        summary = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(matrix_, dim=1)
            do j = 1, SIZE(matrix_, dim=2)
                if (mask_(i, j))  then
                    summary = summary + weight_(i, j)
                    total = total + matrix_(i, j)*weight_(i, j)
                end if
            end do
        end do
        
        matrix = 0.0D0
        do i = 1, SIZE(matrix_, dim=1)
            do j = 1, SIZE(matrix_, dim=2)
                matrix(i, j) = matrix_(i, j)*weight_(i, j) * (summary / total)
            end do
        end do
        
    end subroutine normalize_vector2d_with_mask
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine normalize_vector1d_with_zero (vector, is_zero, weight)
        
        real(KREAL), intent(in out)        :: vector(:)
        logical, intent(in)                :: is_zero
        real(KREAL), intent(in), optional  :: weight(:)
        
        real(KREAL)  :: vector_(SIZE(vector))
        real(KREAL)  :: weight_(SIZE(vector))
        real(KREAL)  :: total, summary
        integer      :: i
        
        weight_ = 1.0D0
        if (PRESENT(weight))  then
            weight_ = weight
        end if
        
        vector_ = vector
        summary = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if ((is_zero) .or. (vector(i) > EPS_ZERO))  then
                summary = summary + weight_(i)
                total = total + vector_(i)*weight_(i)
            end if
        end do
        
        vector = 0.0D0
        do i = 1, SIZE(vector_, dim=1)
            vector(i) = vector_(i)*weight_(i) * (summary / total)
        end do
        
    end subroutine normalize_vector1d_with_zero
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine normalize_vector2d_with_zero (matrix, is_zero, weight)
        
        real(KREAL), intent(in out)        :: matrix(:, :)
        logical, intent(in)                :: is_zero
        real(KREAL), intent(in), optional  :: weight(:, :)
        
        real(KREAL)  :: matrix_(SIZE(matrix,dim=1), SIZE(matrix,dim=2))
        real(KREAL)  :: weight_(SIZE(matrix,dim=1), SIZE(matrix,dim=2))
        real(kind=8)     :: total, summary
        integer          :: i, j
        
        weight_ = 1.0D0
        if (PRESENT(weight))  then
            weight_ = weight
        end if
        
        matrix_ = matrix
        summary = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if ((is_zero) .or. (matrix(i,j) > EPS_ZERO))  then
                    summary = summary + weight_(i, j)
                    total = total + matrix_(i, j)*weight_(i, j)
                end if
            end do
        end do
        
        matrix = 0.0D0
        do i = 1, SIZE(matrix_, dim=1)
            do j = 1, SIZE(matrix_, dim=2)
                matrix(i, j) = matrix_(i, j)*weight_(i, j) * (summary / total)
            end do
        end do
    
    end subroutine normalize_vector2d_with_zero
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine sort_vector (vec_full, vec_reduce, num)
        
        real(KREAL), intent(in)  :: vec_full(:)
        real(KREAL), intent(in out), allocatable :: vec_reduce(:)
        integer, intent(out)  :: num 
        
        type(SetReal)  :: vec_
        real(KREAL)    :: tmp
        integer  :: i, j
        integer  :: i_allocate 
        
        do i = 1, SIZE(vec_full)
            call vec_%add (vec_full(i))
        end do 
        
        num = vec_%size ()
        allocate(vec_reduce(num), stat=i_allocate)
        do i = 1, num
            vec_reduce(i) = vec_%get_item (i)
        end do 
        
        ! increasing 
        do i = 1, num
            do j = i, num 
                if (vec_reduce(j) <= vec_reduce(i))  then
                    tmp = vec_reduce(i)
                    vec_reduce(i) = vec_reduce(j)
                    vec_reduce(j) = tmp 
                end if 
            end do 
        end do 
        
    end subroutine sort_vector
    
end module vector_operation
