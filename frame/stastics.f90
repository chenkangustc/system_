!$
!===================================================================================================
!
!   performs some stastics calculation. 
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    stastics_max_value
!                               stastics_max_id
!                               stastics_min_value
!                               stastics_average_value
!                               stastics_variance
!                               stastics_std_variance
!
!   Public type lists:          No
!
!===================================================================================================
module stastics

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private 
    public  :: stastics_max_value, stastics_max_id, stastics_min_value
    public  :: stastics_average_value, stastics_variance, stastics_std_variance
    
    ! interface statement
    interface  stastics_max_value
        module procedure  stastics_max_value_vector
        module procedure  stastics_max_value_matrix
    end interface  stastics_max_value
    
    interface stastics_max_id
        module procedure  stastics_max_id_vector
        module procedure  stastics_max_id_matrix
    end interface stastics_max_id
    
    interface  stastics_min_value
        module procedure  stastics_min_value_vector
        module procedure  stastics_min_value_matrix
    end interface  stastics_min_value
    
    interface  stastics_average_value
        module procedure  stastics_average_value_vector_with_mask
        module procedure  stastics_average_value_matrix_with_mask
        module procedure  stastics_average_value_vector_with_zero
        module procedure  stastics_average_value_matrix_with_zero
        module procedure  stastics_average_value_vector_with_weightzero
        module procedure  stastics_average_value_matrix_with_weightzero
    end interface  stastics_average_value
    
    interface  stastics_variance
        module procedure  stastics_variance_vector
        module procedure  stastics_variance_matrix
    end interface  stastics_variance
    
    interface stastics_std_variance
        module procedure  stastics_std_variance_vector
        module procedure  stastics_std_variance_matrix
    end interface  stastics_std_variance
    
contains
    !$
    !===============================================================================================
    ! get maxima value of a vector with an optional mask
    !===============================================================================================
    function stastics_max_value_vector (vector, mask)  result(output)
        
        real(KREAL), intent(in)        :: vector(:)
        logical, intent(in), optional  :: mask(:)
        real(kind=8)  :: output
        
        logical  :: is_true(SIZE(vector))
        integer  :: cnt
        integer  :: index_i, i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = - HUGE(0.0E0)                                                  ! "E" while not "D" here to prevent overflow
        index_i = 1
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + 1
                if (vector(i) > output)  then
                    output = vector(i)
                    index_i = i
                end if
            end if
        end do
    
    end function stastics_max_value_vector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function stastics_max_value_matrix (matrix, mask)  result(output)
        
        real(KREAL), intent(in)        :: matrix(:, :)
        logical, intent(in), optional  :: mask(:, :)
        real(kind=8)  :: output
        
        logical  :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        integer  :: cnt
        integer  :: index_i, i
        integer  :: index_j, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = - HUGE(0.0E0)
        index_i = 1
        index_j = 1
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + 1
                    if (matrix(i, j) > output)  then
                        output = matrix(i, j)
                        index_i = i
                        index_j = j
                    end if
                end if
            end do
        end do
    
    end function stastics_max_value_matrix
    
    !$
    !===============================================================================================
    ! get maxima value of a vector with an optional mask
    !===============================================================================================
    subroutine  stastics_max_id_vector (vector, id, mask)
        
        real(KREAL), intent(in)        :: vector(:)
        integer, intent(in out)        :: id
        logical, intent(in), optional  :: mask(:)
        
        real(KREAL)  :: output
        logical  :: is_true(SIZE(vector))
        integer  :: cnt
        integer  :: index_i, i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = - HUGE(0.0E0)
        index_i = 1
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + 1
                if (vector(i) > output)  then
                    output = vector(i)
                    index_i = i
                end if
            end if
        end do
        
        id = index_i
    
    end subroutine stastics_max_id_vector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine  stastics_max_id_matrix (matrix, id, jd, mask)
        
        real(KREAL), intent(in)        :: matrix(:, :)
        integer, intent(in out)        :: id
        integer, intent(in out)        :: jd
        logical, intent(in), optional  :: mask(:, :)
        
        real(KREAL)  :: output
        logical  :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        integer  :: cnt
        integer  :: index_i, i
        integer  :: index_j, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = - HUGE(0.0E0)
        index_i = 1
        index_j = 1
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + 1
                    if (matrix(i, j) > output)  then
                        output = matrix(i, j)
                        index_i = i
                        index_j = j
                    end if
                end if
            end do
        end do
        
        id = index_i
        jd = index_j
    
    end subroutine stastics_max_id_matrix
    
    !$
    !===============================================================================================
    ! get minima value of a vector with an optional mask
    !===============================================================================================
    function stastics_min_value_vector (vector, mask)  result(output)
        
        real(KREAL), intent(in)        :: vector(:)
        logical, intent(in), optional  :: mask(:)
        real(kind=8)  :: output
        
        logical  :: is_true(SIZE(vector))
        integer  :: cnt
        integer  :: index_i, i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = HUGE(0.0E0)
        index_i = 1
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + 1
                if (vector(i) < output)  then
                    output = vector(i)
                    index_i = i
                end if
            end if
        end do
    
    end function stastics_min_value_vector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function stastics_min_value_matrix (matrix, mask)  result(output)
        
        real(KREAL), intent(in)        :: matrix(:, :)
        logical, intent(in), optional  :: mask(:, :)
        real(kind=8)  :: output
        
        logical  :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        integer  :: cnt
        integer  :: index_i, i
        integer  :: index_j, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        output = HUGE(0.0E0)
        index_i = 1
        index_j = 1
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + 1
                    if (matrix(i, j) < output)  then
                        output = matrix(i, j)
                        index_i = i
                        index_j = j
                    end if
                end if
            end do
        end do
    
    end function stastics_min_value_matrix
    
    !$
    !===============================================================================================
    ! get average value of a vector with an optional mask
    !===============================================================================================
    function stastics_average_value_vector_with_mask (vector, weight, mask)  result(output)
        
        real(KREAL), intent(in)            :: vector(:)
        real(KREAL), intent(in), optional  :: weight(:)
        logical, intent(in), optional      :: mask(:)
        real(kind=8)  :: output
        
        real(KREAL)      :: wt(SIZE(vector))
        logical          :: is_true(SIZE(vector))
        real(kind=8)     :: total
        real(kind=8)     :: cnt 
        integer          :: i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        wt = 1.0D0
        if (PRESENT(weight))  then
            wt = weight
        end if
        
        cnt = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + wt(i)
                total = total + vector(i)*wt(i)
            end if
        end do
        
        output = total / cnt
    
    end function stastics_average_value_vector_with_mask
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function stastics_average_value_matrix_with_mask (matrix, weight, mask)  result(output)
        
        real(KREAL), intent(in)            :: matrix(:, :)
        real(KREAL), intent(in), optional  :: weight(:, :)
        logical, intent(in), optional      :: mask(:, :)
        real(kind=8)  :: output
        
        real(KREAL)      :: wt(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        logical          :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        real(kind=8)     :: total
        real(kind=8)     :: cnt 
        integer          :: i, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        wt = 1.0D0
        if (PRESENT(weight))  then
            wt = weight
        end if
        
        cnt = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + wt(i, j)
                    total = total + matrix(i, j)*wt(i, j)
                end if
            end do
        end do
        
        output = total / cnt
    
    end function stastics_average_value_matrix_with_mask
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    function stastics_average_value_vector_with_zero (vector, is_zero)  result(output)
        
        real(KREAL), intent(in)    :: vector(:)
        logical, intent(in)        :: is_zero
        real(kind=8)  :: output
        
        real(kind=8)     :: total
        integer          :: cnt 
        integer          :: i
        
        cnt = 0
        total = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if ((is_zero) .or. (vector(i) > EPS_ZERO) .or. (vector(i) < -EPS_ZERO))  then
                cnt = cnt + 1
                total = total + vector(i)
            end if
        end do
        
        output = total / cnt
    
    end function stastics_average_value_vector_with_zero
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    function stastics_average_value_matrix_with_zero (matrix, is_zero)  result(output)
        
        real(KREAL), intent(in)    :: matrix(:, :)
        logical, intent(in)        :: is_zero
        real(kind=8)  :: output
        
        real(kind=8)     :: total
        integer          :: cnt 
        integer          :: i, j
        
        cnt = 0
        total = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if ((is_zero) .or. (matrix(i,j) > EPS_ZERO) .or. (matrix(i,j) < -EPS_ZERO))  then
                    cnt = cnt + 1
                    total = total + matrix(i, j)
                end if
            end do
        end do
        
        output = total / cnt
    
    end function stastics_average_value_matrix_with_zero
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    function stastics_average_value_vector_with_weightzero (vector, weight, is_zero)  result(output)
        
        real(KREAL), intent(in)    :: vector(:)
        real(KREAL), intent(in)    :: weight(:)
        logical, intent(in)        :: is_zero
        real(kind=8)  :: output
        
        real(kind=8)     :: total
        real(kind=8)     :: cnt 
        integer          :: i
        
        cnt = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if ((is_zero) .or. (vector(i) > EPS_ZERO))  then
                cnt = cnt + weight(i)
                total = total + vector(i)*weight(i)
            end if
        end do
        
        output = total / cnt
    
    end function stastics_average_value_vector_with_weightzero
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    function stastics_average_value_matrix_with_weightzero (matrix, weight, is_zero)  result(output)
        
        real(KREAL), intent(in)    :: matrix(:, :)
        real(KREAL), intent(in)    :: weight(:, :)
        logical, intent(in)        :: is_zero
        real(kind=8)  :: output
        
        real(kind=8)     :: total
        real(kind=8)     :: cnt 
        integer          :: i, j
        
        cnt = 0.0D0
        total = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if ((is_zero) .or. (matrix(i,j) > EPS_ZERO))  then
                    cnt = cnt + weight(i, j)
                    total = total + matrix(i, j)*weight(i, j)
                end if
            end do
        end do
        
        output = total / cnt
    
    end function stastics_average_value_matrix_with_weightzero
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! get variance of a vector with an optional mask
    !===============================================================================================
    function stastics_variance_vector (vector, mask)  result(output)
        
        real(KREAL), intent(in)        :: vector(:)
        logical, intent(in), optional  :: mask(:)
        real(kind=8)  :: output
        
        logical          :: is_true(SIZE(vector))
        real(kind=8)     :: total, square
        integer          :: cnt 
        integer          :: i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        cnt = 0
        total = 0.0D0
        square = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + 1
                total = total + vector(i)
                square = square + vector(i)**2
            end if
        end do
        
        output = (square - cnt * (total/cnt)**2) / (cnt - 1) 
    
    end function stastics_variance_vector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function stastics_variance_matrix (matrix, mask)  result(output)
        
        real(KREAL), intent(in)        :: matrix(:, :)
        logical, intent(in), optional  :: mask(:, :)
        real(kind=8)  :: output
        
        logical          :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        real(kind=8)     :: total, square
        integer          :: cnt 
        integer          :: i, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        cnt = 0
        total = 0.0D0
        square = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + 1
                    total = total + matrix(i, j)
                    square = square + matrix(i, j)**2
                end if
            end do
        end do
        
        output = (square - cnt * (total/cnt)**2) / (cnt - 1) 
    
    end function stastics_variance_matrix
    
    !$
    !===============================================================================================
    ! get standard variance of a vector with an optional mask
    !===============================================================================================
    function stastics_std_variance_vector (vector, mask)  result(output)
        
        real(KREAL), intent(in)        :: vector(:)
        logical, intent(in), optional  :: mask(:)
        real(kind=8)  :: output
        
        logical          :: is_true(SIZE(vector))
        real(kind=8)     :: total, square
        integer          :: cnt 
        integer          :: i
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        cnt = 0
        total = 0.0D0
        square = 0.0D0
        do i = 1, SIZE(vector, dim=1)
            if (is_true(i))  then
                cnt = cnt + 1
                total = total + vector(i)
                square = square + vector(i)**2
            end if
        end do
        
        output = (square - cnt * (total/cnt)**2) / (cnt - 1) 
        output = SQRT(output)
    
    end function stastics_std_variance_vector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function stastics_std_variance_matrix (matrix, mask)  result(output)
        
        real(KREAL), intent(in)        :: matrix(:, :)
        logical, intent(in), optional  :: mask(:, :)
        real(kind=8)  :: output
        
        logical          :: is_true(SIZE(matrix, dim=1), SIZE(matrix, dim=2))
        real(kind=8)     :: total, square
        integer          :: cnt 
        integer          :: i, j
        
        is_true = .TRUE.
        if (PRESENT(mask))  then
            is_true = mask
        end if
        
        cnt = 0
        total = 0.0D0
        square = 0.0D0
        do i = 1, SIZE(matrix, dim=1)
            do j = 1, SIZE(matrix, dim=2)
                if (is_true(i, j))  then
                    cnt = cnt + 1
                    total = total + matrix(i, j)
                    square = square + matrix(i, j)**2
                end if
            end do
        end do
        
        output = (square - cnt * (total/cnt)**2) / (cnt - 1) 
        output = SQRT(output)
    
    end function stastics_std_variance_matrix
    
end module stastics
