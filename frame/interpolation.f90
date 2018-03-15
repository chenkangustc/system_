!$
!===================================================================================================
!
!   module for interpolation(Linear-Linear, Linear-Log, Log-Log, Log-Linear)
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    array_interpolation
!                               two_variables_interpolation
!
!   Public type lists:          No
!
!===================================================================================================
module interpolation
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,       only  : ErrorCollector
    
    implicit none 
    private
    public  :: array_interpolation, two_variables_interpolation
    
    interface  array_interpolation
        module procedure  vector_interpolation
        module procedure  matrix_interpolation
        module procedure  sets_interpolation
    end interface array_interpolation
    
    interface two_variables_interpolation
        module procedure  two_variables_interpolation_point
        module procedure  two_variables_interpolation_1D
        module procedure  two_variables_interpolation_2D
        module procedure  two_variables_interpolation_3D
    end interface two_variables_interpolation
    
    ! interpolation scheme
    integer, parameter  :: LIN_LIN = 1
    integer, parameter  :: LIN_LOG = 2
    integer, parameter  :: LOG_LIN = 3
    integer, parameter  :: LOG_LOG = 4
    
    type(ErrorCollector)  :: a_error
    
contains
    !$
    !===============================================================================================
    ! input is a vector, use index as key, the point may not be index
    !===============================================================================================
    function vector_interpolation (vector, point, scheme)  result(output)
        
        real(KREAL), intent(in)  :: vector(:)
        real(KREAL), intent(in)  :: point
        integer, intent(in)      :: scheme
        real(KREAL)    :: output
        
        ! local variables
        real(KREAL), allocatable  :: key(:)
        real(KREAL), allocatable  :: value(:)
        real(KREAL)  :: x0, x1
        real(KREAL)  :: y0, y1
        real(KREAL)  :: factor
        integer  :: scale, i
        integer  :: i_allocate
        
        scale  = SIZE(vector)
        output = HUGE(0.0E0)
        
        allocate(key(scale), stat=i_allocate)
        allocate(value(scale), stat=i_allocate)
        do i = 1, scale
            key(i) = i
            value(i) = vector(i)
        end do
        
        ! out of range
        if (point <= key(1) )  then
            output = value(1)
            return
        else if (point >= key(scale))  then
            output = value(scale)
            return
        end if
        
        ! deternime interval
        do i = 2, SIZE(key)
            if (key(i) >= point)  then
                x0 = key(i-1);   x1 = key(i)
                y0 = value(i-1); y1 = value(i)
                exit
            end if
        end do
        
        ! perform interpolation
        select case (scheme)
        case (LIN_LIN)
            factor = (point - x0)/(x1 - x0)
            output = (1.0D0 - factor)*y0 + factor*y1

        case (LIN_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = (1.0D0 - factor)*y0 + factor*y1
            
        case (LOG_LIN)
            factor = (point - x0)/(x1 - x0)
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
            
        case (LOG_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
        end select
        
        deallocate(key, stat=i_allocate)
        deallocate(value, stat=i_allocate)
        
    end function vector_interpolation
    
    !$
    !===============================================================================================
    ! input is a matixs(1:2, :), use matrix(1, :) as key
    !===============================================================================================
    function matrix_interpolation (matrix, point, scheme)  result(output)
        
        real(KREAL), intent(in)  :: matrix(:, :)
        real(KREAL), intent(in)  :: point
        integer, intent(in)      :: scheme
        real(KREAL)    :: output
        
        ! local variables
        real(KREAL), allocatable  :: key(:)
        real(KREAL), allocatable  :: value(:)
        real(KREAL)  :: x0, x1
        real(KREAL)  :: y0, y1
        real(KREAL)  :: factor
        integer  :: scale, i
        integer  :: i_allocate
        
        scale  = SIZE(matrix, dim=2)
        output = HUGE(0.0E0)
        
        allocate(key(scale), stat=i_allocate)
        allocate(value(scale), stat=i_allocate)
        key(:) = matrix(1, :)
        value(:) = matrix(2, :)
        
        ! out of range
        if (point <= key(1) )  then
            output = value(1)
            return
        else if (point >= key(scale))  then
            output = value(scale)
            return
        end if
        
        ! deternime interval
        do i = 2, SIZE(key)
            if (key(i) >= point)  then
                x0 = key(i-1);   x1 = key(i)
                y0 = value(i-1); y1 = value(i)
                exit
            end if
        end do
        
        ! perform interpolation
        select case (scheme)
        case (LIN_LIN)
            factor = (point - x0)/(x1 - x0)
            output = (1.0D0 - factor)*y0 + factor*y1
            
        case (LIN_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = (1.0D0 - factor)*y0 + factor*y1
            
        case (LOG_LIN)
            factor = (point - x0)/(x1 - x0)
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
            
        case (LOG_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
        end select

        deallocate(key, stat=i_allocate)
        deallocate(value, stat=i_allocate)
    
    end function matrix_interpolation
    
    !$
    !===============================================================================================
    ! input is a key & value, use key as key
    !===============================================================================================
    function sets_interpolation (input_key, input_value, point, scheme)  result(output)
        
        real(KREAL), intent(in)  :: input_key(:)
        real(KREAL), intent(in)  :: input_value(:)
        real(KREAL), intent(in)  :: point
        integer, intent(in)      :: scheme
        real(KREAL)    :: output
        
        ! local variables
        real(KREAL), allocatable  :: key(:)
        real(KREAL), allocatable  :: value(:)
        real(KREAL)  :: x0, x1
        real(KREAL)  :: y0, y1
        real(KREAL)  :: factor
        integer  :: scale, i
        integer  :: i_allocate
        
        scale  = SIZE(input_key)
        output = HUGE(0.0E0)
        
        allocate(key(scale), stat=i_allocate)
        allocate(value(scale), stat=i_allocate)
        key = input_key
        value = input_value
        
        if (SIZE(key) /= SIZE(value))  then
            call a_error%set (INFO_LIST_FRAMEWORK, 'bad input, array number do not equal')
            call a_error%print (OUTPUT_UNIT)
        end if
        
        ! out of range
        if (point <= key(1) )  then
            output = value(1)
            return
        else if (point >= key(scale))  then
            output = value(scale)
            return
        end if
        
        ! deternime interval
        do i = 2, SIZE(key)
            if (key(i) > point)  then
                x0 = key(i-1);   x1 = key(i)
                y0 = value(i-1); y1 = value(i)
                exit
            end if
        end do
        
        ! perform interpolation
        select case (scheme)
        case (LIN_LIN)
            factor = (point - x0)/(x1 - x0)
            output = (1.0D0 - factor)*y0 + factor*y1
            
        case (LIN_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = (1.0D0 - factor)*y0 + factor*y1
            
        case (LOG_LIN)
            factor = (point - x0)/(x1 - x0)
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
            
        case (LOG_LOG)
            factor = (LOG(point) - LOG(x0))/(LOG(x1) - LOG(x0))
            output = EXP((1.0D0-factor)*LOG(y0) + factor*LOG(y1))
            
        end select
        
        deallocate(key, stat=i_allocate)
        deallocate(value, stat=i_allocate)
    
    end function sets_interpolation
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine two_variables_interpolation_point (p, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        real(KREAL), intent(in out)  :: p
        real(KREAL), intent(in)      :: p00
        real(KREAL), intent(in)      :: p01
        real(KREAL), intent(in)      :: p10
        real(KREAL), intent(in)      :: p11
        real(KREAL), intent(in)      :: x_min
        real(KREAL), intent(in)      :: x
        real(KREAL), intent(in)      :: x_max
        real(KREAL), intent(in)      :: y_min
        real(KREAL), intent(in)      :: y
        real(KREAL), intent(in)      :: y_max
        
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        real(KREAL)  :: p0, p1
        integer, parameter  :: scheme = LIN_LIN
        
        !    (0, 1) ------ (1, 1)    <-  y_max
        !      |             |
        !      |             |
        !      |             |
        !    (0, 0) ------ (1, 0)    <-  y_min
        !      ^             ^
        !     x_min         x_max
        
        key(1) = x_min
        key(2) = x_max
        value(1) = p00
        value(2) = p01
        p0 = sets_interpolation (key, value, x, scheme)
        
        key(1) = x_min
        key(2) = x_max
        value(1) = p10
        value(2) = p11
        p1 = sets_interpolation (key, value, x, scheme)
        
        key(1) = y_min
        key(2) = y_max
        value(1) = p0
        value(2) = p1
        p = sets_interpolation (key, value, y, scheme)
    
    end subroutine two_variables_interpolation_point
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine two_variables_interpolation_1D (p, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        real(KREAL), intent(in out)  :: p(:)
        real(KREAL), intent(in)      :: p00(:)
        real(KREAL), intent(in)      :: p01(:)
        real(KREAL), intent(in)      :: p10(:)
        real(KREAL), intent(in)      :: p11(:)
        real(KREAL), intent(in)      :: x_min
        real(KREAL), intent(in)      :: x
        real(KREAL), intent(in)      :: x_max
        real(KREAL), intent(in)      :: y_min
        real(KREAL), intent(in)      :: y
        real(KREAL), intent(in)      :: y_max
        
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        real(KREAL)  :: p0, p1
        integer, parameter  :: scheme = LIN_LIN
        integer  :: i
        
        do i = 1, SIZE(p)
            key(1) = x_min
            key(2) = x_max
            value(1) = p00(i)
            value(2) = p01(i)
            p0 = sets_interpolation (key, value, x, scheme)
            
            key(1) = x_min
            key(2) = x_max
            value(1) = p10(i)
            value(2) = p11(i)
            p1 = sets_interpolation (key, value, x, scheme)
            
            key(1) = y_min
            key(2) = y_max
            value(1) = p0
            value(2) = p1
            p(i) = sets_interpolation (key, value, y, scheme)
        end do
    
    end subroutine two_variables_interpolation_1D
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine two_variables_interpolation_2D (p, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        real(KREAL), intent(in out)  :: p(:, :)
        real(KREAL), intent(in)      :: p00(:, :)
        real(KREAL), intent(in)      :: p01(:, :)
        real(KREAL), intent(in)      :: p10(:, :)
        real(KREAL), intent(in)      :: p11(:, :)
        real(KREAL), intent(in)      :: x_min
        real(KREAL), intent(in)      :: x
        real(KREAL), intent(in)      :: x_max
        real(KREAL), intent(in)      :: y_min
        real(KREAL), intent(in)      :: y
        real(KREAL), intent(in)      :: y_max
        
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        real(KREAL)  :: p0, p1
        integer, parameter  :: scheme = LIN_LIN
        integer  :: i, j
        
        do i = 1, SIZE(p, dim=1)
            do j = 1, SIZE(p, dim=2)
                key(1) = x_min
                key(2) = x_max
                value(1) = p00(i, j)
                value(2) = p01(i, j)
                p0 = sets_interpolation (key, value, x, scheme)
                
                key(1) = x_min
                key(2) = x_max
                value(1) = p10(i, j)
                value(2) = p11(i, j)
                p1 = sets_interpolation (key, value, x, scheme)
                
                key(1) = y_min
                key(2) = y_max
                value(1) = p0
                value(2) = p1
                p(i, j) = sets_interpolation (key, value, y, scheme)
            end do
        end do
    
    end subroutine two_variables_interpolation_2D
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine two_variables_interpolation_3D (p, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        real(KREAL), intent(in out)  :: p(:, :, :)
        real(KREAL), intent(in)      :: p00(:, :, :)
        real(KREAL), intent(in)      :: p01(:, :, :)
        real(KREAL), intent(in)      :: p10(:, :, :)
        real(KREAL), intent(in)      :: p11(:, :, :)
        real(KREAL), intent(in)      :: x_min
        real(KREAL), intent(in)      :: x
        real(KREAL), intent(in)      :: x_max
        real(KREAL), intent(in)      :: y_min
        real(KREAL), intent(in)      :: y
        real(KREAL), intent(in)      :: y_max
        
        real(KREAL)  :: key(2)
        real(KREAL)  :: value(2)
        real(KREAL)  :: p0, p1
        integer, parameter  :: scheme = LIN_LIN
        integer  :: i, j, k
        
        do i = 1, SIZE(p, dim=1)
            do j = 1, SIZE(p, dim=2)
                do k = 1, SIZE(p, dim=3)
                    key(1) = x_min
                    key(2) = x_max
                    value(1) = p00(i, j, k)
                    value(2) = p01(i, j, k)
                    p0 = sets_interpolation (key, value, x, scheme)
                    
                    key(1) = x_min
                    key(2) = x_max
                    value(1) = p10(i, j, k)
                    value(2) = p11(i, j, k)
                    p1 = sets_interpolation (key, value, x, scheme)
                    
                    key(1) = y_min
                    key(2) = y_max
                    value(1) = p0
                    value(2) = p1
                    p(i, j, k) = sets_interpolation (key, value, y, scheme)
                end do
            end do
        end do
    
    end subroutine two_variables_interpolation_3D
    
end module interpolation
