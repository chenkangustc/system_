!===================================================================================================
!
!   Module for some general special function subroutine or function
!   NOTE: for associated Legendre function, we use Ferrer definition (another one is Hobson)
!         for spherical harmonic funcition, we do not use "Condon Shortley" coefficient
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Legendre
!                               Legendre_associated
!                               Spherical_harmonic
!
!   Public type lists:          No
!
!===================================================================================================
module special_function

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,       only : ErrorCollector
    
    implicit none
    private
    public  :: Legendre, Legendre_associated, Spherical_harmonic
    
    type(ErrorCollector)  :: a_error
    logical, parameter    :: is_Ferrer            = .TRUE.                      ! for associated Legendre
    logical, parameter    :: is_Condon_Shortley   = .FALSE.                     ! for spherical harmonic 
    
contains
    !$
    !===============================================================================================
    ! l-th order Legendre polynomial Pl(x)
    !       -- l is integer >=0
    !       -- x is real [-1, 1]
    !===============================================================================================
    function Legendre (l, x)  result(output)
    
        integer, intent(in)          :: l
        real(KREAL), intent(in)  :: x
        real(KREAL)              :: output
        
        ! local variables
        integer :: il 
        real(KDOUBLE)  :: low_order
        real(KDOUBLE)  :: low_low_order

        ! whether argument is out of range
        if (ABS(x) > 1) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The second argument for Legendre is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        else if (l < 0) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The first argument for Legendre is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        end if
        
        ! recurrence formula for Legendre polynomial
        if (l == 0) then
            output = 1.0
        else if (l == 1) then
            output = x
        else
            low_low_order   = 1.0
            low_order       = x
            do il = 2, l
                output   = ((2*il-1)*x*low_order - (il-1)*low_low_order) / REAL(il, KDOUBLE)
                low_low_order = low_order
                low_order = output
            end do
        end if
            
    end function Legendre
    
    !$
    !===============================================================================================
    ! l-th order associated Legendre polynomial Pl,m(x)
    !       -- l is integer >=0
    !       -- m is intger [-l, l]
    !       -- x is real [-1, 1]
    ! refer to function plgndr_s @Numerical Recipes
    !===============================================================================================
    function Legendre_associated (l, m, x)  result(output)
    
        integer, intent(in)          :: l
        integer, intent(in)          :: m
        real(KREAL), intent(in)  :: x
        real(KREAL)              :: output
        
        ! local variables
        integer  :: il, im
        integer  :: new_m
        real(KDOUBLE)  :: coefficient
        real(KDOUBLE)  :: new_coefficient
        real(KDOUBLE)  :: low_order
        real(KDOUBLE)  :: low_low_order
        real(KDOUBLE)  :: numerator
        real(KDOUBLE)  :: denominator
        
        ! whether argument is out of range
        if (ABS(x) > 1) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The third argument for Legendre_associated is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        else if (l < 0) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The first argument for Legendre_associated is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        else if (m < -l .or. l < m) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The second argument for Legendre_associated is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        end if
        
        output = 1.0
        coefficient = 1.0
        if (m == 0) then
            output = Legendre(l, x)
            
        else if (m > 0)  then
            do im = 1, 2*m, 2
                coefficient = coefficient * REAL(im, KDOUBLE)
            end do
            low_low_order = coefficient * ((SQRT((1.0-x)*(1.0+x))) **m)
            
            if (.NOT. is_Ferrer)  then
                if (MOD(m,2) == 1) low_low_order = -low_low_order
            end if
        
            if (l == m) then
                output = low_low_order
            else 
                low_order = low_low_order*x*(2*m+1)
                if (l == m+1) then
                    output = low_order
                else 
                    do il = m+2, l
                        output = ((2*il-1)*x*low_order - (il-1+m)*low_low_order) / (il-m)
                        low_low_order = low_order
                        low_order = output
                    end do
                end if
            end if
        
        else
            new_m = -m
            do im = 1, 2*new_m, 2
                coefficient = coefficient * REAL(im, KDOUBLE)
            end do
            low_low_order = coefficient * ((SQRT((1.0-x)*(1.0+x))) **new_m)
            
            if (.NOT. is_Ferrer)  then
                if (MOD(m,2) == 1) low_low_order = -low_low_order
            end if
        
            if (l == new_m) then
                output = low_low_order
            else 
                low_order = low_low_order*x*(2*new_m+1)
                if (l == new_m+1) then
                    output = low_order
                else 
                    do il = new_m+2, l
                        output = ((2*il-1)*x*low_order - (il-1+new_m)*low_low_order) / (il-new_m)
                        low_low_order = low_order
                        low_order = output
                    end do
                end if
            end if
            
            new_coefficient = 1.0
            numerator = 1.0
            denominator = 1.0
            do im = 1, l-m
                numerator = numerator * im
            end do
            do im = 1, l+m
                denominator = denominator * im
            end do
            new_coefficient = new_coefficient * numerator / denominator
            if (MOD(new_m,2) == 1) new_coefficient = -new_coefficient
            
            output = new_coefficient * output
        end if
            
    end function Legendre_associated
    
    !$
    !===============================================================================================
    ! l-th order spherical harmonic Rl,m(theta, psi), only the real part
    !       -- l is integer >=0
    !       -- m is intger [-l, l]
    !       -- theta is polar angle in radians
    !       -- psi is azimuthal angle in radians
    !===============================================================================================
    function Spherical_harmonic (l, m, theta, psi)  result(output)
    
        integer, intent(in)          :: l
        integer, intent(in)          :: m
        real(KREAL), intent(in)  :: theta
        real(KREAL), intent(in)  :: psi
        real(KREAL)              :: output
        
        ! local variables
        integer  :: im
        real(KDOUBLE)  :: numerator
        real(KDOUBLE)  :: denominator
        real(KDOUBLE)    :: coefficient
        
        ! whether the argument is out of range
        if (l < 0) then 
            call a_error%set (INFO_LIST_FRAMEWORK, 'The first argument for Spherical harmonic is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        else if (m < -l .or. l < m) then
            call a_error%set (INFO_LIST_FRAMEWORK, 'The second argument for Spherical harmonic is illegal')
            call a_error%print (OUTPUT_UNIT)
            return
        end if
        
        ! the relationship of spherical harmonic and assciated Legendre function
        numerator = 1.0
        denominator = 1.0
        do im = 1, l-m
            numerator = numerator * im
        end do
        do im = 1, l+m
            denominator = denominator * im
        end do
        
        coefficient = (2*l+1) / (4*PI)
        
        output = SQRT(coefficient*numerator/denominator) * Legendre_associated(l,m,COS(theta)) * COS(m*psi)
        
        ! their is a coefficient
        if (is_Condon_Shortley )  then
            output = output * (-1)**m
        end if
    
    end function Spherical_harmonic

end module special_function
