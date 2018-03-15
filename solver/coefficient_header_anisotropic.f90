!$
!===================================================================================================
!
!   prepare for anisotropic scatter treatment, generate expand coefficient 
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!                               
!   Public type lists:          AnisotropicSourceCoefficient
!
!===================================================================================================
module coefficient_header_anisotropic
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
        
    use quadrature_header,      only : QuadratureSet
        
    implicit none
    private
    public  :: AnisotropicSourceCoefficient
    
    ! --------------------------------------------------------------------------
    ! type for source expanse coefficient when anisotropic scatter
    type  :: AnisotropicSourceCoefficient
        real(KREAL), public, allocatable  ::  pl_zero(:, :)
        real(KREAL), public, allocatable  ::  pl_sin(:, :, :)
        real(KREAL), public, allocatable  ::  pl_cos(:, :, :)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_AnisotropicSourceCoefficient
        procedure, public  :: set => Set_AnisotropicSourceCoefficient
        procedure, public  :: clean =>  Free_AnisotropicSourceCoefficient
    end type AnisotropicSourceCoefficient
        
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_AnisotropicSourceCoefficient (this)
        
        class(AnisotropicSourceCoefficient), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%pl_zero(ns%state%scat_order,ns%deduce%direction), stat=i_allocate)
        allocate(this%pl_sin(ns%state%scat_order,ns%state%scat_order,ns%deduce%direction), stat=i_allocate)
        allocate(this%pl_cos(ns%state%scat_order,ns%state%scat_order,ns%deduce%direction), stat=i_allocate)
        
        ! set initial value
        this%pl_zero = REAL_ZERO
        this%pl_sin = REAL_ZERO
        this%pl_cos = REAL_ZERO

        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%pl_zero)
        this%memory = this%memory + REAL_BYTE * SIZE(this%pl_sin)
        this%memory = this%memory + REAL_BYTE * SIZE(this%pl_cos)
        
    end subroutine Allocate_AnisotropicSourceCoefficient
    
    !$
    !===============================================================================================
    ! finalizer for class of AnisotropicSourceCoefficient
    !===============================================================================================
    subroutine Free_AnisotropicSourceCoefficient (this)
        
        class(AnisotropicSourceCoefficient), intent(in out)  :: this
        
        if (allocated(this%pl_zero))            deallocate(this%pl_zero)
        if (allocated(this%pl_sin))             deallocate(this%pl_sin)
        if (allocated(this%pl_cos))             deallocate(this%pl_cos)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_AnisotropicSourceCoefficient
    
    !$
    !===============================================================================================
    ! generate the source expanse coefficient when anisotropic scatter, only once 
    !===============================================================================================
    subroutine Set_AnisotropicSourceCoefficient (this, quad)
        
        class(AnisotropicSourceCoefficient), intent(in out)  :: this
        type(QuadratureSet), intent(in)  :: quad
        
        ! local variables
        integer  :: i, l, k, n, m, is

        real(KDOUBLE)  :: azimuth(ns%deduce%direction)                                          ! azimuth
        real(KDOUBLE)  :: p_lm(ns%deduce%direction, ns%deduce%scat_xs, ns%deduce%scat_xs)       ! associated Legendre polynomial
        real(KDOUBLE)  :: factorial(ns%state%scat_order*2+1)                                    ! factorial times
        real(KDOUBLE)  :: e1, e2, root, pmID
        integer  :: i1, i2, ifm

        ! generate a dimension contains factorial times
        ifm = 2*ns%state%scat_order + 1
        factorial(1) = 1.0
        do i = 2, ifm
            factorial(i) = REAL(i-1) * factorial(i-1)
        end do
        
        ! get associated Legendre polynomial: p_lm(x,l,m), 'm' times 'l' order, based on recurrence formula
        do is = 1, ns%deduce%direction
            e1 = SQRT(1.0-quad%directions(is)%xmu(1)**2 - quad%directions(is)%xmu(2)**2)
            azimuth(is) = ATAN(e1 / quad%directions(is)%xmu(2)) 
            
            if (quad%directions(is)%xmu(2) < 0.0)  then
                azimuth(is) = azimuth(is) + PI
            end if
            if (quad%directions(is)%xmu(3) < 0.0)  then 
                azimuth(is) = -azimuth(is)
            end if
            
            p_lm(is,1,1) = 1.0
            p_lm(is,2,1) = quad%directions(is)%xmu(1)
            p_lm(is,2,2) = SQRT(1.0 - quad%directions(is)%xmu(1)**2)
            if (ns%state%scat_order <= 1)  then
                cycle
            end if
            
            root = 1.0 / p_lm(is,2,2)
            do n = 2, ns%state%scat_order
                e1 = 1.0 - 1.0/REAL(n, KDOUBLE)
                e2 = e1 + 1.0
                p_lm(is,n+1,1) = e2*quad%directions(is)%xmu(1)*p_lm(is,n,1) - e1*p_lm(is,n-1,1)
                do m = 1, n
                    i1 = n - m + 1
                    i2 = n + m - 1
                    e1 = i2
                    e2 = i1
                    p_lm(is,n+1,m+1) = root*(e1*p_lm(is,n,m) - e2*quad%directions(is)%xmu(1)*p_lm(is,n+1,m))
                end do
            end do
        end do 

        ! refer to p59 of the dissertation, just the coefficient
        do is = 1, ns%deduce%direction
            do l = 2, ns%deduce%scat_xs
                pmID = SQRT(2*l - 1.0)
                this%pl_zero(l-1, is) =  pmID * p_lm(is,l,1)      
                do k = 2, l
                    this%pl_cos(k-1, l-1, is) = pmID * SQRT(factorial(l-k+1)/factorial(l+k-1)) * p_lm(is,l,k) * COS((k-1)*azimuth(is))
                    this%pl_sin(k-1, l-1, is) = pmID * SQRT(factorial(l-k+1)/factorial(l+k-1)) * p_lm(is,l,k) * SIN((k-1)*azimuth(is))
                end do
            end do
        end do

    end subroutine Set_AnisotropicSourceCoefficient

end module coefficient_header_anisotropic
