!$
!===================================================================================================
!
!   module for lapack interface, the capital name is the same as Lapack lib.
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    self_DGESV
!                               self_DGTSV
!                               self_DGETRF
!                               self_DGETRI
!
!   Public type lists:          No
!
!===================================================================================================
module lapack_interface
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: self_DGESV, self_DGTSV, self_DGETRF, self_DGETRI
    
    interface  self_DGESV
        module procedure  self_DGESV_vector 
        module procedure  self_DGESV_matrix 
    end interface self_DGESV
    
    interface self_DGTSV
        module procedure  self_DGTSV_vector 
        module procedure  self_DGTSV_matrix 
    end interface self_DGTSV
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGESV_vector  (lhs, rhs, error_info)
        
        real(KREAL), intent(in)       :: lhs(:, :)
        real(KREAL), intent(in out)   :: rhs(:)
        integer, intent(in out), optional :: error_info
        
        ! local variables
        integer  :: info
        integer  :: n, nrhs
        integer  :: lda, ldb
        integer, allocatable  :: ipiv(:)
        double precision, allocatable  :: a(:, :)
        double precision, allocatable  :: b(:)
        
        n    = SIZE(rhs)
        nrhs = 1
        lda  = n
        ldb  = n
        
        allocate(ipiv(n))
        allocate(a(lda, n))
        allocate(b(ldb))
        
        a = lhs
        b = rhs
        
        call DGESV (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
        
        if (info == 0)  then
            rhs = b
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(ipiv)
        deallocate(a)
        deallocate(b)
    
    end subroutine self_DGESV_vector 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGESV_matrix  (lhs, rhs, error_info)
        
        real(KREAL), intent(in)       :: lhs(:, :)
        real(KREAL), intent(in out)   :: rhs(:, :)
        integer, intent(in out), optional :: error_info
        
        ! local variables
        integer  :: info
        integer  :: n, nrhs
        integer  :: lda, ldb
        integer, allocatable  :: ipiv(:)
        double precision, allocatable  :: a(:, :)
        double precision, allocatable  :: b(:, :)
        
        n    = SIZE(rhs, dim=1)
        nrhs = SIZE(rhs, dim=2)
        lda  = n
        ldb  = n
        
        allocate(ipiv(n))
        allocate(a(lda, n))
        allocate(b(ldb, nrhs))
        
        a = lhs
        b = rhs
        
        call DGESV (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
        
        if (info == 0)  then
            rhs = b
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(ipiv)
        deallocate(a)
        deallocate(b)
    
    end subroutine self_DGESV_matrix 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGTSV_vector  (lhs, rhs, error_info)
        
        real(KREAL), intent(in)       :: lhs(:, :)
        real(KREAL), intent(in out)   :: rhs(:)
        integer, intent(in out), optional :: error_info
        
        ! local variables
        integer  :: i
        
        integer  :: info
        integer  :: n, nrhs
        integer  :: ldb
        double precision, allocatable  :: dl(:)
        double precision, allocatable  :: d(:)
        double precision, allocatable  :: du(:)
        double precision, allocatable  :: b(:)
        
        n    = SIZE(rhs)
        nrhs = 1
        ldb  = n
        
        allocate(dl(n-1))
        allocate(d(n))
        allocate(du(n-1))
        allocate(b(ldb))
        
        do i = 1, n
            d(i) = lhs(i, i)
        end do
        do i = 1, n-1
            dl(i) = lhs(i+1, i)
        end do
        do i = 1, n-1
            du(i) = lhs(i, i+1)
        end do
        
        b = rhs
        
        call DGTSV (N, NRHS, DL, D, DU, B, LDB, INFO)
        
        if (info == 0)  then
            rhs = b
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(dl)
        deallocate(d)
        deallocate(du)
        deallocate(b)
    
    end subroutine self_DGTSV_vector 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGTSV_matrix  (lhs, rhs, error_info)
        
        real(KREAL), intent(in)       :: lhs(:, :)
        real(KREAL), intent(in out)   :: rhs(:, :)
        integer, intent(in out), optional :: error_info
        
        ! local variables
        integer  :: i
        
        integer  :: info
        integer  :: n, nrhs
        integer  :: ldb
        double precision, allocatable  :: dl(:)
        double precision, allocatable  :: d(:)
        double precision, allocatable  :: du(:)
        double precision, allocatable  :: b(:, :)
        
        n    = SIZE(rhs, dim=1)
        nrhs = SIZE(rhs, dim=2)
        ldb  = n
        
        allocate(dl(n-1))
        allocate(d(n))
        allocate(du(n-1))
        allocate(b(ldb, nrhs))
        
        do i = 1, n
            d(i) = lhs(i, i)
        end do
        do i = 1, n-1
            dl(i) = lhs(i+1, i)
        end do
        do i = 1, n-1
            du(i) = lhs(i, i+1)
        end do
        
        b = rhs
        
        call DGTSV (N, NRHS, DL, D, DU, B, LDB, INFO)
        
        if (info == 0)  then
            rhs = b
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(dl)
        deallocate(d)
        deallocate(du)
        deallocate(b)
    
    end subroutine self_DGTSV_matrix 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGETRF  (matrix, error_info)
        
        real(KREAL), intent(in out)    :: matrix(:, :)
        integer, intent(in out), optional  :: error_info
        
        ! local variables
        integer  :: info
        integer  :: lda, lwork
        integer  :: n, m
        integer, allocatable  :: ipiv(:)
        double precision, allocatable  :: a(:, :)
        
        m = SIZE(matrix, dim=1)
        n = SIZE(matrix, dim=2)
        lda   = m
        
        allocate(ipiv(MIN(m,n)))
        allocate(a(lda, n))
        
        a = matrix
    
        call DGETRF( M, N, A, LDA, IPIV, INFO )
        
        if (info == 0)  then
            matrix = a
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(ipiv)
        deallocate(a)
    
    end subroutine self_DGETRF 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine self_DGETRI  (matrix, error_info)

        real(KREAL), intent(in out)    :: matrix(:, :)
        integer, intent(in out), optional  :: error_info
        
        ! local variables
        integer  :: info
        integer  :: lda, lwork
        integer  :: n
        integer, allocatable  :: ipiv(:)
        double precision, allocatable  :: a(:, :)
        double precision, allocatable  :: work(:)
        
        n = SIZE(matrix, dim=1)
        lda   = n
        lwork = n
        
        allocate(ipiv(n))
        allocate(a(lda, n))
        allocate(work(ABS(lwork)))
        
        a = matrix
    
        call DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        
        if (info == 0)  then
            matrix = a
        end if
        
        if (PRESENT(error_info))  then
            error_info = info
        end if
        
        deallocate(ipiv)
        deallocate(a)
        deallocate(work)
        
    end subroutine self_DGETRI 
    
end module lapack_interface
