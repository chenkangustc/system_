!$
!===================================================================================================
!
!   class for point kinetics solver
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          PKSolver
!                               PKParameter
!
!===================================================================================================
module pksolver_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use self_vode
    use lapack_interface
    use state_header,           only : TransientState
    
    implicit none
    private
    public  :: PKSolver, PKParameter

    ! --------------------------------------------------------------------------
    ! type for point kinetics solver
    type, private  :: density_vector_tp
        real(KREAL), allocatable   :: vector(:)
    contains
        procedure, public  :: alloc => Allocate_density_vector_tp
        procedure, public  :: clean => Free_density_vector_tp
        generic,   public  :: assignment(=) => Equal_density_vector_tp
        procedure          :: Equal_density_vector_tp
    end type density_vector_tp
    
    type  PKSolver
        real(KREAL), allocatable              :: partial_beta(:)
        real(KREAL), allocatable              :: partial_lambda(:)
        real(KREAL)                           :: beta
        real(KREAL)                           :: generation_time
        real(KREAL)                           :: rho
        real(KREAL)                           :: q
        real(KREAL)                           :: neutron
        real(KREAL), allocatable              :: precursor(:)
        
        integer                               :: dg         = 6
        integer                               :: n_richard  = 12
        real(KREAL)                           :: EPS_       = 1.0D-12
        character(len=MAX_WORD_LEN)           :: method     = 'VODE'
                                              
        real(KREAL)                           :: left
        real(KREAL)                           :: right
        real(KREAL)                           :: pace
        type(density_vector_tp)               :: last
        type(density_vector_tp)               :: current
        type(density_vector_tp), allocatable  :: Rmatrix(:, :)
    contains
        procedure, public  :: alloc => Alloc_PKSolver
        procedure, public  :: clean => Free_PKSolver
        procedure, public  :: advance => Time_advance_PKSolver
        procedure, private  :: rbfd => PKSolver_Perform_rbfd
        procedure, private  :: vode => PKSolver_Perform_vode
    end type PKSolver
    
    ! type for point kinetics parameter
    type  PKParameter
        real(KREAL), allocatable  :: partial_beta(:)                        ! group delayed neutron fraction
        real(KREAL), allocatable  :: partial_lambda(:)                      ! group decay constant
        real(KREAL)               :: beta                                   ! total delayed neutron fraction
        real(KREAL)               :: generation_time                        ! neutron generation time
        real(KREAL)               :: rho                                    ! current reactivity
        real(KREAL)               :: q                                      ! external source intensity
        real(KREAL)               :: neutron
        real(KREAL), allocatable  :: precursor(:)
    contains 
        procedure, public  :: alloc => Allocate_PKParameter
        procedure, public  :: clean => Free_PKParameter
        procedure, public  :: print => Print_PKParameter
    end type PKParameter

    ! --------------------------------------------------------------------------
    private  :: Allocate_density_vector_tp, Free_density_vector_tp, Equal_density_vector_tp
    private  :: Alloc_PKSolver, Free_PKSolver, Time_advance_PKSolver
    private  :: PKSolver_Perform_rbfd, PKSolver_Perform_vode
    private  :: Allocate_PKParameter, Free_PKParameter, Print_PKParameter
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_density_vector_tp (this, dg)
        
        class(density_vector_tp), intent(in out)  :: this
        integer, intent(in)                       :: dg
        
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        allocate(this%vector(dg+1), stat=i_allocate)
        this%vector = REAL_ZERO

    end subroutine Allocate_density_vector_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_density_vector_tp (this)
        
        class(density_vector_tp), intent(in out)  :: this
        
        if (allocated(this%vector))         deallocate(this%vector)
    
    end subroutine Free_density_vector_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Equal_density_vector_tp (left, right)
        
        class(density_vector_tp), intent(in out)  :: left
        type(density_vector_tp), intent(in)       :: right
        
        left%vector = right%vector
    
    end subroutine Equal_density_vector_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Alloc_PKSolver (this, nt, n_richard)
        
        class(PKSolver), intent(in out)   :: this
        type(TransientState), intent(in)  :: nt
        integer, intent(in), optional     :: n_richard
        
        integer  :: i, j
        integer  :: i_allocate
        
        this%dg = nt%state%dg
        if (PRESENT(n_richard))  then
            this%n_richard = n_richard
        end if
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%partial_beta(this%dg), stat=i_allocate)
        allocate(this%partial_lambda(this%dg), stat=i_allocate)
        allocate(this%precursor(this%dg), stat=i_allocate)
        
        call this%last%alloc (this%dg)
        call this%current%alloc (this%dg)
        allocate(this%Rmatrix(this%n_richard, this%n_richard), stat=i_allocate)
        do i = 1, SIZE(this%Rmatrix, dim=1)
            do j = 1, SIZE(this%Rmatrix, dim=2)
                call this%Rmatrix(i, j)%alloc (this%dg)
            end do
        end do
    
    end subroutine Alloc_PKSolver
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_PKSolver (this)
        
        class(PKSolver), intent(in out)  :: this
        integer  :: i, j
        
        if (allocated(this%partial_beta))           deallocate(this%partial_beta)
        if (allocated(this%partial_lambda))         deallocate(this%partial_lambda)
        if (allocated(this%precursor))              deallocate(this%precursor)
        
        call this%last%clean ()
        call this%current%clean ()
        if (allocated(this%Rmatrix))  then
            do i = 1, SIZE(this%Rmatrix, dim=1)
                do j = 1, SIZE(this%Rmatrix, dim=2)
                    call this%Rmatrix(i, j)%clean ()
                end do
            end do
        end if
    
    end subroutine Free_PKSolver
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Time_advance_PKSolver (this, pk, left, right)
        
        class(PKSolver), intent(in out)               :: this
        type(PKParameter), intent(in out)  :: pk
        real(KREAL), intent(in)                       :: left
        real(KREAL), intent(in)                       :: right
        
        integer  :: i, j
        
        this%partial_beta = pk%partial_beta
        this%partial_lambda = pk%partial_lambda
        this%beta = pk%beta
        this%generation_time = pk%generation_time
        this%rho = pk%rho
        this%q = pk%q
        this%neutron = pk%neutron
        this%precursor = pk%precursor
        
        this%left = left
        this%right = right
        this%pace = this%right - this%left
        this%last%vector(1) = this%neutron
        this%last%vector(2: ) = this%precursor
        
        select case (TRIM(this%method))
        case ('RBFD')
            call this%rbfd ()
        case ('VODE')
            call this%vode ()
        end select
        
        this%neutron = this%current%vector(1)
        this%precursor = this%current%vector(2: )
        
        pk%neutron = this%neutron
        pk%precursor = this%precursor
    
    end subroutine Time_advance_PKSolver
    
    !$
    !===============================================================================================
    ! solve point kinetics by RBFD
    !===============================================================================================
    subroutine PKSolver_Perform_rbfd (this)
        
        class(PKSolver), intent(in out)  :: this
        
        real(KREAL)  :: tmp(SIZE(this%last%vector))        
        real(KREAL)  :: length                                              ! step length per back euler
        real(KREAL)  :: abs_error
        integer  :: i, j, k
        
        ! the first point
        length = this%right - this%left
        call Back_Euler (this, this%last, this%left, this%right, length, tmp)
        this%Rmatrix(1, 1)%vector = tmp
        
        ! refinement
        do i = 2, this%n_richard
            length = length / 2.0D0
            call Back_Euler (this, this%last, this%left, this%right, length, tmp)
            this%Rmatrix(i, 1)%vector = tmp
            
            ! extrapolation
            do j = 2, i
                do k = 1, SIZE(this%Rmatrix(1, 1)%vector)
                    this%Rmatrix(i, j)%vector(k) = (2**(j-1) * this%Rmatrix(i, j-1)%vector(k) - this%Rmatrix(i-1,j-1)%vector(k)) / (2**(j-1) - 1)
                end do
            end do
            
            call Get_Richard_error (this%Rmatrix(i, i), this%Rmatrix(i-1, i-1), abs_error)
            if (abs_error < this%EPS_)  then
                this%current%vector = this%Rmatrix(i, i)%vector
                return
            end if
        end do
        
        ! update vector
        this%current%vector = this%Rmatrix(this%n_richard, this%n_richard)%vector
        
    end subroutine PKSolver_Perform_rbfd
    
    !$
    !===============================================================================================
    ! get max relative error of solve vector
    !===============================================================================================
    subroutine Get_Richard_error (first, second, abs_error)
        
        type(density_vector_tp), intent(in)  :: first
        type(density_vector_tp), intent(in)  :: second
        real(KREAL), intent(in out)          :: abs_error
        
        real(KREAL)   :: tmp
        integer       :: i
        
        abs_error = 1.0D0
        do i = 1, SIZE(first%vector)
            tmp = ABS(second%vector(i) - first%vector(i)) / ABS(first%vector(i))
            if (tmp <= abs_error)  then
                abs_error = tmp
            end if
        end do
    
    end subroutine Get_Richard_error
    
    !$
    !===============================================================================================
    ! this is the calculation kernel, different selection
    ! -- evoke DGESV from LAPACK to solve directly
    !===============================================================================================
    subroutine Back_Euler (pk, last, left, right, length, output)
        
        type(PKSolver), intent(in)           :: pk
        type(density_vector_tp), intent(in)  :: last
        real(KREAL), intent(in)              :: left
        real(KREAL), intent(in)              :: right
        real(KREAL), intent(in)              :: length                          ! step length per back euler
        real(KREAL), intent(in out)          :: output(:)
        
        real(KREAL)  :: lhs(SIZE(last%vector), SIZE(last%vector))
        real(KREAL)  :: rhs(SIZE(last%vector))
        
        real(KREAL)  :: old(SIZE(last%vector))
        real(KREAL)  :: new(SIZE(last%vector))
        real(KREAL)  :: tmp(SIZE(last%vector)) 
        real(KREAL)  :: this_time
        integer  :: scale
        integer  :: n_step
        integer  :: i, j, k
        integer  :: i_allocate
        
        ! ----------------------------------------------------------------------
        scale = SIZE(last%vector)
        n_step = NINT((right - left) / length)
        
        do i = 1, n_step
            this_time = left + i*length
            lhs = 0.0D0; rhs = 0.0D0
            tmp = 0.0D0
            
            if (i == 1)  then
                old = last%vector
            else
                old = new
            end if
            
            ! lhs matrix
            do k = 1, scale
                if (k == 1)  then
                    lhs(1, k) = (pk%rho - pk%beta) / pk%generation_time
                else 
                    lhs(1, k) = pk%partial_lambda(k-1)
                end if
            end do
            
            do k = 2, scale 
                lhs(k, 1) = pk%partial_beta(k-1) / pk%generation_time
                lhs(k, k) = - pk%partial_lambda(k-1)
            end do
            
            lhs = - length*lhs
            do k = 1, scale
                lhs(k, k) = 1.0 + lhs(k, k)
            end do
            
            ! rhs matrix
            rhs = old 
            do k = 1, scale 
                if (k == 1)  then
                    rhs(k) = rhs(k) + length*pk%q
                else 
                    rhs(k) = rhs(k)
                end if
            end do
            
            ! get solve for this equations
            call self_DGESV (lhs, rhs)
            new = rhs
            
        end do
        
        output = new
    
    end subroutine Back_Euler
    
    !$
    !===============================================================================================
    ! solve point kinetics by VODE
    ! communicate information by rpar & ipar 
    !===============================================================================================
    subroutine PKSolver_Perform_vode (this)
        
        class(PKSolver), intent(in out)  :: this
    
        ! ----------------------------------------------------------------------
        ! VODE parameter
        integer i, iopar, iopt, iout, ipar, istate, itask, itol, iwork,         &
            &   jsv, leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter,      &
            &   ml, mu, ncfn, neq, nerr, netf, nfe, nfea, nje, nlu, nni, nout,  &
            &   nqu, nst
        double precision atol, dtout, er, erm, ero, hu, rpar, rtol, rwork,      &
            &   t, tout, tout1, y
        dimension y(25), rwork(847), iwork(55), rpar(44)
        dimension rtol(1), atol(1), ipar(1)
        !-----------------------------------------------------------------------
    
        ! VODE solver method
        mf = 21
        
        ! pack into rpar & ipar
        ipar(1) = this%dg
        rpar(1:this%dg) = this%partial_beta
        rpar(this%dg+1:2*this%dg) = this%partial_lambda
        rpar(2*this%dg+1) = this%beta
        rpar(2*this%dg+2) = this%generation_time
        rpar(2*this%dg+3) = this%rho
        rpar(2*this%dg+4) = this%q
        
        itol = 1
        rtol = 1.0D-10
        atol = 0.0D0
        lrw = 847
        liw = 55
        iopt = 0
        
        neq = this%dg + 1
        t = this%left
        y(1:neq) = this%last%vector
        itask = 1
        istate = 1
        tout = this%right
        ero = 0.0D0
        
        call DVODE(VODE_function,neq,y,t,tout,itol,rtol,atol,itask,istate,      &
            &       iopt,rwork,lrw,iwork,liw,VODE_jacobi,mf,rpar,ipar)
        hu = rwork(11)
        nqu = iwork(14)
        er = ABS(y(1))/atol(1)
        ero = MAX(ero,er)
        
        this%current%vector = y(1:neq)
        
    end subroutine PKSolver_Perform_vode
    
    !$
    !===================================================================================================
    ! define function form
    !===================================================================================================
    subroutine VODE_function (neq, t, y, ydot, rpar, ipar)
    
        integer neq, ipar
        double precision t, y, ydot, rpar
        dimension y(neq), ydot(neq), rpar(*), ipar(*)
        
        real(8), allocatable  :: partial_beta(:)
        real(8), allocatable  :: partial_lambda(:)
        real(8)  :: beta, generation_time
        real(8)  :: rho, q
        integer  :: dg, k
        integer  :: i_allocate
        
        ! unpack from ipar & rpar
        dg = ipar(1)
        allocate(partial_beta(dg), stat=i_allocate)
        allocate(partial_lambda(dg), stat=i_allocate)
        partial_beta = rpar(1:dg)
        partial_lambda = rpar(dg+1 :2*dg)
        beta = rpar(2*dg+1)
        generation_time = rpar(2*dg+2) 
        rho = rpar(2*dg+3)
        q = rpar(2*dg+4)
        
        ydot(1) = (rho - beta) / generation_time * y(1) + q
        do k = 2, dg+1
            ydot(1) = ydot(1) + partial_lambda(k-1) * y(k)
        end do
        
        do k = 2, dg+1
            ydot(k) = partial_beta(k-1)/generation_time*y(1) - partial_lambda(k-1)*y(k)
        end do
        
        if (allocated(partial_beta))           deallocate(partial_beta)
        if (allocated(partial_lambda))         deallocate(partial_lambda)
        
    end subroutine VODE_function
          
    !$
    !===================================================================================================
    ! define Jacobi matrix
    !===================================================================================================
    subroutine VODE_jacobi (neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
        
        integer neq, ml, mu, nrowpd, ipar
        double precision t, y, pd, rpar
        dimension y(neq), pd(nrowpd,neq), rpar(*), ipar(*)
        
        real(8), allocatable  :: partial_beta(:)
        real(8), allocatable  :: partial_lambda(:)
        real(8)  :: beta, generation_time
        real(8)  :: rho, q
        integer  :: dg, k
        integer  :: i_allocate
        
        dg = ipar(1)
        allocate(partial_beta(dg), stat=i_allocate)
        allocate(partial_lambda(dg), stat=i_allocate)
        partial_beta = rpar(1:dg)
        partial_lambda = rpar(dg+1 :2*dg)
        beta = rpar(2*dg+1)
        generation_time = rpar(2*dg+2) 
        rho = rpar(2*dg+3)
        q = rpar(2*dg+4)
        
        pd = 0.0D0
        do k = 1, dg+1
            if (k == 1)  then
                pd(1,k) = (rho - beta) / generation_time
            else 
                pd(1,k) = partial_lambda(k-1)
            end if
        end do
        
        do k = 2, dg+1
            pd(k,1) = partial_beta(k-1) / generation_time
            pd(k,k) = - partial_lambda(k-1)
        end do
        
        if (allocated(partial_beta))           deallocate(partial_beta)
        if (allocated(partial_lambda))         deallocate(partial_lambda)
        
    end subroutine VODE_jacobi
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_PKParameter (this, nt)
    
        class(PKParameter), intent(in out)  :: this
        type(TransientState), intent(in)               :: nt
        integer  :: i_allocate
        
        ! check alloated status first
        call this%clean ()
        
        allocate(this%partial_beta(nt%state%dg), stat=i_allocate)
        allocate(this%partial_lambda(nt%state%dg), stat=i_allocate)
        allocate(this%precursor(nt%state%dg), stat=i_allocate)
        
        this%partial_beta   = REAL_ZERO
        this%partial_lambda = REAL_ZERO
        this%precursor      = REAL_ZERO
        
    end subroutine Allocate_PKParameter
    
    !$
    !===============================================================================================
    ! finalizer for class of PKParameter
    !===============================================================================================
    subroutine Free_PKParameter (this)
        
        class(PKParameter), intent(in out)  :: this
        
        if (allocated(this%partial_beta))       deallocate(this%partial_beta)
        if (allocated(this%partial_lambda))     deallocate(this%partial_lambda)
        if (allocated(this%precursor))          deallocate(this%precursor)
        
    end subroutine Free_PKParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_PKParameter (this, unit_)
        
        class(PKParameter), intent(in out)  :: this
        integer, intent(in)  :: unit_
        
        write(unit=unit_, fmt="(1x, A)")  '__________________________________________________________'
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'rho             :', this%rho
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'generation_time :', this%generation_time
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'beta            :', this%beta
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'partial_beta    :', this%partial_beta
        write(unit=unit_, fmt="(1x, A, *(TR3, ES13.6))")  'partial_lambda  :', this%partial_lambda
    
    end subroutine Print_PKParameter
    
end module pksolver_header
