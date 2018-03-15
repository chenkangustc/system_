!$
!===================================================================================================
!
!   module for point kinetics method
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Pk_solver_RBDF
!                               Pk_solver_VODE
!                               Print_step_info
!
!   Public type lists:          No
!
!===================================================================================================
module process_pk
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use self_vode
    use lapack_interface
    use global
    
    use transient_header,           only : AmplitudeFunction
    use coefficient_pk_header
    
    implicit none 
    private
    public  :: Pk_solver_RBDF, Pk_solver_VODE, Print_step_info
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Pk_solver_RBDF (last_vector, ltime, ctime, current_vector)
        
        type(DensityVector), intent(in)     :: last_vector
        real(KREAL), intent(in)         :: ltime
        real(KREAL), intent(in)         :: ctime
        type(DensityVector), intent(in out) :: current_vector
        
        real(KREAL)   :: length                                             ! step length per back euler
        real(KREAL)   :: abs_error
        integer  :: i, j, k
        
        ! generation matrix
        do i = 1, pk_state%n_richard
            do j = 1, pk_state%n_richard
                richard%matrix(i,j)%vector = 0.0
            end do
        end do
        
        ! the first point
        length = ctime - ltime
        richard%matrix(1, 1)%vector = Back_Euler(last_vector, ltime, ctime, length)
        
        ! refinement
        do i = 2, pk_state%n_richard
            length = length / 2.0
            richard%matrix(i, 1)%vector = Back_Euler(last_vector, ltime, ctime, length)
            
            ! extrapolation
            do j = 2, i
                do k = 1, SIZE(richard%matrix(1, 1)%vector)
                    richard%matrix(i, j)%vector(k) = (2**(j-1) * richard%matrix(i, j-1)%vector(k) - richard%matrix(i-1,j-1)%vector(k)) / (2**(j-1) - 1)
                end do
            end do
            
            call Richard_error (richard%matrix(i, i), richard%matrix(i-1, i-1), abs_error)
            if (abs_error < pk_state%tolerance)  then
                current_vector%vector = richard%matrix(i, i)%vector
                return
            end if
        end do
        
        ! update vector
        current_vector%vector = richard%matrix(pk_state%n_richard, pk_state%n_richard)%vector
        
    end subroutine Pk_solver_RBDF

    !$
    !===============================================================================================
    ! get max relative error of solve vector
    !===============================================================================================
    subroutine Richard_error (first, second, abs_error)
        
        type(DensityVector), intent(in)  :: first
        type(DensityVector), intent(in)  :: second
        real(KREAL), intent(in out)  :: abs_error
        
        real(KREAL)   :: tmp_error
        integer  :: i
        
        abs_error = 1.0
        do i = 1, SIZE(first%vector)
            tmp_error = ABS(second%vector(i) - first%vector(i)) / ABS(first%vector(i))
            
            if (tmp_error <= abs_error)  then
                abs_error = tmp_error
            end if
        end do
    
    end subroutine Richard_error
    
    !$
    !===============================================================================================
    ! this is the calculation kernel, different selection
    ! -- ivoke DGESV from LAPACK to solve directly
    ! -- ivoke DGETRI from LAPCK to inverse matrix, then get result
    ! -- perform analytical inverse, then get result
    ! -- perform Pickard iteration to get result, has not been finish
    !===============================================================================================
    function Back_Euler (last_vector, ltime, ctime, length)  result(output)
        
        type(DensityVector), intent(in)  :: last_vector
        real(KREAL), intent(in)    :: ltime
        real(KREAL), intent(in)    :: ctime
        real(KREAL), intent(in)    :: length                                ! step length per back euler
        real(KREAL)   :: output(nt%state%dg+1)
        
        integer  :: scale
        integer  :: n_step
        integer  :: i, j, k
        integer  :: i_allocate
        real(KREAL)  :: lhs(nt%state%dg+1, nt%state%dg+1)
        real(KREAL)  :: rhs(nt%state%dg+1)
        
        real(KREAL)  :: old(nt%state%dg+1)
        real(KREAL)  :: new(nt%state%dg+1)
        real(KREAL)  :: tmp(nt%state%dg+1) 
        
        real(KREAL)  :: rho, source
        real(KREAL)  :: this_time
        
        character(len=MAX_WORD_LEN)  :: RBDF_method
        
        ! ----------------------------------------------------------------------
        scale = SIZE(last_vector%vector)
        n_step = NINT((ctime - ltime) / length)
        
        do i = 1, n_step
            
            this_time = ltime + i*length
            lhs = 0.0
            rhs = 0.0
            tmp = 0.0
            
            if (i == 1)  then
                old = last_vector%vector
            else
                old = new
            end if
            
            ! ------------------------------------------------------------------
            ! generation equations coefficient

            ! NOTE: -- at every coarse point
            ! get current rho & source
            rho = pk_parameter%rho
            source = pk_parameter%source
            
            ! lhs matrix
            do k = 1, scale
                if (k == 1)  then
                    lhs(1,k) = (rho - pk_parameter%beta) / pk_parameter%generation_time
                else 
                    lhs(1,k) = pk_parameter%partial_lambda(k-1)
                end if
            end do
            
            do k = 2, scale 
                lhs(k,1) = pk_parameter%partial_beta(k-1) / pk_parameter%generation_time
                lhs(k,k) = - pk_parameter%partial_lambda(k-1)
            end do
            
            lhs = - length*lhs
            do k = 1, scale
                lhs(k,k) = 1.0 + lhs(k,k)
            end do
            
            ! rhs matrix
            rhs = old 
            do k = 1, scale 
                if (k == 1)  then
                    rhs(k) = rhs(k) + length*source
                else 
                    rhs(k) = rhs(k)
                end if
            end do
            
            ! ------------------------------------------------------------------
            ! get solve for this equations
            RBDF_method = 'DIRECT'
            
            select case (RBDF_method)
            case ('DIRECT')
                call Direct_solve(lhs, rhs, length, new)
                
            case('LU')
                call LU_inverse(lhs, rhs, length, new)
                
            case('ANALYTICAL')
                call Analytical_inverse(lhs, rhs, length, rho, new)
            end select
            
        end do
        
        ! set return 
        output = new
    
    end function Back_Euler
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! solve equations directly
    !===============================================================================================
    subroutine Direct_solve (lhs, rhs, length, solve)
        
        real(KREAL), intent(in out)  :: lhs(:, :)
        real(KREAL), intent(in out)  :: rhs(:)
        real(KREAL), intent(in)  :: length
        real(KREAL), intent(in out)  :: solve(:)
        
        call self_DGESV (lhs, rhs)
        
        solve = rhs
    
    end subroutine Direct_solve
    
    !$
    !===============================================================================================
    ! inverse matrix by LU decomposition
    !===============================================================================================
    subroutine LU_inverse (lhs, rhs, length, solve)
        
        real(KREAL), intent(in out)  :: lhs(:, :)
        real(KREAL), intent(in out)  :: rhs(:)
        real(KREAL), intent(in)  :: length
        real(KREAL), intent(in out)  :: solve(:)
        
        integer  :: j, k
        integer  :: scale
        
        ! perform LU decomposition first
        call self_DGETRF (lhs)
        
        ! then inverse the orignal matrix
        call self_DGETRI (lhs)
        
        ! get result
        solve = 0.0
        do j = 1, scale
            do k = 1, scale
                solve(j) = solve(j) + lhs(j,k) * rhs(k)
            end do
        end do
        
    end subroutine LU_inverse 
    
    !$
    !===============================================================================================
    ! inverse matrix by analytical method
    !===============================================================================================
    subroutine Analytical_inverse (lhs, rhs, length, rho, solve)
        
        real(KREAL), intent(in out)  :: lhs(:, :)
        real(KREAL), intent(in out)  :: rhs(:)
        real(KREAL), intent(in)  :: length
        real(KREAL), intent(in)  :: rho
        real(KREAL), intent(in out)  :: solve(:)
    
        real(KREAL)  :: row(nt%state%dg+1)
        real(KREAL)  :: column(nt%state%dg+1)
        real(KREAL)  :: tmp(nt%state%dg+1, nt%state%dg+1)
        
        real(KREAL)  :: sigma
        integer  :: i, j, k
        integer  :: scale
        
        ! ----------------------------------------------------------------------
        row = 0.0
        column = 0.0
        tmp = 0.0
        sigma = 0.0
        solve = 0.0
        
        scale = SIZE(rhs)
            
        do k = 1, nt%state%dg
            sigma = sigma + 1.0 / (1.0+length*pk_parameter%partial_lambda(k)) * (length*pk_parameter%partial_beta(k)/pk_parameter%generation_time)
        end do
        sigma = sigma + 1.0 - length*(rho/pk_parameter%generation_time)
        
        ! generation tmp
        do k = 1, scale
            if (k == 1)  then
                tmp(k,k) = 0.0
            else
                tmp(k,k) = 1.0 / (1.0+length*pk_parameter%partial_lambda(k-1))
            end if
        end do
        
        do k = 1, scale
            if (k == 1)  then
                row(k) = 1.0
                column(k) = 1.0
            else
                row(k) = 1.0 / (1.0+length*pk_parameter%partial_lambda(k-1)) * (length*pk_parameter%partial_beta(k-1)/pk_parameter%generation_time)
                column(K) = 1.0 / (1.0+length*pk_parameter%partial_lambda(k-1)) * (length*pk_parameter%partial_lambda(k-1))
            end if
        end do
        
        do j = 1, scale
            do k = 1, scale
                lhs(j,k) = tmp(j,k) + (1.0/sigma) * row(j) * column(K)
            end do
        end do
        
        ! get result
        do j = 1, scale
            do k = 1, scale
                solve(j) = solve(j) + lhs(j,k) * rhs(k)
            end do
        end do
        
    end subroutine Analytical_inverse
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! solve point kinetics by VODE
    !===============================================================================================
    subroutine Pk_solver_VODE (last_vector, ltime, ctime, current_vector)
        
        type(DensityVector), intent(in)     :: last_vector
        real(KREAL), intent(in)         :: ltime
        real(KREAL), intent(in)         :: ctime
        type(DensityVector), intent(in out) :: current_vector
    
        ! ----------------------------------------------------------------------
        ! VODE parameter
        integer i, iopar, iopt, iout, ipar, istate, itask, itol, iwork,         &
            &   jsv, leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter,      &
            &   ml, mu, ncfn, neq, nerr, netf, nfe, nfea, nje, nlu, nni, nout,  &
            &   nqu, nst
        double precision atol, dtout, er, erm, ero, hu, rpar, rtol, rwork,      &
            &   t, tout, tout1, y
        dimension y(25), rwork(847), iwork(55), rpar(10)
        dimension rtol(1), atol(1), ipar(1)
    
        !-----------------------------------------------------------------------
        ! VODE solver
        mf = 21
        rpar(1) = pk_parameter%rho
        rpar(2) = pk_parameter%source
        
        itol = 1
        rtol = 1.0d-10
        atol = 0.0d0
        lrw = 847
        liw = 55
        iopt = 0
        
        neq = nt%state%dg + 1
        t = ltime
        y(1:neq) = last_vector%vector
        itask = 1
        istate = 1
        tout = ctime
        ero = 0.0d0
        
        call DVODE(VODE_function,neq,y,t,tout,itol,rtol,atol,itask,istate,      &
            &       iopt,rwork,lrw,iwork,liw,VODE_jacobi,mf,rpar,ipar)
        hu = rwork(11)
        nqu = iwork(14)
        er = ABS(y(1))/atol(1)
        ero = MAX(ero,er)
        
        current_vector%vector = y(1:neq)
    
    end subroutine Pk_solver_VODE
    
    !$
    !===================================================================================================
    ! define function form
    !===================================================================================================
    subroutine VODE_function (neq, t, y, ydot, rpar, ipar)
    
        use constants
        use global_state
        use, intrinsic  :: ISO_FORTRAN_ENV
        
        use global
    
        integer neq, ipar
        double precision t, y, ydot, rpar
        dimension y(neq), ydot(neq), rpar(*), ipar(*)
        
        integer  :: k
        
        ! --------------------------------------------------------------------------
        ! rpar(1)--pk_parameter%rho
        ! rpar(2)--pk_parameter%source
        ydot(1) = (rpar(1) - pk_parameter%beta) / pk_parameter%generation_time * y(1) + rpar(2)
        do k = 2, nt%state%dg+1
            ydot(1) = ydot(1) + pk_parameter%partial_lambda(k-1) * y(k)
        end do
        
        do k = 2, nt%state%dg+1
            ydot(k) = pk_parameter%partial_beta(k-1)/pk_parameter%generation_time*y(1) - pk_parameter%partial_lambda(k-1)*y(k)
        end do
        
    end subroutine VODE_function
          
    !$
    !===================================================================================================
    ! define Jacobi matrix
    !===================================================================================================
    subroutine VODE_jacobi (neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    
        use constants
        use global_state
        use, intrinsic  :: ISO_FORTRAN_ENV
        
        use global
        
        integer neq, ml, mu, nrowpd, ipar
        double precision t, y, pd, rpar
        dimension y(neq), pd(nrowpd,neq), rpar(*), ipar(*)
        
        integer  :: k
        
        ! --------------------------------------------------------------------------
        ! rpar(1)--pk_parameter%rho
        ! rpar(2)--pk_parameter%source
        pd = 0.0
        
        do k = 1, nt%state%dg+1
            if (k == 1)  then
                pd(1,k) = (rpar(1) - pk_parameter%beta) / pk_parameter%generation_time
            else 
                pd(1,k) = pk_parameter%partial_lambda(k-1)
            end if
        end do
        
        do k = 2, nt%state%dg+1
            pd(k,1) = pk_parameter%partial_beta(k-1) / pk_parameter%generation_time
            pd(k,k) = - pk_parameter%partial_lambda(k-1)
        end do
        
    end subroutine VODE_jacobi
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! Print the neutron number and the presursor concentration at every time point. 
    !===============================================================================================
    subroutine Print_step_info (amplitude, tidx, ctime, unit_)
        
        type(AmplitudeFunction), intent(in)  :: amplitude
        integer, intent(in) :: tidx
        real(KREAL), intent(in)  :: ctime
        integer, intent(in) :: unit_
        
        write(unit=unit_, fmt="(1x, I5, 8x, ES11.4, 5x, E13.6, 4x, *(E13.6, 3x))")  tidx, ctime,     &
            &  amplitude%flux, amplitude%precursor(:)
        
    end subroutine Print_step_info
    
    !$
    !===============================================================================================
    ! Print the neutron number and the presursor concentration at every output point. 
    ! see the input parameter of print_frequency.
    !===============================================================================================
    subroutine Print_matrix (matrix, tidx, unit_)
        
        type(DensityVector), intent(in), allocatable     :: matrix(:, :)
        integer, intent(in)  :: tidx
        integer, intent(in)  :: unit_
        
        ! local variables
        integer     :: scale 
        integer     :: i, j
        real(KREAL), allocatable  :: tmp(:)
        
        scale = SIZE(matrix, dim=1)
        allocate(tmp(scale))
        
        write(unit=unit_, fmt="(1x, '_________________________________________', &
            &  '__________________________________________')")
        write(unit=unit_, fmt="(1x, A, 2x, I5)" ) 'The time index is:', tidx
        
        ! only output the neutron density
        do i = 1, scale 
            do j = 1, scale
                tmp(j) = matrix(i,j)%vector(1)
            end do
            write(unit=unit_, fmt="(1x, *(E17.10, 3x))")  tmp(:)
        end do 
        
        write(unit=unit_, fmt="(1x, 2/)" )
        
        deallocate(tmp)
    
    end subroutine Print_matrix
    
end module process_pk
