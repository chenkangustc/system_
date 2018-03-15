!$
!===================================================================================================
!
!   self version: lagrange interpolant to data defined on a ND grid 
!---------------------------------------------------------------------------------------------------
!   public subroutine lists:    No
!
!   public type lists:          No
!
!===================================================================================================
module self_lagrange_interp
    
    implicit none
    public 
    
contains
    !$
    !===============================================================================================
    !  CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
    !  Discussion:
    !    The abscissas are numbered from left to right. The rule is defined on [-1,1].
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) N, the order. 1 <= N.
    !    Output, real ( kind = 8 ) POINTS(N), the abscissas.
    !===============================================================================================
    subroutine cc_compute_points (n, points)

        integer ( kind = 4 ) n
        real ( kind = 8 ) points(n)
        integer ( kind = 4 ) i
        real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
        
        if (n < 1)  then
            write (*, '(a)') ' '
            write (*, '(a)') 'CC_COMPUTE_POINTS - Fatal error!'
            write (*, '(a,i8)') '  Illegal value of N = ', n
            stop
        else if (n == 1)  then
            points(1) = 0.0D+00
        else
            do i = 1, n
                points(i) = cos(real(n-i, kind=8) * pi / real(n-1, kind=8))
            end do
            
            points(1) = -1.0D+00
            if (mod(n, 2) == 1) then
                points((n+1)/2) = 0.0D+00
            end if
            points(n) = +1.0D+00
        end if
        
    end subroutine cc_compute_points
    
    !$
    !===============================================================================================
    !  LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) ND, the number of data points.
    !    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
    !    Input, integer ( kind = 4 ) NI, the number of evaluation points.
    !    Input, real ( kind = 8 ) XI(NI), the evaluation points.
    !    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, of the Jth basis function.
    !===============================================================================================
    subroutine lagrange_basis_1d (nd, xd, ni, xi, lb) 

        integer ( kind = 4 ) nd
        real ( kind = 8 ) xd(nd)
        integer ( kind = 4 ) ni
        real ( kind = 8 ) xi(ni)
        real ( kind = 8 ) lb(ni,nd)
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        
        do i = 1, ni
            do j = 1, nd
                lb(i,j) = product((xi(i) - xd(1:j-1)) / (xd(j) - xd(1:j-1))) * product((xi(i) - xd(j+1:nd)) / (xd(j) - xd(j+1:nd)))
            end do
        end do
        
    end subroutine lagrange_basis_1d
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used in each dimension.
    !    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
    !    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
    !    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
    !===============================================================================================
    subroutine lagrange_interp_nd_grid (m, n_1d, a, b, nd, xd)
        
        integer ( kind = 4 ) m
        integer ( kind = 4 ) n_1d(m)
        real ( kind = 8 ) a(m)
        real ( kind = 8 ) b(m)
        integer ( kind = 4 ) nd
        real ( kind = 8 ) xd(m,nd)
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) n
        real ( kind = 8 ), allocatable :: x_1d(:)
        
        xd(1:m,1:nd) = 0.0D+00
        do i = 1, m
            n = n_1d(i)
            allocate (x_1d(1:n))
            call cc_compute_points (n, x_1d)
            
            x_1d(1:n) = 0.5D+00 * ((1.0D+00 - x_1d(1:n)) * a(i) + (1.0D+00 + x_1d(1:n)) * b(i))
            
            call r8vec_direct_product (i, n, x_1d, m, nd, xd)
            deallocate (x_1d)
        end do
        
    end subroutine lagrange_interp_nd_grid
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule to be used in each dimension.
    !    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
    !    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
    !    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
    !===============================================================================================
    subroutine lagrange_interp_nd_grid2 (m, ind, a, b, nd, xd)
        
        integer ( kind = 4 ) m
        integer ( kind = 4 ) ind(m)
        real ( kind = 8 ) a(m)
        real ( kind = 8 ) b(m)
        integer ( kind = 4 ) nd
        real ( kind = 8 ) xd(m,nd)
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) n
        real ( kind = 8 ), allocatable :: x_1d(:)
        
        xd(1:m,1:nd) = 0.0D+00
        do i = 1, m
            call order_from_level_135 (ind(i), n)
            allocate (x_1d(1:n))
            call cc_compute_points (n, x_1d)
            
            x_1d(1:n) = 0.5D+00 * ((1.0D+00 - x_1d(1:n)) * a(i) + (1.0D+00 + x_1d(1:n)) * b(i))
            
            call r8vec_direct_product (i, n, x_1d, m, nd, xd)
            deallocate (x_1d)
        end do
        
    end subroutine lagrange_interp_nd_grid2
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used in each dimension.
    !    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
    !===============================================================================================
    subroutine lagrange_interp_nd_size (m, n_1d, nd)
        
        integer ( kind = 4 ) m
        integer ( kind = 4 ) n_1d(m)
        integer ( kind = 4 ) nd
        
        nd = product (n_1d(1:m))
        
    end subroutine lagrange_interp_nd_size
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule to be used in each dimension.
    !    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
    !===============================================================================================
    subroutine lagrange_interp_nd_size2 (m, ind, nd)

        integer ( kind = 4 ) m
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ind(m)
        integer ( kind = 4 ) n
        integer ( kind = 4 ) nd
        
        nd = 1
        do i = 1, m
            call order_from_level_135 (ind(i), n)
            nd = nd * n
        end do
        
    end subroutine lagrange_interp_nd_size2
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used in each dimension.
    !    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
    !    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
    !    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
    !    Input, integer ( kind = 4 ) NI, the number of points at which the interpolant is to be evaluated.
    !    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is to be evaluated.
    !    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the points XI.
    !===============================================================================================
    subroutine lagrange_interp_nd_value (m, n_1d, a, b, nd, zd, ni, xi, zi)

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n_1d(m)
        real ( kind = 8 ) a(m)
        real ( kind = 8 ) b(m)
        integer ( kind = 4 ) nd
        real ( kind = 8 ) zd(nd)
        integer ( kind = 4 ) ni
        real ( kind = 8 ) xi(m,ni)
        real ( kind = 8 ) zi(ni)
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) n
        real ( kind = 8 ), allocatable :: value(:)
        real ( kind = 8 ) w(nd)
        real ( kind = 8 ), allocatable :: x_1d(:)
        
        do j = 1, ni
            w(1:nd) = 1.0D+00
            do i = 1, m
                n = n_1d(i)
                allocate (x_1d(1:n))
                allocate (value(1:n))
                call cc_compute_points (n, x_1d)
                
                x_1d(1:n) = 0.5D+00 * ((1.0D+00 - x_1d(1:n)) * a(i) + (1.0D+00 + x_1d(1:n)) * b(i))
                
                call lagrange_basis_1d (n, x_1d, 1, xi(i,j), value)
                call r8vec_direct_product2 (i, n, value, m, nd, w)
                deallocate (value)
                deallocate (x_1d)
            end do
            
            zi(j) = dot_product (w, zd)
        end do
        
    end subroutine lagrange_interp_nd_value
    
    !$
    !===============================================================================================
    !  LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule to be used in each dimension.
    !    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
    !    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
    !    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
    !    Input, integer ( kind = 4 ) NI, the number of points at which the interpolant is to be evaluated.
    !    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is to be evaluated.
    !    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the points XI.
    !===============================================================================================
    subroutine lagrange_interp_nd_value2 (m, ind, a, b, nd, zd, ni, xi, zi)

        integer ( kind = 4 ) m
        integer ( kind = 4 ) ind(m)
        real ( kind = 8 ) a(m)
        real ( kind = 8 ) b(m)
        integer ( kind = 4 ) nd
        real ( kind = 8 ) zd(nd)
        integer ( kind = 4 ) ni
        real ( kind = 8 ) xi(m,ni)
        real ( kind = 8 ) zi(ni)
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) n
        real ( kind = 8 ), allocatable :: value(:)
        real ( kind = 8 ) w(nd)
        real ( kind = 8 ), allocatable :: x_1d(:)
        
        do j = 1, ni
            w(1:nd) = 1.0D+00
            do i = 1, m
                call order_from_level_135 (ind(i), n)
                allocate (x_1d(1:n))
                allocate (value(1:n))
                call cc_compute_points (n, x_1d)
                
                x_1d(1:n) = 0.5D+00 * ((1.0D+00 - x_1d(1:n)) * a(i) + (1.0D+00 + x_1d(1:n)) * b(i))
                
                call lagrange_basis_1d (n, x_1d, 1, xi(i,j), value)
                call r8vec_direct_product2 (i, n, value, m, nd, w)
                deallocate (value)
                deallocate (x_1d)
            end do
            
            zi(j) = dot_product (w, zd)
        end do
        
    end subroutine lagrange_interp_nd_value2
    
    !$
    !===============================================================================================
    !  ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
    !    Clenshaw Curtis rules, and some others, often use the following scheme:
    !    L: 0  1  2  3   4   5
    !    N: 1  3  5  9  17  33 ... 2^L+1
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) L, the level, which should be 0 or greater.
    !    Output, integer ( kind = 4 ) N, the order.
    !===============================================================================================
    subroutine order_from_level_135 (l, n)

        integer ( kind = 4 ) l
        integer ( kind = 4 ) n
        
        if (l < 0)  then
            write ( *, '(a)' ) ''
            write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
            write ( *, '(a)' ) '  Illegal input value of L!'
            stop
        else if (l == 0) then
            n = 1
        else
            n = (2**l) + 1
        end if
        
    end subroutine order_from_level_135
    
    !$
    !===============================================================================================
    !  R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
    !===============================================================================================
    subroutine r8mat_uniform_01 (m, n, seed, r)

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ), parameter :: i4_huge = 2147483647
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        real ( kind = 8 ) r(m,n)
        
        do j = 1, n
            do i = 1, m
                k = seed / 127773
                seed = 16807 * ( seed - k * 127773 ) - k * 2836
                if (seed < 0) then
                    seed = seed + i4_huge
                end if
                r(i,j) = real (seed, kind=8) * 4.656612875D-10
            end do
        end do
        
    end subroutine r8mat_uniform_01

    !$
    !===============================================================================================
    !  R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
    !  Discussion:
    !    An R8VEC is a vector of R8's.
    !    The affine vector L2 norm is defined as:
    !      R8VEC_NORM_AFFINE(V0,V1) = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
    ! --------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) N, the order of the vectors.
    !    Input, real ( kind = 8 ) V0(N), the base vector.
    !    Input, real ( kind = 8 ) V1(N), the vector whose affine norm is desired.
    !    Output, real ( kind = 8 ) R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
    !===============================================================================================
    function r8vec_norm_affine (n, v0, v1)

        integer ( kind = 4 ) n
        
        real ( kind = 8 ) r8vec_norm_affine
        real ( kind = 8 ) v0(n)
        real ( kind = 8 ) v1(n)
        
        r8vec_norm_affine = sqrt(sum((v0(1:n) - v1(1:n))**2))
        
    end function r8vec_norm_affine

    !$
    !===============================================================================================
    !  TIMESTAMP prints the current YMDHMS date as a time stamp.
    !---------------------------------------------------------------------------
    !===============================================================================================
    subroutine timestamp ()

        character ( len = 8 ) ampm
        integer ( kind = 4 ) d
        integer ( kind = 4 ) h
        integer ( kind = 4 ) m
        integer ( kind = 4 ) mm
        character ( len = 9 ), parameter, dimension(12) :: month = [    &
            &   'January  ', 'February ', 'March    ', 'April    ',     &
            &   'May      ', 'June     ', 'July     ', 'August   ',     &
            &   'September', 'October  ', 'November ', 'December ' ]
        integer ( kind = 4 ) n
        integer ( kind = 4 ) s
        integer ( kind = 4 ) values(8)
        integer ( kind = 4 ) y
        
        call date_and_time (values=values)
        
        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)
        
        if (h < 12) then
            ampm = 'AM'
        else if (h == 12) then
            if (n == 0 .and. s == 0) then
                ampm = 'Noon'
            else
                ampm = 'PM'
            end if
        else
            h = h - 12
            if (h < 12) then
                ampm = 'PM'
            else if (h == 12) then
                if (n == 0 .and. s == 0) then
                    ampm = 'Midnight'
                else
                    ampm = 'AM'
                end if
            end if
        end if
        
        write (*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) d, trim (month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)
        
    end subroutine timestamp

    !$
    !===============================================================================================
    !  R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
    !  Discussion:
    !    An R8VEC is a vector of R8's.
    !    To explain what is going on here, suppose we had to construct
    !    a multidimensional quadrature rule as the product of K rules for 1D quadrature.
    !    The product rule will be represented as a list of points and weights.
    !    The J-th item in the product rule will be associated with
    !      item J1 of 1D rule 1,
    !      item J2 of 1D rule 2,
    !      ...,
    !      item JK of 1D rule K.
    !
    !    In particular,
    !      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    !    and
    !      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
    !
    !    So we can construct the quadrature rule if we can properly
    !    distribute the information in the 1D quadrature rules.
    !
    !    This routine carries out that task for the abscissas X.
    !
    !    Another way to do this would be to compute, one by one, the
    !    set of all possible indices (J1,J2,...,JK), and then index
    !    the appropriate information.  An advantage of the method shown
    !    here is that you can process the K-th set of information and
    !    then discard it.
    !
    !  Example:
    !
    !    Rule 1:
    !      Order = 4
    !      X(1:4) = ( 1, 2, 3, 4 )
    !
    !    Rule 2:
    !      Order = 3
    !      X(1:3) = ( 10, 20, 30 )
    !
    !    Rule 3:
    !      Order = 2
    !      X(1:2) = ( 100, 200 )
    !
    !    Product Rule:
    !      Order = 24
    !      X(1:24) =
    !        ( 1, 10, 100 )
    !        ( 2, 10, 100 )
    !        ( 3, 10, 100 )
    !        ( 4, 10, 100 )
    !        ( 1, 20, 100 )
    !        ( 2, 20, 100 )
    !        ( 3, 20, 100 )
    !        ( 4, 20, 100 )
    !        ( 1, 30, 100 )
    !        ( 2, 30, 100 )
    !        ( 3, 30, 100 )
    !        ( 4, 30, 100 )
    !        ( 1, 10, 200 )
    !        ( 2, 10, 200 )
    !        ( 3, 10, 200 )
    !        ( 4, 10, 200 )
    !        ( 1, 20, 200 )
    !        ( 2, 20, 200 )
    !        ( 3, 20, 200 )
    !        ( 4, 20, 200 )
    !        ( 1, 30, 200 )
    !        ( 2, 30, 200 )
    !        ( 3, 30, 200 )
    !        ( 4, 30, 200 )
    !---------------------------------------------------------------------------
    !  Parameters:
    !    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
    !    processed.  The first factor processed must be factor 1!
    !    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
    !    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values for factor FACTOR_INDEX.
    !    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the direct product.
    !    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of
    !    the direct product, which are built up gradually.
    !  Local Parameters:
    !    Local, integer START, the first location of a block of values to set.
    !    Local, integer CONTIG, the number of consecutive values to set.
    !    Local, integer SKIP, the distance from the current value of START
    !    to the next location of a block of values to set.
    !    Local, integer REP, the number of blocks of values to set.
    !===============================================================================================
    subroutine r8vec_direct_product (factor_index, factor_order, factor_value, factor_num, point_num, x)

        integer ( kind = 4 ) factor_num
        integer ( kind = 4 ) factor_order
        integer ( kind = 4 ) point_num
        
        integer ( kind = 4 ), save :: contig
        integer ( kind = 4 ) factor_index
        real ( kind = 8 ) factor_value(factor_order)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ), save :: rep
        integer ( kind = 4 ), save :: skip
        integer ( kind = 4 ) start
        real ( kind = 8 ) x(factor_num,point_num)
        
        if (factor_index == 1) then
            contig = 1
            skip = 1
            rep = point_num
            x(1:factor_num,1:point_num) = 0.0D+00
        end if
        
        rep = rep / factor_order
        skip = skip * factor_order
        
        do j = 1, factor_order
            start = 1 + (j - 1) * contig
            do k = 1, rep
                x(factor_index,start:start+contig-1) = factor_value(j)
                start = start + skip
            end do
        end do
        contig = contig * factor_order
        
    end subroutine r8vec_direct_product

    !$
    !===============================================================================================
    !! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
    !  Discussion:
    !    An R8VEC is a vector of R8's.
    !    To explain what is going on here, suppose we had to construct
    !    a multidimensional quadrature rule as the product of K rules
    !    for 1D quadrature.
    !    The product rule will be represented as a list of points and weights.
    !    The J-th item in the product rule will be associated with
    !      item J1 of 1D rule 1,
    !      item J2 of 1D rule 2,
    !      ...,
    !      item JK of 1D rule K.
    !
    !    In particular,
    !      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    !    and
    !      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
    !
    !    So we can construct the quadrature rule if we can properly
    !    distribute the information in the 1D quadrature rules.
    !
    !    This routine carries out the task involving the weights W.
    !
    !    Another way to do this would be to compute, one by one, the
    !    set of all possible indices (J1,J2,...,JK), and then index
    !    the appropriate information.  An advantage of the method shown
    !    here is that you can process the K-th set of information and
    !    then discard it.
    !
    !  Example:
    !
    !    Rule 1:
    !      Order = 4
    !      W(1:4) = ( 2, 3, 5, 7 )
    !
    !    Rule 2:
    !      Order = 3
    !      W(1:3) = ( 11, 13, 17 )
    !
    !    Rule 3:
    !      Order = 2
    !      W(1:2) = ( 19, 23 )
    !
    !    Product Rule:
    !      Order = 24
    !      W(1:24) =
    !        ( 2 * 11 * 19 )
    !        ( 3 * 11 * 19 )
    !        ( 4 * 11 * 19 )
    !        ( 7 * 11 * 19 )
    !        ( 2 * 13 * 19 )
    !        ( 3 * 13 * 19 )
    !        ( 5 * 13 * 19 )
    !        ( 7 * 13 * 19 )
    !        ( 2 * 17 * 19 )
    !        ( 3 * 17 * 19 )
    !        ( 5 * 17 * 19 )
    !        ( 7 * 17 * 19 )
    !        ( 2 * 11 * 23 )
    !        ( 3 * 11 * 23 )
    !        ( 5 * 11 * 23 )
    !        ( 7 * 11 * 23 )
    !        ( 2 * 13 * 23 )
    !        ( 3 * 13 * 23 )
    !        ( 5 * 13 * 23 )
    !        ( 7 * 13 * 23 )
    !        ( 2 * 17 * 23 )
    !        ( 3 * 17 * 23 )
    !        ( 5 * 17 * 23 )
    !        ( 7 * 17 * 23 )
    !---------------------------------------------------------------------------
    !  Parameters:
    !    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
    !    processed.  The first factor processed must be factor 1!
    !    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
    !    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values for factor FACTOR_INDEX.
    !    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the direct product.
    !    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the direct product, which are built up gradually.
    !  Local Parameters:
    !    Local, integer ( kind = 4 ) START, the first location of a block of values to set.
    !    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values to set.
    !    Local, integer SKIP, the distance from the current value of START
    !    to the next location of a block of values to set.
    !    Local, integer REP, the number of blocks of values to set.
    !===============================================================================================
    subroutine r8vec_direct_product2 (factor_index, factor_order, factor_value, factor_num, point_num, w)

        integer ( kind = 4 ) factor_num
        integer ( kind = 4 ) factor_order
        integer ( kind = 4 ) point_num
        
        integer ( kind = 4 ), save :: contig
        integer ( kind = 4 ) factor_index
        real ( kind = 8 ) factor_value(factor_order)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ), save :: rep
        integer ( kind = 4 ), save :: skip
        integer ( kind = 4 ) start
        real ( kind = 8 ) w(point_num)
        
        if (factor_index == 1) then
            contig = 1
            skip = 1
            rep = point_num
            w(1:point_num) = 1.0D+00
        end if
        
        rep = rep / factor_order
        skip = skip * factor_order
        
        do j = 1, factor_order
            start = 1 + (j - 1) * contig
            do k = 1, rep
                w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
                start = start + skip
            end do
        end do
        contig = contig * factor_order
        
    end subroutine r8vec_direct_product2

end module self_lagrange_interp
