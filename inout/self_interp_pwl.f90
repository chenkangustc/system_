!$
!===================================================================================================
!
!   self version: piecewise linear interpolant to data defined on a 2D grid 
!---------------------------------------------------------------------------------------------------
!   public subroutine lists:    No
!
!   public type lists:          No
!
!===================================================================================================
module self_interp_pwl

    implicit none
    public 

contains
    !$
    !===============================================================================================
    ! Input, integer ( kind = 4 ) NXD, NYD, the number of X and Y data values.
    ! Input, real ( kind = 8 ) XD(NXD), YD(NYD), the sorted X and Y data.
    ! Input, real ( kind = 8 ) ZD(NXD,NYD), the Z data.
    ! Input, integer ( kind = 4 ) NI, the number of interpolation points.
    ! Input, real ( kind = 8 ) XI(NI), YI(NI), the coordinates of the interpolation points.
    ! Output, real ( kind = 8 ) ZI(NI), the value of the interpolant.
    !===============================================================================================
    subroutine pwl_interp_2d (nxd, nyd, xd, yd, zd, ni, xi, yi, zi)

        integer ( kind = 4 ) nxd
        integer ( kind = 4 ) nyd
        integer ( kind = 4 ) ni
        real ( kind = 8 ) xd(nxd)
        real ( kind = 8 ) yd(nyd)
        real ( kind = 8 ) zd(nxd,nyd)
        real ( kind = 8 ) xi(ni)
        real ( kind = 8 ) yi(ni)
        real ( kind = 8 ) zi(ni)
        
        real ( kind = 8 ) alpha
        real ( kind = 8 ) beta
        real ( kind = 8 ) det
        real ( kind = 8 ) dxa
        real ( kind = 8 ) dxb
        real ( kind = 8 ) dxi
        real ( kind = 8 ) dya
        real ( kind = 8 ) dyb
        real ( kind = 8 ) dyi
        real ( kind = 8 ) gamma
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        
        do k = 1, ni
            i = r8vec_bracket5 (nxd, xd, xi(k))
            if (i == -1)  then
                zi(k) = r8_huge ()
                cycle
            end if
            
            j = r8vec_bracket5 (nyd, yd, yi(k))
            if ( j == -1 ) then
                zi(k) = r8_huge ()
                cycle
            end if
            
            if (yi(k) < yd(j+1) + (yd(j) - yd(j+1)) * (xi(i) - xd(i)) / (xd(i+1) - xd(i)))  then
                dxa = xd(i+1) - xd(i)
                dya = yd(j)   - yd(j)
                dxb = xd(i)   - xd(i)
                dyb = yd(j+1) - yd(j)
                dxi = xi(k)   - xd(i)
                dyi = yi(k)   - yd(j)
                
                det = dxa * dyb - dya * dxb
                alpha = ( dxi * dyb - dyi * dxb ) / det
                beta =  ( dxa * dyi - dya * dxi ) / det
                gamma = 1.0D+00 - alpha - beta
                zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)
                
            else
                dxa = xd(i)   - xd(i+1)
                dya = yd(j+1) - yd(j+1)
                dxb = xd(i+1) - xd(i+1)
                dyb = yd(j)   - yd(j+1)
                dxi = xi(k)   - xd(i+1)
                dyi = yi(k)   - yd(j+1)
                
                det = dxa * dyb - dya * dxb
                alpha = ( dxi * dyb - dyi * dxb ) / det
                beta =  ( dxa * dyi - dya * dxi ) / det
                gamma = 1.0D+00 - alpha - beta
                zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)
            end if
        end do
        
    end subroutine pwl_interp_2d
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !  R8VEC_LINSPACE creates a vector of linearly spaced values.
    !  Discussion:
    !    An R8VEC is a vector of R8's.
    !    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
    !    In other words, the interval is divided into N-1 even subintervals,
    !    and the endpoints of intervals are used as the points.
    !---------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
    !    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
    !    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
    !===============================================================================================
    subroutine r8vec_linspace (n, a, b, x)

        integer ( kind = 4 ) n
        real ( kind = 8 ) a
        real ( kind = 8 ) b
        real ( kind = 8 ) x(n)
        integer ( kind = 4 ) i
        
        if (n == 1) then
            x(1) = (a + b) / 2.0D+00
        else
            do i = 1, n
                x(i) = (real(n-i, kind=8) * a + real(i-1, kind=8) * b ) / real(n-1, kind=8)
            end do
        end if
        
    end subroutine r8vec_linspace
        
    !$
    !===============================================================================================
    !  R8VEC3_PRINT prints an R8VEC3.
    !  Discussion:
    !    An R8VEC3 is a dataset consisting of N triples of R8's, stored
    !    as three separate vectors A1, A2, A3.
    !---------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) N, the number of components of the vector.
    !    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), the vectors to be printed.
    !    Input, character ( len = * ) TITLE, a title.
    !===============================================================================================
    subroutine r8vec3_print (n, a1, a2, a3, title)
    
        integer ( kind = 4 ) n
        real ( kind = 8 ) a1(n)
        real ( kind = 8 ) a2(n)
        real ( kind = 8 ) a3(n)
        character ( len = * ) title
        integer ( kind = 4 ) i
        
        write (*, '(a)') ' '
        write (*, '(a)') trim (title)
        write (*, '(a)') ' '
        
        do i = 1, n
            write (*, '(i8,3g14.6)') i, a1(i), a2(i), a3(i)
        end do
        
    end subroutine r8vec3_print

    !$
    !===============================================================================================
    !  R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
    !  Discussion:
    !    An R8VEC is a vector of R8's.
    !    The affine vector L2 norm is defined as:
    !      R8VEC_NORM_AFFINE(V0,V1) = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
    !---------------------------------------------------------------------------
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
        
        r8vec_norm_affine = sqrt(sum((v0(1:n)-v1(1:n))**2))
        
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
    !  R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
    !  Discussion:
    !    We assume XD is sorted.
    !    If XI is contained in the interval [XD(1),XD(N)], then the returned 
    !    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
    !    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
    !    This code implements a version of binary search which is perhaps more
    !    understandable than the usual ones.
    !---------------------------------------------------------------------------
    !    Input, integer ( kind = 4 ) ND, the number of data values.
    !    Input, real ( kind = 8 ) XD(N), the sorted data.
    !    Input, real ( kind = 8 ) XD, the query value.
    !    Output, integer ( kind = 4 ) R8VEC_BRACKET5, the bracket information.
    !===============================================================================================
    function r8vec_bracket5 (nd, xd, xi)

        integer ( kind = 4 ) nd
        real ( kind = 8 ) xd(nd)
        real ( kind = 8 ) xi
        integer ( kind = 4 ) b
        integer ( kind = 4 ) l
        integer ( kind = 4 ) m
        integer ( kind = 4 ) r
        integer ( kind = 4 ) r8vec_bracket5
        
        if (xi < xd(1) .or. xd(nd) < xi) then
            b = -1
        else
            l = 1
            r = nd
            do while (l + 1 < r)
                m = (l + r) / 2
                if (xi < xd(m)) then
                    r = m
                else
                    l = m
                end if
            end do
            b = l
        end if
        
        r8vec_bracket5 = b
        
    end function r8vec_bracket5
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    function r8_huge ()
        real ( kind = 8 ) r8_huge
        r8_huge = 1.0D+30
    end function r8_huge
    
end module self_interp_pwl
