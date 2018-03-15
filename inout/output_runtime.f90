!$
!===================================================================================================
!
!   print run-time information for this code
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Print_memory_size
!
!   Public type lists:          No
!
!===================================================================================================
module output_runtime
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use coefficient_iteration
        
    implicit none 
    private
    public  :: Print_memory_size
    
contains
    !$
    !===============================================================================================
    ! print momery information
    !===============================================================================================
    subroutine Print_memory_size (unit_, case_name)

        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: case_name
        
        ! ----------------------------------------------------------------------
        ! parameters define
        integer, parameter  :: N_ARRAY = 10                                     ! number of array per type
        integer, parameter  :: N_TYPE  = 10                                     ! number of type this program
        
        real(KREAL)  :: array_memory(N_ARRAY)                               ! memory size for single array
        real(KREAL)  :: part_memory(N_TYPE)                                 ! memory size for this part
        real(KREAL)  :: total_memory                                        ! memory size for all program
        
        ! ----------------------------------------------------------------------
        ! header
        write(unit=unit_, fmt="(1x, A)") '(Begin)'
        write(unit=unit_, fmt="(1x, '_________________________________________', &
            &  '__________________________________________')")
        write(unit=unit_, fmt="(1x, A, 1x)", advance='no') 'Case name is:'
        write(unit=unit_, fmt="('[', A, ']')") TRIM(case_name)
        write(unit=unit_, fmt="(/)")
        
        ! ----------------------------------------------------------------------
        ! get basic information
        array_memory = 0.0
        part_memory = 0.0
        total_memory = 0.0
        
        ! ----------------------------------------------------------------------
        ! print geometry @1
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for geometry:'
        write(unit=unit_, fmt=*)  '    '
            
        array_memory( 1) = mesh%memory
        array_memory( 2) = geom%memory
        array_memory( 3) = bound%memory
        
        part_memory(1) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(1), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'mesh  --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(1), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'geom  --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(1), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'bound --', array_memory( 3), 'MB', 100.0*array_memory( 3)/part_memory(1), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print material @2
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for material:'
        write(unit=unit_, fmt=*)  '    '
        
        array_memory( 1) = xsec_inp%memory + xsec%memory + xsec_iter%memory
        array_memory( 2) = param_inp%memory + param%memory
        array_memory( 3) = Q_ext%memory
        
        part_memory(2) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(2), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'xsec                --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(2), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'kinetics parameter  --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(2), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'external source     --', array_memory( 3), 'MB', 100.0*array_memory( 3)/part_memory(2), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print quadrature & sweeping @3
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for quadrature & sweeping:'
        write(unit=unit_, fmt=*)  '    '
        
        array_memory( 1) = quad%memory
        array_memory( 2) = sweep%memory
        
        part_memory(3) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(3), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'quad   --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(3), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'sweep  --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(3), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print iteration parameter @4
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for iteration parameter:'
        write(unit=unit_, fmt=*)  '    '
        
        array_memory( 1) = iter_q%memory
        array_memory( 2) = iter_flux%memory
        array_memory( 3) = flux_scat%memory
        array_memory( 4) = iter_adjoint%memory
        
        part_memory(4) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(4), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'iter_q        --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(4), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'iter_flux     --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(4), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'flux_scat     --', array_memory( 3), 'MB', 100.0*array_memory( 4)/part_memory(4), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'iter_adjoint  --', array_memory( 4), 'MB', 100.0*array_memory( 5)/part_memory(4), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print result container @5
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for result container:'
        write(unit=unit_, fmt=*)  '    '
        
        array_memory( 1) = flux_forward%memory
        array_memory( 2) = flux_adjoint%memory
        
        array_memory( 3) = dist_flux%memory
        array_memory( 4) = dist_power%memory
        array_memory( 5) = dist_fission_rate%memory
        
        if (SIZE(dist_dnps) > 0)  then
            array_memory( 6) = SIZE(dist_dnps) * dist_dnps(1)%memory
        end if
        
        part_memory(5) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(5), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'flux_forward       --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(5), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'flux_adjoint       --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(5), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'dist_flux          --', array_memory( 3), 'MB', 100.0*array_memory( 3)/part_memory(5), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'dist_power         --', array_memory( 4), 'MB', 100.0*array_memory( 4)/part_memory(5), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'dist_fission_rate  --', array_memory( 5), 'MB', 100.0*array_memory( 5)/part_memory(5), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'dist_dnps     --', array_memory( 6), 'MB', 100.0*array_memory( 6)/part_memory(5), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print coefficient for iteration @6
        array_memory = 0.0
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory size for coefficient in iteration:'
        write(unit=unit_, fmt=*)  '    '
        
        array_memory( 1) = coeff_source%memory
        array_memory( 2) = coeff_surface%memory
        array_memory( 3) = coeff_nodal%memory
        
        part_memory(6) = SUM(array_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(6), 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'coeff_source    --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(6), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'coeff_surface   --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(6), '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'coeff_nodal     --', array_memory( 3), 'MB', 100.0*array_memory( 3)/part_memory(6), '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print transient @7
        if (nt%flag%is_transient)  then
            array_memory = 0.0
            write(unit=unit_, fmt=*)  '_______________________________________________________'
            write(unit=unit_, fmt=*)  'memory size for transient:'
            write(unit=unit_, fmt=*)  '    '
            
            array_memory( 2) = shape_last%memory + shape_current%memory + shape_predict%memory
            
            part_memory(7) = SUM(array_memory)
            
            write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'memory for this part is: ', part_memory(7), 'MB'
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'concentration    --', array_memory( 1), 'MB', 100.0*array_memory( 1)/part_memory(7), '%'
            write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'shape function   --', array_memory( 2), 'MB', 100.0*array_memory( 2)/part_memory(7), '%'
            write(unit=unit_, fmt=*)  '    '
        end if
        
        ! ----------------------------------------------------------------------
        ! print summary
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  '_______________________________________________________'
        write(unit=unit_, fmt=*)  'memory SUMMARY:'
        write(unit=unit_, fmt=*)  '    '
        
        total_memory = SUM(part_memory)
        
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A)")  'total memory for this case is: ', total_memory, 'MB'
        write(unit=unit_, fmt=*)  '    '
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'geometry       --', part_memory( 1), 'MB', 100.0*part_memory( 1)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'material       --', part_memory( 2), 'MB', 100.0*part_memory( 2)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'quad & sweep   --', part_memory( 3), 'MB', 100.0*part_memory( 3)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'iteration      --', part_memory( 4), 'MB', 100.0*part_memory( 4)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'container      --', part_memory( 5), 'MB', 100.0*part_memory( 5)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'coefficient    --', part_memory( 6), 'MB', 100.0*part_memory( 6)/total_memory, '%'
        write(unit=unit_, fmt="(1x, A, 1x, ES11.4, 1x, A, 1x, F7.3, A)")  'transient      --', part_memory( 7), 'MB', 100.0*part_memory( 7)/total_memory, '%'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! tail
        write(unit=unit_, fmt="(1x, /)")
        write(unit=unit_, fmt="(1x, '_________________________________________', &
            &  '__________________________________________')")
        write(unit=unit_, fmt="(1x, A, /)")  '(End)'
        
    end subroutine Print_memory_size
    
end module output_runtime
