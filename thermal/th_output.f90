!$
!===================================================================================================
!
!   module for thermal calculation information output
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Print_hotpoint
!                               Print_RBFD_info
!
!   Public type lists:          No
!
!===================================================================================================
module th_output

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    use th_global
    
    implicit none
    private 
    public  :: Print_hotpoint, Print_RBFD_info

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_hotpoint (unit_, tidx, ctime)
        
        integer, intent(in)          :: unit_
        integer, intent(in)          :: tidx
        real(KREAL), intent(in)  :: ctime
        
        write(unit=unit_, fmt="(1x, I4, TR2, ES11.4, TR4, *(F11.4, TR4))") tidx, ctime,          &
                                                    &   hot_point%tfuel_center, hot_point%tclad_inner,          &
                                                    &   hot_point%tclad_surf, hot_point%tcoolant
    
    end subroutine Print_hotpoint
    
    !$
    !===============================================================================================
    ! print information for the consequent RBFD analysis
    !===============================================================================================
    subroutine Print_RBFD_info (unit_)
    
        integer, intent(in)     :: unit_
        
        integer  :: ir, ia
        
        ! power info
        write(unit=unit_, fmt=*)  '___________________________________'
        write(unit=unit_, fmt=*)  'average_linear_power:'
        do ia = 1, nth%na
            write(unit=unit_, fmt="(1x, I4, TR3, *(ES13.6, TR3))")  ia, th_power%avg_linear(ia, :)
        end do
    
        write(unit=unit_, fmt=*)  '___________________________________'
        write(unit=unit_, fmt=*)  'maxima_linear_power:'
        do ia = 1, nth%na
            write(unit=unit_, fmt="(1x, I4, TR3, *(ES13.6, TR3))")  ia, th_power%max_linear(ia, :)
        end do
                
        ! flow info
        write(unit=unit_, fmt=*)  '___________________________________'
        write(unit=unit_, fmt=*)  'assembly flow:'
        write(unit=unit_, fmt="(1x, *(I10, TR6))")     geom_th%geom_type
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  design%assembly_flow
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  design%channel_flowrate
        write(unit=unit_, fmt="(1x, *(L10, TR6))")     design%is_active_channel
    
    end subroutine Print_RBFD_info
    
end module th_output
