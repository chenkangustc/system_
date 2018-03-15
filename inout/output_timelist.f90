!$
!===================================================================================================
!
!   print timelist parameter per time point
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Print_header_timelist
!                               Print_timelist
!                               
!   Public type lists:          No
!
!===================================================================================================
module output_timelist
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use th_global,              only : hot_channel
    use contain_header,         only : TimeListParameter
    
	use imp_assm_global
    implicit none 
    private
    public  :: Print_header_timelist, Print_timelist
    
contains
    !$
    !===============================================================================================
    ! Print head information to the output file and the screen
    !===============================================================================================
    subroutine Print_header_timelist (unit_, casename)
        
        integer, intent(in)             :: unit_
        character(len=*), intent(in)    :: casename
        
        write(unit=unit_, fmt="(1x, A)") '(Begin)'
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        write(unit=unit_, fmt="(1x, A, 1x)", advance='no') 'Case name is:'
        write(unit=unit_, fmt="('[', A, ']')") TRIM(casename)
        
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt="(1X, A, /)") 'Listed parameters per time step:'
        write(unit=unit_, fmt="(1x, A)") &
            & '  step     time        reactivity       beta        power-level    <------   fxy   ------>     fz@pos     <------    Tm    ------>   <----- Tf ----->  <---- Tclad ----> '
        write(unit=unit_, fmt="(1x, A)") &
            & '   No.     (s)            ($)          (pcm)           (%)           nodal@pos    FA@pos                    max     outlet    avg      max      avg     max-in   max-out '
        write(unit=unit_, fmt="(2x, A)") &
            & '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
        
    end subroutine Print_header_timelist
    
    !$
    !===============================================================================================
    ! print integral parameter 
    !===============================================================================================
    subroutine Print_timelist (timelist, tidx, ctime, unit_)
        
        type(TimeListParameter), intent(in) :: timelist
        integer, intent(in)       :: tidx
        real(KREAL), intent(in)   :: ctime
        integer, intent(in)       :: unit_
    
        integer  :: ip		
        real(KREAL)  :: relative_power
        real(KREAL)  :: reactivity_dollar
		
	    !local
		integer i,j
		real(KREAL)  ::max_T_inner,max_T_outer
        
        call time_program%update (silent=.TRUE.)
        
        relative_power = timelist%power / ns%flag%rated_power * 100.0D0
        if (timelist%beta >= EPS_ZERO)  then
            reactivity_dollar = timelist%reactivity / timelist%beta
        else
            reactivity_dollar = 0.0D0 
        end if 
        !calculate the max inner and outer Temperature
		max_T_inner=0.0
		max_T_outer=0.0
		do i=1,ns%state%zone,1
			do j=1,assm1(i)%mesh%Ny,1
				if(max_T_inner<assm1(i)%thermal%Tgs(j)) max_T_inner=assm1(i)%thermal%Tgs(j)
				if(max_T_outer<assm1(i)%thermal%Tsc(j)) max_T_outer=assm1(i)%thermal%Tsc(j)
			enddo
		enddo
        write(unit=unit_, fmt="(1X, (I6, TR1, ES14.5, TR1), 3(ES13.6, TR2), 3(F7.4, TR1, I3, TR2), *(F7.2, TR2))")   &
            &   tidx, ctime,                                                                                       &
            &   reactivity_dollar, timelist%beta*1.0E5, relative_power,                                                         &
            &   timelist%fxy_nodal, timelist%idx_nodal, timelist%fxy_FA, timelist%idx_FA,                                   &
            &   timelist%fz, timelist%idx_z,                                                                    &
            &   timelist%Tm_max, timelist%Tm_outlet, timelist%Tm_avg, timelist%Tf_max, timelist%Tf_avg,         &
            &   max_T_inner, max_T_outer 
			!&   MAXVAL(hot_channel%tclad_inner), MAXVAL(hot_channel%tclad_surf)
        ! to screen
        write(unit=OUTPUT_UNIT, fmt="(1x, 2(A, ES12.5), (A, ES10.3, A))") 't =',  ctime, ' sec., Power =', relative_power, '%;   t-elapsed =', time_program%total, ' sec.'
        
    end subroutine Print_timelist
    
end module output_timelist
