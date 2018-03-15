!$
!===================================================================================================
!
!   perform steady thermal analysis by parallel channel
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_ParallelChannel_steady
!                               Driving_ParallelChannel_transient
!
!   Public type lists:          No
!
!===================================================================================================
module th_parallel_channel

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    use global_state
    use th_global
    use th_output
    
    use gth_parallel_kernel
        
    implicit none
    private 
    public  :: Driving_ParallelChannel_steady, Driving_ParallelChannel_transient
    
contains
    !$
    !===============================================================================================
    ! drivier subroutine of this module
    !===============================================================================================
    subroutine Driving_ParallelChannel_steady (power, fq_core)

        real(KREAL), intent(in)  :: power(nth%na, nth%nr)                       ! power, in Watt
        real(KREAL), intent(in)  :: fq_core(nth%na, nth%nr)                     ! power peak from core calculation
        
        ! set power
        call th_power%set_power (nth, geom_th, geom_assm, power, fq_core)
        
        ! start calculation
        call design%set (a_coolant, nth, geom_th, geom_assm, th_power%max_linear)
        
        ! get flowrate
        if (design%is_search)  then
            call Get_flowrate (design, a_coolant, nth, geom_th, geom_assm, th_power%max_linear, is_init=.TRUE.)
        end if
        call avg_channel%fix_bottom (a_coolant, nth, design)
        call hot_channel%fix_bottom (a_coolant, nth, design)
        
        ! average channel calculation to obtain feedback parameter
        call CalcuThermal (design, nth, geom_th, geom_assm, avg_channel, th_power%avg_linear)
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, avg_channel, th_power%avg_linear)
        call avg_channel%fix_top (nth, design)
        
        ! hot channel calculation to obtain maximum temperature to testify the design limit
        call CalcuThermal (design, nth, geom_th, geom_assm, hot_channel, th_power%max_linear)
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, hot_channel, th_power%max_linear)
        call hot_channel%fix_top (nth, design)
        
        ! get hot channel inforamtion and feedback information
        call hot_point%set (hot_channel, nth)
        
!        call avg_channel%print_avg (FILES%TH_AVERAGE, design, nth, 0, 0.0_KREAL)
!        call Print_hotpoint (FILES%TH_HOT, 0, 0.0_KREAL)
        call Print_RBFD_info (FILES%TH_RBFD)
        
!        call avg_channel%print (nth, geom_th, 144)
!        call avg_channel%print_avg (145, design, nth, 0, 0.0)
!        call avg_channel%print_max (146, design, nth, 0, 0.0)
!        call hot_channel%print (nth, geom_th, 145)
!        call design%print (146)
!        stop(0)
        
!        write(333, fmt="(1x, *(ES12.5, TR3))") 0.0, SUM(design%init_flow), SUM(design%init_flowrate), SUM(design%assembly_flow), SUM(design%channel_flowrate)
        
    end subroutine Driving_ParallelChannel_steady
    
    !$
    !===============================================================================================
    ! drive subroutine of this module
    !===============================================================================================
    subroutine Driving_ParallelChannel_transient (power, fq_core, tidx, ltime, ctime)
                
        real(KREAL), intent(in)  :: power(nth%na, nth%nr)
        real(KREAL), intent(in)  :: fq_core(nth%na, nth%nr)                     ! power peak from core calculation
        integer, intent(in)      :: tidx
        real(KREAL), intent(in)  :: ltime
        real(KREAL), intent(in)  :: ctime
        
        real(KREAL) :: step_length
        integer     :: i 
        integer, save  :: idx = 0 
        
!!!        if (idx == 0)  then
!!!            call th_power%print (300+idx, 300+idx)
!!!            call avg_channel%print (nth, geom_th, 400+idx) 
!!!            do i = 1, 20
!!!                call avg_channel%print_rod (500+idx, nth, ir=i)
!!!            end do 
!!!        end if 
        
        ! update inlet flow & Tm
        call Update_inlet_condition (ctime)
        
        ! set initial
        step_length = ctime - ltime
        call th_power%set_power (nth, geom_th, geom_assm, power, fq_core)
        call avg_channel%fix_bottom (a_coolant, nth, design)
        call hot_channel%fix_bottom (a_coolant, nth, design)
        
        ! average channel calculation to obtain feedback parameter
        call CalcuThermal_transient (design, nth, geom_th, geom_assm, th_power, avg_channel, th_power%avg_linear, step_length)
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, avg_channel, th_power%avg_linear, step_length)
        call avg_channel%fix_top (nth, design)
               
        ! hot channel calculation to obtain maximum temperature to testify the design limit
        call CalcuThermal_transient (design, nth, geom_th, geom_assm, th_power, hot_channel, th_power%max_linear, step_length)
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, hot_channel, th_power%max_linear, step_length)
        call hot_channel%fix_top (nth, design)
        
        ! get hot channel inforamtion and feedback information
        call hot_point%set (hot_channel, nth)
        
!        call Print_hotpoint (FILES%TH_HOT, tidx, ctime)
!        call avg_channel%print_avg (FILES%TH_AVERAGE, design, nth, tidx, ctime)
        
!!!        idx = idx + 1
!!!        if (idx <= 5)  then
!!!            call th_power%print (300+idx, 300+idx)
!!!            call avg_channel%print (nth, geom_th, 400+idx)
!!!            do i = 1, 20 
!!!                call avg_channel%print_rod (500+idx, nth, ir=i)
!!!            end do 
!!!        end if 
!!!        
        
    end subroutine Driving_ParallelChannel_transient
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Update_inlet_condition (ctime)
        
        real(KREAL), intent(in)  :: ctime
        
        real(KREAL)  :: frac
        real(KREAL)  :: natural
        real(KREAL)  :: dt
        real(KREAL)  :: a, b, c, d 
        real(KREAL)  :: ff 
        
        integer  :: i_pert
        integer  :: irbeg, irend 
        integer  :: ir 
        
        if (design%is_blockage .and. ctime>design%block_time)  then
            natural = 1.0D0 - (design%percentage * 0.01D0)
            dt = ctime - design%block_time
            frac = 1.0D0 / (1.0D0 + dt/5.5D0)
            if (frac <= natural)  then
                frac = natural
            end if
            do ir = 1, nth%nr
                design%assembly_flow(ir) = design%init_flow(ir) * frac
                design%channel_flowrate(ir) = design%init_flowrate(ir) * frac
            end do
        end if
        
        if (nt%perturb%is_flow) then
            perturb1: do i_pert = 1, pert_th%n_flow
                associate(pth => pert_th%flow_perts(i_pert))
                dt = ctime - pth%time_start
                natural = pth%natural
                if (SIZE(pth%variables) /= 4)  then
                end if 
                a = pth%variables(1)
                b = pth%variables(2)
                c = pth%variables(3)
                d = pth%variables(4)
                
                if (dt < 0.0D0)  then
                    frac = 1.0D0
                else
                    select case(pth%type)
                    case(1)
                        frac = 1.0D0 / (1.0D0 + dt/a)
                    case default
                    end select
                end if 
                
                if (frac <= natural)  then
                    frac = natural
                end if
                
                if (pth%channelID == 0)  then
                    irbeg = 1
                    irend = nth%nr
                else
                    irbeg = pth%channelID
                    irend = pth%channelID
                end if 
                
                do ir = irbeg, irend 
                    design%assembly_flow(ir) = design%init_flow(ir) * frac
                    design%channel_flowrate(ir) = design%init_flowrate(ir) * frac
                end do
                end associate 
            end do perturb1
        end if
        
        if (nt%perturb%is_Tm)  then
            perturb2: do i_pert = 1, pert_th%n_Tm
                associate(pth => pert_th%Tm_perts(i_pert))
                dt = ctime - pth%time_start
                if (SIZE(pth%variables) /= 4)  then
                end if 
                a = pth%variables(1)
                b = pth%variables(2)
                c = pth%variables(3)
                d = pth%variables(4)
                
                select case(pth%type)
                case(1)
                    if (dt <= 0.0D0)  then
                        design%tcoolin = design%init_tcoolin
                    else
                        design%tcoolin = design%init_tcoolin + a 
                    end if 
                    
                case(2)
                    if (dt <= 0.0D0)  then
                        design%tcoolin = design%init_tcoolin
                    else if (dt <= (pth%time_end - pth%time_start))  then 
                        design%tcoolin = design%init_tcoolin + a * dt
                    else
                        design%tcoolin = design%init_tcoolin + a * (pth%time_end - pth%time_start)
                    end if 
                
                case default
                end select 
                
                end associate
            end do perturb2
        end if 

        write(778, fmt="(1x, L, *(ES12.5, TR3))") nt%perturb%is_flow, ctime, SUM(design%init_flow), SUM(design%init_flowrate), SUM(design%assembly_flow), SUM(design%channel_flowrate), frac, natural
        
    end subroutine Update_inlet_condition
    
end module th_parallel_channel
