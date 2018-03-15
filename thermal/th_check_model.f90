!$
!===================================================================================================
!
!   module for thermal calculation pre process
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_th_check_model 
!
!   Public type lists:          No
!
!===================================================================================================
module th_check_model

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    use th_global
    
    use gth_parallel_kernel
    
    implicit none
    private 
    public  :: Driving_th_check_model
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_th_check_model ()
        
        real(KREAL)  :: power(30,1), fq_core(30,1)
        real(KREAL)  :: dst
        
        integer  :: itasm
        integer  :: ir, ia, i, j 
        integer  :: i_allocate
        
        real(KREAL)  :: tbeg, tend
        real(KREAL)  :: ltime, ctime, step_length
        integer  :: istep, nstep
        
        
        nth%nf = 50
        nth%nc = 5
        nth%n_mesh = nth%nf + nth%nc + 2
        
        nth%nr = 1
        nth%na = 30
        nth%na_start = 1
        nth%na_end = 30
        nth%n_assm_geom = 1
        
        call design%alloc (nth)
        call avg_channel%alloc (nth)
        call hot_channel%alloc (nth)
        
        call geom_th%alloc (nth)
        allocate(geom_assm(nth%n_assm_geom), stat=i_allocate)
        do i = 1, SIZE(geom_assm)
            call geom_assm(i)%alloc (nth)
        end do
        
        call th_power%alloc (nth)
        
        ! pointer allocated
        call Prepare_property_pointer ()
        
        ! ----------------------------------------------------------------------
        ! set initial value
        a_coolant => a_coolant_HLM
        call a_coolant%set (3, option=0)
        
        ! geometry information
        geom_th%height(1:nth%na)    = 3.0E-2                                    ! m
        geom_th%geom_type(1:nth%nr) = 1
        geom_th%cladding_type(:)    = 21 
        geom_th%gap_type(:)         = 1
        geom_th%fuel_type(:)        = 3 
        
        do itasm = 1, nth%n_assm_geom
            call geom_assm(itasm)%alloc (nth)
            geom_assm(itasm)%n_pin       = 90
            geom_assm(itasm)%n_fuelpin   = 90
            geom_assm(itasm)%pitch    = 100.0*1.3406E-2                               ! cm
            geom_assm(itasm)%rod      = 100.0*8.5E-3/2.0                              ! cm
            geom_assm(itasm)%cladth   = 100.0*0.565E-3                                ! cm
            geom_assm(itasm)%bond     = 100.0*0.115E-3                                ! cm
            geom_assm(itasm)%hole     = 100.0*0.90E-3                                 ! cm
            call  geom_assm(itasm)%set (nth)
        end do
        
        ! power input
        fq_core(:,1) = 1.0
        power(:,1) = [ 6.1245,  6.4590,  6.7875,  7.1075,  7.4150,  7.7060,   &
                   &   7.9780,  8.2280,  8.4530,  8.6505,  8.8185,  8.9550,   &
                   &   9.0590,  9.1290,  9.1640,  9.1640,  9.1290,  9.0590,   &
                   &   8.9550,  8.8185,  8.6505,  8.4530,  8.2280,  7.9780,   &
                   &   7.7060,  7.4150,  7.1075,  6.7875,  6.4590,  6.1245 ]
        power = power * 1000.0 * 0.03 * 90 
        
        write(*, *) SUM(power)
        
        call th_power%set_power (nth, geom_th, geom_assm, power, fq_core) 
        
        do ia = 1, nth%na
            write(600, fmt='(1x, I4, *(ES12.5, TR3))') ia, th_power%avg_linear(ia, 1), th_power%max_linear(ia, 1)
        end do 
        
        ! TH input
        design%tmaxcoolout      = 673.0                                         ! k
        design%tmincoolout      = 653.0                                         ! k
        design%tcoolin          = 573.0                                         ! k
        design%tcoolout         = 673.0                                         ! k
        design%is_search        = .FALSE.
        design%is_active_channel = .TRUE.
        
        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------
        ! start calculation
        th_power%max_linear = th_power%avg_linear
        
        call design%set (a_coolant, nth, geom_th, geom_assm, th_power%max_linear)
        
        call Get_flowrate (design, a_coolant, nth, geom_th, geom_assm, th_power%max_linear, is_init=.TRUE.)
        call avg_channel%fix_bottom (a_coolant, nth, design)
        call hot_channel%fix_bottom (a_coolant, nth, design)
        
        do ir = 1, nth%nr
            write(600, fmt='(1x, I4, *(ES12.5, TR3))') ir, design%channel_flowrate(ir), design%assembly_flow(ir)
        end do 
        
        ! average channel calculation to obtain feedback parameter
        call CalcuThermal (design, nth, geom_th, geom_assm, avg_channel, th_power%avg_linear)
        
        do ia = 1, nth%na
            write(600, fmt='(1x, I4, *(ES12.5, TR3))') ia, avg_channel%tcoolant(ia,1), avg_channel%rhocoolant(ia,1), avg_channel%convection(ia,1)
        end do 
        
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, avg_channel, th_power%avg_linear)
        call avg_channel%fix_top (nth, design)
        
        ! hot channel calculation to obtain maximum temperature to testify the design limit
        call CalcuThermal (design, nth, geom_th, geom_assm, hot_channel, th_power%max_linear)
        call RodTemperature (design, nth, geom_th, geom_assm, th_power, hot_channel, th_power%max_linear)
        call hot_channel%fix_top (nth, design)
        
        ! get hot channel inforamtion and feedback information
        call hot_point%set (hot_channel, nth)
        
        
        dst = 0.0
        do ia = 1, nth%na
            dst = dst + 1.0*geom_th%height(ia)
            write(601, fmt = "(1x, I4, TR2, F8.4, TR4, *(F11.4, TR4))")  ia, dst, &
                                                    &   avg_channel%tfuel_center(ia,1), avg_channel%tfuel_surf(ia,1),   &
                                                    &   avg_channel%tclad_surf(ia,1), avg_channel%tcoolant(ia,1)
            dst = dst + 0.0*geom_th%height(ia)
        end do
        
        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------
        ! transient results
        call th_power%print (300)
        call avg_channel%print (nth, geom_th, 400)        
        
        tbeg = 0.0; tend = 10.0; nstep = 50;
        do istep = 1, nstep
            ltime = (istep-1)*(tend-tbeg)/nstep
            ctime = (istep)*(tend-tbeg)/nstep
            step_length = ctime - ltime
            
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
            
            if (istep <= 10)  then
                call th_power%print (300+istep)
                call avg_channel%print (nth, geom_th, 400+istep)
            end if
        end do 
        
        stop('0')
        
    end subroutine Driving_th_check_model
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Prepare_property_pointer ()
        
        integer  :: i_allocate
        
        if (associated(a_coolant_parcs))         deallocate(a_coolant_parcs, stat=i_allocate)
!        if (associated(a_coolant_refprop))       deallocate(a_coolant_refprop, stat=i_allocate)
        if (associated(a_coolant_HLM))           deallocate(a_coolant_HLM, stat=i_allocate)
        if (associated(a_clad_Zr))               deallocate(a_clad_Zr, stat=i_allocate)
        if (associated(a_clad_steels))           deallocate(a_clad_steels, stat=i_allocate)
        if (associated(a_gap_gas))               deallocate(a_gap_gas, stat=i_allocate)
        if (associated(a_fuel_metallic))         deallocate(a_fuel_metallic, stat=i_allocate)
        if (associated(a_fuel_ceramic))          deallocate(a_fuel_ceramic, stat=i_allocate)
        
        allocate(a_coolant_parcs, stat=i_allocate)
!        allocate(a_coolant_refprop, stat=i_allocate)
        allocate(a_coolant_HLM, stat=i_allocate)
        allocate(a_clad_Zr, stat=i_allocate)
        allocate(a_clad_steels, stat=i_allocate)
        allocate(a_gap_gas, stat=i_allocate)
        allocate(a_fuel_metallic, stat=i_allocate)
        allocate(a_fuel_ceramic, stat=i_allocate)
    
    end subroutine Prepare_property_pointer
    
end module th_check_model
