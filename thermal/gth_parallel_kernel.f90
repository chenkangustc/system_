!$
!===================================================================================================
!
!   calculation kernel of parallel model
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Get_flowrate
!                               CalcuThermal
!                               CalcuThermal_transient
!                               RodTemperature
!
!   Public type lists:          No
!
!===================================================================================================
module gth_parallel_kernel

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV

    use lapack_interface
    use stastics
    use exception_header,            only : WarningCollector
    
    use th_global,                   only : a_coolant, a_clad, a_gap, a_fuel, a_clad_Zr, a_clad_steels, a_gap_gas, a_fuel_metallic, a_fuel_ceramic
    use abstract_property_header,    only : CoolantProperty
    use gth_geometry_header,         only : ThermalScale, ThermalGeometry, ThermalAssemblyGeometry
    use gth_thermal_header,          only : ThermalDesign, ThermalChannel
    use gth_power_header,            only : LinearPower
    
    implicit none
    private
    public  :: Get_flowrate, CalcuThermal, CalcuThermal_transient, RodTemperature

    type(WarningCollector)      :: a_warning                                    ! to print warning information when excess define field
    
    real(KREAL), parameter      :: RMS_LIMIT  = 1.0D-6
    integer, parameter          :: ITER_LIMIT = 20

contains
    !$
    !===============================================================================================
    ! get flowrate distribution (pin & assembly)
    !===============================================================================================
    subroutine Get_flowrate (design, a_coolant, nth, geom_th, geom_assm, lnpower, is_init)
        
        type(ThermalDesign), intent(in out)         :: design
        class(CoolantProperty), intent(in), pointer :: a_coolant
        type(ThermalScale), intent(in)              :: nth
        type(ThermalGeometry), intent(in)           :: geom_th
        type(ThermalAssemblyGeometry), intent(in)   :: geom_assm(:)
        real(KREAL), intent(in)                     :: lnpower(:, :)
        logical, intent(in), optional               :: is_init 
        
        ! ----------------------------------------------------------------------
        ! local parameter
        type(ThermalChannel)  :: a_channel
        
        real(KREAL)      :: massflow                                    ! mass flow
        real(KREAL)      :: massflow_density
        real(KREAL)      :: velocity
        real(KREAL)      :: t_out
                         
        real(KREAL)      :: chpower                                     ! channel power density
        integer          :: ia, ir, ityp
        
        logical, allocatable  :: is_satisfied(:)
        real(KREAL)      :: max_clad     = 0.0D0                        ! determine the lower limit
        real(KREAL)      :: max_velocity = 0.0D0                        ! determine the upper limit also
        integer          :: step  = 20
        integer          :: i, j
        integer          :: i_allocate
        
        call a_channel%alloc(nth)
        
        allocate(is_satisfied(nth%nr), stat=i_allocate)
        is_satisfied = .FALSE.
        where (.NOT. design%is_active_channel)
            is_satisfied = .TRUE.
        end where
        
        step = NINT(design%tmaxcoolout - design%tmincoolout)
        search: do i = 1, step + 1
            t_out = design%tmaxcoolout - (design%tmaxcoolout - design%tmincoolout) * (i-1)/step
            
            ! get testing flowrate
            do ir = 1, nth%nr
                ityp = geom_th%geom_type(ir)
                chpower = 0.0D0
                do ia = 1, nth%na
                    chpower = chpower + lnpower(ia, ir)*geom_th%height(ia)
                end do
                
                ! get initial mass flow and mass flow density by outlet temperature
                if (design%is_active_channel(ir) .and. .NOT. is_satisfied(ir))  then
                    massflow = chpower / (a_coolant%get_enthalpy (t_out) - a_coolant%get_enthalpy(design%tcoolin))
                    massflow_density = massflow / geom_assm(ityp)%flowarea
                    
                    design%channel_flowrate(ir) = massflow_density
                    design%assembly_flow(ir) = design%channel_flowrate(ir) * (geom_assm(ityp)%flowarea * geom_assm(ityp)%n_fuelpin)
                end if
            end do
            
            ! testing
            call a_channel%fix_bottom (a_coolant, nth, design)
            call CalcuThermal (design, nth, geom_th, geom_assm, a_channel, lnpower)
            call a_channel%fix_top (nth, design)
            
            if (PRESENT(is_init) .and. is_init)  then
                exit search 
            end if 
            
            do ir = 1, nth%nr
                if (design%is_active_channel(ir) .and. .NOT. is_satisfied(ir))  then
                    max_clad = stastics_max_value (a_channel%tclad_inner(:, ir))
                    max_velocity = stastics_max_value (a_channel%flow_velocity(:, ir))
                    
                    ! selection
                    if (max_velocity > design%max_velocity)  then               ! can not find
                        exit search
                    else if (max_clad < design%tmaxcladsurf)  then              ! this channel ok
                        is_satisfied(ir) = .TRUE.
                    end if
                end if
            end do
            
            if (ALL(is_satisfied))  then
                exit search
            end if
        end do search
        
        if (.NOT. ALL(is_satisfied))  then
!            call a_warning%set (INFO_LIST_FEEDBACK, 'can not search a proper flowrate distribution')
!            call a_warning%print (OUTPUT_UNIT)
!            call a_warning%print (FILES%MAIN)
        end if
        
        if (allocated(is_satisfied))            deallocate(is_satisfied)
        
        ! set initial flow & flowrate
        design%init_flow = design%assembly_flow
        design%init_flowrate = design%channel_flowrate
        
    end subroutine Get_flowrate
    
    !$
    !===============================================================================================
    ! get heat transfered by convection 
    !===============================================================================================
    subroutine CalcuThermal (design, nth, geom_th, geom_assm, a_channel, lnpower)
        
        type(ThermalDesign), intent(in)           :: design
        type(ThermalScale), intent(in)            :: nth        
        type(ThermalGeometry), intent(in)         :: geom_th      
        type(ThermalAssemblyGeometry), intent(in) :: geom_assm(:)        
        type(ThermalChannel), intent(in out)      :: a_channel
        real(KREAL), intent(in)                   :: lnpower(:, :)
        
        ! local parameter
        real(KREAL)      :: sum_hout = 0.0D0
        real(KREAL)      :: sum_flow = 0.0D0
        real(KREAL)      :: hout_avg = 0.0D0
        real(KREAL)      :: tout_avg = 0.0D0
                         
        real(KREAL)      :: massflow                                            ! mass flow
        real(KREAL)      :: t_old = 0.0D0
        real(KREAL)      :: t_new = 0.0D0
                         
        real(KREAL)      :: plusH 
        real(KREAL)      :: rho                                                 ! density
        real(KREAL)      :: Cp                                                  ! heat capacity
        real(KREAL)      :: Cd                                                  ! heat conductivity
        real(KREAL)      :: Nu                                                  ! nusselt number
        real(KREAL)      :: tmid
        integer          :: ia, ir, ityp, niter
        
        ! heat convection between fuel rod and channel
        sum_hout = 0.0; sum_flow = 0.0;
        do ir = 1, nth%nr
            if (design%is_active_channel(ir) )  then
                ityp = geom_th%geom_type(ir)
                
                a_channel%hjunction(0, ir) = a_coolant%get_enthalpy (design%tcoolin)
                massflow = design%channel_flowrate(ir) * geom_assm(ityp)%flowarea
                Cp = a_coolant%get_capacity (design%tcoolin)
                
                associate (pin => geom_assm(ityp)) 
                do ia = nth%na_start, nth%na_end
                    niter = 0
                    iteration: do
                        niter = niter + 1
!                        plusH = lnpower(ia,ir)*geom_th%height(ia) / massflow + a_coolant%get_enthalpy (a_channel%tcoolant(ia-1, ir))
!                        a_channel%tcoolant(ia, ir) = a_coolant%get_temperature (plusH)
                        
                        a_channel%tcoolant(ia, ir) = a_channel%tcoolant(ia-1, ir) + lnpower(ia,ir)*geom_th%height(ia) / (massflow*Cp)
!!!                        tmid = (a_channel%tcoolant(ia, ir) + a_channel%tcoolant(ia-1, ir)) / 2.0
                        tmid = a_channel%tcoolant(ia, ir) 
                        
                        ! heat transfer coefficient
                        rho = a_coolant%get_density (tmid)
                        a_channel%flow_velocity(ia,ir) = massflow / (rho*pin%flowarea)
                        
                        Cp = a_coolant%get_capacity (tmid)
                        Cd = a_coolant%get_conductivity (tmid)
                        Nu = a_coolant%get_nusselt (pin%P2D, a_channel%flow_velocity(ia,ir), pin%Dh, tmid)
                        a_channel%rhocoolant(ia,ir) = rho
                        a_channel%convection(ia, ir) =  Nu * Cd / pin%Dh
                        
                        ! update enthalpy
                        t_new = a_channel%tcoolant(ia, ir)
                        if (niter>=ITER_LIMIT .or. ABS(t_old-t_new)/t_new <= RMS_LIMIT)  then     
                            a_channel%hjunction(ia, ir) = a_coolant%get_enthalpy (a_channel%tcoolant(ia, ir))
                            a_channel%tclad_surf(ia, ir) = a_channel%tcoolant(ia, ir) + lnpower(ia,ir) / (2.0*PI*pin%rod*a_channel%convection(ia, ir))
                            
                            exit iteration
                        end if
                        
                        t_old = t_new
                    end do iteration
                end do
                
                sum_hout = sum_hout + a_channel%hjunction(nth%na_end,ir)*design%channel_flowrate(ir)*pin%flowarea
                sum_flow = sum_flow + design%channel_flowrate(ir)*pin%flowarea
                end associate
            end if
        end do
        
        hout_avg = sum_hout / sum_flow
        tout_avg = a_coolant%get_temperature (hout_avg)
        a_channel%tcoolant_out = tout_avg
        a_channel%rhocoolant_out = a_coolant%get_density (a_channel%tcoolant_out)
        
    end subroutine CalcuThermal
    
    !$
    !===============================================================================================
    ! get heat transfered by convection--in transient
    !===============================================================================================
    subroutine CalcuThermal_transient (design, nth, geom_th, geom_assm, th_power, a_channel, lnpower, step_length)
        
        type(ThermalDesign), intent(in)           :: design
        type(ThermalScale), intent(in)            :: nth        
        type(ThermalGeometry), intent(in)         :: geom_th      
        type(ThermalAssemblyGeometry), intent(in) :: geom_assm(:)  
        type(LinearPower), intent(in)             :: th_power
        type(ThermalChannel), intent(in out)      :: a_channel
        real(KREAL), intent(in)                   :: lnpower(:, :)
        real(KREAL), intent(in)                   :: step_length
        
        ! local parameter
        real(KREAL)      :: sum_hout = 0.0D0
        real(KREAL)      :: sum_flow = 0.0D0
        real(KREAL)      :: hout_avg = 0.0D0
        real(KREAL)      :: tout_avg = 0.0D0
                         
        real(KREAL)      :: massflow                                            ! mass flow
        real(KREAL)      :: gamma_power
        real(KREAL)      :: tmp_t1, tmp_t2
        real(KREAL)      :: coeff_1, coeff_2
        real(KREAL)      :: rhs_1, rhs_2, rhs_3
        real(KREAL)      :: lhs_1, lhs_2, lhs_3
                         
        real(KREAL)      :: rho                                                 ! density
        real(KREAL)      :: Cp                                                  ! heat capacity
        real(KREAL)      :: Cd                                                  ! heat conductivity
        real(KREAL)      :: Nu                                                  ! nusselt number
        real(KREAL)      :: tmid
                         
        real(KREAL)      :: old_h(0:nth%na)                                     ! old value for channel enthalpy
        real(KREAL)      :: old_t(0:nth%na)                                     ! old value for channel bulk temperature
        real(KREAL)      :: new_t(0:nth%na)                                     
        real(KREAL)      :: last_t(0:nth%na)                                    ! coolant temperature for last time point
        real(KREAL)      :: error_iter
        real(KREAL)      :: error_max
        integer          :: ia, ir, ityp, niter
            
        sum_hout = 0.0D0; sum_flow = 0.0D0
        do ir = 1, nth%nr
            if (design%is_active_channel(ir) )  then
                ityp = geom_th%geom_type(ir)
        
                last_t = a_channel%tcoolant(:, ir)
                new_t = a_channel%tcoolant(:, ir)
                niter = 0
                associate (pin => geom_assm(ityp))
                iteration: do 
                    niter = niter + 1
                    old_h = a_channel%hjunction(:, ir)
                    old_t = new_t 
                    new_t(nth%na_start:nth%na_end) = 0.0D0
                    a_channel%hjunction(nth%na_start:nth%na_end, ir)  = 0.0D0
                    a_channel%tcoolant(nth%na_start:nth%na_end, ir)   = 0.0D0
            
                    a_channel%hjunction(0, ir) = a_coolant%get_enthalpy (design%tcoolin)
                    massflow = design%channel_flowrate(ir) * pin%flowarea
                    
                    ! heat transfer between fuel rod and channel   
                    do ia = nth%na_start, nth%na_end
                            
                        ! heat transfer coefficient
!                        tmid = (old_t(ia) + old_t(ia-1)) / 2.0
                        tmid = old_t(ia) 
                        rho = a_coolant%get_density (tmid)
                        a_channel%flow_velocity(ia,ir) = massflow / (rho*pin%flowarea)
                        
                        Cp = a_coolant%get_capacity (tmid)
                        Cd = a_coolant%get_conductivity (tmid)
                        Nu = a_coolant%get_nusselt (pin%P2D, a_channel%flow_velocity(ia,ir), pin%Dh, tmid)
                        a_channel%rhocoolant(ia, ir) = rho
                        a_channel%convection(ia, ir) =  Nu * Cd / pin%Dh
                        
                        ! update temperature, gamma heating
                        gamma_power = lnpower(ia, ir)*th_power%gamma*0.01D0 / (pin%flowarea*rho*Cp)
                        coeff_1 = 2.0D0*PI*pin%rod*a_channel%convection(ia,ir) / (pin%flowarea*rho*Cp)
                        coeff_2 = massflow / (rho*pin%flowarea*geom_th%height(ia))
                        tmp_t1 = a_channel%tclad_surf(ia,ir) - old_t(ia)
                        tmp_t2 = a_coolant%get_temperature (old_h(ia)) - a_coolant%get_temperature (old_h(ia-1))
                        
                        a_channel%tcoolant(ia, ir) = last_t(ia) + step_length * (gamma_power + coeff_1*tmp_t1 - coeff_2*tmp_t2)
                        a_channel%hjunction(ia, ir) = a_coolant%get_enthalpy (a_channel%tcoolant(ia, ir))
                        
                        ! with new 
                        rhs_1 = 2.0D0*PI*pin%rod*a_channel%convection(ia,ir) * a_channel%tclad_surf(ia,ir)
                        rhs_2 = pin%flowarea*rho*Cp * last_t(ia) / step_length
                        rhs_3 = a_channel%flow_velocity(ia,ir)*pin%flowarea*rho*Cp * old_t(ia-1) / geom_th%height(ia)
                        
                        lhs_1 = pin%flowarea*rho*Cp / step_length
                        lhs_2 = a_channel%flow_velocity(ia,ir)*pin%flowarea*rho*Cp / geom_th%height(ia)
                        lhs_3 = 2.0D0*PI*pin%rod*a_channel%convection(ia,ir)
                        
                        new_t(ia) = (rhs_1 + rhs_2 + rhs_3) / (lhs_1 + lhs_2 + lhs_3)
                        a_channel%tcoolant(ia, ir) = new_t(ia)
                        a_channel%hjunction(ia, ir) = a_coolant%get_enthalpy (a_channel%tcoolant(ia, ir))
                    end do
            
                    ! iteration handle
                    error_iter = 0.0D0
                    error_max  = 0.0D0
                    do ia = nth%na_start, nth%na_end
                        error_iter = ABS(a_channel%tcoolant(ia,ir)-old_t(ia)) / old_t(ia)
                        if (error_iter > error_max)  then
                            error_max = error_iter
                        end if
                    end do
                    
                    if (niter>=ITER_LIMIT .or. error_max<RMS_LIMIT)  then
                        exit iteration
                    end if
                end do iteration
                    
                sum_hout = sum_hout + a_channel%hjunction(nth%na_end,ir)*design%channel_flowrate(ir)*pin%flowarea
                sum_flow = sum_flow + design%channel_flowrate(ir)*pin%flowarea
                end associate
            end if
        end do
        
        hout_avg = sum_hout / sum_flow
        tout_avg = a_coolant%get_temperature (hout_avg)
        a_channel%tcoolant_out = tout_avg
        a_channel%rhocoolant_out = a_coolant%get_density (a_channel%tcoolant_out)
        
    end subroutine CalcuThermal_transient
    
    !$
    !===============================================================================================
    ! calculate the temperature distribution of the fuel rod, including fuel,bond and cladding
    !===============================================================================================
    subroutine RodTemperature (design, nth, geom_th, geom_assm, th_power, a_channel, lnpower, step_length)
        
        type(ThermalDesign), intent(in)           :: design
        type(ThermalScale), intent(in)            :: nth        
        type(ThermalGeometry), intent(in)         :: geom_th      
        type(ThermalAssemblyGeometry), intent(in) :: geom_assm(:)  
        type(LinearPower), intent(in)             :: th_power
        type(ThermalChannel), intent(in out)      :: a_channel
        real(KREAL), intent(in)                   :: lnpower(1:nth%na, 1:nth%nr)    ! linear power density
        real(KREAL), intent(in), optional         :: step_length
        
        ! local parameter
        real(KREAL)      :: weight_Pu                                           ! weight fraction for Pu
        real(KREAL)      :: weight_other                                        ! other weight fraction except Pu
                         
        real(KREAL)      :: pellet
        real(KREAL)      :: rho_Cp
        real(KREAL)      :: h_gap, t_gap
        real(KREAL)      :: k_east                                              ! conductivity for east surface
        real(KREAL)      :: k_west                                              ! conductivity for west surface
        real(KREAL)      :: x_1, x_2
        real(KREAL)      :: k_1, k_2
        
        ! for matrix generation
        real(KREAL)      :: none_gamma
        real(KREAL)      :: lhs(nth%n_mesh, nth%n_mesh)
        real(KREAL)      :: rhs(nth%n_mesh)
        real(KREAL)      :: solve(nth%n_mesh)
        real(KREAL)      :: solve_old(nth%n_mesh)
        real(KREAL)      :: error_iter
        real(KREAL)      :: error_max
        
        integer          :: ia, ir, i, j
        integer          :: niter, ityp
        integer          :: i_kind, i_type                                      ! property kind & type
        
        ! ----------------------------------------------------------------------
        ! begin calculation
        none_gamma = 1.0D0 - th_power%gamma*0.01D0
        solve      = 0.0D0
        solve_old  = 0.0D0
        
        do ir = 1, nth%nr
            if (design%is_active_channel(ir) )  then
                ityp = geom_th%geom_type(ir)
                
                ! select cladding
                i_type = MOD(geom_th%cladding_type(ityp), 20)
                i_kind = (geom_th%cladding_type(ityp) - i_type) / 20
                select case(i_kind)
                case(0)
                    a_clad => a_clad_Zr
                case(1)
                    a_clad => a_clad_steels
                case default
                end select
                call a_clad%set (i_type)
                
                ! select gap, only one kind here
                i_type = MOD(geom_th%gap_type(ityp), 20)
                i_kind = (geom_th%gap_type(ityp) - i_type) / 20
                a_gap => a_gap_gas
                call a_gap%set (geom_th%gap_type(ityp))
                
                ! select fuel
                i_type = MOD(geom_th%fuel_type(ityp), 20)
                i_kind = (geom_th%fuel_type(ityp) - i_type) / 20
                select case(i_kind)
                case(0)
                    a_fuel => a_fuel_ceramic
                case(1)
                    a_fuel => a_fuel_metallic
                case default
                end select
                call a_fuel%set (i_type)
                
                do ia = nth%na_start, nth%na_end
                    niter  = 0
                    solve = 800.0D0                                             ! initial guess
                    weight_Pu = design%TRU_weight(ir)
                    weight_other = 1.0D0 - weight_Pu
                    
                    iteration: do
                        niter = niter + 1
                        lhs = 0.0D0
                        rhs = 0.0D0
                        solve_old = solve
                        
                        ! set up matrixs for equation
                        associate (pin => geom_assm(ityp))
                        matrix: do i = 1, nth%n_mesh
                            
                            ! most inner point
                            if (i == 1)  then
                                rhs(i) = -lnpower(ia,ir)*none_gamma / pin%area
                                
                                x_1 = pin%x_surface(i) - pin%x_point(i)
                                x_2 = pin%x_point(i+1) - pin%x_surface(i)
                                k_1 = a_fuel%get_conductivity (solve_old(i))
                                k_2 = a_fuel%get_conductivity (solve_old(i+1))
                                
                                k_east = (pin%x_point(i+1)-pin%x_point(i)) / ((x_1/k_1) + (x_2/k_2))
                                k_west = 0.0D0
                                lhs(i, i+1) = 2.0D0*(pin%x_point(i)+0.5D0*pin%df)*k_east/((pin%x_point(i)+0.25D0*pin%df)*pin%df**2)
                                lhs(i, i)   = - lhs(i, i+1)
                                
                            ! fuel internal point
                            else if (i <= nth%nf)  then
                                rhs(i) = -lnpower(ia,ir)*none_gamma / pin%area
                                
                                x_1 = pin%x_surface(i) - pin%x_point(i)
                                x_2 = pin%x_point(i+1) - pin%x_surface(i)
                                k_1 = a_fuel%get_conductivity (solve_old(i))
                                k_2 = a_fuel%get_conductivity (solve_old(i+1))
                
                                k_east = (pin%x_point(i+1)-pin%x_point(i)) / ((x_1/k_1) + (x_2/k_2))
                                
                                x_1 = pin%x_point(i) - pin%x_surface(i-1)
                                x_2 = pin%x_surface(i-1) - pin%x_point(i-1)
                                k_1 = a_fuel%get_conductivity (solve_old(i))
                                k_2 = a_fuel%get_conductivity (solve_old(i-1))
                                
                                k_west = (pin%x_point(i)-pin%x_point(i-1)) / ((x_1/k_1) + (x_2/k_2))
                                
                                lhs(i, i-1) = (pin%x_point(i)-0.5D0*pin%df)*k_west / (pin%x_point(i)*pin%df**2)
                                lhs(i, i+1) = (pin%x_point(i)+0.5D0*pin%df)*k_east / (pin%x_point(i)*pin%df**2)
                                lhs(i, i)   = - (lhs(i, i-1) + lhs(i, i+1))
                                
                            ! fuel surface point
                            else if (i == nth%nf + 1)  then
                                rhs(i) = -lnpower(ia,ir)*none_gamma / pin%area
                                
                                x_1 = pin%x_point(i) - pin%x_surface(i-1)
                                x_2 = pin%x_surface(i-1) - pin%x_point(i-1)
                                k_1 = a_fuel%get_conductivity (solve_old(i))
                                k_2 = a_fuel%get_conductivity (solve_old(i-1))
                                
                                k_west = (pin%x_point(i)-pin%x_point(i-1)) / ((x_1/k_1) + (x_2/k_2))
                                t_gap  = 0.5D0*solve_old(i)+0.5D0*solve_old(i+1)
                                h_gap  = a_gap%get_transfer (t_gap, pin%pellet, pin%bond, is_inner=.TRUE.)
                                
                                lhs(i, i-1) = (pin%x_point(i)-0.5D0*pin%df)*k_west*2.0D0 / ((pin%x_point(i)-0.25D0*pin%df)*pin%df**2)
                                lhs(i, i+1) = pin%x_point(i)*h_gap*2.0D0 / ((pin%x_point(i)-0.25D0*pin%df)*pin%df)
                                lhs(i, i)   = - (lhs(i, i-1) + lhs(i, i+1))
                                
                            ! clad inner surface point
                            else if (i == nth%nf + 2)  then
                                rhs(i) = 0.0D0
                                
                                x_1 = pin%x_surface(i) - pin%x_point(i)
                                x_2 = pin%x_point(i+1) - pin%x_surface(i)
                                k_1 = a_clad%get_conductivity (solve_old(i))
                                k_2 = a_clad%get_conductivity (solve_old(i+1))
                
                                k_east = (pin%x_point(i+1)-pin%x_point(i)) / ((x_1/k_1) + (x_2/k_2))
                                t_gap  = 0.5D0*solve_old(i-1)+0.5D0*solve_old(i)
                                h_gap  = a_gap%get_transfer (t_gap, pin%pellet, pin%bond, is_inner=.FALSE.)
                                
                                lhs(i, i-1) = 2.0D0*pin%x_point(i)*h_gap / ((pin%x_point(i)+0.25D0*pin%dc)*pin%dc)
                                lhs(i, i+1) = 2.0D0*(pin%x_point(i)+0.5D0*pin%dc)*k_east / ((pin%x_point(i)+0.25D0*pin%dc)*pin%dc**2)
                                lhs(i, i)   = - (lhs(i, i-1) + lhs(i, i+1))
                                
                            ! clad internal point
                            else if (i <= nth%nf + nth%nc + 1)  then
                                rhs(i) = 0.0D0
                                
                                x_1 = pin%x_surface(i) - pin%x_point(i)
                                x_2 = pin%x_point(i+1) - pin%x_surface(i)
                                k_1 = a_clad%get_conductivity (solve_old(i))
                                k_2 = a_clad%get_conductivity (solve_old(i+1))
                
                                k_east = (pin%x_point(i+1)-pin%x_point(i)) / ((x_1/k_1) + (x_2/k_2))
                                
                                x_1 = pin%x_point(i) - pin%x_surface(i-1)
                                x_2 = pin%x_surface(i-1) - pin%x_point(i-1)
                                k_1 = a_clad%get_conductivity (solve_old(i))
                                k_2 = a_clad%get_conductivity (solve_old(i-1))
                                
                                k_west = (pin%x_point(i)-pin%x_point(i-1)) / ((x_1/k_1) + (x_2/k_2))
                                
                                lhs(i, i-1) = (pin%x_point(i)-0.5D0*pin%dc)*k_west / (pin%x_point(i)*pin%dc**2)
                                lhs(i, i+1) = (pin%x_point(i)+0.5D0*pin%dc)*k_east / (pin%x_point(i)*pin%dc**2)
                                lhs(i, i)   = - (lhs(i, i-1) + lhs(i, i+1))
                                
                            ! clad outer surface point
                            else 
                                rhs(i) = -2.0D0*a_channel%convection(ia,ir)*a_channel%tcoolant(ia,ir)*pin%x_point(i) / ((pin%x_point(i)-0.25D0*pin%dc)*pin%dc)
                                
                                x_1 = pin%x_point(i) - pin%x_surface(i-1)
                                x_2 = pin%x_surface(i-1) - pin%x_point(i-1)
                                k_1 = a_clad%get_conductivity (solve_old(i))
                                k_2 = a_clad%get_conductivity (solve_old(i-1))
                                
                                k_west = (pin%x_point(i)-pin%x_point(i-1)) / ((x_1/k_1) + (x_2/k_2))
                                
                                lhs(i, i-1) = 2.0D0*(pin%x_point(i)-0.5D0*pin%dc)*k_west / ((pin%x_point(i)-0.25D0*pin%dc)*pin%dc**2)
                                lhs(i, i)   = -lhs(i, i-1) - 2.0D0*a_channel%convection(ia,ir)*pin%x_point(i)/((pin%x_point(i)-0.25D0*pin%dc)*pin%dc)
                                
                            end if
                        end do matrix
                        end associate
                        
                        ! eleminate rho_Cp
                        do i = 1, nth%n_mesh
                            if (i <= nth%nf+1)  then
                                rho_Cp = a_fuel%get_density (solve_old(i)) * a_fuel%get_capacity (solve_old(i))
                            else
                                rho_Cp = a_clad%get_density (solve_old(i)) * a_clad%get_capacity (solve_old(i))
                            end if
                            
                            rhs(i) = rhs(i) / rho_Cp
                            lhs(i, :) = lhs(i, :) / rho_Cp
                        end do
                        
                        ! NOTE: update matrix for transient
                        if (PRESENT(step_length))  then
                            do i = 1, SIZE(lhs, dim=1)
                                do j = 1, SIZE(lhs, dim=2)
                                    if (i == j)  then
                                        lhs(i, j) = 1.0D0 - step_length * lhs(i, j)
                                    else
                                        lhs(i, j) = - step_length * lhs(i, j)
                                    end if
                                end do
                            end do
                            
                            do i = 1, SIZE(rhs, dim=1)
                                rhs(i) = a_channel%trod(i, ia, ir) - step_length * rhs(i)
                            end do
                        end if
                        
                        ! ----------------------------------------------------------
                        ! solve equation by LAPACK routine
                        call self_DGESV (lhs, rhs)
                        solve = rhs
                        
                        ! iteration update
                        error_iter = 0.0D0
                        error_max  = 0.0D0
                        do i = 1, nth%n_mesh
                            error_iter = ABS(solve(i)-solve_old(i)) / solve_old(i)
                            if (error_iter > error_max)  then
                                error_max = error_iter
                            end if
                        end do
                        
                        if (niter>=ITER_LIMIT .or. error_max<RMS_LIMIT)  then
                            exit iteration
                        end if
                        
                    end do iteration
                    
                    ! update information
                    call a_channel%update (geom_assm(ityp), nth, design, solve, ia, ir)
                                        
                end do
            end if
        end do
        
    end subroutine RodTemperature
        
end module gth_parallel_kernel