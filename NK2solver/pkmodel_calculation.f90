!$
!===================================================================================================
!
!   point-kinetics model calculatoin
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module pkmodel_calculation

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    use th_global
    
    use stastics,                   only : stastics_max_value
    use timestep_header,            only : TimeStepInfo
    use output_timelist,            only : Print_timelist
    use driver_post_process,        only : Run_post_process
    use TH2NK_interface_self,       only : Perform_TH_self
    use th_pre_process,             only : Driving_th_pre_process
    use th_post_process,            only : Driving_th_post_process
    
    implicit none
    private
    public  :: pkmodel_prepare, pkmodel_run
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine pkmodel_prepare ()

        call pk_param%alloc (nt)
        call pk_solver%alloc (nt)
        call timelist%alloc ()
        
        nt%flag%is_transient = .TRUE.
        
    end subroutine pkmodel_prepare

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine pkmodel_run ()
        
        type(TimeStepInfo)  :: step_Macro, step_Medial
        integer             :: i_Macro, i_Medial
        integer             :: id_Macro, id_Medial
        real(KREAL)         :: ctime, ltime
        real(KREAL)         :: tr_double, tl_double
        integer             :: id_Total
        real(KREAL)         :: rho 
        
        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: toutlet
        real(KREAL)  :: last 
        real(KREAL)  :: current 
        real(KREAL)  :: hactive
        integer      :: ia 
        
        real(KREAL)  :: T_fuel, T_coolant, Rho_coolant
        
        pk_rho%is_feedback = ns%feedback%is_feedback 
        pk_rho%initial = 0.0D0
        
        pk_param%beta = SUM(pk_param%partial_beta)
        pk_param%neutron = 1.0D0
        pk_param%precursor = pk_param%partial_beta / pk_param%partial_lambda / pk_param%generation_time 
        
        id_Total = 0
        id_Macro = 0
        ctime = time_step%get_start ()
        
        ! for feedback 
        if (ns%feedback%is_feedback)  then
            call Driving_th_pre_process ()
            
            transient_flag = .FALSE.
            Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
            max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
            last = 0.0; current = 0.0;
            
            hactive = SUM(geom%height(nth%na_start: nth%na_end))
            power = 0.0 
            do ia = nth%na_start, nth%na_end
                power(:, ia) = geom%height(ia)/hactive * ns%flag%power_level * pk_param%neutron
            end do 
            call Perform_TH_self(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
            
            T_fuel = 0.0; T_coolant = 0.0; Rho_coolant = 0.0;
            do ia = nth%na_start, nth%na_end
                T_fuel = T_fuel + geom%height(ia)/hactive * Tfuel(1,ia)
                T_coolant = T_coolant + geom%height(ia)/hactive * Tcoolant(1,ia)
                Rho_coolant = Rho_coolant + geom%height(ia)/hactive * Rhocoolant(1,ia)
            end do 
            call pk_rho%set_init_fdbk (T_fuel, T_coolant, Rho_coolant)
        end if 
        
        call pkmodel_print (id_Total, ctime, Tfuel, Tcoolant, T_fuel, T_coolant, max_Tfuel, toutlet)
        
        macro: do i_Macro = 1, time_step%get_section_count ()
            id_Macro = id_Macro + 1
            step_Macro = time_step%info (id_Macro, is_step=.FALSE.)    
            
            ! ------------------------------------------------------------------
            ! cycle for micro time step
            id_Medial = 0
            medial: do i_Medial = 1, time_step%step_per_section(i_Macro)
                id_Medial = id_Medial + 1
                id_Total = id_Total + 1
                step_Medial = time_step%info (id_Total, is_step=.TRUE.)  
                
                ctime = step_Medial%right
                tr_double = ctime
                tl_double = step_Medial%left
                
                if (ns%feedback%is_feedback)  then
                    transient_flag = .TRUE.
                    Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
                    max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
                    last = tl_double; current = tr_double;
                    
                    hactive = SUM(geom%height(nth%na_start: nth%na_end))
                    power = 0.0 
                    do ia = nth%na_start, nth%na_end
                        power(:, ia) = geom%height(ia)/hactive * ns%flag%power_level * pk_param%neutron
                    end do 
                    call Perform_TH_self(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
                    
                    T_fuel = 0.0; T_coolant = 0.0; Rho_coolant = 0.0;
                    do ia = nth%na_start, nth%na_end
                        T_fuel = T_fuel + geom%height(ia)/hactive * Tfuel(1,ia)
                        T_coolant = T_coolant + geom%height(ia)/hactive * Tcoolant(1,ia)
                        Rho_coolant = Rho_coolant + geom%height(ia)/hactive * Rhocoolant(1,ia)
                    end do 
                    call pk_rho%set_fdbk (T_fuel, T_coolant, Rho_coolant)
                end if 
                
                call pk_rho%get (tr_double, pk_param%rho)
                call pk_solver%advance (pk_param, left=tl_double, right=tr_double)
                call pkmodel_print (id_Total, ctime, Tfuel, Tcoolant, T_fuel, T_coolant, max_Tfuel, toutlet)
                
            end do medial
        end do macro
            
        call Run_post_process ()
        if (ns%feedback%is_feedback)  then
            call Driving_th_post_process ()
        end if 
            
    end subroutine pkmodel_run

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine pkmodel_print (tidx, tcurrent, Tfuel, Tcoolant, T_fuel, T_coolant, max_Tfuel, toutlet)
        
        integer, intent(in)      :: tidx
        real(KREAL), intent(in)  :: tcurrent 
        real(KREAL), intent(in)  :: Tfuel(:, :)
        real(KREAL), intent(in)  :: Tcoolant(:, :)
        real(KREAL), intent(in)  :: T_fuel
        real(KREAL), intent(in)  :: T_coolant
        real(KREAL), intent(in)  :: max_Tfuel
        real(KREAL), intent(in)  :: toutlet
        
!        ns%flag%rated_power = 1.0
!        timelist%initial_power = 1.0
        timelist%power = pk_param%neutron * ns%flag%power_level 
        timelist%precursor = pk_param%precursor
        timelist%reactivity = pk_param%rho 
        timelist%beta = pk_param%beta 
        
        if (ns%feedback%is_feedback)  then
            timelist%Tm_max = stastics_max_value(Tcoolant)
            timelist%Tm_outlet = toutlet
            timelist%Tm_avg = T_coolant
            timelist%Tf_max = max_Tfuel
            timelist%Tf_avg = T_fuel
        end if 
        
        call Print_timelist (timelist, tidx, tcurrent, FILES%TIMELIST)
    
    end subroutine pkmodel_print

end module pkmodel_calculation 

