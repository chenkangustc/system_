!$
!===================================================================================================
!
!   module for feedback interface
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Check_feedback_steady
!                               Check_feedback_transient
!                               Check_xsec_transient
!
!   Public type lists:          No
!
!===================================================================================================
module feedback
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
	    
    use stastics
    use input_xsec,                 only : read_xsec_unknown, read_xsec_LRA
    use TH2NK_interface_self,       only : Perform_TH_self
	
	!IMPC
    use imp_assm_global
	use TH2NK_interface_imp,        only : Perform_TH_imp
	
    implicit none 
    private
    public  :: Check_feedback_steady, Check_feedback_steady2, Check_feedback_transient, Check_xsec_transient
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_feedback_steady (is_pass)
        
        logical, intent(in out)  :: is_pass
        
        ! local variables
        real(KREAL)  :: power_density(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_core(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_zone(ns%state%zone)

        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq(ns%state%zone, ns%state%layer)
        integer  :: ia, iz
        
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: last 
        real(KREAL)  :: current 
        
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        last = 0.0; current = 0.0;
        
        call dist_power%zone_layer (mesh, geom, .FALSE., power_density)
        call dist_power%fq_zone_layer (mesh, geom, fq_core)
        call dist_power%fq_zone (mesh, geom, fq_zone)
        
        power = 0.0
        fq = 0.0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                power(iz, ia) = power_density(iz, ia) * geom%zone_area(iz) * geom%height(ia)
				!if(ns%feedback%is_feedback .and. ns%feedback%is_inner) imp_pow(iz,ia)=power(iz,ia)
                fq(iz, ia) = fq_core (iz, ia)
            end do
        end do
		
      ! open(1,file='.\output\powDistribution.txt')
      ! write(1,100) power
	  ! 100 Format(F15.5)
      ! close(1) 
	  ! read(*,*)
		 
        is_pass = .FALSE.
        if (ns%feedback%is_model)  then
            select case(ns%feedback%model_name)
            case ('LRA')
                Tfuel = self_lra%tf0
                max_Tfuel = stastics_max_value(Tfuel)
                call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, mask=mat_info%mask_core)
                is_pass = .TRUE.
                return
            end select
        end if
        
        if (ns%feedback%is_inner)  then
            !call Perform_TH_self(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
            call Perform_TH_imp(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
			call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, hot_Tm=max_Tcoolant, hot_Rho_m=min_Rhocoolant, out_Tm=toutlet, mask=mat_info%mask_core)
            call self_fdbk%check (is_pass, .FALSE.)
            call self_fdbk%relax ()
            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param, cr_bank, iter_count%kcritical)
        end if
        
    end subroutine Check_feedback_steady
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_feedback_steady2 (is_pass)
        
        logical, intent(in out)  :: is_pass
        
        ! local variables
        real(KREAL)  :: power_density(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_core(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_zone(ns%state%zone)

        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq(ns%state%zone, ns%state%layer)
        integer  :: ia, iz
        
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: last 
        real(KREAL)  :: current 
        
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        last = 0.0; current = 0.0;
        
        call dist_power%zone_layer (mesh, geom, .FALSE., power_density)
        call dist_power%fq_zone_layer (mesh, geom, fq_core)
        call dist_power%fq_zone (mesh, geom, fq_zone)
        
        power = 0.0
        fq = 0.0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                power(iz, ia) = power_density(iz, ia) * geom%zone_area(iz) * geom%height(ia)
                fq(iz, ia) = fq_core (iz, ia)
            end do
        end do
        
        is_pass = .FALSE.
        if (ns%feedback%is_model)  then
            select case(ns%feedback%model_name)
            case ('LRA')
                Tfuel = self_lra%tf0
                max_Tfuel = stastics_max_value(Tfuel)
                call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, mask=mat_info%mask_core)
                is_pass = .TRUE.
                return
            end select
        end if
        
        if (ns%feedback%is_inner)  then
            call Perform_TH_self(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
            call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, hot_Tm=max_Tcoolant, hot_Rho_m=min_Rhocoolant, out_Tm=toutlet, mask=mat_info%mask_core)
            call self_fdbk%check (is_pass, .FALSE.)
!            call self_fdbk%relax ()
            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param, cr_bank, iter_count%kcritical)
        end if
        
    end subroutine Check_feedback_steady2
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_feedback_transient (is_pass, t_left, t_right)
        
        logical, intent(in out)  :: is_pass
        real(KREAL), intent(in)  :: t_left
        real(KREAL), intent(in)  :: t_right
        
        ! local variables
        real(KREAL)  :: power_density(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_core(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq_zone(ns%state%zone)

        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq(ns%state%zone, ns%state%layer)
        integer  :: ia, iz
        real(KREAL)  :: ltime
        
        logical  :: transient_flag = .TRUE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer)  
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer)  
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer)  
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel  
        real(KREAL)  :: max_Tcoolant  
        real(KREAL)  :: min_Rhocoolant  
        real(KREAL)  :: last  
        real(KREAL)  :: current  
        
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        
        last = t_left
        current = t_right
        
        call dist_power%zone_layer (mesh, geom, .FALSE., power_density)
        call dist_power%fq_zone_layer (mesh, geom, fq_core)
        call dist_power%fq_zone (mesh, geom, fq_zone)
        
        power = 0.0
        fq = 0.0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                power(iz, ia) = power_density(iz, ia) * geom%zone_area(iz) * geom%height(ia)
				!if(ns%feedback%is_feedback .and. ns%feedback%is_inner) imp_pow(iz,ia)=power(iz,ia)
                fq(iz, ia) = fq_core(iz, ia)
            end do
        end do
        
        is_pass = .FALSE.
        if (ns%feedback%is_model)  then
            select case(ns%feedback%model_name)
            case ('LRA')
                call self_lra%update_tf (ns, mesh, geom, mat_info, xsec, flux_forward, current-last, Tfuel)
                max_Tfuel = stastics_max_value(Tfuel)
                call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, mask=mat_info%mask_core)
!                call read_xsec_LRA (self_lra, mat_info, xsec_inp, xsec, param, cr_bank)
!                call cr_bank%map (xsec, param)
                is_pass = .TRUE.
                return 
            end select
        end if
        
        if (ns%feedback%is_inner)  then
            !call Perform_TH_self(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
            call Perform_TH_imp(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
			call self_fdbk%update (geom, Tfuel, Tcoolant, Rhocoolant, hot_Tf=max_Tfuel, hot_Tm=max_Tcoolant, hot_Rho_m=min_Rhocoolant, out_Tm=toutlet, mask=mat_info%mask_core)
            call self_fdbk%check (is_pass, .FALSE.)
!            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param, cr_bank)
!            call cr_bank%map (xsec, param)
        end if
        is_pass = .TRUE.
    
    end subroutine Check_feedback_transient
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_xsec_transient ()
                
        logical  :: is_pass
                
        is_pass = .FALSE.
        if (ns%feedback%is_model)  then
            select case(ns%feedback%model_name)
            case ('LRA')
                call read_xsec_LRA (self_lra, mat_info, xsec_inp, xsec, param, cr_bank)
                call cr_bank%map (xsec, param)
                is_pass = .TRUE.
                return 
            end select
        end if
        
        if (ns%feedback%is_inner)  then
            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param, cr_bank, iter_count%kcritical)
        end if
        is_pass = .TRUE.
    
    end subroutine Check_xsec_transient
    
end module feedback
