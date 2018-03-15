!$
!===================================================================================================
!
!   module for performing perturbation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Perform_perturbation
!
!   Public type lists:          No
!
!===================================================================================================
module perturbation
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use material_header,    only : cross_section_info_tp, kinetics_parameter_info_tp
    
    implicit none 
    private
    public  :: Perform_perturbation, Move_CRWorth
    
    ! convert percentage to acutal value
    real(KREAL), parameter  :: PERCENTAGE_ = 1.0D-2
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Perform_perturbation (tidx, ctime)
    
        integer, intent(in)  :: tidx
        real(KREAL), intent(in)  :: ctime
        
        if (nt%perturb%is_source) then
            call Perturb_external_source(ctime)
        end if
        
        if (nt%perturb%is_xsec) then
            call Perturb_cross_section(ctime)
        end if 
        
        if (nt%perturb%is_CR_move) then
            call Perturb_control_rod(ctime)
        end if
        
        call pert_cr%set_trip (timelist, ctime)
        if (pert_cr%is_trip)  then 
            call Trip_control_rod (ctime)
        end if 
        
!        write(111, fmt='(1x, L, L, *(ES12.5, TR2))') pert_cr%trip%is_active, pert_cr%is_trip, &
!            &   ctime, pert_cr%trip%tstart, timelist%power, ns%flag%rated_power
        
        if (nt%perturb%is_mat) then
            call Perturb_material(ctime)
        end if
        
        call mat_info%set_mask (xsec)
        
    end subroutine Perform_perturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Move_CRWorth (ipert, istep)
        
        integer, intent(in)  :: ipert
        integer, intent(out)  :: istep 
        
        real(KREAL)  :: step(cr_bank%state%n_bank)
        
        istep = worth_cr%getistep(ipert)
        step = cr_bank%init_step
        step(worth_cr%bank) = istep
            
        call cr_bank%move (mat_info, step)
        call cr_bank%homo (xsec_inp, param_inp, geom)
        call cr_bank%map (xsec, param)
        call cr_bank%print_conf(geom, FILES%MAIN)
        
        call mat_info%set_mask (xsec)
    
    end subroutine Move_CRWorth

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! Perform perturbation by external source change directely
    !===============================================================================================
    subroutine Perturb_external_source (ctime)
        
        real(KREAL), intent(in)  :: ctime 
        
        ! local variables
        real(KREAL)  :: Fall, Fpart 
        integer  :: i_pert, ig, ia, iz, ie
        
        ! ----------------------------------------------------------------------
        ! reset Q_ext by initial condition
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                ie = Q_ext%adding(iz, ia)
                ! '0' means no external source
                if (ie == 0)  then
                    Q_ext%matrixs(iz, ia)%intensity = 0.0D0
                else
                    Q_ext%matrixs(iz, ia)%intensity = Q_ext%kinds(ie)%intensity
                end if
            end do
        end do
        
        ! perform perturbation
        do i_pert = 1, pert_q%n_pert
            associate (q => pert_q%perts(i_pert))
            Fall = q%percent*PERCENTAGE_
            Fpart = q%percent*PERCENTAGE_*(ctime-q%time_start)/(q%time_end-q%time_start)
            
            if (ctime <= q%time_start) then
                
            else if (ctime <= q%time_end) then
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (Q_ext%adding(iz,ia) == q%kindID) then
                            do ig = q%ng_start, q%ng_end
                                select case(TRIM(q%type))
                                case ('STEP')
                                    Q_ext%matrixs(iz,ia)%intensity(ig) = Q_ext%matrixs(iz,ia)%intensity(ig) + Q_ext%kinds(q%kindID)%intensity(ig)*Fall
                                case ('RAMP')
                                    Q_ext%matrixs(iz,ia)%intensity(ig) = Q_ext%matrixs(iz,ia)%intensity(ig) + Q_ext%kinds(q%kindID)%intensity(ig)*Fpart
                                case default
                                end select
                            end do
                        end if 
                    end do
                end do
                
            else
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (Q_ext%adding(iz,ia) == q%kindID) then
                            do ig = q%ng_start, q%ng_end
                                select case(TRIM(q%type))
                                case ('STEP', 'RAMP')
                                    Q_ext%matrixs(iz,ia)%intensity(ig) = Q_ext%matrixs(iz,ia)%intensity(ig) + Q_ext%kinds(q%kindID)%intensity(ig)*Fall
                                case default
                                end select
                            end do
                        end if
                    end do
                end do
            end if
            end associate
        end do
        
    end subroutine Perturb_external_source
    
    !$
    !===============================================================================================
    ! Perform perturbation by cross section change
    !===============================================================================================
    subroutine Perturb_cross_section (ctime)
        
        real(KREAL), intent(in)  :: ctime 
        
        ! local varialbles
        real(KREAL)  :: Fall, Fpart         
        integer  :: i_pert, ig, ia, iz, iig, im

        ! ----------------------------------------------------------------------
        ! reset manual perterted xsec by initial condition
        do i_pert = 1, pert_xsec%n_pert
        associate (xsp => pert_xsec%perts(i_pert))
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    if (mat_info%loading(iz,ia) == xsp%matID) then
                        xsec%matrixs(iz,ia) = xsec_init%matrixs(iz,ia)
                    end if 
                end do 
            end do 
        end associate
        end do
        
        ! perform perturbation 
        perturb: do i_pert = 1, pert_xsec%n_pert
            associate (xsp => pert_xsec%perts(i_pert))
            Fall = xsp%percent*PERCENTAGE_
            Fpart = xsp%percent*PERCENTAGE_*(ctime-xsp%time_start)/(xsp%time_end-xsp%time_start)            
            
            ! before perturbation
            if (ctime <= xsp%time_start) then
            
            ! during perturbation
            else if (ctime <= xsp%time_end) then
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (mat_info%loading(iz,ia) == xsp%matID) then
                            do ig = xsp%ng_start, xsp%ng_end
                                select case(TRIM(xsp%xs_type))
                                
                                ! for total cross section
                                case ('SIGMA_T')
                                    select case(TRIM(xsp%type))
                                    case ('STEP')
                                        xsec%matrixs(iz,ia)%sigma_t(ig) = xsec%matrixs(iz,ia)%sigma_t(ig) + xsec_inp%mats(xsp%matID)%sigma_t(ig)*Fall
                                    case ('RAMP')
                                        xsec%matrixs(iz,ia)%sigma_t(ig) = xsec%matrixs(iz,ia)%sigma_t(ig) + xsec_inp%mats(xsp%matID)%sigma_t(ig)*Fpart
                                    end select
                                
                                ! for nu*sigma_f 
                                case ('SIGMA_F_NU')
                                    select case(TRIM(xsp%type))
                                    case ('STEP')
                                        xsec%matrixs(iz,ia)%sigma_f_nu(ig) = xsec%matrixs(iz,ia)%sigma_f_nu(ig) + xsec_inp%mats(xsp%matID)%sigma_f_nu(ig)*Fall
                                        xsec%matrixs(iz,ia)%sigma_f_kappa(ig) = xsec%matrixs(iz,ia)%sigma_f_kappa(ig) + xsec_inp%mats(xsp%matID)%sigma_f_kappa(ig)*Fall
                                    case ('RAMP')
                                        xsec%matrixs(iz,ia)%sigma_f_nu(ig) = xsec%matrixs(iz,ia)%sigma_f_nu(ig) + xsec_inp%mats(xsp%matID)%sigma_f_nu(ig)*Fpart
                                        xsec%matrixs(iz,ia)%sigma_f_kappa(ig) = xsec%matrixs(iz,ia)%sigma_f_kappa(ig) + xsec_inp%mats(xsp%matID)%sigma_f_kappa(ig)*Fpart
                                    end select
                                
                                ! for scatter cross section
                                case ('SIGMA_S')
                                    do iig = xsp%to_ng_start, xsp%to_ng_end
                                        select case(TRIM(xsp%type))
                                        case ('STEP')
                                            xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) = xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) + xsec_inp%mats(xsp%matID)%sigma_s(ig, iig, xsp%scat_index)*Fall
                                        case ('RAMP')
                                            xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) = xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) + xsec_inp%mats(xsp%matID)%sigma_s(ig, iig, xsp%scat_index)*Fpart
                                        end select
                                    end do
            
                                end select
                            end do
                        end if
                    end do
                end do
                
            ! after perturbation
            else 
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (mat_info%loading(iz,ia) == xsp%matID) then
                            do ig = xsp%ng_start, xsp%ng_end
                                select case(TRIM(xsp%xs_type))
                                
                                ! for total cross section
                                case ('SIGMA_T')
                                    select case(TRIM(xsp%type))
                                    case ('STEP', 'RAMP')
                                        xsec%matrixs(iz,ia)%sigma_t(ig) = xsec%matrixs(iz,ia)%sigma_t(ig) + xsec_inp%mats(xsp%matID)%sigma_t(ig)*Fall
                                    end select
                                
                                ! for nu*sigma_f 
                                case ('SIGMA_F_NU')
                                    select case(TRIM(xsp%type))
                                    case ('STEP', 'RAMP')
                                        xsec%matrixs(iz,ia)%sigma_f_nu(ig) = xsec%matrixs(iz,ia)%sigma_f_nu(ig) + xsec_inp%mats(xsp%matID)%sigma_f_nu(ig)*Fall
                                        xsec%matrixs(iz,ia)%sigma_f_kappa(ig) = xsec%matrixs(iz,ia)%sigma_f_kappa(ig) + xsec_inp%mats(xsp%matID)%sigma_f_kappa(ig)*Fall
                                    end select
                                
                                ! for scatter cross section
                                case ('SIGMA_S')
                                    do iig = xsp%to_ng_start, xsp%to_ng_end
                                        select case(TRIM(xsp%type))
                                        case ('STEP', 'RAMP')
                                            xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) = xsec%matrixs(iz,ia)%sigma_s(ig, iig, xsp%scat_index) + xsec_inp%mats(xsp%matID)%sigma_s(ig, iig, xsp%scat_index)*Fall
                                        end select
                                    end do
                                    
                                end select
                            end do
                        end if
                    end do
                end do
            
            end if 
            end associate
        end do perturb
    
    end subroutine Perturb_cross_section
    
    !$
    !===============================================================================================
    ! Perform perturbation by control rod movement
    !===============================================================================================
    subroutine Perturb_control_rod (ctime)
        
        real(KREAL), intent(in)  :: ctime
        
        ! local variables
        real(KREAL)  :: step(cr_bank%state%n_bank)
        real(KREAL)  :: move
        integer  :: i_pert, i_rod, i_bank
        integer  :: ig, ia, iz
    
        ! perform perturbation
        step = cr_bank%init_step
        do i_pert = 1, pert_cr%n_pert
            associate (cr => pert_cr%perts(i_pert))
            if (ctime <= cr%time_start) then
!                select case (TRIM(cr%type))
!                case ('STEP')
!                    step(cr%bankID) = step(cr%bankID)
!                case ('RAMP')
!                    step(cr%bankID) = step(cr%bankID)
!                end select
                
            else if (ctime <= cr%time_end) then
                select case (TRIM(cr%type))
                case ('STEP')
                    step(cr%bankID) = cr%step_end
                case ('RAMP')
                    move = (cr%step_end-cr_bank%init_step(cr%bankID)) * (ctime-cr%time_start) / (cr%time_end-cr%time_start)
                    step(cr%bankID) = step(cr%bankID) + move
                end select
            
            else
                step(cr%bankID) = cr%step_end
            end if 
            end associate
        end do
            
        ! move CR bank
        call cr_bank%move (mat_info, step)
        call cr_bank%homo (xsec_inp, param_inp, geom)
        call cr_bank%map (xsec, param)
        
!        call cr_bank%print_conf(geom, FILES%MAIN)
        
    end subroutine Perturb_control_rod
    
    !$
    !===============================================================================================
    ! Perform perturbation by control rod trip 
    !===============================================================================================
    subroutine Trip_control_rod (t_right)
        
!        type(TransientState), intent(in)         :: nt
!        type(Geometry), intent(in)               :: geom
!        type(MaterialInfo), intent(in out)       :: mat_info
!        type(CorePerturbation), intent(in)       :: pert
!        type(CrossSection), intent(in out)       :: xsec
!        type(KineticsParameter), intent(in out)  :: param
!        type(ControlRodBank), intent(in out)     :: cr_bank
        real(KREAL), intent(in)                  :: t_right
        
        real(KREAL)  :: step(cr_bank%state%n_bank)
        real(KREAL)  :: move
        integer      :: i_pert, i_rod, i_bank
        integer      :: ig, ia, iz
        
        if ((.NOT. pert_cr%trip%is_active) .OR. (t_right <= pert_cr%trip%tstart + pert_cr%trip%tdelay))  then
            return 
        end if
        
        ! move CR bank to "0"
        step = 0 
        call cr_bank%move (mat_info, step)
        call cr_bank%homo (xsec_inp, param_inp, geom)
        call cr_bank%map (xsec, param)
            
    end subroutine Trip_control_rod
    
    !$
    !===============================================================================================
    ! Perform perturbation by material change
    !===============================================================================================
    subroutine Perturb_material (ctime)
        
        real(KREAL), intent(in)  :: ctime
        
        type(cross_section_info_tp)  :: tmp1, tmp2, tmp3, tmp11, tmp12, tmp22
        real(KREAL)  :: FL, FR
        integer  :: i_pert, ig, ia, iz, iig, im
        
        ! ----------------------------------------------------------------------
        ! reset manual perterted xsec by initial condition
        do i_pert = 1, pert_xsec%n_pert
        associate (xsp => pert_xsec%perts(i_pert))
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    if (mat_info%loading(iz,ia) == xsp%matID) then
                        xsec%matrixs(iz,ia) = xsec_init%matrixs(iz,ia)
                    end if 
                end do 
            end do 
        end associate
        end do
        
        ! perform perturbation 
        perturb: do i_pert = 1, pert_mat%n_pert
            associate (matp => pert_mat%perts(i_pert))
            ! before perturbation
            if (ctime <= matp%time_start) then
            
            ! during perturbation
            else if (ctime <= matp%time_end) then
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (mat_info%loading(iz,ia) == matp%mat_beg) then
                            tmp1 = xsec_inp%mats(matp%mat_beg) .MULTI. (matp%per_beg*PERCENTAGE_)           ! left by xsec-1
                            tmp2 = xsec_inp%mats(matp%mat_beg) .MULTI. (1.0D0 - matp%per_beg*PERCENTAGE_)   ! sub by xsec-1
                            tmp3 = xsec_inp%mats(matp%mat_end) .MULTI. (matp%per_end*PERCENTAGE_)           ! add by xsec-2
                            select case(TRIM(matp%type))
                            case('STEP')
                                xsec%matrixs(iz,ia) = tmp1 .ADD. tmp3
                            case('RAMP')
                                ! = LL*(t2-t)/(t2-t1) + RR*(t-t1)/(t2-t1)
                                FL = (ctime-matp%time_start) / (matp%time_end-matp%time_start)
                                FR = (matp%time_end-ctime) / (matp%time_end-matp%time_start)
                                tmp11 = tmp1 
                                tmp12 = tmp2 .MULTI. FR
                                tmp22 = tmp3 .MULTI. FL
                                xsec%matrixs(iz,ia) = tmp11 .ADD. tmp12 .ADD. tmp22
                            case default
                            end select
                        end if
                    end do
                end do
                
            ! after perturbation
            else 
                do ia = 1, ns%state%layer
                    do iz = 1, ns%state%zone
                        if (mat_info%loading(iz,ia) == matp%mat_beg) then
                            tmp1 = xsec_inp%mats(matp%mat_beg) .MULTI. (matp%per_beg*PERCENTAGE_)           ! left by xsec-1
                            tmp2 = xsec_inp%mats(matp%mat_beg) .MULTI. (1.0D0 - matp%per_beg*PERCENTAGE_)   ! sub by xsec-1
                            tmp3 = xsec_inp%mats(matp%mat_end) .MULTI. (matp%per_end*PERCENTAGE_)           ! add by xsec-2
                            select case(TRIM(matp%type))
                            case('STEP')
                                xsec%matrixs(iz,ia) = tmp1 .ADD. tmp3
                            case('RAMP')
                                xsec%matrixs(iz,ia) = tmp1 .ADD. tmp3
                            case default
                            end select
                        end if
                    end do
                end do
            
            end if 
            end associate
        end do perturb
        
        call tmp1%clean ()
        call tmp2%clean ()
        call tmp3%clean ()
        call tmp11%clean ()
        call tmp12%clean ()
        call tmp22%clean ()
        
    end subroutine Perturb_material
    
end module perturbation
