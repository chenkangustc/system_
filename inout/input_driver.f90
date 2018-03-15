!$
!===================================================================================================
!
!   driver for input parameter from files
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_input
!
!   Public type lists:          No
!
!===================================================================================================
module input_driver
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use th_global
    use coefficient_iteration
    
    use files_dirs
    use input_keyword
    use input_plain,                only : Driving_plain_scan, Driving_plain_read
    use input_xml,                  only : Driving_xml_scan, Driving_xml_read
    
    use output_main
    use output_runtime
    
    use pkmodel_calculation
    
    implicit none 
    private
    public  :: Driving_input
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_input (file_main)
        
        character(len=*), intent(in)  :: file_main
        
        ! local variables
        character(len=MAX_WORD_LEN)  :: extension
        integer  :: io_error
        integer  :: i, j, k
        
        ! ----------------------------------------------------------------------
        ! set key words
        call Set_section_keyword ()
        call Set_card_keyword ()
        call Set_option_keyword ()
        
        call Set_default_parameter ()
        
        ! ----------------------------------------------------------------------
        call Get_file_extension (file_main, extension)
        
        ! scan input file and get parameter for memory allocate
        if (TRIM(ADJUSTL(extension)) /= EXTEND_XML)  then
            call Driving_plain_scan ()
        
            ! rewind the input file pointor
            rewind(unit=FILES%CASENAME, iostat=io_error)
        else 
!            call Driving_xml_scan (file_main)
        end if
                
        ! pk-model case, 
        if (nt%flag%is_pkmodel)  then
            call pkmodel_prepare ()
            call Allocate_after_scanning ()
            call Driving_plain_read ()
            call Set_deduce_parameter ()
            
            call pkmodel_run ()
            stop ('0')
            
        else
            call Set_deduce_parameter ()
            
            ! allocate memory for input directly-determinate variables
            call Allocate_after_scanning ()
            
            ! ----------------------------------------------------------------------
            ! get the section key word from input file
            if (TRIM(ADJUSTL(extension)) /= EXTEND_XML)  then
                call Driving_plain_read ()
            else
!                call Driving_xml_read (file_main)
            end if
            
            ! allocatable some not-input directely determinate variables
            call Allocate_after_reading ()
            
            ! ----------------------------------------------------------------------        
            ! output memory information
            call Print_memory_size (FILES%MEMORY, ns%flag%case_title)
            
            call Print_case_parsing (FILES%MAIN)
        end if 
        
        call det%phead (FILES%DET, TRIM(ns%flag%case_title))
        
    end subroutine Driving_input
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! set default control parameter
    ! left blank for this set in header file
    !===============================================================================================
    subroutine Set_default_parameter ()
            
        criteria_eigen%error_type          = 2
        criteria_eigen%max_inner           = 10
        criteria_eigen%max_outer           = 500
        criteria_eigen%error_inner_flux    = 5.0E-6
        criteria_eigen%error_outer_flux    = 1.0E-5
        criteria_eigen%error_fission_rate  = 1.0E-5
        criteria_eigen%error_eigen         = 1.0E-5
        
        criteria_fsp%error_type            = 2
        criteria_fsp%max_inner             = 13
        criteria_fsp%max_outer             = 600
        criteria_fsp%error_inner_flux      = 1.0E-6
        criteria_fsp%error_outer_flux      = 5.0E-6
        criteria_fsp%error_fission_rate    = 5.0E-6
        criteria_fsp%error_eigen           = 5.0E-6
            
    end subroutine Set_default_parameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_deduce_parameter ()
    
        ! number of matrix = scatter order + 1(P0)
        ns%deduce%scat_xs = ns%state%scat_order + 1
        ns%deduce%nodal_total = ns%state%nodal * ns%state%layer
        
        if (ns%flag%is_60degree )  then
            ns%deduce%direction = 12 * ns%state%sn                              ! 12 directions per latitude
        else
            ns%deduce%direction = ns%state%sn * (ns%state%sn + 2)
        end if
    
    end subroutine Set_deduce_parameter
    
    !$
    !===============================================================================================
    ! allocatable some input directely determinate variables
    !===============================================================================================
    subroutine Allocate_after_scanning ()
        
        integer  :: i
        integer  :: i_allocate
    
        ! for geometry information
        call geom%alloc ()
        call mesh%alloc ()
        call bound%alloc ()
        call mesh_vtk%alloc ()

        ! for material information
        call mat_info%alloc ()
        call xsec_inp%alloc ()
        call param_inp%alloc ()
        call Q_ext%alloc ()        
        
        ! for quadrature and sweep information
        call quad%alloc ()
        call sweep%alloc ()
        
        ! for perturbation 
        call pert_xsec%alloc ()
        call pert_q%alloc ()
        call pert_cr%alloc ()
        call pert_mat%alloc ()
        call pert_case%alloc ()
        call pert_th%alloc ()
        
        ! for feedback parameter
        call self_fdbk%alloc (ns)
        call det%alloc (ns)
        
        ! for thermal calculation
        call nth%set (ns%state%zone, ns%state%layer, ns%state%layer_top, ns%state%layer_bottom)
        
        call design%alloc (nth)
        call geom_th%alloc (nth)
        call th_power%alloc (nth)

        allocate(geom_assm(nth%n_assm_geom), stat=i_allocate)
        do i = 1, SIZE(geom_assm)
            call geom_assm(i)%alloc (nth)
        end do
    
    end subroutine Allocate_after_scanning
    
    !$
    !===============================================================================================
    ! allocatable some none input directely determinate variables
    !===============================================================================================
    subroutine Allocate_after_reading ()
    
        integer  :: i_allocate
        integer  :: i, ia, iz, im, ie
        
        ! for xsec information per nodal per layer
        call xsec%alloc ()
        call xsec_init%alloc ()
        call xsec_iter%alloc ()
        
        call param%alloc ()
        
        ! for iteration process
        call iter_q%alloc ()
        call iter_flux%alloc ()
        call iter_flux%alloc_bound (bound)
        call iter_count%alloc ()

        call iter_adjoint%alloc ()
        
        call flux_scat%alloc ()
        
        ! for iteration
        call coeff_surface%alloc ()
        call coeff_nodal%alloc ()
        call coeff_source%alloc ()
        
        ! for SOR 
        call accele_SOR%alloc ()
        
        ! for result container        
        call timelist%alloc ()

        call flux_forward%alloc (is_angular=.TRUE.)
        call flux_adjoint%alloc (is_angular=.TRUE.)
        
        call dist_flux%alloc ()
        call dist_power%alloc ()
        call dist_fission_rate%alloc ()
        
        ! ----------------------------------------------------------------------
        ! for perturbation
        if (nt%flag%is_perturb )  then
            call xsec_unpert%alloc ()
            
            call flux_forward_unpert%alloc (is_angular=.TRUE.)
            call flux_adjoint_pt%alloc (is_angular=.TRUE.)
            call flux_adjoint_gpt%alloc (is_angular=.TRUE.)
        end if
        
        ! ----------------------------------------------------------------------
        ! for transient 
        if (nt%flag%is_transient )  then
            call dnp_solver%alloc ()
            call pk_param%alloc (nt)
            call pk_solver%alloc (nt)
            
            allocate(dist_dnps(nt%state%dg), stat=i_allocate)
            do i = 1, SIZE(dist_dnps)
                call dist_dnps(i)%alloc ()
            end do

            call shape_last%alloc ()
            call shape_current%alloc ()
            call shape_predict%alloc ()

            call amplitude%alloc ()
            call pk_parameter%alloc ()
        end if
        
        ! ----------------------------------------------------------------------
        ! set initial value for nodal-information by material-information
        
        ! store cross section information per zone per layer
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = mat_info%loading(iz, ia)
                xsec%matrixs(iz,ia) = xsec_inp%mats(im)
            end do
        end do
    
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = mat_info%loading(iz, ia)
                param%matrixs(iz,ia) = param_inp%mats(im)
            end do
        end do
        
        ! for feedback parameter
        call self_lra%alloc (ns, mat_info, xsec_inp)
        
        ! homogeneous for CR bank
        if (nt%flag%is_CR_rod)  then 
            call cr_bank%set_xsec (geom, mat_info, xsec, param)
            call cr_bank%move (mat_info, cr_bank%init_step)
            call cr_bank%homo (xsec_inp, param_inp, geom)
            call cr_bank%map (xsec, param)
        end if 
        
        call mat_info%set_mask (xsec) 
!        call cr_bank%print_conf(geom, 150)
        
        ! generate matrix distribution
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                ie = Q_ext%adding(iz, ia)
                
                ! '0' means no external source
                if (ie == 0)  then
                    Q_ext%matrixs(iz, ia)%intensity = 0.0
                else
                    Q_ext%matrixs(iz, ia)%intensity = Q_ext%kinds(ie)%intensity
                end if
            end do
        end do
        
!!!        call pert_th%print(790)
    
    end subroutine Allocate_after_reading
    
end module input_driver
