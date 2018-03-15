!$
!===================================================================================================
!
!   print header and tail for output file or standard output unit
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Print_header_main
!                               Print_tail_main
!                               Print_case_parsing
!                               Print_iteration_step_eigen
!                               Print_iteration_step_FSP
!                               Print_iteration_title
!                               Print_iteration_summary
!
!   Public type lists:          No
!
!===================================================================================================
module output_main
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
        
    implicit none
    private
    public  :: Print_header_main, Print_tail_main, Print_case_parsing
    public  :: Print_iteration_title, Print_iteration_summary
    public  :: Print_iteration_step_eigen, Print_iteration_step_FSP

contains
    !$
    !===============================================================================================
    ! Print head information to the output file and the screen
    !===============================================================================================
    subroutine Print_header_main (unit_, case_name)
        
        ! intent parameters
        integer, intent(in)             :: unit_
        character(len=*), intent(in)    :: case_name
        
        ! print the code name and the beginning time
        write(unit=unit_, fmt="(1x, A)") '(Begin)'
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        
        call current_date%print (unit_)
        
!         write(unit=unit_, fmt="(/, *(1x, A, /), /)")                                          &
!             &  '                                =============                                ' ,    &
!             &  '                                  D A I S Y                                  ' ,    &
!             &  '                                =============                                '
        write(unit=unit_, fmt="(/, *(1x, A, /), /)")                                          &
             &  '                               ===================                               ' ,    &
             &  '                                  IMPC-transient                                  ' ,    &
             &  '                               ===================                                '
        
        ! print the copyright and basic information
!        write(unit=unit_, fmt="(1x, 12x, '   ______________________________________________________')")
!        write(unit=unit_, fmt="(1x, 12x, '   under Developping by NECP laboratory')")
!        write(unit=unit_, fmt="(1x, 12x, '                     at Xi''an Jiaotong University')")
        write(unit=unit_, fmt="(1x, 12x, '   ------------------------------------------------------')")
        write(unit=unit_, fmt="(1x, 12x, 3x, A, A, I1, '.', I1, I1)") TRIM(CODENAME), ': ver ', VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
        write(unit=unit_, fmt="(1x, 12x, 3x, A)") TRIM(COPYRIGHT)
        write(unit=unit_, fmt="(1x, 12x, 3x, A)") TRIM(LABORATORY)
		write(unit=unit_, fmt="(1x, 12x, 3x, A)") TRIM(IMP)
        write(unit=unit_, fmt="(1x, 12x, '   ------------------------------------------------------')")
    
        ! print the case title
        write(unit=unit_, fmt="(/)")
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        write(unit=unit_, fmt="(1x, A, 1x)", advance='no') 'Case name is:'
        write(unit=unit_, fmt="(1x, A)") TRIM(case_name)
        write(unit=unit_, fmt="(/)")
        
    end subroutine Print_header_main
    
    !$
    !===============================================================================================
    !    Print tail information to the output file and the screen
    !===============================================================================================
    subroutine Print_tail_main (unit_, case_name)
        
        ! intent parameters
        integer, intent(in)           :: unit_
        character(len=*), intent(in)  :: case_name
    
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        
!        call time_program%stop ()
        call cputime_program%stop ()
        
        write(unit=unit_, fmt="(1x, A, 1x, '[', A, ']', /)")  'End of calculation of case:', TRIM(case_name)
        
        call current_date%print (unit_)
        
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        
        ! write the end mark
        write(unit=unit_, fmt="(1x, A, /)")  '(End)'
        
    end subroutine Print_tail_main
    
    !$
    !===============================================================================================
    ! parse the case information based on the input
    !===============================================================================================
    subroutine Print_case_parsing (unit_)
        
        ! intent parameters
        integer, intent(in)     :: unit_
        
        integer, parameter  :: REVERSE_ORDER  = -1
        integer  :: ia, iz, it, i_pert
        real(KREAL)   :: distance
        
        ! ----------------------------------------------------------------------
        ! parsing title
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        write(unit=unit_, fmt="(2X, A)") 'Parsing the input files:'
        
        ! ----------------------------------------------------------------------
        ! problem state
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Problem state:')")
        
        write(unit=unit_, fmt="(2X, 'Total energy group number is                :', TR5, I6)")  ns%state%ng
        write(unit=unit_, fmt="(2X, 'Energy group number with upscatter is       :', TR5, I6)")  ns%state%upscatter
        
        write(unit=unit_, fmt="(/)")  
        write(unit=unit_, fmt="(2X, 'Nodal point number per layer is             :', TR5, I6)")  ns%state%point
        write(unit=unit_, fmt="(2X, 'Nodal number per layer is                   :', TR5, I6)")  ns%state%nodal
        write(unit=unit_, fmt="(2X, 'Total axial layer divided is                :', TR5, I6)")  ns%state%layer
        write(unit=unit_, fmt="(2X, 'Number of boundary segment is               :', TR5, I6)")  ns%state%segment
        
        write(unit=unit_, fmt="(/)")                                                                                
        write(unit=unit_, fmt="(2X, 'Number of material zone per layer is        :', TR5, I6)")  ns%state%zone
        write(unit=unit_, fmt="(2X, 'Total material number is                    :', TR5, I6)")  ns%state%mat
        write(unit=unit_, fmt="(2X, 'Total external source number is             :', TR5, I6)")  ns%state%source
        write(unit=unit_, fmt="(2X, 'Subscript if SN is                          :', TR5, I6)")  ns%state%sn
        write(unit=unit_, fmt="(2X, 'Anisotropic scatter order is                :', TR5, I6)")  ns%state%scat_order
        
        write(unit=unit_, fmt="(/)")                                                                                
        write(unit=unit_, fmt="(2X, 'Number of scatter matrix is                 :', TR5, I6)")  ns%deduce%scat_xs
        write(unit=unit_, fmt="(2X, 'Total number of nodal is                    :', TR5, I6)")  ns%deduce%nodal_total
        write(unit=unit_, fmt="(2X, 'Total number of ordinate is                 :', TR5, I6)")  ns%deduce%direction
        
        write(unit=unit_, fmt="(2x, 'Number of delay neutron precursor group is  :', TR5, I6)")  nt%state%dg
        
        ! ----------------------------------------------------------------------
        ! problem flag
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Problem flag:')")
        
        write(unit=unit_, fmt="(2x, 'The case name is                            :', TR5, A)") TRIM(ns%flag%case_title)
        write(unit=unit_, fmt="(2X, 'The initial power level is                  :', TR5, ES11.4, ' Watt')")  ns%flag%power_level
        
        if (ns%flag%is_eigen ) then
            write(unit=unit_, fmt="(2X, 'This case is eigenvalue problem             ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'This case is eigenvalue problem             ?', TR5, A)")  '[No]'
        end if
        
        if (ns%flag%is_bevel_edge ) then
            write(unit=unit_, fmt="(2X, 'This case has bevel edge                    ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'This case has bevel edge                    ?', TR5, A)")  '[No]'
        end if
        
        if (ns%flag%is_60degree ) then
            write(unit=unit_, fmt="(2X, 'This case use 60 degree symmetry            ?', TR5, A)")  '[Yes]'
            write(unit=unit_, fmt="(2X, 'The number of 60 degree sector is           :', TR5, I4)") ns%flag%n_theta / 60
        else
            write(unit=unit_, fmt="(2X, 'This case use 60 degree symmetry            ?', TR5, A)")  '[No]'
        end if
        
        if (ns%flag%is_square ) then
            write(unit=unit_, fmt="(2X, 'The assembly is regular square              ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'The assembly is regular square              ?', TR5, A)")  '[No]'
        end if
        
        if (ns%flag%is_hexagonal ) then
            write(unit=unit_, fmt="(2X, 'The assembly is regular hexagonal           ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'The assembly is regular hexagonal           ?', TR5, A)")  '[No]'
        end if
        
        if (nt%flag%is_transient ) then
            write(unit=unit_, fmt="(2X, 'Perform transient simulation                ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'Perform transient simulation                ?', TR5, A)")  '[No]'
        end if
        
        if (nt%flag%is_adjoint ) then
            write(unit=unit_, fmt="(2X, 'Perform adjoint simulation                  ?', TR5, A)")  '[Yes]'
            select case(nt%flag%adjoint_type)
            case ('HOMOGENEOUS')
                write(unit=unit_, fmt="(2X, 'The adjoint type is                         :', TR5, A)")  '[Homogeneous]'
            case ('ONE')
                write(unit=unit_, fmt="(2X, 'The adjoint type is                         :', TR5, A)")  '[One]'
            case ('FISSION_TIMES')
                write(unit=unit_, fmt="(2X, 'The adjoint type is                         :', TR5, A)")  '[sigma_f]'
            case ('FISSION_NEUTRON')
                write(unit=unit_, fmt="(2X, 'The adjoint type is                         :', TR5, A)")  '[nu * sigma_f]'
            case ('FISSION_ENERGY')
                write(unit=unit_, fmt="(2X, 'The adjoint type is                         :', TR5, A)")  '[kappa * sigma_f]'
            end select
        else
            write(unit=unit_, fmt="(2X, 'Perform adjoint simulation                  ?', TR5, A)")  '[No]'
        end if
        
        if (nt%flag%is_CR_rod ) then
            write(unit=unit_, fmt="(2X, 'This case has control rod                   ?', TR5, A)")  '[Yes]'
            if (nt%flag%is_decusping ) then
                write(unit=unit_, fmt="(2X, 'Perform decusping for control rod           ?', TR5, A)")  '[Yes]'
            else
                write(unit=unit_, fmt="(2X, 'Perform decusping for control rod           ?', TR5, A)")  '[No]'
            end if
        else
            write(unit=unit_, fmt="(2X, 'This case has control rod                   ?', TR5, A)")  '[No]'
        end if
        
        ! ----------------------------------------------------------------------
        ! method selection
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Method selection:')")
        
        if (ns%method%is_LW ) then
            write(unit=unit_, fmt="(2X, 'This case use iteration extrapolation       ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'This case use iteration extrapolation       ?', TR5, A)")  '[No]'
        end if
        
        if (ns%method%is_upscatter_cycle ) then
            write(unit=unit_, fmt="(2X, 'This case perform upscatter cycle           ?', TR5, A)")  '[Yes]'
            write(unit=unit_, fmt="(2X, 'Number of upscatter cycle per iteration     :', TR5, I4)")  ns%method%n_upscatter_cycle
        else
            write(unit=unit_, fmt="(2X, 'This case perform upscatter cycle           ?', TR5, A)")  '[No]'
        end if

        if (nt%flag%is_transient ) then
            select case(TRIM(nt%method%scheme)) 
            case ('THETA')
                write(unit=unit_, fmt="(2X, 'The transient method is                     :', TR5, A)")  '[THETA]'
                if (nt%method%is_extrapolation ) then
                    write(unit=unit_, fmt="(2X, 'Perform transient flux extrapoltion         ?', TR5, A)")  '[Yes]'
                else 
                    write(unit=unit_, fmt="(2X, 'Perform transient flux extrapoltion         ?', TR5, A)")  '[No]'
                end if
                
            case ('PCQS')
                write(unit=unit_, fmt="(2X, 'The transient method is                     :', TR5, A)")  '[PCQS]'
                if (nt%method%is_extrapolation ) then
                    write(unit=unit_, fmt="(2X, 'Perform transient flux extrapoltion         ?', TR5, A)")  '[Yes]'
                else 
                    write(unit=unit_, fmt="(2X, 'Perform transient flux extrapoltion         ?', TR5, A)")  '[No]'
                end if
                
            case ('PK')
                write(unit=unit_, fmt="(2X, 'The transient method is                     :', TR5, A)")  '[PK]'
                
            case ('SCM')
                write(unit=unit_, fmt="(2X, 'The transient method is                     :', TR5, A)")  '[SCM]'
                
            case ('GIQS')
                write(unit=unit_, fmt="(2X, 'The transient method is                     :', TR5, A)")  '[GIQS]'
            
            end select
        end if
        
        ! ----------------------------------------------------------------------
        ! iteration convergence criteria
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        
        write(unit=unit_, fmt="(   2x, 'eigenvalue convergence criteria:')")
        write(unit=unit_, fmt="(2X, 'Max iteration number for inner iteration    :', TR11, I5)") criteria_eigen%max_inner
        write(unit=unit_, fmt="(2X, 'Max iteration number for outer iteration    :', TR11, I5)") criteria_eigen%max_outer        
        write(unit=unit_, fmt="(2X, 'Relative flux error for inner iteration     :', TR5, ES11.4)") criteria_eigen%error_inner_flux
        write(unit=unit_, fmt="(2X, 'Relative flux error for outer iteration     :', TR5, ES11.4)") criteria_eigen%error_outer_flux
        write(unit=unit_, fmt="(2X, 'Relative flux error for fission rate        :', TR5, ES11.4)") criteria_eigen%error_fission_rate
        write(unit=unit_, fmt="(2X, 'Relative flux error for eigenvalue          :', TR5, ES11.4)") criteria_eigen%error_eigen
        
        write(unit=unit_, fmt="(   2x, 'fsp convergence criteria:')")
        write(unit=unit_, fmt="(2X, 'Max iteration number for inner iteration    :', TR11, I5)") criteria_fsp%max_inner
        write(unit=unit_, fmt="(2X, 'Max iteration number for outer iteration    :', TR11, I5)") criteria_fsp%max_outer        
        write(unit=unit_, fmt="(2X, 'Relative flux error for inner iteration     :', TR5, ES11.4)") criteria_fsp%error_inner_flux
        write(unit=unit_, fmt="(2X, 'Relative flux error for outer iteration     :', TR5, ES11.4)") criteria_fsp%error_outer_flux
        write(unit=unit_, fmt="(2X, 'Relative flux error for fission rate        :', TR5, ES11.4)") criteria_fsp%error_fission_rate
        write(unit=unit_, fmt="(2X, 'Relative flux error for eigenvalue          :', TR5, ES11.4)") criteria_fsp%error_eigen
        
        ! ----------------------------------------------------------------------
        ! output selection
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Steady output selection:')")
        
        if (ns%output%is_HDF5 ) then
            write(unit=unit_, fmt="(2X, 'Output HDF5 format distribution             ?', TR5, A)")  '[Yes]'
        else 
            write(unit=unit_, fmt="(2X, 'Output HDF5 format distribution             ?', TR5, A)")  '[No]'
        end if
        
        if (ns%output%is_vtk ) then
            write(unit=unit_, fmt="(2X, 'Output vtk file for visualization           ?', TR5, A)")  '[Yes]'
        else
            write(unit=unit_, fmt="(2X, 'Output vtk file for visualization           ?', TR5, A)")  '[No]'
        end if
        
        ! ----------------------------------------------------------------------
        ! feedback control
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Feedback control:')")
        
        if (ns%feedback%is_feedback ) then
            write(unit=unit_, fmt="(2X, 'Transient state consider feedback           ?', TR5, A)")  '[Yes]'
            if (ns%feedback%is_model ) then
                write(unit=unit_, fmt="(2X, 'Feedback is considered by                   :', TR5, A)")  '[Model TH]'
                write(unit=unit_, fmt="(2x, 'Feedback model name is                      :', TR5, A)")  ns%feedback%model_name
            end if 
            if (ns%feedback%is_inner ) then
                write(unit=unit_, fmt="(2X, 'Feedback is considered by                   :', TR5, A)")  '[Inner TH]'
            end if 
            if (ns%feedback%is_coupled ) then
                write(unit=unit_, fmt="(2X, 'Feedback is considered by                   :', TR5, A)")  '[Couple TH]'
            end if
        else
            write(unit=unit_, fmt="(2X, 'Transient state consider feedback           ?', TR5, A)")  '[No]'
        end if
        
        ! ----------------------------------------------------------------------
        ! material configure per zone per layer
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2X, 'Material configuration (Top --> Bottom):')") 
        write(unit=unit_, fmt="(2X, 'Layer index:', TR4, 'Height(cm):', TR4, 'Distance to bottom(cm):', TR4, 'Material ID per zone:')")
        write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
        ! output by layer reversing
        distance = SUM(geom%height)
        do ia = ns%state%layer, 1, REVERSE_ORDER
            distance = distance - geom%height(ia) / 2.0
            write(unit=unit_, fmt="(2X, I4, TR11, F8.4, TR8, F8.4, TR18, *(I4, TR2))") ia, geom%height(ia), distance, (mat_info%loading(iz,ia), iz=1, ns%state%zone)
            distance = distance - geom%height(ia) / 2.0
        end do
        write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
        
        ! ----------------------------------------------------------------------
        ! external source configure per zone per layer
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2X, 'External source configuration (Top --> Bottom):')") 
        write(unit=unit_, fmt="(2X, 'Layer index:', TR4, 'Height(cm):', TR4, 'Distance to bottom(cm):', TR4, 'Source ID per zone:')")
        write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
        ! output by layer reversing
        distance = SUM(geom%height)
        do ia = ns%state%layer, 1, REVERSE_ORDER
            distance = distance - geom%height(ia) / 2.0
            write(unit=unit_, fmt="(2X, I4, TR11, F8.4, TR8, F8.4, TR18, *(I4, TR2))") ia, geom%height(ia), distance, (Q_ext%adding(iz,ia), iz=1, ns%state%zone)
            distance = distance - geom%height(ia) / 2.0
        end do
        write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
        
        ! ----------------------------------------------------------------------
        ! control rod information
        if (nt%flag%is_CR_rod ) then
            write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
            write(unit=unit_, fmt="(   2X, 'Control rod configuration information:')") 
            
            write(unit=unit_, fmt="(2X, 'Number of control rod is                    :', TR5, I3)")  cr_bank%state%n_rod
            write(unit=unit_, fmt="(2X, 'Number of control rod bank is               :', TR5, I3)")  cr_bank%state%n_bank
            write(unit=unit_, fmt="(2X, 'Number of upper reflector layer is          :', TR5, I3)")  cr_bank%state%n_top
            write(unit=unit_, fmt="(2X, 'Number of lower reflector layer is          :', TR5, I3)")  cr_bank%state%n_bottom
            
            write(unit=unit_, fmt="(/)")
            write(unit=unit_, fmt="(2X, 'Step size for control rod is                :', TR5, F8.4)")  cr_bank%state%step_size
            write(unit=unit_, fmt="(2X, 'Minimum step for control rod is             :', TR5, F6.1)")  cr_bank%state%min_step
            write(unit=unit_, fmt="(2X, 'Maximum step for control rod is             :', TR5, F6.1)")  cr_bank%state%max_step
            
            ! 
            ! output by layer reversing
            write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
            write(unit=unit_, fmt="(   2X, 'Control rod configuration:')") 
            
            write(unit=unit_, fmt="(/)")
            write(unit=unit_, fmt="(2X, TR15, 'Item:', TR21, 'Value per control rod:')")
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
            
            ! rod IDs
            write(unit=unit_, fmt="(2X, A20, TR10)", advance='no')  'Rod index --'
            do iz = 1, cr_bank%state%n_rod
                if (iz < cr_bank%state%n_rod) then
                    write(unit=unit_, fmt="(TR8, I4)", advance='no')  iz
                else 
                    write(unit=unit_, fmt="(TR8, I4)")  iz
                end if
            end do
            
            ! bank IDs
            write(unit=unit_, fmt="(2X, A20, TR10)", advance='no')  'Bank index --'
            do iz = 1, cr_bank%state%n_rod
                if (iz < cr_bank%state%n_rod) then
                    write(unit=unit_, fmt="(TR8, I4)", advance='no')  cr_bank%rods(iz)%bank
                else 
                    write(unit=unit_, fmt="(TR8, I4)")  cr_bank%rods(iz)%bank
                end if
            end do
            
            ! control bank zone
            write(unit=unit_, fmt="(2X, A20, TR10)", advance='no')  'Zone index --'
            do iz = 1, cr_bank%state%n_rod
                if (iz < cr_bank%state%n_rod) then
                    write(unit=unit_, fmt="(TR8, I4)", advance='no')  cr_bank%rods(iz)%zone
                else 
                    write(unit=unit_, fmt="(TR8, I4)")  cr_bank%rods(iz)%zone
                end if
            end do
            
            ! half insert
            write(unit=unit_, fmt="(2X, A20, TR10)", advance='no')  'Half insert --'
            do iz = 1, cr_bank%state%n_rod
                if (iz < cr_bank%state%n_rod) then
                    write(unit=unit_, fmt="(TR10, L2)", advance='no')  cr_bank%rods(iz)%is_gray
                else 
                    write(unit=unit_, fmt="(TR10, L2)")  cr_bank%rods(iz)%is_gray
                end if
            end do
            
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------')")
            
            ! output by layer reversing
            call cr_bank%print_conf (geom, unit_)
            
        end if
        
        ! ----------------------------------------------------------------------
        ! boundary condition
        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
        write(unit=unit_, fmt="(   2x, 'Boundary condition:')")
        
        write(unit=unit_, fmt="(/, 2X, 'Boundary condition for axial:')")
        if (ABS(bound%axial(1)-bound%VACCUM) < EPS_ZERO)  then
            write(unit=unit_, fmt="(2X, 'Lower axial boundary is                     :', TR5, A)")   '[Vaccum]'
        else if (ABS(bound%axial(1)-bound%REFLECT) < EPS_ZERO)  then
            write(unit=unit_, fmt="(2X, 'Lower axial boundary is                     :', TR5, A)")   '[Reflect]'
        end if 
        
        if (ABS(bound%axial(2)-bound%VACCUM) < EPS_ZERO)  then
            write(unit=unit_, fmt="(2X, 'Upper axial boundary is                     :', TR5, A)")   '[Vaccum]'
        else if (ABS(bound%axial(2)-bound%REFLECT) < EPS_ZERO)  then
            write(unit=unit_, fmt="(2X, 'Upper axial boundary is                     :', TR5, A)")   '[Reflect]'
        end if 
        
        write(unit=unit_, fmt="(/, 2X, 'Boundary condition for radial:')")
        write(unit=unit_, fmt="(2X, 'Radial boundary is splited into             :', TR5, I4)")  ns%state%segment
        write(unit=unit_, fmt="(2X, 'Boundary ID:', TR4, 'Boundary condition:')")
        write(unit=unit_, fmt="(2X, '------------------------------------')")
        do ia = 1, ns%state%segment
            if (ABS(bound%radial(ia)-bound%VACCUM) < EPS_ZERO)  then
                write(unit=unit_, fmt="(2X, I4, TR12, A)")  ia, '[Vaccum]'
            else if (ABS(bound%radial(ia)-bound%REFLECT) < EPS_ZERO)  then 
                write(unit=unit_, fmt="(2X, I4, TR12, A)")  ia, '[Reflect]'
            end if 
        end do
        write(unit=unit_, fmt="(2X, '------------------------------------')")
        
!        ! ----------------------------------------------------------------------  
!        ! time section and step
!        write(unit=unit_, fmt="(/, 2X, '_____________________________________________')") 
!        write(unit=unit_, fmt="(   2x, 'Transient TimeStep:')")
!        
!        write(unit=unit_, fmt="(2X, 'Start time of transient is                  :', TR5, F8.4, ' second')")  time_step%get_start ()
!        write(unit=unit_, fmt="(2X, 'Final time of transient is                  :', TR5, F8.4, ' second')")  time_step%get_end ()
!        write(unit=unit_, fmt="(2X, 'Time section number is                      :', TR8, I5)")  time_step%get_section_count ()
!        write(unit=unit_, fmt="(2X, 'TimeStep number is                         :', TR8, I5)")  time_step%get_step_count ()
!        
!        write(unit=unit_, fmt="(/, 2X, 'Section ID:', TR4, 'Start time:', TR4, 'End time:', TR4, 'Step number:', TR4, 'Step length:')")
!        write(unit=unit_, fmt="(2X, '----------------------------------------------------------------------')")
!        do it = 1, time_step%get_section_count ()
!            if (it == 1) then
!                write (unit=unit_, fmt="(2X, I5, TR9, F8.4, TR7, F8.4, TR5, I5, TR11, F8.4)")      &
!                    &  it, time_step%start_time, time_step%sections(it)%point, time_step%step_per_section(it), time_step%sections(it)%length
!            else
!                write (unit=unit_, fmt="(2X, I5, TR9, F8.4, TR7, F8.4, TR5, I5, TR11, F8.4)")      &
!                    &  it, time_step%sections(it-1)%point, time_step%sections(it)%point, time_step%step_per_section(it), time_step%sections(it)%length
!            end if
!        end do
!        write(unit=unit_, fmt="(2X, '----------------------------------------------------------------------')")
            
        ! ------------------------------------------------------------------
        ! perturbation details
        ! for cross section perturbation directly
        if (nt%perturb%is_xsec ) then
            write(unit=unit_, fmt="(/, 2x, '_____________________________________________')") 
            write(unit=unit_, fmt="(   2X, 'Perturb by Cross section directly:', /)")
            
            write(unit=unit_, fmt="(2X, 'Number of cross section perturbation is     :', TR5, I3)") pert_xsec%n_pert
            write(unit=unit_, fmt="(2X, 'Perturb ID:', TR4, 'Matrial ID:', TR4, 'Energy group start:', TR4, 'Energy group end:', TR4,   &
                &   'Time start:', TR4, 'Time end:', TR4, 'Cross section type:', TR4, 'Perturb type:', TR4, 'Perturb perentage:')")
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '------------------------------------------------------------------')")
            do i_pert = 1, pert_xsec%n_pert
                associate (xsp => pert_xsec%perts(i_pert))
                    select case(TRIM(xsp%xs_type))
                    
                    ! for total cross section
                    case ('SIGMA_T')
                        select case(TRIM(xsp%type))
                        case ('STEP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'TOTAL XS', 'STEP', xsp%percent, ' %'
                        case ('RAMP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'TOTAL XS', 'RAMP', xsp%percent, ' %'
                        end select
                    
                    ! for nu multi fission cross section
                    case ('SIGMA_F_NU')
                        select case(TRIM(xsp%type))
                        case ('STEP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'FISSION XS', 'STEP', xsp%percent, ' %'
                        case ('RAMP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'FISSION XS', 'RAMP', xsp%percent, ' %'
                        end select
                        
                    ! for scatter cross section
                    case ('SIGMA_S')
                        select case(TRIM(xsp%type))
                        case ('STEP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'SCATTER XS', 'STEP', xsp%percent, ' %'
                        case ('RAMP')
                            write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A12, TR11, A, TR10, F7.2, A)")     &
                                &  i_pert, xsp%matID, xsp%ng_start, xsp%ng_end, xsp%time_start, xsp%time_end, 'SCATTER XS', 'RAMP', xsp%percent, ' %'
                        end select
                        
                    end select
                    
                end associate
            end do
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '------------------------------------------------------------------')")
        end if
        
        ! for external source perturbation
        if (nt%perturb%is_source ) then
            write(unit=unit_, fmt="(/, 2x, '_____________________________________________')") 
            write(unit=unit_, fmt="(   2X, 'Perturb by external source:', /)")
            
            write(unit=unit_, fmt="(2X, 'Number of external source perturbation is   :', TR5, I3)") pert_q%n_pert
            write(unit=unit_, fmt="(2X, 'Perturb ID:', TR4, 'Matrial ID:', TR4, 'Energy group start:', TR4, 'Energy group end:', TR4,   &
                &   'Time start:', TR4, 'Time end:', TR4, 'Perturb type:', TR4, 'Perturb perentage:')")
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '------------------------------------------')")
            do i_pert = 1, pert_q%n_pert
                associate (source => pert_q%perts(i_pert))

                    select case(TRIM(source%type))
                    case ('step', 'STEP')
                        write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A, TR11, F7.2, A)")     &
                            &  i_pert, source%kindID, source%ng_start, source%ng_end, source%time_start, source%time_end, 'STEP', source%percent, ' %'
                    case ('ramp', 'RAMP')
                        write(unit=unit_, fmt="(2X, I3, TR9, I4, TR11, I4, TR20, I4, TR18, F8.4, TR7, F8.4, TR7, A, TR11, F7.2, A)")     &
                            &  i_pert, source%kindID, source%ng_start, source%ng_end, source%time_start, source%time_end, 'RAMP', source%percent, ' %'
                    end select
                    
                end associate      
            end do
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '------------------------------------------')")
        end if
        
        ! for control rod moving
        if (nt%perturb%is_CR_move) then
        
            write(unit=unit_, fmt="(/, 2x, '_____________________________________________')") 
            write(unit=unit_, fmt="(   2X, 'Perturb by Control rod moving:', /)")
            
            write(unit=unit_, fmt="(2X, 'Number of control bank moving is   :', TR5, I3)") pert_cr%n_pert
            write(unit=unit_, fmt="(2X, 'Perturb ID:', TR4, 'CR bank ID:', TR4, 'Time start:', TR8, 'Time end:', TR10,   &
                &  'Postion start:', TR4, 'Postion end:')")
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '----------------------------------------------------')")
            do i_pert = 1, pert_cr%n_pert
                associate (cr => pert_cr%perts(i_pert))
                    write(unit=unit_, fmt="(2X, I3, TR9, I4, TR12, F8.4, TR12, F8.4, TR10, F7.3, TR12, F7.3)")   &
                        &  i_pert, cr%bankID, cr%time_start, cr%time_end, cr%step_start, cr%step_end
                end associate
            end do
            write(unit=unit_, fmt="(2X, '------------------------------------------------------------------------------------------------',   &
                &  '----------------------------------------------------')")
        end if
        
    end subroutine Print_case_parsing
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! print iteration result for eigenvalue problem
    !===============================================================================================
    subroutine Print_iteration_step_eigen (unit_, iter_count, error_eigen, error_flux, error_inner)
    
        type(IterationCounter), intent(in)  :: iter_count
        real(KREAL), intent(in)  :: error_eigen
        real(KREAL), intent(in)  :: error_flux
        real(KREAL), intent(in)  :: error_inner
        integer, intent(in)            :: unit_
        
        write(unit=unit_, fmt="(6x,I4,10X,I4,7X,F10.7,4x,ES11.4,3X,ES11.4,3X,ES11.4)")    &
            &   iter_count%out, iter_count%in, iter_count%eigenvalue, error_eigen, error_flux, error_inner
    
    end subroutine Print_iteration_step_eigen
    
    !$
    !===============================================================================================
    ! print iteration result for FSP problem
    !===============================================================================================
    subroutine Print_iteration_step_FSP (unit_, iter_count, error_eigen, error_flux, error_inner)
    
        type(IterationCounter), intent(in)  :: iter_count
        real(KREAL), intent(in)  :: error_eigen
        real(KREAL), intent(in)  :: error_flux
        real(KREAL), intent(in)  :: error_inner
        integer, intent(in)            :: unit_
        
        write(unit=unit_, fmt="(6x,I4,10X,I4,7X,F10.7,4x,ES11.4,3X,ES11.4,3X,ES11.4)")    &
            &   iter_count%out, iter_count%in, iter_count%ks, error_eigen, error_flux, error_inner
    
    end subroutine Print_iteration_step_FSP

    !$
    !===============================================================================================
    ! print title information of iteration
    !===============================================================================================
    subroutine Print_iteration_title (unit_, is_eigen, is_adjoint, is_transient, tidx, ctime)
        
        ! intent parameters
        integer, intent(in)  :: unit_
        logical, intent(in)  :: is_eigen
        logical, intent(in)  :: is_adjoint
        logical, intent(in), optional          :: is_transient
        integer, intent(in), optional          :: tidx
        real(KREAL), intent(in), optional  :: ctime
        
        ! for eigenvalue problem
        if (is_eigen ) then 
            if (.NOT. is_adjoint) then
!                write(unit=unit_, fmt="(/)")
                write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
                if (present(is_transient) .and. is_transient )  then
                    write(unit=unit_, fmt="(1X, 'TimeStep', I5, ' @ ', ES11.4, ' s.')") tidx, ctime
                end if
                write(unit=unit_, fmt="(1X, A)")  'Forward, Eigenvalue calculation:'
                write(unit=unit_, fmt="(1X, A)")  '----------------------------------'
                write(unit=unit_, fmt="(4x,'Outer Index:',5x,'Inner #:',6x      &
                    & ,'K_eff:',6x,'ERROR_k',7x,'ERROR_out',5x,'ERROR_in')")
                write(unit=unit_, fmt="(2x,  '----------------------------------------', &
                    &  '------------------------------------------')")
            else 
!                write(unit=unit_, fmt="(/)")
                write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
                if (present(is_transient) .and. is_transient )  then
                    write(unit=unit_, fmt="(1X, 'TimeStep', I5, ' @ ', ES11.4, ' s.')") tidx, ctime
                end if
                write(unit=unit_, fmt="(1X, A)")  'Adjoint, Eigenvalue calculation:'
                write(unit=unit_, fmt="(1X, A)")  '----------------------------------'
                write(unit=unit_, fmt="(4x,'Outer Index:',5x,'Inner #:',6x      &
                    & ,'K_eff:',6x,'ERROR_k',7x,'ERROR_out',5x,'ERROR_in')")
                write(unit=unit_, fmt="(2x,  '----------------------------------------', &
                    &  '------------------------------------------')")
            end if 
            
        ! for FSP problem
        else 
            if (.NOT. is_adjoint) then
!                write(unit=unit_, fmt="(/)")
                write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
                if (present(is_transient) .and. is_transient )  then
                    write(unit=unit_, fmt="(1X, 'TimeStep', I5, ' @ ', ES11.4, ' s.')") tidx, ctime
                end if
                write(unit=unit_, fmt="(1X, A)")  'Forward, Fixed-Source calculation:'
                write(unit=unit_, fmt="(1X, A)")  '----------------------------------'
                write(unit=unit_, fmt="(4x,'Outer Index:',5x,'Inner #:',6x      &
                    & ,'K_s:  ',6x,'ERROR_k',7x,'ERROR_out',5x,'ERROR_in')")
                write(unit=unit_, fmt="(2x,  '----------------------------------------', &
                    &  '------------------------------------------')")
            else 
!                write(unit=unit_, fmt="(/)")
                write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
                if (present(is_transient) .and. is_transient )  then
                    write(unit=unit_, fmt="(1X, 'TimeStep', I5, ' @ ', ES11.4, ' s.')") tidx, ctime
                end if
                write(unit=unit_, fmt="(1X, A)")  'Adjoint, Fixed-Source calculation:'
                write(unit=unit_, fmt="(1X, A)")  '----------------------------------'
                write(unit=unit_, fmt="(4x,'Outer Index:',5x,'Inner #:',6x      &
                    & ,'K_s:  ',6x,'ERROR_k',7x,'ERROR_out',5x,'ERROR_in')")
                write(unit=unit_, fmt="(2x,  '----------------------------------------', &
                    &  '------------------------------------------')")
            end if
        end if
        
    end subroutine Print_iteration_title
    
    !$
    !===============================================================================================
    ! print summay information per iteration
    !===============================================================================================
    subroutine Print_iteration_summary (unit_, is_eigen, is_adjoint)

        integer, intent(in)  :: unit_
        logical, intent(in)  :: is_eigen
        logical, intent(in)  :: is_adjoint
        
        ! loop index
        integer  :: ia, ir, im, ig, id, is, iz, il
        
        ! local variables
        real(KREAL)  :: flux_volume(ns%state%ng, ns%state%mat)              ! flux multi volume per matetrial ID per energy group
        real(KREAL)  :: volume_per_mat(ns%state%mat)                        ! volume per material ID
        
        ! ----------------------------------------------------------------------
        if(is_eigen ) then
            write(unit=unit_, fmt="(1X, /)") 
            write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
            write(unit=unit_, fmt="(1X, 'Eigenvalue is               :', F10.7)") iter_count%eigenvalue
            write(unit=unit_, fmt="(1X, 'Number of out iteration is  :', I4)") iter_count%out
            write(unit=unit_, fmt="(1X, /)") 
            
        else
            write(unit=unit_, fmt="(1X, /)") 
            write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
            write(unit=unit_, fmt="(1X, 'Number of out iteration is   :', I4)") iter_count%out
        end if
        
        ! ----------------------------------------------------------------------
        flux_volume     = 0.0
        volume_per_mat  = 0.0
        
        ! obtain volume sum, volume*flux sum per material zone
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = (ia-1) * ns%state%nodal + ir
                im = mat_info%loading(iz, ia)
                do ig = 1, ns%state%ng
                    if (.NOT. is_adjoint)  then
                        flux_volume(ig, im) = flux_volume(ig,im) + flux_forward%ngs(ig)%scalar(ir,ia)*geom%area(ir)*geom%height(ia)
                    else 
                        flux_volume(ig, im) = flux_volume(ig,im) + flux_adjoint%ngs(ig)%scalar(ir,ia)*geom%area(ir)*geom%height(ia)
                    end if
                end do
                
                volume_per_mat(im) = volume_per_mat(im) + geom%area(ir)*geom%height(ia)
            end do
        end do
        
        ! flux per material zone after nomarlize, volume of material ID as weighting function
        do im = 1, ns%state%mat
            if (ABS(volume_per_mat(im)) > EPS_ZERO) then
                do ig = 1, ns%state%ng
                    flux_volume(ig,im) = flux_volume(ig,im) / volume_per_mat(im)
                end do
            end if
        end do
        
        ! ----------------------------------------------------------------------
        ! output normalize flux per material zone
        if (ns%state%mat <= 50)  then
            write(unit=unit_, fmt="(1X, 'The nomorlize coefficient is: ', ES11.4)") timelist%normal_factor
            write(unit=unit_, fmt="(1X, /)") 
            write(unit=unit_, fmt="(1X, 'The flux per matrial ID per energy group:')")   
            write(unit=unit_, fmt="(1X, ' material ID:', TR9)", advance='no')
            do im = 1, ns%state%mat
                write(unit=unit_, fmt="(I4, TR7)", advance='no') im
            end do
            write(unit=unit_, fmt="()")
            do ig = 1, ns%state%ng
                write(unit=unit_, fmt="(1X, ' Energy group=', I3, TR2, *(ES11.4))") ig, (flux_volume(ig,im), im=1, ns%state%mat)
            end do
            
            ! output volume per material zone
            write(unit=unit_, fmt="(1X, /)") 
            write(unit=unit_, fmt="(1X, 'The volume per material zone:')")
            write(unit=unit_, fmt="(1X, ' material ID:', TR9)", advance='no')
            do im = 1, ns%state%mat
                write(unit=unit_, fmt="(I4, TR7)", advance='no') im
            end do
            write(unit=unit_, fmt="()")
            write(unit=unit_, fmt="(1X, TR19, *(ES11.4))") (volume_per_mat(im), im=1, ns%state%mat)
        end if
        
    end subroutine Print_iteration_summary
    
end module output_main
