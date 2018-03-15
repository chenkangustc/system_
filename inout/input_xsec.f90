!$
!===================================================================================================
!
!   cross section library format translate:
!       1-interpolation in HDF5, 2-fitting by Lilac, 3-fitting by NEACRP, 4-LRA
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    xsec_read_known
!                               read_xsec_unknown
!                               read_xsec_LRA
!                               
!
!   Public type lists:          No
!
!===================================================================================================
module input_xsec
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector
    use interpolation,              only : two_variables_interpolation
    use vector_operation,           only : sort_vector
    use self_lagrange_interp
    use string 
    use files_dirs
    use hdf5_interface
    use HDF5
    use Link
    
    use material_header,            only : Material, CrossSection, CrossSectionInput, cross_section_info_tp,        &
                                        &   KineticsParameter, KineticsParameterInput, kinetics_parameter_info_tp
    use CRbank_header,              only : ControlRodBank
    use link_header,                only : LinkParameter
    use LRAmodel_header,            only : LRAmodel
    use feedback_header,            only : FeedbackParameter

    implicit none 
    public
!    public  :: xsec_read_known, read_xsec_unknown, read_xsec_LRA
    
    interface  xsec_read_known
        module procedure  Import_to_XSEC
        module procedure  HDF5_to_XSEC
        module procedure  Plain_to_XSEC_by_file
        module procedure  Plain_to_PARAM_by_file
        module procedure  Plain_to_XSEC_by_unit
        module procedure  Plain_to_PARAM_by_unit
    end interface xsec_read_known
    
    interface  read_xsec_unknown
        module procedure  read_xsec_unknown_core
        module procedure  read_xsec_unknown_CR
        module procedure  read_xsec_unknown_all
    end interface read_xsec_unknown
    
    interface  read_xsec_LRA
        module procedure  read_xsec_LRA_core
        module procedure  read_xsec_LRA_CR
        module procedure  read_xsec_LRA_all
    end interface read_xsec_LRA
    
    ! --------------------------------------------------------------------------
    ! store lilac information into memory for fitting, while not read from files
    ! type for fitting terms information
    type, private  :: NumTerms_Type
        integer, public, allocatable              :: Option(:)
        integer, public, allocatable              :: XOrder(:)
        integer, public, allocatable              :: IndexXvars(:)
        real(8), public, allocatable              :: XvarsInTerm(:)
        real(8), public, allocatable              :: AV(:)
        real(8), public, allocatable              :: STD(:)
        
        integer, public, allocatable              :: OrderMatrix(:, :)
        real(8), public, allocatable              :: xVector(:, :)
        integer, public                           :: NumXinTerm
        integer, public                           :: NumCoefs
    end type NumTerms_Type
    
    ! type for fitting information
    type, private  :: lilac_parameter_tp 
        integer, public                           :: ID
        character(len=MAX_WORD_LEN), public       :: file
        integer, public                           :: NumXvars
        integer, public                           :: NumTYvars
        type(Link_Type), public                   :: a_link
        integer, public, allocatable              :: TOption(:)
        type(NumTerms_Type), public, allocatable  :: Terms(:)
    end type lilac_parameter_tp
    
    type  LilacParameter
        type(lilac_parameter_tp), public, allocatable   :: mats(:, :)
        logical, public                                 :: is_define = .FALSE.
    contains
        procedure, public  :: clean => Free_LilacParameter
    end type LilacParameter
    
    ! --------------------------------------------------------------------------
    ! store NEACRP information into memory for fitting, while not read from files
    type, private  :: neacrp_parameter_tp
        integer, public                         :: ID
        character(len=MAX_WORD_LEN), public     :: file
        logical, public                         :: is_CR
        real(KREAL), public, allocatable        :: xsec0(:)
        real(KREAL), public, allocatable        :: partial(:)
        real(KREAL), public, allocatable        :: delataCR(:)
    end type neacrp_parameter_tp
    
    type  NEACRPParameter
        type(neacrp_parameter_tp), public, allocatable  :: mats(:)
        logical, public                                 :: is_define = .FALSE.
    contains
        procedure, public  :: clean => Free_NEACRPParameter
    end type NEACRPParameter
    
    ! --------------------------------------------------------------------------
    ! type for variable information in HDF5 interpolation
    type, private  ::  state_variables_tp
        real(KREAL)  :: t_fuel
        real(KREAL)  :: rho_coolant
    contains
        procedure, private :: Equal_state_variables_tp
        generic,   public  :: assignment (=) => Equal_state_variables_tp
        procedure, public  :: in => Is_in_field
    end type state_variables_tp
    
    type, private  ::  hdf5_state_point_tp
        integer, public                           :: ID
        character(len=MAX_WORD_LEN), public       :: file
        integer, public                           :: n_case
        integer, public                           :: flag                       ! flag for material type
        real(KREAL), allocatable, public          :: t_fuel(:)                  ! fuel temperature in K
        real(KREAL), allocatable, public          :: rho_coolant(:)             ! coolant density in g/cm^3
        
        integer, public                           :: n_xTf 
        integer, public                           :: n_yDc
        integer, allocatable, public              :: xyto1D(:, :)
        real(KREAL), allocatable, public          :: xTf(:)
        real(KREAL), allocatable, public          :: yDc(:)
        
        type(cross_section_info_tp), allocatable, public       :: xsecs(:)
        type(kinetics_parameter_info_tp), allocatable, public  :: params(:)
    end type hdf5_state_point_tp
    
    type  HDF5Parameter
        type(hdf5_state_point_tp), public, allocatable  :: mats(:)
        logical, public                                 :: is_define = .FALSE.
    contains
        procedure, public  :: clean => Free_HDF5Parameter
        procedure, public  :: print => Print_HDF5Parameter
    end type HDF5Parameter
        
    type(cross_section_info_tp)  :: xsec_00
    type(cross_section_info_tp)  :: xsec_01
    type(cross_section_info_tp)  :: xsec_10
    type(cross_section_info_tp)  :: xsec_11
    type(kinetics_parameter_info_tp)  :: param_00
    type(kinetics_parameter_info_tp)  :: param_01
    type(kinetics_parameter_info_tp)  :: param_10
    type(kinetics_parameter_info_tp)  :: param_11
    
    type(cross_section_info_tp)  :: xsec_a
    type(cross_section_info_tp)  :: xsec_b
    type(cross_section_info_tp)  :: xsec_c
    type(kinetics_parameter_info_tp)  :: param_a
    type(kinetics_parameter_info_tp)  :: param_b
    type(kinetics_parameter_info_tp)  :: param_c 
        
    type(NEACRPParameter)  :: neacrp_fit 
    type(LilacParameter)   :: lilac_fit
    type(HDF5Parameter)    :: hdf5_interp
    real(8), parameter  :: NECP_RealZero  = 1.0D-18
    real(8), parameter  :: NECP_RealEqual = 1.0D-8 
    
    integer, parameter  :: LILAC_FITTING  = 1
    integer, parameter  :: NEACRP_FITTING = 2
    
    integer, parameter  :: INTERP_LINNEAR   = 1 
    integer, parameter  :: INTERP_PICEWISL  = 2
    integer, parameter  :: INTERP_LAGRANGE  = 3 
    
    ! --------------------------------------------------------------------------
    type(ErrorCollector)  :: a_error 
    
contains
    !$
    !===============================================================================================
    ! read xsec from burn-simulator output, for all nodals
    !===============================================================================================
    subroutine Import_to_XSEC (dir, file_in, burnstep, a_xsec, a_param)
        
        character(len=*), intent(in)  :: dir
        character(len=*), intent(in)  :: file_in
        integer, intent(in)           :: burnstep
        type(CrossSectionInput), intent(in out)                 :: a_xsec
        type(KineticsParameterInput), intent(in out), optional  :: a_param
        
        ! local
        character(len=MAX_WORD_LEN)  :: file_from
        
        character(len=MAX_WORD_LEN)  :: groupname
        character(len=MAX_WORD_LEN)  :: dataname
        character(len=MAX_WORD_LEN)  :: tmp, tmp_ia, tmp_iz
        integer(HID_T)  :: file
        integer(HID_T)  :: burngroup, nodalgroup

        integer  :: hdferr
        integer  :: rank
        integer  :: dim_scale(1:10)
        
        integer  :: i_allocate
        integer  :: io_error
        integer  :: ia, iz, im
        
        ! get file name
        file_from = TRIM(ADJUSTL(dir)) // TRIM(ADJUSTL(file_in))
        
        ! initialize fortran interface
        call h5open_f(hdferr)

        ! create a new file using the default properties
        call hdf5_file_open(TRIM(ADJUSTL(file_from)), file, 'r')
        
        write(unit=tmp, fmt=*) burnstep
        groupname = 'step_' // TRIM(ADJUSTL(tmp))
        call hdf5_open_group(file, groupname, burngroup)
        call hdf5_get_attr(file, groupname, 'layer', ns%state%layer)
        call hdf5_get_attr(file, groupname, 'zone', ns%state%zone)
        call hdf5_get_attr(file, groupname, 'ng', ns%state%ng)
        
        im = 0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = im + 1
                
                write(unit=tmp_ia, fmt=*) ia
                write(unit=tmp_iz, fmt=*) iz
                groupname = 'layer_' // TRIM(ADJUSTL(tmp_ia)) // '_zone_' // TRIM(ADJUSTL(tmp_iz))
                call hdf5_open_group(burngroup, groupname, nodalgroup)
                
                ! for cross section
                dataname = 'chi_steady'
                rank = 1
                dim_scale(1) = ns%state%ng
                call hdf5_read_data(nodalgroup, dataname, a_xsec%mats(im)%chi_steady, dim_scale(1:rank))
                
                dataname = 'sigma_t'
                rank = 1
                dim_scale(1) = ns%state%ng
                call hdf5_read_data(nodalgroup, dataname, a_xsec%mats(im)%sigma_t, dim_scale(1:rank))

                dataname = 'sigma_f_nu'
                rank = 1
                dim_scale(1) = ns%state%ng
                call hdf5_read_data(nodalgroup, dataname, a_xsec%mats(im)%sigma_f_nu, dim_scale(1:rank))

                dataname = 'sigma_f_kappa'
                rank = 1
                dim_scale(1) = ns%state%ng
                call hdf5_read_data(nodalgroup, dataname, a_xsec%mats(im)%sigma_f_kappa, dim_scale(1:rank))
                ! NOTE:
                a_xsec%mats(im)%sigma_f_kappa = a_xsec%mats(im)%sigma_f_kappa * 1.0E6

                dataname = 'sigma_s'
                rank = 3
                dim_scale(1) = ns%state%ng
                dim_scale(2) = ns%state%ng
                dim_scale(3) = ns%deduce%scat_xs
                call hdf5_read_data(nodalgroup, dataname, a_xsec%mats(im)%sigma_s, dim_scale(1:rank))
                
                if (present(a_param))  then
                    ! for kinetics paramters
                    dataname = 'chi_delay'
                    rank = 2
                    dim_scale(1) = nt%state%dg
                    dim_scale(2) = ns%state%ng
                    call hdf5_read_data(nodalgroup, dataname, a_param%mats(im)%chi_delay, dim_scale(1:rank))
                    
                    dataname = 'velocity'
                    rank = 1
                    dim_scale(1) = ns%state%ng
                    call hdf5_read_data(nodalgroup, dataname, a_param%mats(im)%velocity, dim_scale(1:rank))
                    ! NOTE:
                    a_param%mats(im)%velocity = 1383168.0 / a_param%mats(im)%velocity
                    
                    dataname = 'sigma_bvf'
                    rank = 2
                    dim_scale(1) = nt%state%dg
                    dim_scale(2) = ns%state%ng
                    call hdf5_read_data(nodalgroup, dataname, a_param%mats(im)%sigma_bvf, dim_scale(1:rank))
                    
                    dataname = 'sigma_bvl'
                    rank = 2
                    dim_scale(1) = nt%state%dg
                    dim_scale(2) = ns%state%ng
                    call hdf5_read_data(nodalgroup, dataname, a_param%mats(im)%sigma_bvl, dim_scale(1:rank))
                end if
                
                call hdf5_close_group(nodalgroup)
            end do
        end do
        
        ! close this case
        call hdf5_close_group(burngroup)
        
        ! exit hdf5 environment
        call hdf5_file_close(file)
        call h5close_f(hdferr)
        
    end subroutine Import_to_XSEC
    
    !$
    !===============================================================================================
    ! read xsec from lattice output, for one material
    !===============================================================================================
    subroutine HDF5_to_XSEC (dir, file_in, case_id, burn_id, a_xsec, a_param)
    
        character(len=*), intent(in)  :: dir
        character(len=*), intent(in)  :: file_in
        integer, intent(in)           :: case_id
        integer, intent(in)           :: burn_id
        type(cross_section_info_tp), intent(in out)                 :: a_xsec
        type(kinetics_parameter_info_tp), intent(in out), optional  :: a_param
        
        character(len=MAX_WORD_LEN)   :: file_from
        
        character(len=MAX_WORD_LEN) :: groupname
        character(len=MAX_WORD_LEN) :: subname
        character(len=MAX_WORD_LEN) :: tmp
        integer(HID_T) :: file, group
        integer(HID_T) :: casegroup, burngroup, datagroup
        
        integer  :: hdferr
        integer  :: rank
        integer  :: dim_scale(1:10)
        
        integer  :: i_allocate
        integer  :: io_error
        integer  :: i_case, i_burnup, i_xs, i_index
        integer  :: ig
        integer  :: ip
        
        type  state_point_tp
            real(KREAL)  :: burnup                                              ! current burnup value
            real(KREAL)  :: t_fuel                                              ! fuel temperature in K
            real(KREAL)  :: rho_coolant                                         ! coolant density in g/cm^3
            integer          :: n_case                                          ! number of case
            integer          :: n_burnup                                        ! number of burn point
            integer          :: flag                                            ! flag for material type
        end type state_point_tp
        type(state_point_tp)  :: a_stat
            
        ! type for xsec status define
        type  xsec_status_tp
            integer                     :: ng                                   ! energy group 
            integer                     :: dg                                   ! delay neutron group
            integer                     :: scat_xs                              ! number of scatter matrix
        end type xsec_status_tp
        type(xsec_status_tp)            :: local_ns
        
        real(KREAL), allocatable  :: sigma_a(:)
        real(KREAL), allocatable  :: sigma_t_add(:)                             ! this sigma_t comes from adding sigma_s and sigma_a
        
        real(KREAL), allocatable  :: sigma_t(:)
        real(KREAL), allocatable  :: sigma_s(:, :, :)
        real(KREAL), allocatable  :: sigma_f_nu(:)
        real(KREAL), allocatable  :: sigma_f_kappa(:)
        real(KREAL), allocatable  :: chi_steady(:)
        
        real(KREAL), allocatable  :: sigma_bvf(:, :)
        real(KREAL), allocatable  :: sigma_bvl(:, :)
        real(KREAL), allocatable  :: chi_delay(:, :)
        real(KREAL), allocatable  :: velocity(:)
        real(KREAL), allocatable  :: beta(:)
        real(KREAL), allocatable  :: lambda(:)
        
        integer                       :: n_pin
        real(KREAL), allocatable  :: pow_shape(:)
        
        ! ----------------------------------------------------------------------
        local_ns%ng = ns%state%ng
        local_ns%dg = nt%state%dg
        local_ns%scat_xs = ns%deduce%scat_xs
        
        allocate(sigma_a(local_ns%ng), stat=i_allocate)
        allocate(sigma_t_add(local_ns%ng), stat=i_allocate)
        
        allocate(sigma_t(local_ns%ng), stat=i_allocate)
        allocate(sigma_s(local_ns%ng, local_ns%ng, local_ns%scat_xs), stat=i_allocate)
        allocate(sigma_f_nu(local_ns%ng), stat=i_allocate)
        allocate(sigma_f_kappa(local_ns%ng), stat=i_allocate)
        allocate(chi_steady(local_ns%ng), stat=i_allocate)
        
        allocate(sigma_bvf(local_ns%dg, local_ns%ng), stat=i_allocate)
        allocate(sigma_bvl(local_ns%dg, local_ns%ng), stat=i_allocate)
        allocate(chi_delay(local_ns%dg, local_ns%ng), stat=i_allocate)
        allocate(velocity(local_ns%ng), stat=i_allocate)
        allocate(beta(local_ns%dg), stat=i_allocate)
        allocate(lambda(local_ns%dg), stat=i_allocate)
            
        ! get file name
        file_from = TRIM(ADJUSTL(dir)) // TRIM(ADJUSTL(file_in))
        
        ! initialize fortran interface
        call h5open_f(hdferr)

        ! create a new file using the default properties
        call hdf5_file_open(TRIM(ADJUSTL(file_from)), file, 'r')
        
        ! ----------------------------------------------------------------------
        ! start transfer
        ! open the top group in the file
        call Get_file_basename(TRIM(file_in), groupname)
        tmp = ''
        call hdf5_open_group(file, groupname, group)
        call hdf5_get_attr(file, groupname, 'NumLIST', a_stat%n_case)
        call hdf5_get_attr(file, groupname, 'Fuel_Flag', a_stat%flag)
        
        ! find case
        groupname = 'CASE' // TRIM(Int_to_string(case_id, digit=4))
        call hdf5_open_group(group, groupname, casegroup)
        call hdf5_get_attr(group, groupname, 'NBnPT', a_stat%n_burnup)
        call hdf5_get_attr(group, groupname, 'TF', a_stat%t_fuel)
        call hdf5_get_attr(group, groupname, 'DC', a_stat%rho_coolant)
        
        ! find burnup point
        groupname = 'BUPNT' // TRIM(Int_to_string(burn_id, digit=3))
        call hdf5_open_group(casegroup, groupname, burngroup)
        call hdf5_get_attr(casegroup, groupname, 'BU', a_stat%burnup)
        
        ! xsec
        groupname = 'MACCOM'
        call hdf5_open_group(burngroup, groupname, datagroup)
        
        rank = 1
        dim_scale(1) = local_ns%ng
        call hdf5_read_data(datagroup, 'TOTAL', sigma_t, dim_scale(1:rank))
        
        rank = 3
        dim_scale(1) = local_ns%ng
        dim_scale(2) = local_ns%ng
        dim_scale(3) = local_ns%scat_xs
        call hdf5_read_data(datagroup, 'Scatt', sigma_s, dim_scale(1:rank))
        
        rank = 1
        dim_scale(1) = local_ns%ng
        call hdf5_read_data(datagroup, 'NuFiss', sigma_f_nu, dim_scale(1:rank))
        
        rank = 1
        dim_scale(1) = local_ns%ng
        call hdf5_read_data(datagroup, 'KpFiss', sigma_f_kappa, dim_scale(1:rank))
        
        rank = 1
        dim_scale(1) = local_ns%ng
        call hdf5_read_data(datagroup, 'CHI', chi_steady, dim_scale(1:rank))
                
        rank = 1
        dim_scale(1) = local_ns%ng
        call hdf5_read_data(datagroup, 'ABSORB', sigma_a, dim_scale(1:rank))
        
        ! NOTE: get total xsec, should be Legendre expansion
        sigma_t_add = 0.0
        do ig = 1, local_ns%ng
!            do i_xs = 1, local_ns%scat_xs
                i_xs = 1
                sigma_t_add(ig) = sigma_a(ig) + SUM(sigma_s(:, ig, i_xs))
!            end do
        end do
                
        call hdf5_close_group(datagroup)
        
        if (PRESENT(a_param))  then 
            ! kinetics parameter
            groupname = 'KTPARAM'
            call hdf5_open_group(burngroup, groupname, datagroup)
                
            select case(TRIM(nt%flag%kinetics_type))
            case ('NORMAL')
                rank = 1
                dim_scale(1) = local_ns%dg
                call hdf5_read_data(datagroup, 'beta', beta, dim_scale(1:rank))
                
                rank = 1
                dim_scale(1) = local_ns%dg
                call hdf5_read_data(datagroup, 'lambda', lambda, dim_scale(1:rank))
                
                rank = 2
                dim_scale(1) = local_ns%dg
                dim_scale(2) = local_ns%ng
                call hdf5_read_data(datagroup, 'CHID', chi_delay, dim_scale(1:rank))
                
                rank = 1
                dim_scale(1) = local_ns%ng
                call hdf5_read_data(datagroup, 'SIGMAV', velocity, dim_scale(1:rank))
                
            case ('SRAC')
                rank = 2
                dim_scale(1) = local_ns%dg
                dim_scale(2) = local_ns%ng
                call hdf5_read_data(datagroup, 'BVFSIG', sigma_bvf, dim_scale(1:rank))
                
                rank = 2
                dim_scale(1) = local_ns%dg
                dim_scale(2) = local_ns%ng
                call hdf5_read_data(datagroup, 'BVLSIG', sigma_bvl, dim_scale(1:rank))
                
                rank = 2
                dim_scale(1) = local_ns%dg
                dim_scale(2) = local_ns%ng
                call hdf5_read_data(datagroup, 'CHID', chi_delay, dim_scale(1:rank))
                
                rank = 1
                dim_scale(1) = local_ns%ng
                call hdf5_read_data(datagroup, 'SIGMAV', velocity, dim_scale(1:rank))
                velocity = 1383168.0 / velocity                                     ! see SRAC manual
            
            end select 
            call hdf5_close_group(datagroup)
            
            ! power shape
            if (a_stat%flag == 4) then
                groupname = 'POWSHAPE'
                call hdf5_open_group(burngroup, groupname, datagroup)
                call hdf5_get_attr(burngroup, groupname, 'NUMPIN', n_pin)
                
                if (allocated(pow_shape))  then
                    deallocate(pow_shape, stat=i_allocate)
                end if
                allocate(pow_shape(n_pin), stat=i_allocate)
                
                rank = 1
                dim_scale(1) = n_pin
                call hdf5_read_data(datagroup, 'PIN', pow_shape, dim_scale(1:rank))
            
                call hdf5_close_group(datagroup)
            end if
        end if 
        
        call hdf5_close_group(burngroup)
        
        ! ----------------------------------------------------------------------       
        ! transit out, NOTE scatter xsec
        a_xsec%chi_steady     = chi_steady
!        a_xsec%sigma_t        = sigma_t
        a_xsec%sigma_t        = sigma_t_add                                     ! use this one 
        a_xsec%sigma_f_nu     = sigma_f_nu
        a_xsec%sigma_f_kappa  = sigma_f_kappa
        do i_xs = 1, local_ns%scat_xs
            do ig = 1, local_ns%ng
                a_xsec%sigma_s(ig, :, i_xs)= sigma_s(:, ig, i_xs)
            end do
        end do
        
        if (present(a_param))  then
            a_param%chi_delay = chi_delay
            a_param%velocity = velocity
            a_param%sigma_bvf = sigma_bvf
            a_param%sigma_bvl = sigma_bvl
            a_param%beta = beta
            a_param%lambda = lambda
        end if
        
        call hdf5_close_group(casegroup)
        call hdf5_close_group(group)
        
        ! ----------------------------------------------------------------------
        ! close file
        call hdf5_file_close(file)
        call h5close_f(hdferr)
        
        if (allocated(sigma_a))         deallocate(sigma_a, stat=i_allocate)
        if (allocated(sigma_t_add))     deallocate(sigma_t_add, stat=i_allocate)

        if (allocated(sigma_t))         deallocate(sigma_t, stat=i_allocate)
        if (allocated(sigma_s))         deallocate(sigma_s, stat=i_allocate)
        if (allocated(sigma_f_nu))      deallocate(sigma_f_nu, stat=i_allocate)
        if (allocated(sigma_f_kappa))   deallocate(sigma_f_kappa, stat=i_allocate)
        if (allocated(chi_steady))      deallocate(chi_steady, stat=i_allocate)
           
        if (allocated(sigma_bvf))       deallocate(sigma_bvf, stat=i_allocate)
        if (allocated(sigma_bvl))       deallocate(sigma_bvl, stat=i_allocate)
        if (allocated(chi_delay))       deallocate(chi_delay, stat=i_allocate)
        if (allocated(velocity))        deallocate(velocity, stat=i_allocate)
        if (allocated(beta))            deallocate(beta, stat=i_allocate)
        if (allocated(lambda))          deallocate(lambda, stat=i_allocate)
        
        if (allocated(pow_shape))       deallocate(pow_shape, stat=i_allocate)
        
    end subroutine HDF5_to_XSEC
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Plain_to_XSEC_by_file (dir, file_in, a_xsec)
    
        character(len=*), intent(in)                 :: dir
        character(len=*), intent(in)                 :: file_in
        type(cross_section_info_tp), intent(in out)  :: a_xsec
        
        ! local variables
        type(cross_section_info_tp)  :: tmp_xsec
        character(len=MAX_WORD_LEN)  :: file_from
        character(len=MAX_WORD_LEN)  :: tmp_line
        integer  :: unit_in
        integer  :: io_error
        integer  :: ig, iig, ic, ip
        
        call tmp_xsec%alloc ()
        
        file_from = TRIM(ADJUSTL(dir)) // TRIM(ADJUSTL(file_in))
        open(newunit=unit_in, file=file_from, status='old', action='read', iostat=io_error)
        
        ! read steady xsec
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%chi_steady(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_t(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_nu(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_kappa(ig), ig=1, ns%state%ng)
        
        do ic = 1, ns%deduce%scat_xs
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_s(ig,iig,ic), iig=1, ns%state%ng)
            end do
        end do
        
        close (unit=unit_in, status='keep', iostat=io_error)
        
        a_xsec = tmp_xsec
        
        call tmp_xsec%clean ()
    
    end subroutine Plain_to_XSEC_by_file
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Plain_to_PARAM_by_file (dir, file_in, a_param)
    
        character(len=*), intent(in)  :: dir
        character(len=*), intent(in)  :: file_in
        type(kinetics_parameter_info_tp), intent(in out)  :: a_param
        
        ! local variables
        type(cross_section_info_tp)       :: tmp_xsec
        type(kinetics_parameter_info_tp)  :: tmp_param
        character(len=MAX_WORD_LEN)  :: file_from
        character(len=MAX_WORD_LEN)  :: tmp_line
        integer  :: unit_in
        integer  :: io_error
        integer  :: ig, iig, ic, ip
        
        call tmp_xsec%alloc ()
        call tmp_param%alloc ()
        
        file_from = TRIM(ADJUSTL(dir)) // TRIM(ADJUSTL(file_in))
        open(newunit=unit_in, file=file_from, status='old', action='read', iostat=io_error)
        
        ! read steady xsec
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%chi_steady(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_t(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_nu(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_kappa(ig), ig=1, ns%state%ng)
        
        do ic = 1, ns%deduce%scat_xs
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_s(ig,iig,ic), iig=1, ns%state%ng)
            end do
        end do
        
        ! read kinetics
        ! skip token
        read(unit=unit_in, fmt=*, iostat=io_error)  tmp_line
        
        select case(TRIM(nt%flag%kinetics_type))
        case ('NORMAL')
            do ip = 1, nt%state%dg
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%chi_delay(ip, ig), ig=1, ns%state%ng)
            end do
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%velocity(ig), ig=1,ns%state%ng)
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%beta(ip), ip=1,nt%state%dg)
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%lambda(ip), ip=1,nt%state%dg)

        case ('SRAC')
            do ip = 1, nt%state%dg
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%chi_delay(ip, ig), ig=1, ns%state%ng)
            end do
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%velocity(ig), ig=1,ns%state%ng)
            
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%sigma_bvf(ip,ig), ip=1,nt%state%dg)
            end do
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%sigma_bvl(ip,ig), ip=1,nt%state%dg)
            end do
        end select
        
        close (unit=unit_in, status='keep', iostat=io_error)
        
        a_param = tmp_param
        
        call tmp_xsec%clean ()
        call tmp_param%clean ()
    
    end subroutine Plain_to_PARAM_by_file
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Plain_to_XSEC_by_unit (unit_in, a_xsec)
    
        integer, intent(in)                          :: unit_in
        type(cross_section_info_tp), intent(in out)  :: a_xsec
        
        ! local variables
        type(cross_section_info_tp)  :: tmp_xsec
        integer  :: ig, iig, ic
        integer  :: io_error
        
        call tmp_xsec%alloc ()
        
        ! read steady xsec
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%chi_steady(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_t(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_nu(ig), ig=1, ns%state%ng)
        read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_f_kappa(ig), ig=1, ns%state%ng)
        
        do ic = 1, ns%deduce%scat_xs
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_xsec%sigma_s(ig,iig,ic), iig=1, ns%state%ng)
            end do
        end do
        
        a_xsec = tmp_xsec
        
        call tmp_xsec%clean ()
        
    end subroutine Plain_to_XSEC_by_unit
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Plain_to_PARAM_by_unit (unit_in, a_param)
    
        integer, intent(in)                               :: unit_in
        type(kinetics_parameter_info_tp), intent(in out)  :: a_param
        
        ! local variables
        type(kinetics_parameter_info_tp)  :: tmp_param
        integer  :: ig, iig, ip
        integer  :: io_error
        
        call tmp_param%alloc ()
        
        ! read kinetics parameter
        select case(TRIM(nt%flag%kinetics_type))
        case ('NORMAL')
            do ip = 1, nt%state%dg
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%chi_delay(ip, ig), ig=1, ns%state%ng)
            end do
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%velocity(ig), ig=1,ns%state%ng)
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%beta(ip), ip=1,nt%state%dg)
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%lambda(ip), ip=1,nt%state%dg)

        case ('SRAC')
            do ip = 1, nt%state%dg
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%chi_delay(ip, ig), ig=1, ns%state%ng)
            end do
            read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%velocity(ig), ig=1,ns%state%ng)
            
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%sigma_bvf(ip,ig), ip=1,nt%state%dg)
            end do
            do ig = 1, ns%state%ng
                read(unit=unit_in, fmt=*, iostat=io_error)  (tmp_param%sigma_bvl(ip,ig), ip=1,nt%state%dg)
            end do
        end select
        
        a_param = tmp_param
        
        call tmp_param%clean ()
    
    end subroutine Plain_to_PARAM_by_unit
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine read_xsec_unknown_all (self_fdbk, self_link, mat_info, xsec, param, cr_bank, kcritical)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(CrossSection), intent(in out)       :: xsec
        type(KineticsParameter), intent(in out)  :: param
        type(ControlRodBank), intent(in out)     :: cr_bank
        real(KREAL), intent(in)                  :: kcritical 
        
        if (nt%flag%is_CR_rod)  then
            if (TRIM(ns%flag%link_type) == 'FITTING')  then
                call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                call cr_bank%map (xsec, param)
                
            else if (TRIM(ns%flag%link_type) == 'INTERPOLATION')  then
                call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                call cr_bank%map (xsec, param)
                    
            else if (TRIM(ns%flag%link_type) == 'NEACRP')  then
                call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                call cr_bank%map (xsec, param)
            end if 
            
        else
            if (TRIM(ns%flag%link_type) == 'FITTING')  then
                call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                
            else if (TRIM(ns%flag%link_type) == 'INTERPOLATION')  then
                call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                    
            else if (TRIM(ns%flag%link_type) == 'NEACRP')  then
                call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, xsec, param)
                call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
            end if 
        end if 
        
        ! fix the reserved keff 
        call xsec%critical (kcritical)
        if (nt%flag%is_CR_rod)  then
            call cr_bank%critical (kcritical)
        end if
        
    end subroutine read_xsec_unknown_all
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine read_xsec_unknown_core (self_fdbk, self_link, mat_info, xsec, param)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(CrossSection), intent(in out)       :: xsec
        type(KineticsParameter), intent(in out)  :: param
        
        if (TRIM(ns%flag%link_type) == 'FITTING')  then
            call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, xsec, param)
            
        else if (TRIM(ns%flag%link_type) == 'INTERPOLATION')  then
            call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, xsec, param)
                
        else if (TRIM(ns%flag%link_type) == 'NEACRP')  then
            call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, xsec, param)
        end if  
        
    end subroutine read_xsec_unknown_core
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine read_xsec_unknown_CR (self_fdbk, self_link, mat_info, cr_bank)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(ControlRodBank), intent(in out)     :: cr_bank
        
        if (TRIM(ns%flag%link_type) == 'FITTING')  then
            call read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
            
        else if (TRIM(ns%flag%link_type) == 'INTERPOLATION')  then
            call read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
                
        else if (TRIM(ns%flag%link_type) == 'NEACRP')  then
            call read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, cr_bank=cr_bank)
        end if 
        
    end subroutine read_xsec_unknown_CR
        
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! get xsec by lilac fitting
    !===============================================================================================
    subroutine read_xsec_unknown_Lilac_new (self_fdbk, self_link, mat_info, xsec, param, cr_bank)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(CrossSection), intent(in out), optional       :: xsec
        type(KineticsParameter), intent(in out), optional  :: param
        type(ControlRodBank), intent(in out), optional     :: cr_bank
    
        ! local variables
        type(cross_section_info_tp)       :: xsec1, xsec2
        type(kinetics_parameter_info_tp)  :: param1, param2
        real(KREAL)                       :: Llen, Ulen
        real(KREAL)                       :: Lwt, Uwt
        integer                           :: Lmat, Umat
        
        character(len=MAX_WORD_LEN)  :: extent_ = ''
        integer  :: i_rod
        integer  :: ia, iz, im, ipoint
        integer  :: i, j
        integer  :: nx, ny
        integer  :: unit_
        integer  :: cnt_xsec, cnt_param
        integer  :: tmp_int(2)
        integer  :: i_allocate, io_error
        
        real(8), allocatable     :: Xvars(:)
        real(8), allocatable     :: Yvars(:)
        
        Lwt = 1.0; Uwt = 1.0;
        
        ! set lilac_fit for the first time
        if (.NOT. lilac_fit%is_define)  then
            allocate(lilac_fit%mats(ns%state%mat, self_link%npoint), stat=i_allocate)
            do im = 1, ns%state%mat
                do ipoint = 1, self_link%npoint
                    extent_ = '.' // TRIM(ADJUSTL(Int_to_string(ipoint, digit=3)))
                    lilac_fit%mats(im, ipoint)%ID = mat_info%libs(im)%ID
                    lilac_fit%mats(im, ipoint)%file = TRIM(ADJUSTL(DIR_LINK)) // TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(extent_))
                    
                    ! get 'NumTYvars' from link file
                    open(newunit=unit_, file=lilac_fit%mats(im, ipoint)%file, status='old', action='read', iostat=io_error)
                    read(unit=unit_, fmt=*)  
                    read(unit=unit_, fmt=*)  
                    read(unit=unit_, fmt=*)  
                    read(unit=unit_, fmt=*)  tmp_int(1), tmp_int(2), ny
                    close(unit=unit_, status='keep', iostat=io_error)
                    nx = self_link%n_parameter
                    lilac_fit%mats(im, ipoint)%NumXvars = nx
                    lilac_fit%mats(im, ipoint)%NumTYvars = ny
                    
                    allocate(Xvars(nx), stat=i_allocate)
                    Xvars = 0.0
                    call pre_Link_Interpolate (lilac_fit%mats(im, ipoint), Xvars, nx)
                    if (allocated(Xvars))       deallocate(Xvars)
                end do 
            end do
            lilac_fit%is_define = .TRUE.
        end if
        
        if (PRESENT(xsec) .AND. PRESENT(param))  then
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    im = mat_info%loading(iz, ia)
                    ipoint = self_link%get_point(self_fdbk%Bu(iz,ia))
                    
                    nx = lilac_fit%mats(im, ipoint)%NumXvars
                    ny = lilac_fit%mats(im, ipoint)%NumTYvars
                    allocate(Xvars(nx), stat=i_allocate)
                    allocate(Yvars(ny), stat=i_allocate)
                        
                    call self_fdbk%get_state(self_link, iz, ia, Xvars)
                    call do_Link_Interpolate (lilac_fit%mats(im, ipoint), Xvars, nx, Yvars, ny)
                    
                    xsec1 = xsec%matrixs(iz, ia)
                    cnt_xsec = xsec1%count ()
                    call xsec1%assign (Yvars(: cnt_xsec))
                    
                    param1 = param%matrixs(iz, ia)
                    cnt_param = param1%count ()
                    call param1%assign (Yvars(cnt_xsec+1: cnt_xsec+cnt_param))
                    
                    if (allocated(Xvars))       deallocate(Xvars)
                    if (allocated(Yvars))       deallocate(Yvars)
                            
                    xsec%matrixs(iz, ia) = xsec1
                    param%matrixs(iz, ia) = param1
                end do
            end do
        end if
        
        ! feedback for CR xsec
        if (PRESENT(cr_bank) .AND. allocated(cr_bank%rods))  then
            call xsec1%alloc ()
            call param1%alloc ()
            call xsec2%alloc ()
            call param2%alloc ()
                    
            do i_rod = 1, cr_bank%state%n_rod
                iz = cr_bank%rods(i_rod)%zone
                do ia = cr_bank%state%n_bottom+1, cr_bank%state%na-cr_bank%state%n_top
                    ipoint = self_link%get_point(self_fdbk%Bu(iz,ia))
                    
                    ! lower
                    Lmat = cr_bank%rods(i_rod)%Lmat(ia)
                    Llen = cr_bank%rods(i_rod)%Llen(ia)
                    nx = lilac_fit%mats(Lmat, ipoint)%NumXvars
                    ny = lilac_fit%mats(Lmat, ipoint)%NumTYvars
                    allocate(Xvars(nx), stat=i_allocate)
                    allocate(Yvars(ny), stat=i_allocate)
                    
                    call self_fdbk%get_state(self_link, iz, ia, Xvars)
                    call do_Link_Interpolate (lilac_fit%mats(Lmat, ipoint), Xvars, nx, Yvars, ny)
                    
                    cnt_xsec = xsec1%count ()
                    cnt_param = param1%count ()
                    call xsec1%assign (Yvars(: cnt_xsec))
                    call param1%assign (Yvars(cnt_xsec+1: cnt_xsec+cnt_param))
                    
                    ! upper 
                    Umat = cr_bank%rods(i_rod)%Umat(ia)
                    Ulen = cr_bank%rods(i_rod)%Ulen(ia)
                    call self_fdbk%get_state(self_link, iz, ia, Xvars)
                    call do_Link_Interpolate (lilac_fit%mats(Umat, ipoint), Xvars, nx, Yvars, ny)
                    
                    cnt_xsec = xsec2%count ()
                    cnt_param = param2%count ()
                    call xsec2%assign (Yvars(: cnt_xsec))
                    call param2%assign (Yvars(cnt_xsec+1: cnt_xsec+cnt_param))
                    
                    if (allocated(Xvars))       deallocate(Xvars)
                    if (allocated(Yvars))       deallocate(Yvars)
                    
                    ! 
                    call cr_bank%xsecs(i_rod, ia)%weight ([xsec1, xsec2], [Llen, Ulen], [Lwt, Uwt])
                    call cr_bank%params(i_rod, ia)%weight ([param1, param2], [Llen, Ulen], [Lwt, Uwt])
                end do
            end do
                    
            call xsec1%clean ()
            call param1%clean ()
            call xsec2%clean ()
            call param2%clean ()
        end if
                
    end subroutine read_xsec_unknown_Lilac_new
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_LilacParameter (this)
        
        class(LilacParameter), intent(in out)  :: this
        integer  :: im, ipoint, it
        
        do im = 1, SIZE(this%mats, dim=1)
            do ipoint = 1, SIZE(this%mats, dim=2)
                associate (imat => this%mats(im, ipoint))
                call post_Link_Interpolate (imat%a_link)
                
                if (allocated(imat%TOption))        deallocate(imat%TOption)
                if (allocated(imat%Terms))  then
                    do it = 1, SIZE(imat%Terms)
                        if (allocated(imat%Terms(it)%Option))       deallocate(imat%Terms(it)%Option)
                        if (allocated(imat%Terms(it)%XOrder))       deallocate(imat%Terms(it)%XOrder)
                        if (allocated(imat%Terms(it)%IndexXvars))   deallocate(imat%Terms(it)%IndexXvars)
                        if (allocated(imat%Terms(it)%XvarsInTerm))  deallocate(imat%Terms(it)%XvarsInTerm)
                        if (allocated(imat%Terms(it)%AV))           deallocate(imat%Terms(it)%AV)
                        if (allocated(imat%Terms(it)%STD))          deallocate(imat%Terms(it)%STD)
                        if (allocated(imat%Terms(it)%OrderMatrix))  deallocate(imat%Terms(it)%OrderMatrix)
                        if (allocated(imat%Terms(it)%xVector))      deallocate(imat%Terms(it)%xVector)
                    end do
                    deallocate(imat%Terms)
                end if
                end associate
            end do 
        end do
        
        this%is_define = .FALSE.
        
    end subroutine Free_LilacParameter
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! Interpolate/Get values of fitted variable Y at given X
    ! <step 1>: read information into memory
    !===============================================================================================
    subroutine pre_Link_Interpolate (imat, Xvars, NumXvars)
    
        type(lilac_parameter_tp), intent(in out)   :: imat
        integer, intent(in)                      :: NumXvars                    ! Lengths of X and Y arrays
        real(8), intent(in)                      :: Xvars(NumXvars)             ! X values of the state

        integer  :: NumTerms, NumTXvars, NumTYRead
        integer  :: NumXinTerm, NumStates, NumCoefs
        integer  :: i, j, k
        integer  :: io_error
        integer  :: unit_
        
        open(newunit=unit_, file=imat%file, status='old', action='read', iostat=io_error)
        read(unit_, *)
        read(unit_, *)
        read(unit_, *)
        
        ! Read in term information
        read (unit_, *) NumTerms, NumTXvars, NumTYRead
        call Link_Define (imat%a_link, NumTYRead, NumTXvars, NumTerms)
        read (unit_, *) (imat%a_link%NameXvars(i), i=1, imat%a_link%NumTXvars)
        read (unit_, *) (imat%a_link%NameYvars(i), i=1, imat%a_link%NumTYvars)
        allocate(imat%TOption(NumTerms))
        allocate(imat%Terms(NumTerms))
        do i = 1, NumTerms
           read(unit_, *) imat%a_link%NameTerm(i), imat%a_link%Parentheses(i), imat%a_link%Operator(i), imat%TOption(i)
        end do
        read(unit_, *)
        read(unit_, *)
        read(unit_, *)
        
        do i = 1, imat%a_link%NumTerms
           associate(iterm => imat%Terms(i))
           
           ! Read in coefficients in each term
           read(unit_, *) NumXinTerm
           iterm%NumXinTerm = NumXinTerm
           
           allocate(iterm%Option(NumXinTerm), iterm%XOrder(NumXinTerm), iterm%IndexXvars(NumXinTerm), iterm%XvarsInTerm(NumXinTerm))
           read(unit_, *) (iterm%Option(j), j=1, NumXinTerm)
           do j = 1, NumXinTerm
              if (iterm%Option(j) /= 1) then
                 write (*,*) '[Link] FATAL ERROR: Sorry, Least Square only, Option must be 1'
              end if
           end do
           read(unit_, *) (iterm%XOrder(j), j=1, NumXinTerm)
           read(unit_, *) (iterm%IndexXvars(j), j=1, NumXinTerm)
           read(unit_, *)
           allocate (iterm%AV(NumXinTerm), iterm%STD(NumXinTerm))
           do j = 1, NumXinTerm
              read(unit_, *) iterm%AV(j), iterm%STD(j)
           end do
           
           ! Store information in User_Link
           NumCoefs = 1
           do j = 1, NumXinTerm
              NumCoefs = NumCoefs * (iterm%XOrder(j)+1)
           end do
           iterm%NumCoefs = NumCoefs
           
           NumStates = 1
           call YXdiscrete_Define (imat%a_link%YX, NumXinTerm, NumStates, NumCoefs)
           call LinearSystem_Define (imat%a_link%LSystem, NumCoefs)
           call YXdiscrete_SetXOptions (imat%a_link%YX, iterm%XOrder, NumXinTerm, iterm%Option)
           read(unit_, *)
           read(unit_, *)
           allocate(iterm%OrderMatrix(NumXinTerm, NumCoefs))
           read(unit_, *) ((iterm%OrderMatrix(k,j), j=1, NumCoefs), k=1, NumXinTerm)
           call YXdiscrete_SetOrderMatrix (imat%a_link%YX, iterm%OrderMatrix, NumXinTerm, NumCoefs)
           
           allocate (iterm%xVector(NumTYRead, NumCoefs))
           do k = 1, NumTYRead
              read(unit_, *)
              read(unit_, *)
              read(unit_, *) (iterm%xVector(k, j) ,j=1, NumCoefs)
           end do
           read(unit_, *)
           
           end associate
        end do
        close (unit_)
       
    end subroutine pre_Link_Interpolate
   
    !$
    !===============================================================================================
    ! Interpolate/Get values of fitted variable Y at given X
    ! <step 2>: do fitting
    !===============================================================================================
    subroutine do_Link_Interpolate (imat, Xvars, NumXvars, Yvars, NumTYvars)
    
        type(lilac_parameter_tp), intent(in out)   :: imat
        integer, intent(in)                      :: NumXvars,NumTYvars          ! Lengths of X and Y arrays
        real(8), intent(in)                      :: Xvars(NumXvars)             ! X values of the state
        real(8), intent(out)                     :: Yvars(NumTYvars)            ! Fitted values of Y to pass out

        integer  :: i, j, k
        integer, allocatable :: XOrders(:)
        
        Yvars = 0.0
        do i = 1, imat%a_link%NumTerms
           associate(iterm => imat%Terms(i))
           
           ! Get XBar from X
           do j = 1, iterm%NumXinTerm
              if (ABS(iterm%STD(j)-0.0D0) > NECP_RealZero)  then
                 iterm%XvarsInTerm(j) = (Xvars(iterm%IndexXvars(j))-iterm%AV(j)) / iterm%STD(j)
              else
                 iterm%XvarsInTerm(j) = Xvars(iterm%IndexXvars(j))-iterm%AV(j)
              end if
           end do
           
           ! Firstly : get term values from X
           allocate(XOrders(iterm%NumXinTerm))
           XOrders = iterm%XOrder
           do k = 1, NumTYvars
              call Link_YfromXC (XOrders, iterm%xVector(k,:), iterm%XvarsInTerm, imat%a_link%Terms(k,i), iterm%NumCoefs, iterm%NumXinTerm)
              if (imat%TOption(i) /= 0) then
                 imat%a_link%Terms(k,i) = Xvars(imat%TOption(i)) * imat%a_link%Terms(k,i)
              end if
           end do
           deallocate(XOrders)
           
           end associate
        end do
        
        ! Secondly : get Y values from terms
        do i = 1, NumTYvars
           call Link_YfromTerms (imat%a_link, Yvars(i), i)
        end do
        
    end subroutine do_Link_Interpolate
    
    !$
    !===============================================================================================
    ! Interpolate/Get values of fitted variable Y at given X
    ! <step 3>: free memory
    !===============================================================================================
    subroutine post_Link_Interpolate (User_Link)
    
        type(Link_Type), intent(in out)  :: User_Link
        
        call  YXdiscrete_Void (User_Link%YX)
        call  LinearSystem_Void (User_Link%LSystem)
        call  Link_Void (User_Link)
        
    end subroutine post_Link_Interpolate
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! get xsec by NEA/CRP benchmark fitting
    !===============================================================================================
    subroutine read_xsec_unknown_NEACRP (self_fdbk, self_link, mat_info, xsec, param, cr_bank)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(CrossSection), intent(in out), optional       :: xsec
        type(KineticsParameter), intent(in out), optional  :: param
        type(ControlRodBank), intent(in out), optional     :: cr_bank
    
        ! local variables
        type(cross_section_info_tp)  :: xsec1, xsec2
        real(KREAL)                  :: Llen, Ulen
        real(KREAL)                  :: Lwt, Uwt
        integer                      :: Lmat, Umat
        
        integer, parameter  :: N_BASE    = 9                                    ! number of base xsec
        integer, parameter  :: N_PARTIAL = 40                                   ! number of partial
        integer  :: i_rod
        integer  :: ia, iz, im
        integer  :: i, j
        integer  :: nx, ny
        integer  :: unit_
        integer  :: i_allocate, io_error
        
        real(8), allocatable     :: Xvars(:)
        real(8), allocatable     :: Yvars(:)
        
        Lwt = 1.0; Uwt = 1.0;
        
        ! set neacrp_fit for the first time
        if (.NOT. neacrp_fit%is_define)  then
            allocate(neacrp_fit%mats(ns%state%mat), stat=i_allocate)
            do im = 1, ns%state%mat
                associate (imat => neacrp_fit%mats(im))
                
                imat%ID = mat_info%libs(im)%ID
                imat%file = TRIM(ADJUSTL(DIR_LINK)) // TRIM(ADJUSTL(mat_info%libs(im)%file))
                imat%is_CR = mat_info%libs(im)%is_CR
                
                allocate(imat%xsec0(N_BASE), stat=i_allocate)
                allocate(imat%partial(N_PARTIAL), stat=i_allocate)
                call Open_file (imat%file, is_input=.TRUE., file_unit=unit_)
                read(unit=unit_, fmt=*) imat%xsec0
                read(unit=unit_, fmt=*) imat%partial
                call Close_file (unit_, is_keep=.TRUE.)   
                
                end associate
            end do
            neacrp_fit%is_define = .TRUE.
        end if
        
        if (PRESENT(xsec) .AND. PRESENT(param))  then
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    im = mat_info%loading(iz, ia)
                    associate (imat => neacrp_fit%mats(im))
                    
                    xsec1 = xsec%matrixs(iz, ia)
                    nx = self_link%n_parameter
                    allocate(Xvars(nx), stat=i_allocate)
                    call self_fdbk%get_state(self_link, iz, ia, Xvars)
                    
                    ny = xsec1%count ()
                    allocate(Yvars(ny), stat=i_allocate)
                    
                    call NEACRP2DAISY (Xvars, imat%xsec0, imat%partial, imat%delataCR, Yvars)
                    call xsec1%assign (Yvars)
                    
                    if (allocated(Xvars))       deallocate(Xvars)
                    if (allocated(Yvars))       deallocate(Yvars)
                    
                    xsec%matrixs(iz, ia) = xsec1
                    end associate
                end do
            end do
        end if 
        
        ! feedback for CR xsec
        if (PRESENT(cr_bank) .AND. allocated(cr_bank%rods))  then
            call xsec1%alloc ()
            call xsec2%alloc ()
                    
            do i_rod = 1, cr_bank%state%n_rod
                iz = cr_bank%rods(i_rod)%zone
                do ia = cr_bank%state%n_bottom+1, cr_bank%state%na-cr_bank%state%n_top
                    ! lower
                    Lmat = cr_bank%rods(i_rod)%Lmat(ia)
                    Llen = cr_bank%rods(i_rod)%Llen(ia)
                    associate (imat => neacrp_fit%mats(Lmat))
                    
                    nx = self_link%n_parameter
                    allocate(Xvars(nx), stat=i_allocate)
                    call self_fdbk%get_state(self_link, iz, ia, Xvars)
                    
                    ny = xsec1%count ()
                    allocate(Yvars(ny), stat=i_allocate)
                    
                    call NEACRP2DAISY (Xvars, imat%xsec0, imat%partial, imat%delataCR, Yvars)
                    call xsec1%assign (Yvars)
                    end associate
                    
                    ! upper 
                    Umat = cr_bank%rods(i_rod)%Umat(ia)
                    Ulen = cr_bank%rods(i_rod)%Ulen(ia)
                    associate (imat => neacrp_fit%mats(Umat))
                    
                    ny = xsec2%count ()
                    call NEACRP2DAISY (Xvars, imat%xsec0, imat%partial, imat%delataCR, Yvars)
                    call xsec2%assign (Yvars)
                    
                    if (allocated(Xvars))       deallocate(Xvars)
                    if (allocated(Yvars))       deallocate(Yvars)
                    end associate
                    
                    ! 
                    call cr_bank%xsecs(i_rod, ia)%weight ([xsec1, xsec2], [Llen, Ulen], [Lwt, Uwt])
                end do
            end do
            
            call xsec1%clean ()
            call xsec2%clean ()
        end if
                
    end subroutine read_xsec_unknown_NEACRP
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine NEACRP2DAISY (xin, xsec0, partial, deltaCR, yin)
    
        real(8), intent(in)  :: xin(:)
        real(KREAL), intent(in)  :: xsec0(:)
        real(KREAL), intent(in)  :: partial(:)
        real(KREAL), intent(in), allocatable  :: deltaCR(:)
        real(8), intent(in out)  :: yin(:)
    
        real(KREAL)  :: tmp(SIZE(xsec0))
        real(KREAL)  :: cb, tm, rhom, tf
        real(KREAL)  :: tr1, a1, vf1, f1, s12
        real(KREAL)  :: tr2, a2, vf2, f2
        
        integer, parameter  :: NG_ = 2
        integer  :: ibeg, iend
        
        cb = xin(1) - partial(10)
        tm = xin(2) - CKELVIN - partial(20)
        rhom = xin(3)/1000.0 - partial(30)
        tf = SQRT(xin(4)) - SQRT(partial(40))
        
        tmp = xsec0 + partial(1:9)*cb + partial(11:19)*tm + partial(21:29)*rhom + partial(31:39)*tf
        
!        if (allocated(deltaCR))  then
!            tmp = tmp + deltaCR
!        end if
        tr1 = tmp(1); a1 = tmp(2); vf1 = tmp(3); f1 = tmp(4); s12 = tmp(5);
        tr2 = tmp(6); a2 = tmp(7); vf2 = tmp(8); f2 = tmp(9); 
        
        ibeg = 1
!        iend = NG_ * 4
!        yin(ibeg: iend) = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]              ! adf
        
!        ibeg = iend + 1                                                         ! chi_steady
        iend = NG_
        yin(ibeg: iend) = [1.0, 0.0]
        
        ibeg = iend + 1                                                         ! sigma_t
        iend = iend + NG_
        yin(ibeg: iend) = [tr1, tr2]
        
        ibeg = iend + 1                                                         ! sigma_f_nu
        iend = iend + NG_
        yin(ibeg: iend) = [vf1, vf2]
        
        ibeg = iend + 1                                                         ! sigma_f_kappa
        iend = iend + NG_
        yin(ibeg: iend) = [f1*3.213E-11, f2*3.206E-11]
        
        ibeg = iend + 1                                                         ! sigma_s
        iend = iend + NG_ * 2
        yin(ibeg: iend) = [tr1-a1-s12, s12, 0.0_KREAL, tr2-a2]

    end subroutine NEACRP2DAISY
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_NEACRPParameter (this)
        
        class(NEACRPParameter), intent(in out)  :: this
        integer  :: im, it
        
        do im = 1, SIZE(this%mats)
            if (allocated(this%mats(im)%xsec0))          deallocate(this%mats(im)%xsec0)
            if (allocated(this%mats(im)%partial))        deallocate(this%mats(im)%partial)
            if (allocated(this%mats(im)%delataCR))       deallocate(this%mats(im)%delataCR)
        end do
        
        this%is_define = .TRUE.
        
    end subroutine Free_NEACRPParameter    
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine read_xsec_unknown_interpolation (self_fdbk, self_link, mat_info, xsec, param, cr_bank)
        
        type(FeedbackParameter), intent(in)      :: self_fdbk
        type(LinkParameter), intent(in)          :: self_link
        type(Material), intent(in out)           :: mat_info
        type(CrossSection), intent(in out), optional       :: xsec
        type(KineticsParameter), intent(in out), optional  :: param
        type(ControlRodBank), intent(in out), optional     :: cr_bank
        
        ! local variables
        type(cross_section_info_tp)       :: xsec1, xsec2
        type(kinetics_parameter_info_tp)  :: param1, param2
        real(KREAL)                  :: Llen, Ulen
        real(KREAL)                  :: Lwt, Uwt
        integer                      :: Lmat, Umat
        
        character(len=MAX_WORD_LEN)  :: file_in
        integer  :: i_rod
        integer  :: ia, iz, im
        integer  :: i, j
        integer  :: nx, ny
        integer  :: unit_
        integer  :: i_allocate, io_error
        
        character(len=MAX_WORD_LEN)  :: file_from, dir
        character(len=MAX_WORD_LEN)  :: groupname
        character(len=MAX_WORD_LEN)  :: tmp
        integer(HID_T) :: file, group, casegroup
        
        integer  :: interp_method = 2 
        
        logical  :: log1, log2 
        integer  :: case_id
        integer  :: hdferr
        
        ! set hdf5_interp for the first time
        if (.NOT. hdf5_interp%is_define)  then
            allocate(hdf5_interp%mats(ns%state%mat), stat=i_allocate)
            
            do im = 1, ns%state%mat
            associate (mt => hdf5_interp%mats(im))
                mt%ID = mat_info%libs(im)%ID
                file_in = TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(EXTEND_H5))
                mt%file = TRIM(ADJUSTL(DIR_LINK)) // TRIM(ADJUSTL(file_in))
                
                call h5open_f(hdferr)
                call hdf5_file_open(TRIM(ADJUSTL(mt%file)), file, 'r')
                
                ! open the top group in the file
                call Get_file_basename(TRIM(file_in), groupname)
                call hdf5_open_group(file, groupname, group)
                call hdf5_get_attr(file, groupname, 'NumLIST', mt%n_case)
                call hdf5_get_attr(file, groupname, 'Fuel_Flag', mt%flag)
                
                allocate(mt%t_fuel(mt%n_case), stat=i_allocate)
                allocate(mt%rho_coolant(mt%n_case), stat=i_allocate)
                allocate(mt%xsecs(mt%n_case), stat=i_allocate)
                allocate(mt%params(mt%n_case), stat=i_allocate)
                
                do case_id = 1, mt%n_case
                    groupname = 'CASE' // TRIM(Int_to_string(case_id, digit=4))
                    call hdf5_open_group(group, groupname, casegroup)
                    call hdf5_get_attr(group, groupname, 'TF', mt%t_fuel(case_id))
                    call hdf5_get_attr(group, groupname, 'DC', mt%rho_coolant(case_id))
                    call hdf5_close_group(casegroup)
                end do 
                
                call hdf5_close_group(group)
                
                ! close file
                call hdf5_file_close(file)
                call h5close_f(hdferr)
                
                ! 
                call sort_vector (mt%t_fuel, mt%xTf, mt%n_xTf)
                call sort_vector (mt%rho_coolant, mt%yDc, mt%n_yDc)
                
                if (mt%n_case /= mt%n_xTf*mt%n_yDc)  then
                    call a_error%set (INFO_LIST_INPUT, 'xsec point number is not correct: '//TRIM(ADJUSTL(mt%file)))  
                    call a_error%print (FILES%MAIN)
                end if 
                
                allocate(mt%xyto1D(mt%n_xTf, mt%n_yDc), stat=i_allocate)
                do i = 1, mt%n_xTf 
                    do j = 1, mt%n_yDc
                        do case_id = 1, mt%n_case
                            log1 = ABS(mt%xTf(i) - mt%t_fuel(case_id)) <= EPS_ZERO
                            log2 = ABS(mt%yDc(j) - mt%rho_coolant(case_id)) <= EPS_ZERO
                            if (log1 .and. log2)  then
                                exit 
                            end if 
                        end do 
                        
                        mt%xyto1D(i, j) = case_id
                    end do 
                end do 
                
                do case_id = 1, mt%n_case
                    call mt%xsecs(case_id)%alloc ()
                    call mt%params(case_id)%alloc ()
                    call HDF5_to_XSEC (DIR_LINK, file_in, case_id, 1, mt%xsecs(case_id), mt%params(case_id))
                end do 
                
            end associate 
            end do
            
            hdf5_interp%is_define = .TRUE.
        end if
        
        !
        call xsec_00%alloc (); call param_00%alloc ();
        call xsec_01%alloc (); call param_01%alloc ();
        call xsec_10%alloc (); call param_10%alloc ();
        call xsec_11%alloc (); call param_11%alloc ();
        call xsec_a%alloc (); call param_a%alloc ();
        call xsec_b%alloc (); call param_b%alloc ();
        call xsec_c%alloc (); call param_c%alloc ();
        
        if (PRESENT(xsec) .AND. PRESENT(param))  then
            call xsec1%alloc ()
            call param1%alloc ()
            
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    im = mat_info%loading(iz, ia)
                    file_in = TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(EXTEND_H5))
                    
                    select case(interp_method)
                    case (INTERP_LINNEAR)
                        call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, im, xsec1, param1)
                    case (INTERP_PICEWISL)
                        call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, im, xsec1, param1)
                    case (INTERP_LAGRANGE)
                        call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, im, xsec1, param1)
                    case default 
                    end select 
                    xsec%matrixs(iz, ia) = xsec1
                    param%matrixs(iz, ia) = param1
                end do
            end do
            
            call xsec1%clean ()
            call param1%clean ()
        end if 
        
        ! feedback for CR xsec
        Lwt = 1.0; Uwt = 1.0;
        if (PRESENT(cr_bank) .AND. allocated(cr_bank%rods))  then
            call xsec1%alloc ()
            call xsec2%alloc ()
            call param1%alloc ()
            call param2%alloc ()
                    
            do i_rod = 1, cr_bank%state%n_rod
                iz = cr_bank%rods(i_rod)%zone
                do ia = cr_bank%state%n_bottom+1, cr_bank%state%na-cr_bank%state%n_top
                    ! lower
                    Lmat = cr_bank%rods(i_rod)%Lmat(ia)
                    Llen = cr_bank%rods(i_rod)%Llen(ia)
                    file_in = TRIM(ADJUSTL(mat_info%libs(Lmat)%file)) // TRIM(ADJUSTL(EXTEND_H5))
                    select case(interp_method)
                    case (INTERP_LINNEAR)
                        call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, Lmat, xsec1, param1)
                    case (INTERP_PICEWISL)
                        call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, Lmat, xsec1, param1)
                    case (INTERP_LAGRANGE)
                        call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, Lmat, xsec1, param1)
                    case default 
                    end select
                    
                    ! upper 
                    Umat = cr_bank%rods(i_rod)%Umat(ia)
                    Ulen = cr_bank%rods(i_rod)%Ulen(ia)
                    file_in = TRIM(ADJUSTL(mat_info%libs(Umat)%file)) // TRIM(ADJUSTL(EXTEND_H5))
                    select case(interp_method)
                    case (INTERP_LINNEAR)
                        call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, Umat, xsec2, param2)
                    case (INTERP_PICEWISL)
                        call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, Umat, xsec2, param2)
                    case (INTERP_LAGRANGE)
                        call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, Umat, xsec2, param2)
                    case default 
                    end select 
                    
                    ! 
                    call cr_bank%xsecs(i_rod, ia)%weight ([xsec1, xsec2], [Llen, Ulen], [Lwt, Uwt])
                    call cr_bank%params(i_rod, ia)%weight ([param1, param2], [Llen, Ulen], [Lwt, Uwt])
                end do
            end do
            
            call xsec1%clean ()
            call xsec2%clean ()
            call param1%clean ()
            call param2%clean ()
        end if
        
        call xsec_00%clean (); call param_00%clean ();
        call xsec_01%clean (); call param_01%clean ();
        call xsec_10%clean (); call param_10%clean ();
        call xsec_11%clean (); call param_11%clean ();
        call xsec_a%clean (); call param_a%clean ();
        call xsec_b%clean (); call param_b%clean ();
        call xsec_c%clean (); call param_c%clean ();
            
    end subroutine read_xsec_unknown_interpolation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
        
        type(FeedbackParameter), intent(in)     :: self_fdbk
        type(LinkParameter), intent(in)         :: self_link
        character(len=*)                        :: file_in
        integer, intent(in)                     :: iz
        integer, intent(in)                     :: ia
        integer, intent(in)                     :: im 
        type(cross_section_info_tp), intent(in out)       :: a_xsec
        type(kinetics_parameter_info_tp), intent(in out)  :: a_param
        
        type(state_variables_tp)  :: value_refer
        type(state_variables_tp)  :: value_current
        type(state_variables_tp)  :: value_00
        type(state_variables_tp)  :: value_01
        type(state_variables_tp)  :: value_10
        type(state_variables_tp)  :: value_11
        
        integer  :: index_00
        integer  :: index_01
        integer  :: index_10
        integer  :: index_11
        
        real(KREAL)  :: x_min, x, x_max
        real(KREAL)  :: y_min, y, y_max
        
        real(8)  :: Xvars(2)
        integer  :: nx 
        integer  :: case_id
          
        nx = self_link%n_parameter
        call self_fdbk%get_state(self_link, iz, ia, Xvars)
        value_refer%t_fuel = Xvars(1)
        value_refer%rho_coolant = Xvars(2) 
        
        if (value_refer%t_fuel >= self_link%info(1)%max)  then
            value_refer%t_fuel = self_link%info(1)%max
        end if 
        if (value_refer%t_fuel <= self_link%info(1)%min)  then
            value_refer%t_fuel = self_link%info(1)%min
        end if 
        if (value_refer%rho_coolant >= self_link%info(2)%max)  then
            value_refer%rho_coolant = self_link%info(2)%max
        end if 
        if (value_refer%rho_coolant <= self_link%info(2)%min)  then
            value_refer%rho_coolant = self_link%info(2)%min
        end if 

        index_00 = -1
        value_00%t_fuel = self_link%info(1)%min
        value_00%rho_coolant = self_link%info(2)%min 

        index_01 = -1
        value_01%t_fuel = self_link%info(1)%min
        value_01%rho_coolant = self_link%info(2)%max 

        index_10 = -1
        value_10%t_fuel = self_link%info(1)%max
        value_10%rho_coolant = self_link%info(2)%min 

        index_11 = -1
        value_11%t_fuel = self_link%info(1)%max
        value_11%rho_coolant = self_link%info(2)%max 
        
        do case_id = 1, hdf5_interp%mats(im)%n_case
            value_current%t_fuel = hdf5_interp%mats(im)%t_fuel(case_id)
            value_current%rho_coolant = hdf5_interp%mats(im)%rho_coolant(case_id)
            
            if (value_current%in (value_refer, value_00))  then
                value_00 = value_current
                index_00 = case_id
            end if 
            if (value_current%in (value_refer, value_01))  then
                value_01 = value_current
                index_01 = case_id
            end if 
            if (value_current%in (value_refer, value_10))  then
                value_10 = value_current
                index_10 = case_id
            end if 
            if (value_current%in (value_refer, value_11))  then
                value_11 = value_current
                index_11 = case_id
            end if
        end do

        if ((index_00 < 0) .or. (index_01 < 0) .or. (index_10 < 0) .or. (index_11 < 0)) then 
            write(OUTPUT_UNIT, fmt="(1x, *(A, I6))") 'iz = ', iz, 'ia = ', ia
            write(OUTPUT_UNIT, fmt="(1x, *(A, ES12.5))") 'Tf = ', value_refer%t_fuel, ', Rhom = ', value_refer%rho_coolant
            write(OUTPUT_UNIT, fmt="(1x, *(I3, TR2))") index_00, index_01, index_10, index_11
            stop('0')
        end if 
        
        ! get xsec and interpolation
        xsec_00 = hdf5_interp%mats(im)%xsecs(index_00)
        xsec_01 = hdf5_interp%mats(im)%xsecs(index_01)
        xsec_10 = hdf5_interp%mats(im)%xsecs(index_10)
        xsec_11 = hdf5_interp%mats(im)%xsecs(index_11)
        
        param_00 = hdf5_interp%mats(im)%params(index_00)
        param_01 = hdf5_interp%mats(im)%params(index_01)
        param_10 = hdf5_interp%mats(im)%params(index_10)
        param_11 = hdf5_interp%mats(im)%params(index_11)
        
        x_min = value_00%t_fuel
        x_max = value_11%t_fuel
        x = value_refer%t_fuel

        y_min = value_00%rho_coolant
        y_max = value_11%rho_coolant
        y = value_refer%rho_coolant
        
        if (.TRUE.)  then
            call a_xsec%interpolate (xsec_00, xsec_01, xsec_10, xsec_11, x_min, x, x_max, y_min, y, y_max)
            call a_param%interpolate (param_00, param_01, param_10, param_11, x_min, x, x_max, y_min, y, y_max)
            
        else
            call two_variables_interpolation (a_xsec%chi_steady, xsec_00%chi_steady, xsec_01%chi_steady, xsec_10%chi_steady, xsec_11%chi_steady,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%sigma_t, xsec_00%sigma_t, xsec_01%sigma_t, xsec_10%sigma_t, xsec_11%sigma_t,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%sigma_f_nu, xsec_00%sigma_f_nu, xsec_01%sigma_f_nu, xsec_10%sigma_f_nu, xsec_11%sigma_f_nu,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%sigma_f_kappa, xsec_00%sigma_f_kappa, xsec_01%sigma_f_kappa, xsec_10%sigma_f_kappa, xsec_11%sigma_f_kappa,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%sigma_s, xsec_00%sigma_s, xsec_01%sigma_s, xsec_10%sigma_s, xsec_11%sigma_s,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%sigma_f, xsec_00%sigma_f, xsec_01%sigma_f, xsec_10%sigma_f, xsec_11%sigma_f,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%nu, xsec_00%nu, xsec_01%nu, xsec_10%nu, xsec_11%nu,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_xsec%kappa, xsec_00%kappa, xsec_01%kappa, xsec_10%kappa, xsec_11%kappa,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
        
            call two_variables_interpolation (a_param%chi_delay, param_00%chi_delay, param_01%chi_delay, param_10%chi_delay, param_11%chi_delay,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_param%velocity, param_00%velocity, param_01%velocity, param_10%velocity, param_11%velocity,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_param%beta, param_00%beta, param_01%beta, param_10%beta, param_11%beta,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_param%lambda, param_00%lambda, param_01%lambda, param_10%lambda, param_11%lambda,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_param%sigma_bvf, param_00%sigma_bvf, param_01%sigma_bvf, param_10%sigma_bvf, param_11%sigma_bvf,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
            call two_variables_interpolation (a_param%sigma_bvl, param_00%sigma_bvl, param_01%sigma_bvl, param_10%sigma_bvl, param_11%sigma_bvl,  &
                                            &   x_min, x, x_max, y_min, y, y_max)
        end if 
        
    end subroutine xsec_interpolation_lin
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
        
        type(FeedbackParameter), intent(in)     :: self_fdbk
        type(LinkParameter), intent(in)         :: self_link
        character(len=*)                        :: file_in
        integer, intent(in)                     :: iz
        integer, intent(in)                     :: ia
        integer, intent(in)                     :: im 
        type(cross_section_info_tp), intent(in out)       :: a_xsec
        type(kinetics_parameter_info_tp), intent(in out)  :: a_param
        
        real(8)  :: Xvars(2)
        integer  :: nx 
        real(KREAL) :: xi
        real(KREAL) :: yi
        
        real(KREAL) :: alpha
        real(KREAL) :: beta
        real(KREAL) :: det
        real(KREAL) :: dxa
        real(KREAL) :: dxb
        real(KREAL) :: dxi
        real(KREAL) :: dya
        real(KREAL) :: dyb
        real(KREAL) :: dyi
        real(KREAL) :: gamma
        integer :: i, j
        
        associate (mt => hdf5_interp%mats(im), xd => hdf5_interp%mats(im)%xTf, ydin => hdf5_interp%mats(im)%yDc)
            nx = self_link%n_parameter
            call self_fdbk%get_state(self_link, iz, ia, Xvars)
            xi = Xvars(1)
            yi = Xvars(2)
            
            if (xi >= MAXVAL(xd))   xi = MAXVAL(xd)
            if (yi >= MAXVAL(ydin))   yi = MAXVAL(ydin)
            if (xi <= MINVAL(xd))   xi = MINVAL(xd)
            if (yi <= MINVAL(ydin))   yi = MINVAL(ydin)
            
            i = r8vec_bracket5 (mt%n_xTf, xd, xi)
            j = r8vec_bracket5 (mt%n_yDc, ydin, yi)
            
            if (yi < ydin(j+1) - (ydin(j+1) - ydin(j)) * (xi - xd(i)) / (xd(i+1) - xd(i)))  then
                dxa = xd(i+1) - xd(i)
                dya = ydin(j)   - ydin(j)
                dxb = xd(i)   - xd(i)
                dyb = ydin(j+1) - ydin(j)
                dxi = xi - xd(i)
                dyi = yi - ydin(j)
                det = dxa * dyb - dya * dxb
                alpha = (dxi * dyb - dyi * dxb) / det
                beta =  (dxa * dyi - dya * dxi) / det
                gamma = 1.0D+00 - alpha - beta
!!!!                zi = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)
                
                xsec_a = mt%xsecs(mt%xyto1D(i+1,j)) .MULTI. alpha
                xsec_b = mt%xsecs(mt%xyto1D(i,j+1)) .MULTI. beta
                xsec_c = mt%xsecs(mt%xyto1D(i,j)) .MULTI. gamma
                param_a = mt%params(mt%xyto1D(i+1,j)) .MULTI. alpha
                param_b = mt%params(mt%xyto1D(i,j+1)) .MULTI. beta
                param_c = mt%params(mt%xyto1D(i,j)) .MULTI. gamma
            else
                dxa = xd(i)   - xd(i+1)
                dya = ydin(j+1) - ydin(j+1)
                dxb = xd(i+1) - xd(i+1)
                dyb = ydin(j)   - ydin(j+1)
                dxi = xi - xd(i+1)
                dyi = yi - ydin(j+1)
                det = dxa * dyb - dya * dxb
                alpha = (dxi * dyb - dyi * dxb) / det
                beta =  (dxa * dyi - dya * dxi) / det
                gamma = 1.0D+00 - alpha - beta
!!!!                zi = alpha * zd(i,j+1) + beta * zd(i+1,j) + gamma * zd(i+1,j+1)
                
                xsec_a = mt%xsecs(mt%xyto1D(i,j+1)) .MULTI. alpha
                xsec_b = mt%xsecs(mt%xyto1D(i+1,j)) .MULTI. beta
                xsec_c = mt%xsecs(mt%xyto1D(i+1,j+1)) .MULTI. gamma
                param_a = mt%params(mt%xyto1D(i,j+1)) .MULTI. alpha
                param_b = mt%params(mt%xyto1D(i+1,j)) .MULTI. beta
                param_c = mt%params(mt%xyto1D(i+1,j+1)) .MULTI. gamma
            end if
            a_xsec = xsec_a .ADD. xsec_b .ADD. xsec_c
            a_param = param_a .ADD. param_b .ADD. param_c
        end associate 
        
    end subroutine xsec_interpolation_pwl
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function r8vec_bracket5 (nd, xd, xi)  result (idx)

        integer, intent(in) :: nd
        real(KREAL), intent(in) :: xd(nd)
        real(KREAL), intent(in) :: xi
        integer :: idx
        
        integer :: b
        integer :: l
        integer :: m
        integer :: r
        
        if (xi < xd(1) .or. xd(nd) < xi) then
            b = -1
        else
            l = 1
            r = nd
            do while (l + 1 < r)
                m = (l + r) / 2
                if (xi < xd(m)) then
                    r = m
                else
                    l = m
                end if
            end do
            b = l
        end if
        
        idx = b
        
    end function r8vec_bracket5    

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
        
        type(FeedbackParameter), intent(in)     :: self_fdbk
        type(LinkParameter), intent(in)         :: self_link
        character(len=*)                        :: file_in
        integer, intent(in)                     :: iz
        integer, intent(in)                     :: ia
        integer, intent(in)                     :: im 
        type(cross_section_info_tp), intent(in out)       :: a_xsec
        type(kinetics_parameter_info_tp), intent(in out)  :: a_param
        
        real(8)  :: Xvars(2)
        real(8)  :: xx, yy 
        real(8), allocatable :: a(:)
        real(8), allocatable :: b(:)
        integer(4) :: i, j
        integer(4) :: m, n
        integer(4) :: nd
        integer(4) :: ni
        integer(4), allocatable :: n_1d(:)
        real(8), allocatable :: x_1d(:)
        real(8), allocatable :: value(:)
        real(8), allocatable :: xd(:,:)
        real(8), allocatable :: xi(:,:)
        real(8), allocatable :: zd(:)
        real(8), allocatable :: ze(:)
        real(8), allocatable :: zi(:)
        real(8), allocatable :: w(:)
        integer  :: ii, jj
        integer  :: i_allocate
        
        associate (mt => hdf5_interp%mats(im), xdin => hdf5_interp%mats(im)%xTf, ydin => hdf5_interp%mats(im)%yDc)
            m = self_link%n_parameter
            allocate (n_1d(m), stat=i_allocate)
            allocate (a(m), stat=i_allocate)
            allocate (b(m), stat=i_allocate)
            n_1d(1: 2) = [mt%n_xTf, mt%n_yDc]
            a(1: 2) = [MINVAL(xdin), MINVAL(ydin)]
            b(1: 2) = [MAXVAL(xdin), MAXVAL(ydin)]
            nd = PRODUCT(n_1d)
            
            allocate(xd(m, nd), stat=i_allocate)
            allocate(zd(nd), stat=i_allocate)
            allocate(w(nd), stat=i_allocate)
            
            ! lagrange_interp_nd_grid
            xd = 0.0D0
            do i = 1, m
                n = n_1d(i)
                allocate(x_1d(n), stat=i_allocate)
                if (i == 1) then
                    x_1d = xdin
                else
                    x_1d = ydin
                end if 
!                x_1d(1:n) = 0.5D0 * ((1.0D0 - x_1d(1:n)) * a(i) + (1.0D0 + x_1d(1:n)) * b(i))
                call r8vec_direct_product (i, n, x_1d, m, nd, xd)
                if (allocated(x_1d))        deallocate(x_1d)
            end do 
            
            ni = 1 
            allocate(xi(m, ni), stat=i_allocate)
            allocate(ze(ni), stat=i_allocate)
            allocate(zi(ni), stat=i_allocate)
            
            call self_fdbk%get_state(self_link, iz, ia, Xvars)
            xx = Xvars(1)
            yy = Xvars(2)
            if (xx >= MAXVAL(xdin))   xx = MAXVAL(xdin)
            if (yy >= MAXVAL(ydin))   yy = MAXVAL(ydin)
            if (xx <= MINVAL(xdin))   xx = MINVAL(xdin)
            if (yy <= MINVAL(ydin))   yy = MINVAL(ydin)
            xi(:, 1) = [xx, yy]
            
            ! lagrange_interp_nd_value
            do j = 1, ni
                w(1:nd) = 1.0D0
                do i = 1, m
                    n = n_1d(i)
                    allocate(x_1d(n), stat=i_allocate)
                    allocate (value(n), stat=i_allocate)
                    if (i == 1) then
                        x_1d = xdin
                    else
                        x_1d = ydin
                    end if 
!                    x_1d(1:n) = 0.5D+00 * ((1.0D+00 - x_1d(1:n)) * a(i) + (1.0D+00 + x_1d(1:n)) * b(i))
                    call lagrange_basis_1d (n, x_1d, 1, xi(i,j), value)
                    call r8vec_direct_product2 (i, n, value, m, nd, w)
                    
                    if (allocated(x_1d))        deallocate (x_1d)
                    if (allocated(value))       deallocate (value)
                end do
                
!                zi(j) = dot_product (w, zd)
                
                a_xsec = mt%xsecs(1) .MULTI. w(1)
                a_param = mt%params(1) .MULTI. w(1)
                do i = 2, nd 
                    xsec_a = mt%xsecs(i) .MULTI. w(i)
                    param_a = mt%params(i) .MULTI. w(i)
                    a_xsec = a_xsec .ADD. xsec_a
                    a_param = a_param .ADD. param_a
                end do 
            end do
        end associate 
        
        if (allocated(a))        deallocate(a)
        if (allocated(b))        deallocate(b)
        if (allocated(n_1d))     deallocate(n_1d)
        if (allocated(xd))       deallocate(xd)
        if (allocated(xi))       deallocate(xi)
        if (allocated(zd))       deallocate(zd)
        if (allocated(ze))       deallocate(ze)
        if (allocated(zi))       deallocate(zi)
        if (allocated(w))        deallocate(w)
        
    end subroutine xsec_interpolation_lag
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Equal_state_variables_tp (left, right)
        
        class(state_variables_tp), intent(in out)  :: left
        type(state_variables_tp), intent(in)       :: right
        
        left%t_fuel = right%t_fuel
        left%rho_coolant = right%rho_coolant
    
    end subroutine Equal_state_variables_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_HDF5Parameter (this)
        
        class(HDF5Parameter), intent(in out)  :: this
        integer  :: im, ipoint
        
        do im = 1, SIZE(this%mats)
            if (allocated(this%mats(im)%t_fuel))        deallocate(this%mats(im)%t_fuel)
            if (allocated(this%mats(im)%rho_coolant))   deallocate(this%mats(im)%rho_coolant)
            if (allocated(this%mats(im)%xTf))           deallocate(this%mats(im)%xTf)
            if (allocated(this%mats(im)%yDc))           deallocate(this%mats(im)%yDc)
            if (allocated(this%mats(im)%xyto1D))        deallocate(this%mats(im)%xyto1D)
            
            if (allocated(this%mats(im)%xsecs))  then
                do ipoint = 1, SIZE(this%mats(im)%xsecs)
                    call this%mats(im)%xsecs(ipoint)%clean ()
                end do 
                deallocate(this%mats(im)%xsecs)
            end if 
            
            if (allocated(this%mats(im)%params))  then
                do ipoint = 1, SIZE(this%mats(im)%params)
                    call this%mats(im)%params(ipoint)%clean ()
                end do 
                deallocate(this%mats(im)%params)
            end if 
        end do
        
        this%is_define = .FALSE.
    
    end subroutine Free_HDF5Parameter   

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_HDF5Parameter (this, unit_)
        
        class(HDF5Parameter), intent(in)  :: this
        integer, intent(in)               :: unit_ 
        
        integer  :: im, j, k
        
        write(unit=unit_, fmt=*)  this%is_define, SIZE(this%mats)
        do im = 1, SIZE(this%mats)
            associate(mat => this%mats(im))
            write(unit=unit_, fmt='(1x, I4, TR3, A, I4, I4)')  mat%ID, TRIM(mat%file), mat%n_case, mat%flag
            write(unit=unit_, fmt='(1x, *(ES12.5, TR3))')  mat%t_fuel
            write(unit=unit_, fmt='(1x, *(ES12.5, TR3))')  mat%rho_coolant
            end associate
        end do 
    
    end subroutine Print_HDF5Parameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Is_in_field (this, p1, p2) result(is_true) 
        
        class(state_variables_tp), intent(in)  :: this
        type(state_variables_tp), intent(in)   :: p1
        type(state_variables_tp), intent(in)   :: p2
        logical  :: is_true
        
        is_true = .FALSE.
        
        if (((this%t_fuel-p1%t_fuel) * (this%t_fuel-p2%t_fuel) <= 0.0)  .AND.  &
            &   ((this%rho_coolant-p1%rho_coolant) * (this%rho_coolant-p2%rho_coolant) <= 0.0))  then
            is_true = .TRUE.
        end if
    
    end function Is_in_field
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine read_xsec_LRA_all (self_lra, mat_info, xsec_inp, xsec, param, cr_bank)
        
        type(LRAmodel), intent(in out)           :: self_lra
        type(Material), intent(in out)           :: mat_info
        type(CrossSectionInput), intent(in)      :: xsec_inp
        type(CrossSection), intent(in out)       :: xsec
        type(KineticsParameter), intent(in out)  :: param
        type(ControlRodBank), intent(in out)     :: cr_bank
        
        call read_xsec_LRA_core (self_lra, mat_info, xsec_inp, xsec, param)
        call read_xsec_LRA_CR (self_lra, mat_info, xsec_inp, cr_bank)
        
    end subroutine read_xsec_LRA_all
    
    !$
    !===============================================================================================
    ! get xsec by LRA benchmark 
    !===============================================================================================
    subroutine read_xsec_LRA_core (self_lra, mat_info, xsec_inp, xsec, param)
        
        type(LRAmodel), intent(in out)           :: self_lra
        type(Material), intent(in out)           :: mat_info
        type(CrossSectionInput), intent(in)      :: xsec_inp
        type(CrossSection), intent(in out)       :: xsec
        type(KineticsParameter), intent(in out)  :: param
    
        ! local variables
        type(cross_section_info_tp)  :: xsec1
        
        integer  :: ia, iz, im
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = mat_info%loading(iz, ia)
                xsec1 = xsec%matrixs(iz, ia)
                
                call self_lra%update_xsec (ns, iz, ia, im, xsec1)
                xsec%matrixs(iz, ia) = xsec1
            end do
        end do
                
    end subroutine read_xsec_LRA_core
    
    !$
    !===============================================================================================
    ! get xsec by LRA benchmark 
    !===============================================================================================
    subroutine read_xsec_LRA_CR (self_lra, mat_info, xsec_inp, cr_bank)
        
        type(LRAmodel), intent(in out)           :: self_lra
        type(Material), intent(in out)           :: mat_info
        type(CrossSectionInput), intent(in)      :: xsec_inp
        type(ControlRodBank), intent(in out)     :: cr_bank
    
        ! local variables
        type(cross_section_info_tp)  :: xsec1, xsec2
        real(KREAL)                  :: Llen, Ulen
        real(KREAL)                  :: Lwt, Uwt
        integer                      :: Lmat, Umat
        
        integer  :: i_rod
        integer  :: ia, iz, im
        
        Lwt = 1.0; Uwt = 1.0;
        
        if (allocated(cr_bank%rods))  then
            do i_rod = 1, cr_bank%state%n_rod
                iz = cr_bank%rods(i_rod)%zone
                do ia = cr_bank%state%n_bottom+1, cr_bank%state%na-cr_bank%state%n_top
                    ! lower
                    Lmat = cr_bank%rods(i_rod)%Lmat(ia)
                    Llen = cr_bank%rods(i_rod)%Llen(ia)
                    
                    xsec1 = xsec_inp%mats(Lmat)
                    call self_lra%update_xsec (ns, iz, ia, Lmat, xsec1)
                    
                    ! upper 
                    Umat = cr_bank%rods(i_rod)%Umat(ia)
                    Ulen = cr_bank%rods(i_rod)%Ulen(ia)
                    
                    xsec2 = xsec_inp%mats(Umat)
                    call self_lra%update_xsec (ns, iz, ia, Umat, xsec2)
                    
                    ! 
                    call cr_bank%xsecs(i_rod, ia)%weight ([xsec1, xsec2], [Llen, Ulen], [Lwt, Uwt])
                end do
            end do
        end if
                
    end subroutine read_xsec_LRA_CR
    
end module input_xsec
