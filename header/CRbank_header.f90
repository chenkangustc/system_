!$
!===================================================================================================
!
!   class for control rod information
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          ControlRodBank
!
!===================================================================================================
module CRbank_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,               only : WarningCollector
    use state_header,                   only : SteadyState
    use geometry_header,                only : Geometry
    use material_header,                only : Material, CrossSection, KineticsParameter,       &
                                        &   CrossSectionInput, KineticsParameterInput,          &
                                        &   cross_section_info_tp, kinetics_parameter_info_tp
    
    implicit none 
    private
    public  :: ControlRodBank
    
    ! --------------------------------------------------------------------------
    ! for control rod illustration:
    ! --------------------------------------------------------------------------
    !          Nodal:        CRrod:
    !           ____          ____                                              
    !          | 07 |        |ref |  <--upper-reflector                         
    !          |____|        |____|                                             
    !          | 06 |        |M1+ |                                  
    !          |____|        |CR  |                                             
    !          | 05 |        |____|  
    !          |____|        |M2+ |  
    !          | 04 |        |CR  |                                             
    !          |____|        |    |  <--CRrod-mat                               
    !          | 03 |        |    |                                             
    !          |____|        |____|  
    !          | 02 |        |M3  |                                             
    !          |____|        |____|                                  
    !          | 01 |        |ref |                                         
    !          |____|        |____|  <--lower-reflector                          
    ! --------------------------------------------------------------------------
    ! Note:
    !      consider CR nodal as less as possible
    ! --------------------------------------------------------------------------
    
    ! --------------------------------------------------------------------------
    ! type for control rod information per rod
    type, private  :: control_rod_info_tp
        integer, public                                  :: bank                ! control rod bank index
        integer, public                                  :: zone                ! material zone for control rod
        integer, public                                  :: upper_mat           ! uppper following material 
        integer, public                                  :: lower_mat           ! lower following material
        logical, public                                  :: is_gray             ! is part insert nodal ?
        character(len=MAX_WORD_LEN), public, allocatable :: flags(:)            ! flags for axial material type
        real(KREAL), public                              :: black_beg = -1.0
        real(KREAL), public                              :: black_end = -1.0
        real(KREAL), public                              :: hblack = -HUGE(0.0) ! heigh of black material
        integer, public, allocatable                     :: Umat(:)             ! material ID for the upper half of nodal
        integer, public, allocatable                     :: Lmat(:)             ! material ID for the lower half of nodal
        real(KREAL), public, allocatable                 :: Ulen(:)             ! height of the upper half
        real(KREAL), public, allocatable                 :: Llen(:)             ! height of the lower half
        real(KREAL), public, allocatable                 :: height(:)           ! height of the total nodal
    contains
        procedure, public  :: alloc => Alloc_control_rod_info_tp
        procedure, public  :: clean => Free_control_rod_info_tp
    end type control_rod_info_tp
    
    ! type for control state information
    type, private  :: control_rod_state_tp
        integer, public                         :: na
        integer, public                         :: n_top                        ! layer of upper reflector
        integer, public                         :: n_bottom                     ! layer of lower reflector
        integer, public                         :: n_rod                        ! number of control rods
        integer, public                         :: n_bank                       ! number of control rod banks
        real(KREAL), public                     :: cr_gap                       ! CR gap between lower CR and upper reflector
        real(KREAL), public                     :: max_step                     ! max step
        real(KREAL), public                     :: min_step                     ! min step
        real(KREAL), public                     :: step_size                    ! move step size in (cm)
        real(KREAL), public                     :: h_upper                      ! height of upper reflector, in cm
        real(KREAL), public                     :: h_active                     ! height of core active, in cm
        real(KREAL), public                     :: h_lower                      ! height of lower reflector, in cm
    end type control_rod_state_tp
    
    ! type for axial nodal symbol
    type, private  :: control_rod_symbol_tp
        character(len=MAX_WORD_LEN), public  :: WHITE     = 'WHITE'             ! normal fuel nodal
        character(len=MAX_WORD_LEN), public  :: GRAY      = 'GRAY'              ! half insert nodal
        character(len=MAX_WORD_LEN), public  :: BLACK     = 'BLACK'             ! control rod nodal
        character(len=MAX_WORD_LEN), public  :: REFLECTOR = 'REFLECTOR'         ! reflector nodal
    end type control_rod_symbol_tp
    
    ! --------------------------------------------------------------------------
    ! type for control rod bank information
    type  ControlRodBank 
        type(control_rod_symbol_tp), public                    :: symbol
        type(control_rod_state_tp), public                     :: state
        type(control_rod_info_tp), public, allocatable         :: rods(:)       ! information per CR rod
        type(cross_section_info_tp), public, allocatable       :: xsecs(:, :)   ! CR xsec
        type(kinetics_parameter_info_tp), public, allocatable  :: params(:, :)  ! CR kinetic
        real(KREAL), public, allocatable                       :: init_step(:)  ! initial rod step per bank
        real(KREAL), public, allocatable                       :: last_step(:)  ! last rod step per bank
        logical, public                        :: is_search    = .FALSE.        ! is critical CR search 
        real(KREAL), public                    :: search_value = 1.0D0          ! search target
        real(KREAL), public                    :: search_limit = 3.0D-5         ! search target
        integer, public                        :: bank_idx     = 1              ! search bank 
    contains
        procedure, public  :: alloc => Allocate_ControlRodBank
        procedure, public  :: clean => Free_ControlRodBank
        procedure, public  :: print_xsec => Print_ControlRodBank_xsec
        procedure, public  :: print_conf => Print_ControlRodBank
        procedure, public  :: critical => Adjust_ControlRodBank_xsec
        procedure, public  :: set_hblack => Set_ControlRodBank_hblack
        procedure, public  :: init => Set_ControlRodBank_init_position
        procedure, public  :: black => Set_ControlRodBank_black
        procedure, public  :: set => Set_ControlRodBank
        procedure, public  :: set_xsec => Set_ControlRodBank_xsec
        procedure, public  :: move => Move_ControlRodBank
        procedure, public  :: homo => Homogeneous_ControlRodBank
        procedure, public  :: map => Map_CR2Core
        generic, public    :: assignment (=) => Equal_ControlRodBank
        procedure, private :: Equal_ControlRodBank
        procedure, public  :: rod_search => Search_ControlRodBank
    end type ControlRodBank
    
    ! --------------------------------------------------------------------------
    type(WarningCollector)  :: a_warning
    
contains   
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_control_rod_info_tp (this, ns)
        
        class(control_rod_info_tp), intent(in out)  :: this
        type(SteadyState), intent(in)               :: ns
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%flags(ns%state%layer), stat=i_allocate)
        allocate(this%Umat(ns%state%layer), stat=i_allocate)
        allocate(this%Lmat(ns%state%layer), stat=i_allocate)
        allocate(this%Ulen(ns%state%layer), stat=i_allocate)
        allocate(this%Llen(ns%state%layer), stat=i_allocate)
        allocate(this%height(ns%state%layer), stat=i_allocate)
        
        this%flags = CHAR_NULL
        this%Umat = INT_ONE
        this%Lmat = INT_ONE
        this%Ulen = REAL_ONE
        this%Llen = REAL_ZERO
        this%height = REAL_ONE
    
    end subroutine Alloc_control_rod_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of control_rod_info_tp
    !===============================================================================================
    subroutine Free_control_rod_info_tp (this)
        
        class(control_rod_info_tp), intent(in out)  :: this
        
        if (allocated(this%flags))          deallocate(this%flags)
        if (allocated(this%Umat))           deallocate(this%Umat)
        if (allocated(this%Lmat))           deallocate(this%Lmat)
        if (allocated(this%Ulen))           deallocate(this%Ulen)
        if (allocated(this%Llen))           deallocate(this%Llen)
        if (allocated(this%height))         deallocate(this%height)
    
    end subroutine Free_control_rod_info_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_ControlRodBank (this, ns, geom)
        
        class(ControlRodBank), intent(in out)  :: this
        type(SteadyState), intent(in)          :: ns
        type(Geometry), intent(in)             :: geom
        
        real(KREAL)  :: h_active
        integer  :: ia, i_rod
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        this%state%na = ns%state%layer
        h_active = 0.0D0
        do ia = this%state%n_bottom + 1, this%state%na - this%state%n_top
            h_active = h_active + geom%height(ia)
        end do
        
        ia = this%state%na
        i_rod = this%state%n_rod
        allocate(this%xsecs(i_rod, ia), stat=i_allocate)
        allocate(this%params(i_rod, ia), stat=i_allocate)
        
        do i_rod = 1, SIZE(this%xsecs, dim=1)
            do ia = 1, SIZE(this%xsecs, dim=2)
                call this%xsecs(i_rod, ia)%alloc ()
            end do
        end do
    
        do i_rod = 1, SIZE(this%params, dim=1)
            do ia = 1, SIZE(this%params, dim=2)
                call this%params(i_rod, ia)%alloc ()
            end do
        end do
        
        allocate(this%rods(this%state%n_rod), stat=i_allocate)
        do i_rod = 1, SIZE(this%rods)
            call this%rods(i_rod)%alloc (ns)
        end do
        
        ! set bank step
        allocate(this%init_step(this%state%n_bank), stat=i_allocate)
        allocate(this%last_step(this%state%n_bank), stat=i_allocate)
        
    end subroutine Allocate_ControlRodBank
    
    !$
    !===============================================================================================
    ! finalizer for class of ControlRodBank
    !===============================================================================================
    subroutine Free_ControlRodBank (this)
        
        class(ControlRodBank), intent(in out)  :: this
        integer  :: i, j
                
        if (allocated(this%xsecs))  then
            do i = 1, SIZE(this%xsecs, dim=1)
                do j = 1, SIZE(this%xsecs, dim=2)
                    call this%xsecs(i, j)%clean ()
                end do
            end do
            ! free it self
            deallocate(this%xsecs)
        end if
        
        if (allocated(this%params))   then
            do i = 1, SIZE(this%params, dim=1)
                do j = 1, SIZE(this%params, dim=2) 
                    call this%params(i, j)%clean ()
                end do
            end do
            ! free itself
            deallocate(this%params)
        end if
        
        if (allocated(this%rods))  then
            do i = 1, SIZE(this%rods)
                call this%rods(i)%clean ()
            end do
            ! free it self
            deallocate(this%rods)
        end if
        
        ! free bank
        if (allocated(this%init_step))      deallocate(this%init_step)
        if (allocated(this%last_step))      deallocate(this%last_step)
    
    end subroutine Free_ControlRodBank
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_ControlRodBank_xsec (this, ns, unit_)
        
        class(ControlRodBank), intent(in out)  :: this
        type(SteadyState), intent(in)          :: ns
        integer, intent(in)                    :: unit_
        
        integer  :: ia, i_rod, ig, ic, ip
        
        do i_rod = 1, SIZE(this%xsecs, dim=1)
            do ia = 1, SIZE(this%xsecs, dim=2)
                write(unit=unit_, fmt=*) TRIM(CHAR_SSUBMARK)
                write(unit=unit_, fmt=*) 'cr rod index is    :', i_rod
                write(unit=unit_, fmt=*) 'cr axial index is  :', ia
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%xsecs(i_rod,ia)%chi_steady(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%xsecs(i_rod,ia)%sigma_t(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%xsecs(i_rod,ia)%sigma_f_nu(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%xsecs(i_rod,ia)%sigma_f_kappa(:)
                do ic = 1, (ns%state%scat_order + 1)
                    do ig = 1, ns%state%ng
                        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%xsecs(i_rod,ia)%sigma_s(ig, :, ic)
                    end do
                end do
            end do
        end do
        
        do i_rod = 1, SIZE(this%params, dim=1)
            do ia = 1, SIZE(this%params, dim=2)
                write(unit=unit_, fmt=*) TRIM(CHAR_SSUBMARK)
                write(unit=unit_, fmt=*) 'cr rod index is    :', i_rod
                write(unit=unit_, fmt=*) 'cr axial index is  :', ia
                do ip = 1, SIZE(this%params(i_rod,ia)%chi_delay, dim=1)
                    write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%params(i_rod,ia)%chi_delay(ip, :)
                end do
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%params(i_rod,ia)%velocity(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%params(i_rod,ia)%beta(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%params(i_rod,ia)%lambda(:)
            end do
        end do
    
    end subroutine Print_ControlRodBank_xsec
    
    !$
    !===============================================================================================
    ! print cr rod material configuration in nodal
    !===============================================================================================
    subroutine Print_ControlRodBank (this, geom, unit_)
        
        class(ControlRodBank), intent(in out)  :: this
        type(Geometry), intent(in)             :: geom
        integer, intent(in)                    :: unit_
        
        real(KREAL)  :: dis, hmove
        integer      :: i_rod, i_bank, ia
        
        ! nodal info
        do i_rod = 1, MIN(10, this%state%n_rod)
            associate (rod => this%rods(i_rod))
            i_bank = rod%bank
            hmove = (this%last_step(i_bank) - this%init_step(i_bank)) * this%state%step_size
            
            write(unit=unit_, fmt="(/, 2X, A)")  TRIM(CHAR_SSUBMARK) 
            write(unit=unit_, fmt="(   2X, 'Control rod matrial in nodal (Top --> Bottom): < ', I2, '@', I2, ' >')") i_rod, i_bank
            write(unit=unit_, fmt="(   2X, 'Move info: ', F8.2, ' --> ', F8.2, ':   ', F8.4, ' cm')") this%init_step(i_bank), this%last_step(i_bank), hmove
            write(unit=unit_, fmt="(2X, 'Axial ID:', TR4, 'Height(cm):', TR4, 'To bottom(cm):', TR4, 'Mat flags:', TR4, 'len-PART(cm)', TR4, 'frac-PART(%)', TR4, 'mat-PART')")
            write(unit=unit_, fmt="(2X, '----------------------------------------------------------------------------------------------------')")
            
            dis = SUM(geom%height)
            do ia = this%state%na, 1, REVERSE_ORDER
                dis = dis - geom%height(ia) / 2.0D0
                write(unit=unit_, fmt="(2X, I4, TR11, F8.4, TR8, F8.4, TR7, A10)", advance="no")  ia, geom%height(ia), dis, rod%flags(ia)
                write(unit=unit_, fmt="(TR1,      'U:', TR2, F8.4, TR7, F8.4, TR7, I4)") rod%Ulen(ia), rod%Ulen(ia)/rod%height(ia), rod%Umat(ia)
                write(unit=unit_, fmt="(2X, TR57, 'L:', TR2, F8.4, TR7, F8.4, TR7, I4)") rod%Llen(ia), rod%Llen(ia)/rod%height(ia), rod%Lmat(ia)
                dis = dis - geom%height(ia) / 2.0D0
            end do
            write(unit=unit_, fmt="(2X, '----------------------------------------------------------------------------------------------------')")
            write(unit=unit_, fmt="(/)")
            
            end associate
        end do
        
    end subroutine Print_ControlRodBank
    
    !$
    !===============================================================================================
    ! adjust value of <nu>
    !===============================================================================================
    subroutine Adjust_ControlRodBank_xsec (this, k_eff)
        
        class(ControlRodBank), intent(in out)  :: this
        real(KREAL), intent(in)            :: k_eff
        
        integer  :: ia, i_rod
        
        do ia = 1, SIZE(this%xsecs, dim=2)
            do i_rod = 1, SIZE(this%xsecs, dim=1)
                this%xsecs(i_rod,ia)%nu = this%xsecs(i_rod,ia)%nu / k_eff
                this%xsecs(i_rod,ia)%sigma_f_nu = this%xsecs(i_rod,ia)%sigma_f_nu / k_eff
            end do
        end do
    
    end subroutine Adjust_ControlRodBank_xsec
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! set lower following material for CR
    !===============================================================================================
    subroutine Set_ControlRodBank_lower_mat (this, lower_mat)
        
        class(ControlRodBank), intent(in out)  :: this
        integer, intent(in)  :: lower_mat(:)
        
        integer  :: i
        
        do i = 1, this%state%n_rod
            this%rods(i)%lower_mat = lower_mat(i)
        end do
    
    end subroutine Set_ControlRodBank_lower_mat
    
    !$
    !===============================================================================================
    ! set upper following material for CR
    !===============================================================================================
    subroutine Set_ControlRodBank_upper_mat (this, upper_mat)
        
        class(ControlRodBank), intent(in out)  :: this
        integer, intent(in)  :: upper_mat(:)
        
        integer  :: i
        
        do i = 1, this%state%n_rod
            this%rods(i)%upper_mat = upper_mat(i)
        end do
    
    end subroutine Set_ControlRodBank_upper_mat
    
    !$
    !===============================================================================================
    ! set black material height
    !===============================================================================================
    subroutine Set_ControlRodBank_hblack (this, hblack)
        
        class(ControlRodBank), intent(in out)  :: this
        real(KREAL), intent(in)  :: hblack(:)
        
        integer  :: i
        
        do i = 1, this%state%n_rod
            this%rods(i)%hblack = hblack(i)
        end do
    
    end subroutine Set_ControlRodBank_hblack

    !$
    !===============================================================================================
    ! set initial step information for control rod bank
    !===============================================================================================
    subroutine Set_ControlRodBank_init_position (this, init_step)
        
        class(ControlRodBank), intent(in out)  :: this
        real(KREAL), intent(in)  :: init_step(:)
        
        integer  :: i
        
        do i = 1, SIZE(init_step)
            this%init_step(i) = init_step(i)
            this%last_step(i) = init_step(i)
        end do
    
    end subroutine Set_ControlRodBank_init_position
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! set half insert information
    !===============================================================================================
    subroutine Set_ControlRodBank_black (this, mat_info)
        
        class(ControlRodBank), intent(in out)  :: this
        type(Material), intent(in)         :: mat_info
        
        real(KREAL)  :: blackbottom, blacktop
        real(KREAL)  :: hbottom, htop
        integer  :: Lmat, Umat
        integer  :: i_rod, ia, im
        
        ! update BLACK flag, minus value means no BLACK
        do i_rod = 1, this%state%n_rod
            this%rods(i_rod)%black_beg = -HUGE(0.0)
            this%rods(i_rod)%black_end = -HUGE(0.0)
            hbottom = REAL_ZERO
            htop = SUM(this%rods(i_rod)%height)
            
            ! upstairs
            upstairs: do ia = 1, this%state%na
                im = this%rods(i_rod)%Lmat(ia)
                if (mat_info%libs(im)%is_CR)  then
                    this%rods(i_rod)%black_beg = hbottom
                    exit upstairs
                end if
                hbottom = hbottom + this%rods(i_rod)%Llen(ia)
                
                im = this%rods(i_rod)%Umat(ia)
                if (mat_info%libs(im)%is_CR)  then
                    this%rods(i_rod)%black_beg = hbottom
                    exit upstairs
                end if
                hbottom = hbottom + this%rods(i_rod)%Ulen(ia)
            end do upstairs
            
            ! downstairs
            downstairs: do ia = this%state%na, 1, REVERSE_ORDER
                im = this%rods(i_rod)%Umat(ia)
                if (mat_info%libs(im)%is_CR)  then
                    this%rods(i_rod)%black_end = htop
                    exit downstairs
                end if
                htop = htop - this%rods(i_rod)%Ulen(ia)
                
                im = this%rods(i_rod)%Lmat(ia)
                if (mat_info%libs(im)%is_CR)  then
                    this%rods(i_rod)%black_end = htop
                    exit downstairs
                end if
                htop = htop - this%rods(i_rod)%Llen(ia)
            end do downstairs
        end do
        
        ! absorption mat height
        do i_rod = 1, this%state%n_rod
            blacktop = 0.0; blackbottom = 0.0
            if (this%rods(i_rod)%black_end > 0)  then
                blacktop = this%rods(i_rod)%black_end
            end if
            if (this%rods(i_rod)%black_beg > 0)  then
                blackbottom = this%rods(i_rod)%black_beg
            end if
            
            if ((this%rods(i_rod)%hblack < 0.0) .and. (blacktop-blackbottom > 0.0))  then
                this%rods(i_rod)%hblack = blacktop - blackbottom
            end if
        end do
        
        ! set fuel flag
        do i_rod = 1, this%state%n_rod
            this%rods(i_rod)%is_gray = .FALSE.
            
            do ia = this%state%n_bottom+1, this%state%na-this%state%n_top
                Lmat = this%rods(i_rod)%Lmat(ia)
                Umat = this%rods(i_rod)%Umat(ia)
                
                if (ALL([mat_info%libs(Lmat)%is_CR, mat_info%libs(Umat)%is_CR]))  then
                    this%rods(i_rod)%flags(ia) = this%symbol%BLACK
                else if (ANY([mat_info%libs(Lmat)%is_CR, mat_info%libs(Umat)%is_CR]))  then
                    this%rods(i_rod)%flags(ia) = this%symbol%GRAY
                    this%rods(i_rod)%is_gray = .TRUE.
                else
                    this%rods(i_rod)%flags(ia) = this%symbol%WHITE
                end if
            end do
        end do
        
        ! set reflector flag
        if (this%state%n_bottom /= 0)  then
            do i_rod = 1, this%state%n_rod
                do ia = 1, this%state%n_bottom
                    this%rods(i_rod)%flags(ia) = this%symbol%REFLECTOR
                end do
            end do
        end if
        
        if (this%state%n_top /= 0)  then
            do i_rod = 1, this%state%n_rod
                do ia = this%state%na-this%state%n_top+1, this%state%na
                    this%rods(i_rod)%flags(ia) = this%symbol%REFLECTOR
                end do
            end do
        end if
        
    end subroutine Set_ControlRodBank_black
    
    !$
    !===============================================================================================
    ! set control rod configuration & mesh cross section
    !===============================================================================================
    subroutine Set_ControlRodBank (this, geom, mat_info, xsec, param, cr_conf)
        
        class(ControlRodBank), intent(in out)   :: this
        type(Geometry), intent(in)              :: geom 
        type(Material), intent(in)              :: mat_info
        type(CrossSection), intent(in)          :: xsec
        type(KineticsParameter), intent(in)     :: param
        integer, intent(in)                     :: cr_conf(:)
        
        integer          :: i_rod, ia
        integer          :: iz
        integer          :: im
        
        ! get control rod material information
        i_rod = 0
        do iz = 1, SIZE(cr_conf)
            if (cr_conf(iz) /= 0)  then                                         ! none-zero means it is control rod
                i_rod = i_rod + 1
                this%rods(i_rod)%bank = cr_conf(iz)
                this%rods(i_rod)%zone = iz
                this%rods(i_rod)%height = geom%height                           ! NOTE:---
            end if
        end do
        
        ! get initial fine mesh information
        this%state%h_lower = SUM(geom%height( :this%state%n_bottom))
        this%state%h_upper = SUM(geom%height(this%state%na-this%state%n_top+1: ))
        this%state%h_active = SUM(geom%height) - this%state%h_lower - this%state%h_upper
            
    end subroutine Set_ControlRodBank
    
    !$
    !===============================================================================================
    ! set xsec & kinetics parameters to control rod zone, initial state
    !===============================================================================================
    subroutine Set_ControlRodBank_xsec (this, geom, mat_info, xsec, param)
        
        class(ControlRodBank), intent(in out)   :: this
        type(Geometry), intent(in)              :: geom 
        type(Material), intent(in)              :: mat_info
        type(CrossSection), intent(in)          :: xsec
        type(KineticsParameter), intent(in)     :: param
        
        integer          :: i_rod, ia
        integer          :: iz, im
        
        do i_rod = 1, this%state%n_rod
            iz = this%rods(i_rod)%zone
            do ia = 1, this%state%na
                im = mat_info%loading(iz, ia)
                
                this%rods(i_rod)%Umat(ia) = im
                this%rods(i_rod)%Lmat(ia) = im
                this%rods(i_rod)%Ulen(ia) = geom%height(ia)
                this%rods(i_rod)%Llen(ia) = REAL_ZERO
                this%rods(i_rod)%height(ia) = geom%height(ia)
            end do
        end do
            
        call this%black (mat_info)
    
    end subroutine Set_ControlRodBank_xsec
    
    !$
    !===============================================================================================
    ! update mat info by new step
    !===============================================================================================
    subroutine Move_ControlRodBank (this, mat_info, step)
        
        class(ControlRodBank), intent(in out)  :: this
        type(Material), intent(in)             :: mat_info
        real(KREAL), intent(in), optional      :: step(:)
        
        real(KREAL)  :: new_step(SIZE(step))
        real(KREAL)  :: hmove
        real(KREAL)  :: hl, hu, dis
        real(KREAL)  :: Uface, Lface, hCR 
        integer  :: i_bank                                                      ! count for control rod bank
        integer  :: i_rod                                                       ! count for control rod number
        integer  :: ia, iz, im, im_add, im_sub 
        
        ! if not exists, move to initial conf 
        if (PRESENT(step))  then
            new_step = step 
        else
            new_step = this%last_step 
        end if 
        
        do i_bank = 1, this%state%n_bank
            new_step(i_bank) = step(i_bank)
            
            ! check current step, warning
            if (this%state%min_step < 0.0)  then
            else
                if (new_step(i_bank) < this%state%min_step)  then
                end if
            end if
            if (this%state%max_step < 0.0)  then
            else 
                if (new_step(i_bank) > this%state%max_step)  then
                end if
            end if
        end do
        
        ! update mat info 
        do i_rod = 1, this%state%n_rod 
            iz = this%rods(i_rod)%zone
            i_bank = this%rods(i_rod)%bank
            hCR = new_step(i_bank)*this%state%step_size + this%state%cr_gap + this%state%h_lower
            
            do ia = 1, this%state%na
                im = mat_info%loading(iz, ia)
                im_add = mat_info%libs(im)%ID_add
                im_sub = mat_info%libs(im)%ID_sub
                Uface = SUM(this%rods(i_rod)%height( :ia))
                Lface = SUM(this%rods(i_rod)%height( :ia-1))
                
                if (hCR <= Lface)  then
                    this%rods(i_rod)%Umat(ia) = im_add
                    this%rods(i_rod)%Lmat(ia) = im_add
                    this%rods(i_rod)%Ulen(ia) = Uface - Lface
                    this%rods(i_rod)%Llen(ia) = 0.0 
                else if ((Lface < hCR) .AND. (hCR < Uface))  then
                    this%rods(i_rod)%Umat(ia) = im_add
                    this%rods(i_rod)%Lmat(ia) = im_sub
                    this%rods(i_rod)%Ulen(ia) = Uface - hCR 
                    this%rods(i_rod)%Llen(ia) = hCR - Lface 
                else if (Uface <= hCR)  then
                    this%rods(i_rod)%Umat(ia) = im_sub
                    this%rods(i_rod)%Lmat(ia) = im_sub
                    this%rods(i_rod)%Ulen(ia) = 0.0 
                    this%rods(i_rod)%Llen(ia) = Uface - Lface 
                end if
            end do 
        end do 
        
        call this%black (mat_info)
        
        ! update step position
        this%last_step = new_step
        
    end subroutine Move_ControlRodBank
    
    !$
    !===============================================================================================
    ! get nodal cross section for CR from fine mesh 
    !===============================================================================================
    subroutine Homogeneous_ControlRodBank (this, xsec_inp, param_inp, geom)
        
        class(ControlRodBank), intent(in out)   :: this
        type(CrossSectionInput), intent(in)     :: xsec_inp                         
        type(KineticsParameterInput), intent(in):: param_inp                        
        type(Geometry), intent(in)              :: geom
        
        integer  :: i_rod, iz
        integer  :: ia                                                          ! index for axial nodal
        type(cross_section_info_tp)       :: Lxsec, Uxsec
        type(kinetics_parameter_info_tp)  :: Lparam, Uparam
        real(KREAL)                       :: Llen, Ulen
        real(KREAL)                       :: Lwt, Uwt
        integer                           :: Lmat, Umat
        
        Lwt = 1.0; Uwt = 1.0;
        
        ! update nodal cross section
        do i_rod = 1, this%state%n_rod
            iz = this%rods(i_rod)%zone
            
            do ia = this%state%n_bottom+1, this%state%na-this%state%n_top
                Lmat = this%rods(i_rod)%Lmat(ia)
                Umat = this%rods(i_rod)%Umat(ia)
                Llen = this%rods(i_rod)%Llen(ia)
                Ulen = this%rods(i_rod)%Ulen(ia)
                Lxsec = xsec_inp%mats(Lmat)
                Uxsec = xsec_inp%mats(Umat)
                Lparam = param_inp%mats(Lmat)
                Uparam = param_inp%mats(Umat)
                
                call this%xsecs(i_rod, ia)%weight ([Lxsec, Uxsec], [Llen, Ulen], [Lwt, Uwt])
                call this%params(i_rod, ia)%weight ([Lparam, Uparam], [Llen, Ulen], [Lwt, Uwt])
            end do
        end do
    
    end subroutine Homogeneous_ControlRodBank
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Map_CR2Core (this, xsec, param)
        
        class(ControlRodBank), intent(in out)   :: this
        type(CrossSection), intent(in out)      :: xsec                         ! NOTE: also intent out
        type(KineticsParameter), intent(in out) :: param                        ! NOTE: also intent out
        
        integer  :: i_rod, iz
        integer  :: ia
        
        do i_rod = 1, this%state%n_rod
            iz = this%rods(i_rod)%zone
            do ia = this%state%n_bottom+1, this%state%na-this%state%n_top
                xsec%matrixs(iz, ia) = this%xsecs(i_rod, ia)
                param%matrixs(iz, ia) = this%params(i_rod, ia)
            end do
        end do
    
    end subroutine Map_CR2Core
    
    !$
    !===============================================================================================
    ! assignment implement for ControlRodBank
    !===============================================================================================
    subroutine Equal_ControlRodBank (left, right)
        
        class(ControlRodBank), intent(in out)  :: left
        type(ControlRodBank), intent(in)       :: right
        
        integer  :: ni, nj
        integer  :: i, j
        integer  :: i_allocate
        
        call left%clean ()
        
        left%state%na = right%state%na
        left%state%n_top = right%state%n_top
        left%state%n_bottom = right%state%n_bottom
        left%state%n_rod = right%state%n_rod
        left%state%n_bank = right%state%n_bank
        left%state%max_step = right%state%max_step
        left%state%min_step = right%state%min_step
        left%state%step_size = right%state%step_size
        left%state%h_upper = right%state%h_upper
        left%state%h_active = right%state%h_active
        left%state%h_lower = right%state%h_lower
        
        if (allocated(right%rods))  then
            ni = SIZE(right%rods)
            allocate(left%rods(ni), stat=i_allocate)
            do i = 1, ni
                left%rods(i)%bank = right%rods(i)%bank
                left%rods(i)%zone = right%rods(i)%zone
                left%rods(i)%upper_mat = right%rods(i)%upper_mat
                left%rods(i)%lower_mat = right%rods(i)%lower_mat
                left%rods(i)%is_gray = right%rods(i)%is_gray
                left%rods(i)%black_beg = right%rods(i)%black_beg
                left%rods(i)%black_end = right%rods(i)%black_end
                left%rods(i)%hblack = right%rods(i)%hblack
                
                if (allocated(right%rods(i)%flags))  then
                    nj = SIZE(right%rods(i)%flags)
                    allocate(left%rods(i)%flags(nj), stat=i_allocate)
                    allocate(left%rods(i)%Umat(nj), stat=i_allocate)
                    allocate(left%rods(i)%Lmat(nj), stat=i_allocate)
                    allocate(left%rods(i)%Ulen(nj), stat=i_allocate)
                    allocate(left%rods(i)%Llen(nj), stat=i_allocate)
                    allocate(left%rods(i)%height(nj), stat=i_allocate)
                    left%rods(i)%flags = right%rods(i)%flags
                    left%rods(i)%Umat = right%rods(i)%Umat
                    left%rods(i)%Lmat = right%rods(i)%Lmat
                    left%rods(i)%Ulen = right%rods(i)%Ulen
                    left%rods(i)%Llen = right%rods(i)%Llen
                    left%rods(i)%height = right%rods(i)%height
                end if
            end do
        end if
        
        if (allocated(right%xsecs))  then
            ni = SIZE(right%xsecs, dim=1)
            nj = SIZE(right%xsecs, dim=2)
            allocate(left%xsecs(ni, nj), stat=i_allocate)
            allocate(left%params(ni, nj), stat=i_allocate)
            do i = 1, ni
                do j = 1, nj
                    left%xsecs(i, j) = right%xsecs(i, j)
                    left%params(i, j) = right%params(i, j)
                end do
            end do
        end if
        
        if (allocated(right%init_step))  then
            ni = SIZE(right%init_step, dim=1)
            do i = 1, ni
                allocate(left%init_step(ni), stat=i_allocate)
                allocate(left%last_step(ni), stat=i_allocate)
                left%init_step = right%init_step
                left%last_step = right%init_step
            end do
        end if
    
    end subroutine Equal_ControlRodBank
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Search_ControlRodBank (this, k_eff, is_pass)
        
        class(ControlRodBank), intent(in out)  :: this
        real(KREAL), intent(in)      :: k_eff
        logical, intent(in out)      :: is_pass
        
        integer, save          :: cnt = 0
        logical, save          :: is_damp = .FALSE.
        real(KREAL), save      :: pos_old, keff_old
        real(KREAL)            :: pos_new, keff_new
        real(KREAL)            :: error
        real(KREAL)            :: pos_LL, keff_LL
        
        is_pass = .FALSE.
        if (ABS(k_eff - this%search_value) <= this%search_limit)  then
            is_pass = .TRUE.
            return
        end if
        
        associate (rod => this%init_step(this%bank_idx))
            if ((rod < this%state%min_step) .or. (this%state%max_step < rod))  then
                call a_warning%set (INFO_LIST_CRBANK, 'CR postion exceed')
                call a_warning%print (FILES%MAIN)
                call a_warning%print (OUTPUT_UNIT)
                is_pass = .TRUE.
                return
            end if 
            
            cnt = cnt + 1
            ! save old value 
            if (cnt == 1)  then
                pos_LL = rod
                keff_LL = k_eff
            else 
                pos_LL = pos_old
                keff_LL = keff_old
            end if 
            
            if (cnt == 1)  then
                error = (k_eff - this%search_value) 
                pos_old = rod
                keff_old = k_eff
                if (error < 0.0)  then
                    rod = NINT(0.5 * (rod + this%state%max_step))
                else
                    rod = NINT(0.5 * (rod + this%state%min_step))
                end if 
            else
                pos_new = rod
                keff_new = k_eff
                if ((cnt >= 5) .OR. ((ABS(k_eff-this%search_value)) <= 0.0001) .OR. ((ABS(keff_new - keff_old)) <= 0.0001))  then
                    is_damp = .TRUE.
                end if
                
                if (NINT(pos_new-pos_old) == 0)  then
                    call a_warning%set (INFO_LIST_CRBANK, 'Minial step is searched')
                    call a_warning%print (FILES%MAIN)
                    call a_warning%print (OUTPUT_UNIT)
                    is_pass = .TRUE.
                    return
                end if 
                
                error = NINT((this%search_value-keff_new) / (keff_new-keff_old) * (pos_new-pos_old))
                rod = pos_new + error
                pos_old = pos_new
                keff_old = keff_new
            end if
            
            write(FILES%MAIN, fmt="(1x, A)")  TRIM(CHAR_SUBMARK)
            write(FILES%MAIN, fmt="(1x, A)")  "Rod Search Results:"
            write(FILES%MAIN, fmt="(1x, *(A, TR3))")  "  No. ",  "LL-pos ", "keff     ", "L-pos ", "keff    ", "==>", "this pos       delta-pos" 
            write(FILES%MAIN, fmt="(1x, I4, TR2, 2(F9.1, TR3, F9.6), TR6, F9.1)")  cnt, pos_LL, keff_LL, pos_old, keff_old, rod
            
            if (.TRUE.)  then 
                write(OUTPUT_UNIT, fmt="(1x, A)")  TRIM(CHAR_SUBMARK)
                write(OUTPUT_UNIT, fmt="(1x, A)")  "Rod Search Results:"
                write(OUTPUT_UNIT, fmt="(1x, *(A, TR3))")  "  No. ",  "LL-pos ", "keff     ", "L-pos ", "keff    ", "==>", "this pos       delta-pos" 
                write(OUTPUT_UNIT, fmt="(1x, I4, TR2, 2(F9.1, TR3, F9.6), TR6, F9.1, TR5, F10.1)")  cnt, pos_LL, keff_LL, pos_old, keff_old, rod, (rod-pos_old)
            end if
            
        end associate
        
        this%last_step = this%init_step
        
    end subroutine Search_ControlRodBank
    
end module CRbank_header
