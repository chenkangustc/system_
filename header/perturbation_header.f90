!$
!===================================================================================================
!
!   class for transient perturbation from cross section, external source and control rod
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          XsecPerturbation
!                               SourcePerturbation
!                               ControlRodPerturbation
!                               MaterialPerturbation
!                               CasePerturbation
!
!===================================================================================================
module perturbation_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use files_dirs
    
    use material_header,                only : Material, CrossSection, cross_section_info_tp
    use CRbank_header,                  only : ControlRodBank
    use contain_header,                 only : TimeListParameter
    use input_xsec,                     only : xsec_read_known
    
    implicit none
    private
    public  :: XsecPerturbation, SourcePerturbation, ControlRodPerturbation, MaterialPerturbation, &
        &   CasePerturbation, THPerturbation
    
    ! --------------------------------------------------------------------------
    ! type for cross section perturbation per line
    type, private  :: xsec_perturbation_tp
        character(len=MAX_WORD_LEN), public     :: xs_type       = 'SIGMA_S'    ! which kind of cross section to be perturbed
        integer, public                         :: matID         = 0            ! material ID for perturbation
        integer, public                         :: ng_start      = 1            ! starting energy group with perturbation 
        integer, public                         :: ng_end        = 1            ! end energy group with perturbation
        real(KREAL), public                     :: time_start    = 0.0          ! perturbation start time point
        real(KREAL), public                     :: time_end      = 0.0          ! perturbation end time point
        real(KREAL), public                     :: percent       = 0.0          ! perturbation quantity (%)
        character(len=MAX_WORD_LEN), public     :: type          = 'RAMP'       ! how to perturb
        integer, public                         :: to_ng_start   = 1            ! the following three line use for scatter xs
        integer, public                         :: to_ng_end     = 1
        integer, public                         :: scat_index    = 0
    end type xsec_perturbation_tp
    
    ! type for external source perturbation per line
    type, private  :: source_perturbation_tp
        integer, public                         :: kindID        = 0            ! source kind ID for perturbation
        integer, public                         :: ng_start      = 1            ! starting energy group with perturbation
        integer, public                         :: ng_end        = 1            ! end energy group with perturbation
        real(KREAL), public                     :: time_start    = 0.0          ! perturbation start time point
        real(KREAL), public                     :: time_end      = 0.0          ! perturbation end time point
        real(KREAL), public                     :: percent       = 0.0          ! perturbation quantity (%)
        character(len=MAX_WORD_LEN), public     :: type          = 'RAMP'       ! how to perturb
    end type source_perturbation_tp
    
    ! type for control rod perturbation per line
    type, private  :: control_rod_move_info_tp
        integer, public                         :: bankID        = 0            ! control rod bank index
        real(KREAL), public                     :: step_start    = 0.0          ! position step when starting move
        real(KREAL), public                     :: step_end      = 0.0          ! position step when ending move
        real(KREAL), public                     :: time_start    = 0.0          ! end energy group with perturbation
        real(KREAL), public                     :: time_end      = 0.0          ! perturbation start time point
        character(len=MAX_WORD_LEN), public     :: type          = 'RAMP'       ! how to perturb
    end type control_rod_move_info_tp
    
    ! type for control rod trip information 
    type, private  :: control_rod_trip_info_tp
        real(KREAL), public                     :: pstart       = 0.35          ! power level active CR trip 
        real(KREAL), public                     :: tstart       = 0.0           ! time when power approach 
        real(KREAL), public                     :: tdelay       = 0.6           ! time delay since power approach
        logical, public                         :: is_active    = .FALSE.       ! is trip active ? 
    end type control_rod_trip_info_tp
    
    ! type for material perburbation per line
    type, private  :: mat_perturbation_info_tp
        integer, public                         :: mat_beg      = 1             ! initial material id
        real(KREAL), public                     :: per_beg      = 100.0         ! initial material percentage
        integer, public                         :: mat_end      = 2             ! final material id
        real(KREAL), public                     :: per_end      = 50.0          ! final material percentage
        real(KREAL), public                     :: time_start   = 0.0           ! perburbation start time point
        real(KREAL), public                     :: time_end     = 0.0           ! perburbation end time point
        character(len=MAX_WORD_LEN), public     :: type         = 'RAMP'        ! how to perturb
    end type mat_perturbation_info_tp
    
    ! type for case perturbation per one
    type, private  :: case_xsec_info_tp
        integer, public                         :: ID                           ! material ID
        character(len=MAX_WORD_LEN), public     :: name                         ! material filename
        integer, public                         :: case_ID                      ! case matrix ID in material file
        integer, public                         :: burn_ID                      ! burn ID in material file
    end type case_xsec_info_tp
    
    type, private  :: case_perturbation_tp
        integer, public                               :: caseID                 ! ID for perturbation index
        logical, public                               :: is_new_configure
        integer, public, allocatable                  :: loading(:, :)
        type(case_xsec_info_tp), public, allocatable  :: mats(:)
    contains
        procedure, public  :: alloc => Allocate_case_perturbation_tp
        procedure, public  :: clean => Free_case_perturbation_tp
    end type case_perturbation_tp
    
    ! type for thermal-hydraulic per one 
    type, private  :: flow_perturbation_info_tp
        integer, public                         :: channelID            ! TH channel index, "0" for all
        integer, public                         :: type                 ! function type
        real(KREAL)                             :: natural      
        real(KREAL)                             :: time_start
        real(KREAL), allocatable                :: variables(:)         ! parameters to define perturbation function
    end type flow_perturbation_info_tp
    type, private  :: Tm_perturbation_info_tp
        integer, public                         :: type                 ! function type  
        real(KREAL)                             :: time_start
        real(KREAL)                             :: time_end 
        real(KREAL), allocatable                :: variables(:)         ! parameters to define perturbation function 
    end type Tm_perturbation_info_tp
    
    type  THPerturbation
        integer, public                                       :: n_flow     = 0
        integer, public                                       :: var_flow   = 0
        integer, public                                       :: n_Tm       = 0
        integer, public                                       :: var_Tm     = 0 
        type(flow_perturbation_info_tp), public, allocatable  :: flow_perts(:) 
        type(Tm_perturbation_info_tp), public, allocatable    :: Tm_perts(:) 
    contains
        procedure, public  :: alloc => Allocate_THPerturbation 
        procedure, public  :: clean => Free_THPerturbation
        procedure, public  :: print => Print_THPerturbation
    end type THPerturbation
    
    ! --------------------------------------------------------------------------
    ! type for cross section perturbation total
    type  XsecPerturbation
        integer, public                                     :: n_pert  = 0      ! number of perturbation
        type(xsec_perturbation_tp), public, allocatable     :: perts(:)         ! hold for perturbation
    contains
        procedure, public  :: alloc => Allocate_XsecPerturbation
        procedure, public  :: clean => Free_XsecPerturbation
    end type XsecPerturbation
    
    ! type for external source perturbation total
    type  SourcePerturbation
        integer, public                                     :: n_pert  = 0      ! number of perturbation
        type(source_perturbation_tp), public, allocatable   :: perts(:)         ! hold for perturbation
    contains
        procedure, public  :: alloc => Allocate_SourcePerturbation
        procedure, public  :: clean =>  Free_SourcePerturbation
    end type SourcePerturbation
    
    ! type for control rod perturbation total
    type  ControlRodPerturbation
        integer, public                                      :: n_pert  = 0      ! number of control rod moving
        type(control_rod_move_info_tp), public, allocatable  :: perts(:)
        logical, public                                      :: is_trip = .FALSE.
        type(control_rod_trip_info_tp), public               :: trip             ! trip infomation 
    contains
        procedure, public  :: alloc => Allocate_ControlRodPerturbation
        procedure, public  :: clean => Free_ControlRodPerturbation
        procedure, public  :: init => Init_ControlRodPerturbation
        procedure, public  :: set_trip => Set_CorePerturbation_trip 
    end type ControlRodPerturbation
    
    ! type for material perturbation total
    type  MaterialPerturbation
        integer, public                                      :: n_pert  = 0      
        type(mat_perturbation_info_tp), public, allocatable  :: perts(:)
    contains
        procedure, public  :: alloc => Allocate_MaterialPerturbation
        procedure, public  :: clean => Free_MaterialPerturbation
    end type MaterialPerturbation
    
    ! type for state case pertubation
    type  CasePerturbation
        integer, public                                      :: n_pert  = 0      ! number of perturbation
        type(case_perturbation_tp), public, allocatable      :: perts(:)
    contains
        procedure, public  :: alloc => Allocate_CasePerturbation
        procedure, public  :: clean => Free_CasePerturbation
        procedure, public  :: set_conf => Set_CasePerturbation_configure
        procedure, public  :: fix_conf => Fix_CasePerturbation_configure
        procedure, public  :: set_xsec => Set_CasePerturbation_xsec
        procedure, public  :: get_xsec => Get_CasePerturbation_xsec
    end type CasePerturbation
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_XsecPerturbation (this)
        
        class(XsecPerturbation), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%perts(this%n_pert), stat=i_allocate)
        
        ! do not explicit initialization, use default value
        
    end subroutine Allocate_XsecPerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of XsecPerturbation
    !===============================================================================================
    subroutine Free_XsecPerturbation (this)
        
        class(XsecPerturbation), intent(in out)  :: this
        
        if (allocated(this%perts))           deallocate(this%perts)
    
    end subroutine Free_XsecPerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_SourcePerturbation (this)
        
        class(SourcePerturbation), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%perts(this%n_pert), stat=i_allocate)
        
        ! do noe explicit initializatioin, use default value
    
    end subroutine Allocate_SourcePerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of SourcePerturbation
    !===============================================================================================
    subroutine Free_SourcePerturbation (this)
        
        class(SourcePerturbation), intent(in out)  :: this
        
        if (allocated(this%perts))           deallocate(this%perts)
    
    end subroutine Free_SourcePerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_ControlRodPerturbation (this)
        
        class(ControlRodPerturbation), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%perts(this%n_pert), stat=i_allocate)
    
    end subroutine Allocate_ControlRodPerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of ControlRodPerturbation
    !===============================================================================================
    subroutine Free_ControlRodPerturbation (this)
        
        class(ControlRodPerturbation), intent(in out)  :: this
        
        if (allocated(this%perts))              deallocate(this%perts)
    
    end subroutine Free_ControlRodPerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Init_ControlRodPerturbation (this, cr_bank)
        
        class(ControlRodPerturbation), intent(in out)  :: this
        type(ControlRodBank), intent(in)               :: cr_bank 
        
        integer  :: ipert, ibank
        
        do ipert = 1, SIZE(this%perts)
            ibank = this%perts(ipert)%bankID
            this%perts(ibank)%step_start = cr_bank%init_step(ibank)
        end do 
    
    end subroutine Init_ControlRodPerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CorePerturbation_trip (this, timelist, tin)
        
        class(ControlRodPerturbation), intent(in out)  :: this
        type(TimeListParameter), intent(in)       :: timelist
        real(KREAL), intent(in)                   :: tin
        
        real(KREAL)  :: power_frac
        
        power_frac = timelist%power / ns%flag%rated_power * 100.0D0
        
        ! get trip by power in percentage 
        if ((.NOT. this%trip%is_active) .AND. (power_frac >= this%trip%pstart * 100.0D0))  then
            this%trip%is_active = .TRUE.
            this%trip%tstart = tin 
        end if 
    
    end subroutine Set_CorePerturbation_trip
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_MaterialPerturbation (this)
        
        class(MaterialPerturbation), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%perts(this%n_pert), stat=i_allocate)
        
    end subroutine Allocate_MaterialPerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of MaterialPerturbation
    !===============================================================================================
    subroutine Free_MaterialPerturbation (this)
        
        class(MaterialPerturbation), intent(in out)  :: this
        
        if (allocated(this%perts))          deallocate(this%perts)
    
    end subroutine Free_MaterialPerturbation
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_THPerturbation (this)
        
        class(THPerturbation), intent(in out)  :: this
        integer :: i, i_allocate
        
        allocate(this%flow_perts(this%n_flow), stat=i_allocate)
        do i = 1, SIZE(this%flow_perts)
            allocate(this%flow_perts(i)%variables(this%var_flow), stat=i_allocate)
        end do
        
        allocate(this%Tm_perts(this%n_Tm), stat=i_allocate)
        do i = 1, SIZE(this%Tm_perts)
            allocate(this%Tm_perts(i)%variables(this%var_Tm), stat=i_allocate)
        end do 
        
    end subroutine Allocate_THPerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of THPerturbation
    !===============================================================================================
    subroutine Free_THPerturbation (this)
        
        class(THPerturbation), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%flow_perts))  then
            do i = 1, SIZE(this%flow_perts)
                if (allocated(this%flow_perts(i)%variables))  then
                    deallocate(this%flow_perts(i)%variables)
                end if 
            end do 
            deallocate(this%flow_perts)
        end if 
    
        if (allocated(this%Tm_perts))  then
            do i = 1, SIZE(this%Tm_perts)
                if (allocated(this%Tm_perts(i)%variables))  then
                    deallocate(this%Tm_perts(i)%variables)
                end if 
            end do 
            deallocate(this%Tm_perts)
        end if 
    
    end subroutine Free_THPerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_THPerturbation (this, unit_)
        
        class(THPerturbation), intent(in)  :: this
        integer, intent(in)                :: unit_ 
        
        integer :: i_pert
        
        write(unit=unit_, fmt='(1x, A, I4)') 'n_flow    = ', this%n_flow
        write(unit=unit_, fmt='(1x, A, I4)') 'var_flow  = ', this%var_flow
        write(unit=unit_, fmt='(1x, A, I4)') 'n_Tm      = ', this%n_Tm
        write(unit=unit_, fmt='(1x, A, I4)') 'var_Tm    = ', this%var_Tm
        
        if (this%n_flow > 0)  then
            do i_pert = 1, this%n_flow
            associate(pth => this%flow_perts(i_pert))
                write(unit=unit_, fmt='(1x, A, I4)') '------------FLOW : ', i_pert
                write(unit=unit_, fmt='(1x, A, I4)') 'pth%channelID  = ', pth%channelID
                write(unit=unit_, fmt='(1x, A, I4)') 'pth%type       = ', pth%type
                write(unit=unit_, fmt='(1x, A, ES12.5)') 'pth%natural    = ', pth%natural
                write(unit=unit_, fmt='(1x, A, ES12.5)') 'pth%time_start = ', pth%time_start
                write(unit=unit_, fmt='(1x, A, ES12.5)') 'pth%variables  = ', pth%variables
            end associate 
            end do 
        end if 
        
        if (this%n_Tm > 0)  then
            do i_pert = 1, this%n_Tm
            associate(pth => this%Tm_perts(i_pert))
                write(unit=unit_, fmt='(1x, A, I4)') '--------------TM : ', i_pert
                write(unit=unit_, fmt='(1x, A, I4)') 'pth%type       = ', pth%type
                write(unit=unit_, fmt='(1x, A, ES12.5)') 'pth%time_start = ', pth%time_start
                write(unit=unit_, fmt='(1x, A, ES12.5)') 'pth%variables  = ', pth%variables
            end associate 
            end do 
        end if 
        
    end subroutine Print_THPerturbation
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_case_perturbation_tp (this)
        
        class(case_perturbation_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check status first
        call this%clean ()
        
        allocate(this%loading(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%mats(ns%state%mat), stat=i_allocate)
        
        this%is_new_configure = .FALSE.
        this%loading = 1
    
    end subroutine Allocate_case_perturbation_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of case_perturbation_tp
    !===============================================================================================
    subroutine Free_case_perturbation_tp (this)
    
        class(case_perturbation_tp), intent(in out)  :: this
        
        if (allocated(this%loading))            deallocate(this%loading)       
        if (allocated(this%mats))               deallocate(this%mats)
        
    end subroutine Free_case_perturbation_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_CasePerturbation (this)
        
        class(CasePerturbation), intent(in out)  :: this
        integer  :: i_allocate
        integer  :: i
        
        ! check status first
        call this%clean ()
        
        allocate(this%perts(this%n_pert), stat=i_allocate)
        
        do i = 1, this%n_pert
            call this%perts(i)%alloc ()
        end do
    
    end subroutine Allocate_CasePerturbation
    
    !$
    !===============================================================================================
    ! finalizer for class of CasePerturbation
    !===============================================================================================
    subroutine Free_CasePerturbation (this)
        
        class(CasePerturbation), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%perts))  then
            do i = LBOUND(this%perts, dim=1), UBOUND(this%perts, dim=1)
                call this%perts(i)%clean ()
            end do
            
            ! free itself
            deallocate(this%perts)
        end if
    
    end subroutine Free_CasePerturbation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CasePerturbation_configure (this, i_pert, assign, ID, configure)
        
        class(CasePerturbation), intent(in out)  :: this
        integer, intent(in)  :: i_pert
        integer, intent(in)  :: assign(:)
        integer, intent(in)  :: ID(:)
        integer, intent(in)  :: configure(:, :)
        
        integer  :: ia, i
        
        do ia = 1, ns%state%layer
            plane: do i = 1, SIZE(ID, dim=1)                                    ! number of plane
                if (assign(ia) == ID(i)) then
                    this%perts(i_pert)%loading(:, ia) = configure(:, ID(i))
                    exit plane
                end if
            end do plane
        end do
    
    end subroutine Set_CasePerturbation_configure
    
    !$
    !===============================================================================================
    ! set configure which is the same as initial state
    !===============================================================================================
    subroutine Fix_CasePerturbation_configure (this, mat_info)
        
        class(CasePerturbation), intent(in out)  :: this
        type(Material), intent(in)               :: mat_info
        
        integer  :: i
        
        do i = 1, this%n_pert
            if (.NOT. this%perts(i)%is_new_configure)  then
                this%perts(i)%loading = mat_info%loading
            end if
        end do
    
    end subroutine Fix_CasePerturbation_configure
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CasePerturbation_xsec (this, i_pert, name, ID)
        
        class(CasePerturbation), intent(in out)  :: this
        integer, intent(in)            :: i_pert
        character(len=*), intent(in)   :: name
        integer, intent(in)            :: ID(3)
        
        integer  :: i_mat
        
        i_mat = ID(1)
        
        this%perts(i_pert)%mats(i_mat)%ID = ID(1)
        this%perts(i_pert)%mats(i_mat)%name = name
        this%perts(i_pert)%mats(i_mat)%case_ID = ID(2)
        this%perts(i_pert)%mats(i_mat)%burn_ID = ID(3)
    
    end subroutine Set_CasePerturbation_xsec
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_CasePerturbation_xsec (this, i_pert, xsec)
        
        class(CasePerturbation), intent(in out) :: this
        integer, intent(in)                 :: i_pert
        type(CrossSection), intent(in out)  :: xsec
        
        ! local variables
        type(cross_section_info_tp)  :: a_xsec(ns%state%mat)
        character(len=MAX_WORD_LEN)  :: dir
        character(len=MAX_WORD_LEN)  :: file_name
        character(len=MAX_WORD_LEN)  :: extension
        integer  :: case_ID, burn_ID

        integer  :: im, iz, ia, ig
        
        do im = 1, ns%state%mat
            dir = DIR_XSEC
            file_name = this%perts(i_pert)%mats(im)%name
            case_ID = this%perts(i_pert)%mats(im)%case_ID
            burn_ID = this%perts(i_pert)%mats(im)%burn_ID
            
            call a_xsec(im)%alloc ()
            call Get_file_extension (file_name, extension)
            
            if (TRIM(extension) == EXTEND_H5)  then
                call xsec_read_known (dir, file_name, case_ID, burn_ID, a_xsec(im))
            else
                call xsec_read_known (dir, file_name, a_xsec(im))
            end if
        end do
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = this%perts(i_pert)%loading(iz, ia)
                xsec%matrixs(iz, ia) = a_xsec(im)
            end do
        end do
        
        do im = 1, ns%state%mat
            call a_xsec(im)%clean ()
        end do
    
    end subroutine Get_CasePerturbation_xsec
    
end module perturbation_header
