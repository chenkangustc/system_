!$
!===================================================================================================
!
!   class for material and cross section
!    -- material: this term used to descrpt the input
!    -- cross section: this is really used in core code
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          CrossSection
!                               CrossSectionInput
!                               KineticsParameter
!                               KineticsParameterInput
!                               Material
!                               ExternalSource
!
!===================================================================================================
module material_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV

    use interpolation,              only : two_variables_interpolation
    use vector_operation,           only : vector_to_matrix
    use exception_header,           only : WarningCollector
    
    use geometry_header,            only : Meshing, Geometry
    
    implicit none
    private 
    public  :: CrossSection, CrossSectionInput
    public  :: KineticsParameter, KineticsParameterInput
    public  :: Material
    public  :: ExternalSource
    
    public  :: cross_section_info_tp, kinetics_parameter_info_tp
    
    type(WarningCollector)  :: a_warning

    ! --------------------------------------------------------------------------
    ! NOTE:
    !   "cross_section_info_tp" & "kinetics_parameter_info_tp" is public
    !   but is should be used only for Control Rod cross section in another file
    !   and this point should be fixed in the future
    ! --------------------------------------------------------------------------
    ! type for steady cross section
    type  :: cross_section_info_tp
        real(KREAL), public, allocatable  :: chi_steady(:)                  ! steady fission spectrum (1)
        real(KREAL), public, allocatable  :: sigma_t(:)                     ! macroscopic total xs (cm^-1)
        real(KREAL), public, allocatable  :: sigma_f_nu(:)                  ! macroscopic fission xs * nu (cm^-1)
        real(KREAL), public, allocatable  :: sigma_f_kappa(:)               ! macroscopic fission xs * kappa (J/cm)
        real(KREAL), public, allocatable  :: sigma_s(:, :, :)               ! macroscopic scatter xs (cm^-1)
        
        ! backup only, under some situation may needed
        real(KREAL), public, allocatable  :: sigma_f(:)                     ! macroscopic fission xs (cm^-1)
        real(KREAL), public, allocatable  :: nu(:)                          ! neutron number release per fission (1)
        real(KREAL), public, allocatable  :: kappa(:)                       ! energy release per fission (J)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_cross_section_info_tp
        procedure, public  :: clean => Free_cross_section_info_tp
        procedure, public  :: print => Print_cross_section_info_tp
        procedure, public  :: set => Set_cross_section_info_tp
        procedure, public  :: count => Count_cross_section_info_tp
        procedure, public  :: assign => Assign_cross_section_info_tp
        procedure, public  :: interpolate => Interpolate_cross_section_info_tp
        procedure, public  :: zero => Set_cross_section_info_tp_zero
        procedure, public  :: weight => Weight_cross_section_info_tp
        generic,   public  :: assignment (=) => Equal_cross_section_info_tp
        generic,   public  :: operator (.ADD.) => Add_cross_section_info_tp
        generic,   public  :: operator (.SUB.) => Sub_cross_section_info_tp
        generic,   public  :: operator (.MULTI.) => Multi_cross_section_info_tp
        procedure, private :: Equal_cross_section_info_tp
        procedure, private :: Add_cross_section_info_tp
        procedure, private :: Sub_cross_section_info_tp
        procedure, private :: Multi_cross_section_info_tp
    end type cross_section_info_tp
    
    ! type for transient kinetics parameter
    type  :: kinetics_parameter_info_tp
        real(KREAL), public, allocatable  :: chi_delay(:, :)                ! delay fission spectrum (1)
        real(KREAL), public, allocatable  :: velocity(:)                    ! neutron velocity (cm/s)
        real(KREAL), public, allocatable  :: beta(:)                        ! fraction of delay neutron (1)
        real(KREAL), public, allocatable  :: lambda(:)                      ! decay constant (s^-1)
        
        ! used for SRAC kinetics parameter format
        real(KREAL), public, allocatable  :: sigma_bvf(:, :)                ! see SRAC manual
        real(KREAL), public, allocatable  :: sigma_bvl(:, :)                ! see SRAC manual
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_kinetics_parameter_info_tp
        procedure, public  :: clean => Free_kinetics_parameter_info_tp
        procedure, public  :: print => Print_kinetics_parameter_info_tp 
        procedure, public  :: count => Count_kinetics_parameter_info_tp
        procedure, public  :: assign => Assign_kinetics_parameter_info_tp
        procedure, public  :: interpolate => Interpolate_kinetics_parameter_info_tp
        procedure, public  :: zero => Set_kinetics_parameter_info_tp_zero
        procedure, public  :: weight => Weight_kinetics_parameter_info_tp
        generic,   public  :: assignment (=) => Equal_kinetics_parameter_info_tp
        generic,   public  :: operator (.ADD.) => Add_kinetics_parameter_info_tp
        generic,   public  :: operator (.SUB.) => Sub_kinetics_parameter_info_tp
        generic,   public  :: operator (.MULTI.) => Multi_kinetics_parameter_info_tp
        procedure, private :: Equal_kinetics_parameter_info_tp
        procedure, private :: Add_kinetics_parameter_info_tp
        procedure, private :: Sub_kinetics_parameter_info_tp
        procedure, private :: Multi_kinetics_parameter_info_tp
    end type kinetics_parameter_info_tp
    
    ! type for material informaton
    type, private  :: material_info_tp
        integer, public                         :: ID                           ! material identification
        integer, public                         :: ID_add                       ! material identification with CR
        integer, public                         :: ID_sub                       ! material identification without CR        
        character(len=MAX_WORD_LEN), public     :: name                         ! material name
        character(len=MAX_WORD_LEN), public     :: file                         ! file nam which contain xsec information
        logical, public                         :: is_CR    = .FALSE.           ! is control rod material ?
        logical, public                         :: is_fuel  = .TRUE.            ! is fuel material or is fission material
    end type material_info_tp
    
    ! type for external source information
    type, private  :: external_source_info_tp
        real(KREAL), public, allocatable  :: intensity(:)                       ! external source intensity (cm^-3*s^-1)
        real(KREAL), public  :: memory = REAL_ZERO        
    contains
        procedure, private  :: alloc => Allocate_external_source_info_tp
        procedure, private  :: clean => Free_external_source_info_tp
    end type external_source_info_tp
    
    ! --------------------------------------------------------------------------
    ! type for cross section for input
    type CrossSectionInput 
        type(cross_section_info_tp), public, allocatable   :: mats(:)          ! per material
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_CrossSectionInput
        procedure, public  :: clean =>  Free_CrossSectionInput
        procedure, public  :: critical => Adjust_mats_critical
    end type CrossSectionInput
    
    ! type for cross section of the whole energy group
    type CrossSection
        type(cross_section_info_tp), public, allocatable   :: matrixs(:, :)    ! per zone per layer
        logical, public          :: is_critical = .FALSE.
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_CrossSection
        procedure, public  :: critical => Adjust_matrixs_critical
        procedure, public  :: print => Print_CrossSection
        procedure, public  :: clean =>  Free_CrossSection
        
        procedure, public  :: transpose_scat => Transpose_scatter_xsec
        procedure, public  :: is_fission => Is_contain_fission_material
        procedure, public  :: zero => Set_CrossSection_zero
    end type CrossSection
     
    ! type for kinetics parameter for input
    type KineticsParameterInput
        type(kinetics_parameter_info_tp), public, allocatable  :: mats(:)      ! per material
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_KineticsParameterInput
        procedure, public  :: clean => Free_KineticsParameterInput
    end type KineticsParameterInput
    
    ! type for kinetics parameter of the whole delay neutron group
    type KineticsParameter
        type(kinetics_parameter_info_tp), public, allocatable  :: matrixs(:, :)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_KineticsParameter
        procedure, public  :: clean => Free_KineticsParameter
        procedure, public  :: print => Print_KineticsParameter
        procedure, public  :: zero => Set_KineticsParameter_zero
    end type KineticsParameter
    
    ! type for materials information and configuration
    type Material
        integer, public, allocatable                       :: loading(:, :)     ! material loading ID per zone per layer ?
        type(material_info_tp), public, allocatable        :: libs(:)           ! material information library
        logical, public, allocatable                       :: mask_core(:, :)
        logical, public, allocatable                       :: mask_axial(:)
        logical, public, allocatable                       :: mask_radial(:)
    contains
        procedure, public  :: alloc => Allocate_Material
        procedure, public  :: clean => Free_Material
        procedure, public  :: set_row => Set_Material_by_row
        procedure, public  :: set_col => Set_Material_by_column
        procedure, public  :: set_mask => Set_Material_mask
        procedure, public  :: fix => Fix_Material_import
        procedure, public  :: get_name => Get_material_name
        procedure, public  :: get_file => Get_material_file
        procedure, public  :: get_ID => Get_material_ID
        procedure, public  :: is_fuel => Is_materil_fuel
        procedure, public  :: is_CR => Is_material_CR
    end type Material
    
    ! type for external source
    type  ExternalSource
        integer, public, allocatable                        :: adding(:, :)                         ! source loading ID per zone per layer
        type(external_source_info_tp), public, allocatable  :: kinds(:)                             ! source per kind, for external source input
        type(external_source_info_tp), public, allocatable  :: matrixs(:,:)                         ! source for current status
        
        logical, public                                     :: is_adjust                            ! is adjust to power level
        real(KREAL), public                                 :: iter_normal                          ! normal factor when iteration        
        type(external_source_info_tp), public, allocatable  :: iter_scalar(:, :)                    ! source used iteration, NOTE: per nodal with direction
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_ExternalSource
        procedure, public  :: set_row => Set_ExternalSource_by_row
        procedure, public  :: set_col => Set_ExternalSource_by_column
        procedure, public  :: clean =>  Free_ExternalSource
        procedure, public  :: adjust => Adjust_ExternalSource
        procedure, public  :: normal => Normal_ExternalSource
        procedure, public  :: strength => Get_ExternalSource_strength
        procedure, public  :: print => Print_ExternalSource
    end type ExternalSource
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_cross_section_info_tp (this)
        
        class(cross_section_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%chi_steady(ns%state%ng), stat=i_allocate)
        allocate(this%sigma_t(ns%state%ng), stat=i_allocate)
        allocate(this%sigma_f_nu(ns%state%ng), stat=i_allocate)
        allocate(this%sigma_f_kappa(ns%state%ng), stat=i_allocate)
        allocate(this%sigma_s(ns%state%ng, ns%state%ng, ns%deduce%scat_xs), stat=i_allocate)
        
        allocate(this%sigma_f(ns%state%ng), stat=i_allocate)
        allocate(this%nu(ns%state%ng), stat=i_allocate)
        allocate(this%kappa(ns%state%ng), stat=i_allocate)
        
        this%chi_steady     = REAL_ZERO
        this%sigma_t        = REAL_ZERO
        this%sigma_f_nu     = REAL_ZERO
        this%sigma_f_kappa  = REAL_ZERO
        this%sigma_s        = REAL_ZERO

        this%sigma_f        = REAL_ZERO
        this%nu             = REAL_ZERO
        this%kappa          = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%chi_steady)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_t)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_f_nu)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_f_kappa)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_s)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_f)
        this%memory = this%memory + REAL_BYTE * SIZE(this%nu)
        this%memory = this%memory + REAL_BYTE * SIZE(this%kappa)
        
    end subroutine Allocate_cross_section_info_tp

    !$
    !===============================================================================================
    ! finalizer for class of cross_section_info_tp
    !===============================================================================================
    subroutine Free_cross_section_info_tp (this)
        
        class(cross_section_info_tp), intent(in out)  :: this
        
        if (allocated(this%chi_steady))         deallocate(this%chi_steady)
        if (allocated(this%sigma_t))            deallocate(this%sigma_t)
        if (allocated(this%sigma_f_nu))         deallocate(this%sigma_f_nu)
        if (allocated(this%sigma_f_kappa))      deallocate(this%sigma_f_kappa)
        if (allocated(this%sigma_s))            deallocate(this%sigma_s)

        if (allocated(this%sigma_f))            deallocate(this%sigma_f)
        if (allocated(this%nu))                 deallocate(this%nu)
        if (allocated(this%kappa))              deallocate(this%kappa)
        
        this%memory = REAL_ZERO
        
    end subroutine Free_cross_section_info_tp
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Print_cross_section_info_tp (this, unit_)
    
        class(cross_section_info_tp), intent(in)  :: this
        integer, intent(in)  :: unit_
        
        integer  :: ic, ig
        
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%chi_steady(:)
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%sigma_t(:)
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%sigma_f_nu(:)
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%sigma_f_kappa(:)
        
        do ic = 1, ns%deduce%scat_xs
            do ig = 1, ns%state%ng
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%sigma_s(ig, :, ic)
            end do
        end do
    
    end subroutine Print_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! Set sigma_f & sigma_f_kappa
    !===============================================================================================
    subroutine Set_cross_section_info_tp (this, is_nu, is_kappa)
        
        class(cross_section_info_tp), intent(in out)  :: this
        logical, intent(in)  :: is_nu
        logical, intent(in)  :: is_kappa
        integer  :: i
        
        this%nu     = 2.43
        this%kappa  = 200.0 * 1.6E-13                                           ! 200 MeV
        
        do i = 1, SIZE(this%sigma_f_nu)
            if (is_nu )  then
                this%sigma_f(i) = this%sigma_f_nu(i) / this%nu(i)
            end if

            if (is_kappa )  then
                this%sigma_f_kappa(i) = this%sigma_f_nu(i) * this%kappa(i) / this%nu(i)
            end if
        end do
            
    end subroutine Set_cross_section_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Count_cross_section_info_tp (this)  result(counting)
        
        class(cross_section_info_tp), intent(in out)  :: this
        integer  :: counting
        
        counting = 0
        counting = counting + SIZE(this%chi_steady)
        counting = counting + SIZE(this%sigma_t)
        counting = counting + SIZE(this%sigma_f_nu)
        counting = counting + SIZE(this%sigma_f_kappa)
        counting = counting + SIZE(this%sigma_s)
    
    end function Count_cross_section_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Assign_cross_section_info_tp (this, values)
    
        class(cross_section_info_tp), intent(in out)  :: this
        real(8), intent(in)  :: values(1:)                                      ! real(8) for consisten with Lilac
        
        real(KREAL)  :: tmp_value(SIZE(values))
        integer  :: i_start, i_end
        
        if (SIZE(values) /= this%count ())  then
        
        end if
        tmp_value = values
        
        i_start = 1
        i_end = SIZE(this%chi_steady)
        call vector_to_matrix(this%chi_steady, tmp_value(i_start: i_end))
        
        i_start = i_end + 1
        i_end = i_end + SIZE(this%sigma_t)
        call vector_to_matrix(this%sigma_t, tmp_value(i_start: i_end))
        
        i_start = i_end + 1
        i_end = i_end + SIZE(this%sigma_f_nu)
        call vector_to_matrix(this%sigma_f_nu, tmp_value(i_start: i_end))
        
        i_start = i_end + 1
        i_end = i_end + SIZE(this%sigma_f_kappa)
        call vector_to_matrix(this%sigma_f_kappa, tmp_value(i_start: i_end))
        
        i_start = i_end + 1
        i_end = i_end + SIZE(this%sigma_s)
        call vector_to_matrix(this%sigma_s, tmp_value(i_start: i_end))

    end subroutine Assign_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! interpolation by two variables
    !===============================================================================================
    subroutine Interpolate_cross_section_info_tp (this, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        class(cross_section_info_tp), intent(in out)  :: this
        type(cross_section_info_tp), intent(in)       :: p00
        type(cross_section_info_tp), intent(in)       :: p01
        type(cross_section_info_tp), intent(in)       :: p10
        type(cross_section_info_tp), intent(in)       :: p11
        real(KREAL), intent(in)  :: x_min
        real(KREAL), intent(in)  :: x
        real(KREAL), intent(in)  :: x_max
        real(KREAL), intent(in)  :: y_min
        real(KREAL), intent(in)  :: y
        real(KREAL), intent(in)  :: y_max
        
        call two_variables_interpolation (this%chi_steady, p00%chi_steady, p01%chi_steady, p10%chi_steady, p11%chi_steady,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_t, p00%sigma_t, p01%sigma_t, p10%sigma_t, p11%sigma_t,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_f_nu, p00%sigma_f_nu, p01%sigma_f_nu, p10%sigma_f_nu, p11%sigma_f_nu,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_f_kappa, p00%sigma_f_kappa, p01%sigma_f_kappa, p10%sigma_f_kappa, p11%sigma_f_kappa,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_s, p00%sigma_s, p01%sigma_s, p10%sigma_s, p11%sigma_s,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_f, p00%sigma_f, p01%sigma_f, p10%sigma_f, p11%sigma_f,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%nu, p00%nu, p01%nu, p10%nu, p11%nu,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%kappa, p00%kappa, p01%kappa, p10%kappa, p11%kappa,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
    
    end subroutine Interpolate_cross_section_info_tp
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Set_cross_section_info_tp_zero (this)
        
        class(cross_section_info_tp), intent(in out)  :: this
        
        this%chi_steady     = REAL_ZERO
        this%sigma_t        = REAL_ZERO
        this%sigma_f_nu     = REAL_ZERO
        this%sigma_f_kappa  = REAL_ZERO
        this%sigma_s        = REAL_ZERO
        
        this%sigma_f        = REAL_ZERO
        this%nu             = REAL_ZERO
        this%kappa          = REAL_ZERO
    
    end subroutine Set_cross_section_info_tp_zero
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Weight_cross_section_info_tp (this, xsec, volume, weight)
    
        class(cross_section_info_tp), intent(in out)  :: this
        type(cross_section_info_tp), intent(in)       :: xsec(:)
        real(KREAL), intent(in)                   :: volume(:)
        real(KREAL), intent(in)                   :: weight(:)
        
        real(KREAL)  :: denominator
        integer          :: i_start, i_end, i
        
        i_start = LBOUND(weight, dim=1)
        i_end = UBOUND(weight, dim=1)
        
        call this%zero ()
        denominator = 0.0
        
        do i = i_start, i_end
            this%chi_steady     = this%chi_steady + xsec(i)%chi_steady * volume(i) * weight(i)
            this%sigma_t        = this%sigma_t + xsec(i)%sigma_t * volume(i) * weight(i)
            this%sigma_f_nu     = this%sigma_f_nu + xsec(i)%sigma_f_nu * volume(i) * weight(i)
            this%sigma_f_kappa  = this%sigma_f_kappa + xsec(i)%sigma_f_kappa * volume(i) * weight(i)
            this%sigma_s        = this%sigma_s + xsec(i)%sigma_s * volume(i) * weight(i)
            
            this%sigma_f        = this%sigma_f + xsec(i)%sigma_f * volume(i) * weight(i)
            this%nu             = this%nu + xsec(i)%nu * volume(i) * weight(i)
            this%kappa          = this%kappa + xsec(i)%kappa * volume(i) * weight(i)
            
            denominator = denominator + volume(i) * weight(i)
        end do
        
        this%chi_steady     = this%chi_steady / denominator
        this%sigma_t        = this%sigma_t / denominator
        this%sigma_f_nu     = this%sigma_f_nu / denominator
        this%sigma_f_kappa  = this%sigma_f_kappa / denominator
        this%sigma_s        = this%sigma_s / denominator
        
        this%sigma_f        = this%sigma_f / denominator
        this%nu             = this%nu / denominator
        this%kappa          = this%kappa / denominator
    
    end subroutine Weight_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! assignment implement for cross_section_info_tp
    !===============================================================================================
    subroutine Equal_cross_section_info_tp (left, right)
        
        class(cross_section_info_tp), intent(in out)  :: left
        type(cross_section_info_tp), intent(in)       :: right
        
        ! allocate lhs first
        call left%alloc ()
        
        left%chi_steady     = right%chi_steady
        left%sigma_t        = right%sigma_t
        left%sigma_f_nu     = right%sigma_f_nu
        left%sigma_f_kappa  = right%sigma_f_kappa
        left%sigma_s        = right%sigma_s
       
        left%sigma_f        = right%sigma_f
        left%nu             = right%nu
        left%kappa          = right%kappa
    
    end subroutine Equal_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! (+) operator
    !===============================================================================================
    function Add_cross_section_info_tp (inp1, inp2)  result(outp)
        
        class(cross_section_info_tp), intent(in)  :: inp1
        type(cross_section_info_tp), intent(in)  :: inp2
        type(cross_section_info_tp)  :: outp
        
        ! allocate result first
        call outp%alloc ()
        
        outp%chi_steady     = inp1%chi_steady    + inp2%chi_steady
        outp%sigma_t        = inp1%sigma_t       + inp2%sigma_t
        outp%sigma_f_nu     = inp1%sigma_f_nu    + inp2%sigma_f_nu
        outp%sigma_f_kappa  = inp1%sigma_f_kappa + inp2%sigma_f_kappa
        outp%sigma_s        = inp1%sigma_s       + inp2%sigma_s
                                                 
        outp%sigma_f        = inp1%sigma_f       + inp2%sigma_f
        outp%nu             = inp1%nu            + inp2%nu
        outp%kappa          = inp1%kappa         + inp2%kappa
        
    end function Add_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! (-) operator
    !===============================================================================================
    function Sub_cross_section_info_tp (inp1, inp2)  result(outp)
        
        class(cross_section_info_tp), intent(in)  :: inp1
        type(cross_section_info_tp), intent(in)  :: inp2
        type(cross_section_info_tp)  :: outp
        
        ! allocate result first
        call outp%alloc ()
        
        outp%chi_steady     = inp1%chi_steady    - inp2%chi_steady
        outp%sigma_t        = inp1%sigma_t       - inp2%sigma_t
        outp%sigma_f_nu     = inp1%sigma_f_nu    - inp2%sigma_f_nu
        outp%sigma_f_kappa  = inp1%sigma_f_kappa - inp2%sigma_f_kappa
        outp%sigma_s        = inp1%sigma_s       - inp2%sigma_s
                                                 
        outp%sigma_f        = inp1%sigma_f       - inp2%sigma_f
        outp%nu             = inp1%nu            - inp2%nu
        outp%kappa          = inp1%kappa         - inp2%kappa
        
    end function Sub_cross_section_info_tp
    
    !$
    !===============================================================================================
    ! (x) operator
    !===============================================================================================
    function Multi_cross_section_info_tp (inp1, factor)  result(outp)
        
        class(cross_section_info_tp), intent(in)  :: inp1
        real(KREAL), intent(in) :: factor
        type(cross_section_info_tp)  :: outp
        
        ! allocate result fist
        call outp%alloc ()
        
        outp%chi_steady     = inp1%chi_steady    * factor
        outp%sigma_t        = inp1%sigma_t       * factor
        outp%sigma_f_nu     = inp1%sigma_f_nu    * factor
        outp%sigma_f_kappa  = inp1%sigma_f_kappa * factor
        outp%sigma_s        = inp1%sigma_s       * factor
                                                 
        outp%sigma_f        = inp1%sigma_f       * factor
        outp%nu             = inp1%nu            * factor
        outp%kappa          = inp1%kappa         * factor
        
    end function Multi_cross_section_info_tp
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Allocate_kinetics_parameter_info_tp (this)
        
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%chi_delay(nt%state%dg, ns%state%ng), stat=i_allocate)
        allocate(this%velocity(ns%state%ng), stat=i_allocate)                   ! velocity for different energy group
        allocate(this%beta(nt%state%dg), stat=i_allocate)
        allocate(this%lambda(nt%state%dg), stat=i_allocate)
        
        allocate(this%sigma_bvf(nt%state%dg, ns%state%ng), stat=i_allocate)
        allocate(this%sigma_bvl(nt%state%dg, ns%state%ng), stat=i_allocate)
        
        this%chi_delay = REAL_ZERO
        this%velocity  = REAL_ZERO
        this%beta      = REAL_ZERO
        this%lambda    = REAL_ZERO
        
        this%sigma_bvf = REAL_ZERO
        this%sigma_bvl = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%chi_delay)
        this%memory = this%memory + REAL_BYTE * SIZE(this%velocity)
        this%memory = this%memory + REAL_BYTE * SIZE(this%beta)
        this%memory = this%memory + REAL_BYTE * SIZE(this%lambda)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_bvf)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_bvl)
    
    end subroutine Allocate_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of kinetics_parameter_info_tp
    !===============================================================================================
    subroutine Free_kinetics_parameter_info_tp (this)
        
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        
        if (allocated(this%chi_delay))          deallocate(this%chi_delay)
        if (allocated(this%velocity))           deallocate(this%velocity)
        if (allocated(this%beta))               deallocate(this%beta)
        if (allocated(this%lambda))             deallocate(this%lambda)
        
        if (allocated(this%sigma_bvf))          deallocate(this%sigma_bvf)
        if (allocated(this%sigma_bvl))          deallocate(this%sigma_bvl)
        
        this%memory = REAL_ZERO
        
    end subroutine Free_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_kinetics_parameter_info_tp (this, unit_)
    
        class(kinetics_parameter_info_tp), intent(in)  :: this
        integer, intent(in) :: unit_
        
        integer  :: ia, iz, ip
        
        do ip = 1, nt%state%dg
            write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%chi_delay(ip, :)
        end do
        
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%velocity(:)
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%beta
        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%lambda
    
    end subroutine Print_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Count_kinetics_parameter_info_tp (this)  result(counting)
        
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        integer  :: counting
        
        counting = 0
        counting = counting + SIZE(this%chi_delay)
        counting = counting + SIZE(this%velocity)
        
        select case(TRIM(nt%flag%kinetics_type))
        case ('NORMAL')
            counting = counting + SIZE(this%beta)
            counting = counting + SIZE(this%lambda)
            
        case ('SRAC')
            counting = counting + SIZE(this%sigma_bvf)
            counting = counting + SIZE(this%sigma_bvl)
        end select
        
    end function Count_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Assign_kinetics_parameter_info_tp (this, values)
    
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        real(8), intent(in)  :: values(1:)                                      ! real(8) for consisten with Lilac
        
        real(KREAL)  :: tmp_value(SIZE(values))
        integer  :: i_start, i_end 
        
        if (SIZE(values) /= this%count ())  then
        
        end if
        tmp_value = values
        
        i_start = 1
        i_end = SIZE(this%chi_delay)
        call vector_to_matrix(this%chi_delay, tmp_value(i_start: i_end))
        
        i_start = i_end + 1
        i_end = i_end + SIZE(this%velocity)
        call vector_to_matrix(this%velocity, tmp_value(i_start: i_end))
        
        select case(TRIM(nt%flag%kinetics_type))
        case ('NORMAL')
            i_start = i_end + 1
            i_end = i_end + SIZE(this%beta)
            call vector_to_matrix(this%beta, tmp_value(i_start: i_end))
            
            i_start = i_end + 1
            i_end = i_end + SIZE(this%lambda)
            call vector_to_matrix(this%lambda, tmp_value(i_start: i_end))
        
        case ('SRAC')
            i_start = i_end + 1
            i_end = i_end + SIZE(this%sigma_bvf)
            call vector_to_matrix(this%sigma_bvf, tmp_value(i_start: i_end))
            
            i_start = i_end + 1
            i_end = i_end + SIZE(this%sigma_bvl)
            call vector_to_matrix(this%sigma_bvl, tmp_value(i_start: i_end))
        end select
        
    end subroutine Assign_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! (+) operator
    !===============================================================================================
    function Add_kinetics_parameter_info_tp (inp1, inp2)  result(outp)
        
        class(kinetics_parameter_info_tp), intent(in)  :: inp1
        type(kinetics_parameter_info_tp), intent(in)  :: inp2
        type(kinetics_parameter_info_tp)  :: outp
        
        ! allocate result first
        call outp%alloc ()
        
        outp%chi_delay  = inp1%chi_delay + inp2%chi_delay
        outp%velocity   = inp1%velocity  + inp2%velocity
        outp%beta       = inp1%beta      + inp2%beta
        outp%lambda     = inp1%lambda    + inp2%lambda
                                                 
        outp%sigma_bvf  = inp1%sigma_bvf + inp2%sigma_bvf
        outp%sigma_bvl  = inp1%sigma_bvl + inp2%sigma_bvl
        
    end function Add_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! (-) operator
    !===============================================================================================
    function Sub_kinetics_parameter_info_tp (inp1, inp2)  result(outp)
        
        class(kinetics_parameter_info_tp), intent(in)  :: inp1
        type(kinetics_parameter_info_tp), intent(in)  :: inp2
        type(kinetics_parameter_info_tp)  :: outp
        
        ! allocate result first
        call outp%alloc ()
        
        outp%chi_delay  = inp1%chi_delay - inp2%chi_delay
        outp%velocity   = inp1%velocity  - inp2%velocity
        outp%beta       = inp1%beta      - inp2%beta
        outp%lambda     = inp1%lambda    - inp2%lambda
                                        
        outp%sigma_bvf  = inp1%sigma_bvf - inp2%sigma_bvf
        outp%sigma_bvl  = inp1%sigma_bvl - inp2%sigma_bvl
        
    end function Sub_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! (x) operator
    !===============================================================================================
    function Multi_kinetics_parameter_info_tp (inp1, factor)  result(outp)
        
        class(kinetics_parameter_info_tp), intent(in)  :: inp1
        real(KREAL), intent(in)  :: factor
        type(kinetics_parameter_info_tp)  :: outp
        
        ! allocate result first
        call outp%alloc ()
        
        outp%chi_delay  = inp1%chi_delay * factor
        outp%velocity   = inp1%velocity  * factor
        outp%beta       = inp1%beta      * factor
        outp%lambda     = inp1%lambda    * factor
        
        outp%sigma_bvf  = inp1%sigma_bvf * factor
        outp%sigma_bvl  = inp1%sigma_bvl * factor
        
    end function Multi_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! interpolation by two variables
    !===============================================================================================
    subroutine Interpolate_kinetics_parameter_info_tp (this, p00, p01, p10, p11, x_min, x, x_max, y_min, y, y_max)
        
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        type(kinetics_parameter_info_tp), intent(in)       :: p00
        type(kinetics_parameter_info_tp), intent(in)       :: p01
        type(kinetics_parameter_info_tp), intent(in)       :: p10
        type(kinetics_parameter_info_tp), intent(in)       :: p11
        real(KREAL), intent(in)  :: x_min
        real(KREAL), intent(in)  :: x
        real(KREAL), intent(in)  :: x_max
        real(KREAL), intent(in)  :: y_min
        real(KREAL), intent(in)  :: y
        real(KREAL), intent(in)  :: y_max
        
        call two_variables_interpolation (this%chi_delay, p00%chi_delay, p01%chi_delay, p10%chi_delay, p11%chi_delay,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%velocity, p00%velocity, p01%velocity, p10%velocity, p11%velocity,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%beta, p00%beta, p01%beta, p10%beta, p11%beta,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%lambda, p00%lambda, p01%lambda, p10%lambda, p11%lambda,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_bvf, p00%sigma_bvf, p01%sigma_bvf, p10%sigma_bvf, p11%sigma_bvf,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
        call two_variables_interpolation (this%sigma_bvl, p00%sigma_bvl, p01%sigma_bvl, p10%sigma_bvl, p11%sigma_bvl,  &
                                        &   x_min, x, x_max, y_min, y, y_max)
    
    end subroutine Interpolate_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Set_kinetics_parameter_info_tp_zero (this)
    
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        
        this%chi_delay = REAL_ZERO
        this%velocity  = REAL_ZERO
        this%beta      = REAL_ZERO
        this%lambda    = REAL_ZERO
        
        this%sigma_bvf = REAL_ZERO
        this%sigma_bvl = REAL_ZERO
    
    end subroutine Set_kinetics_parameter_info_tp_zero
    
    !$
    !===============================================================================================
    !  
    !===============================================================================================
    subroutine Weight_kinetics_parameter_info_tp (this, param, volume, weight)
    
        class(kinetics_parameter_info_tp), intent(in out)  :: this
        type(kinetics_parameter_info_tp), intent(in)       :: param(:)
        real(KREAL), intent(in)                        :: volume(:)
        real(KREAL), intent(in)                        :: weight(:)
        
        real(KREAL)  :: denominator
        integer          :: i_start, i_end, i
        
        i_start = LBOUND(weight, dim=1)
        i_end = UBOUND(weight, dim=1)
        
        call this%zero ()    
        denominator = 0.0
        
        do i = i_start, i_end
            this%chi_delay = this%chi_delay + param(i)%chi_delay * volume(i) * weight(i)
            this%velocity  = this%velocity  + param(i)%velocity * volume(i) * weight(i)
            this%beta      = this%beta      + param(i)%beta * volume(i) * weight(i)
            this%lambda    = this%lambda    + param(i)%lambda * volume(i) * weight(i)
            
            this%sigma_bvf = this%sigma_bvf + param(i)%sigma_bvf * volume(i) * weight(i)
            this%sigma_bvl = this%sigma_bvl + param(i)%sigma_bvl * volume(i) * weight(i)
            
            denominator = denominator + volume(i) * weight(i)
        end do
        
        this%chi_delay = this%chi_delay / denominator
        this%velocity  = this%velocity / denominator
        this%beta      = this%beta / denominator
        this%lambda    = this%lambda / denominator
        
        this%sigma_bvf = this%sigma_bvf / denominator
        this%sigma_bvl = this%sigma_bvl / denominator
    
    end subroutine Weight_kinetics_parameter_info_tp
    
    !$
    !===============================================================================================
    ! assigement implement for kinetics_parameter_info_tp
    !===============================================================================================
    subroutine Equal_kinetics_parameter_info_tp (left, right)
        
        class(kinetics_parameter_info_tp), intent(in out)  :: left
        type(kinetics_parameter_info_tp), intent(in)       :: right
        
        ! allocate lhs first
        call left%alloc ()
        
        left%chi_delay  = right%chi_delay
        left%velocity   = right%velocity
        left%beta       = right%beta
        left%lambda     = right%lambda

        left%sigma_bvf  = right%sigma_bvf
        left%sigma_bvl  = right%sigma_bvl
        
    end subroutine Equal_kinetics_parameter_info_tp
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_CrossSectionInput (this)
        
        class(CrossSectionInput), intent(in out)  :: this
        integer  :: i_allocate
        integer  :: i
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%mats(ns%state%mat), stat=i_allocate)
        do i = 1, SIZE(this%mats, dim=1)
            call this%mats(i)%alloc ()
        end do
        
        this%memory = REAL_ZERO
        this%memory = this%mats(1)%memory * SIZE(this%mats)
    
    end subroutine Allocate_CrossSectionInput
    
    !$
    !===============================================================================================
    ! finalizer for class of CrossSectionInput
    !===============================================================================================
    subroutine Free_CrossSectionInput (this)
        
        class(CrossSectionInput), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%mats))  then
            do i = 1, SIZE(this%mats, dim=1)
                call this%mats(i)%clean ()
            end do
            ! free it self
            deallocate(this%mats)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_CrossSectionInput
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_CrossSection (this)
        
        class(CrossSection), intent(in out)  :: this
        integer  :: i_allocate
        integer  :: i, j
        
        ! check allocated status first
        call this%clean ()
        
        ! allocate itself first
        allocate(this%matrixs(ns%state%zone, ns%state%layer), stat=i_allocate)
        
        ! allocate nested type object
        do i = 1, SIZE(this%matrixs, dim=1)
            do j = 1, SIZE(this%matrixs, dim=2)
                call this%matrixs(i, j)%alloc ()
            end do
        end do
        
        this%memory = REAL_ZERO
        this%memory = this%matrixs(1, 1)%memory * SIZE(this%matrixs)
    
    end subroutine Allocate_CrossSection
    
    !$
    !===============================================================================================
    ! finalizer for class of CrossSection
    !===============================================================================================
    subroutine Free_CrossSection (this)
        
        class(CrossSection), intent(in out)  :: this
        integer  :: i, j
        
        if (allocated(this%matrixs))  then
            ! free the nested type objects
            do i = 1, SIZE(this%matrixs, dim=1)
                do j = 1, SIZE(this%matrixs, dim=2)
                    call this%matrixs(i, j)%clean ()
                end do
            end do
            
            ! free itself
            deallocate(this%matrixs)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_CrossSection
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_CrossSection_zero (this)
    
        class(CrossSection), intent(in out)  :: this
        integer  :: i, j
        
        do i = 1, SIZE(this%matrixs, dim=1)
            do j = 1, SIZE(this%matrixs, dim=2)
                call this%matrixs(i, j)%zero ()
            end do
        end do

    end subroutine Set_CrossSection_zero
    
    !$
    !===============================================================================================
    ! adjust value of nu, leading system to rigorous critical 
    !===============================================================================================
    subroutine Adjust_mats_critical (this, k_eff)
        
        class(CrossSectionInput), intent(in out)  :: this
        real(KREAL), intent(in)  :: k_eff                                   ! k_eff when initial iteration
        
        integer  :: im
        
        do im = 1, ns%state%mat
            this%mats(im)%nu = this%mats(im)%nu / k_eff
            this%mats(im)%sigma_f_nu = this%mats(im)%sigma_f_nu / k_eff
        end do
    
    end subroutine Adjust_mats_critical    
    
    !$
    !===============================================================================================
    ! adjust value of nu, leading system to rigorous critical 
    !===============================================================================================
    subroutine Adjust_matrixs_critical (this, k_eff)
        
        class(CrossSection), intent(in out)  :: this
        real(KREAL), intent(in)  :: k_eff                                   ! k_eff when initial iteration
        
        integer  :: iz, ia
        
        do iz = 1, ns%state%zone
            do ia = 1, ns%state%layer
                this%matrixs(iz, ia)%nu = this%matrixs(iz, ia)%nu / k_eff
                this%matrixs(iz, ia)%sigma_f_nu = this%matrixs(iz, ia)%sigma_f_nu / k_eff
            end do
        end do    
        
        this%is_critical = .TRUE.
    
    end subroutine Adjust_matrixs_critical
    
    !$
    !===============================================================================================
    ! transpose scatter xesc for adjoint calculation
    !===============================================================================================
    subroutine Transpose_scatter_xsec (this)
        
        class(CrossSection), intent(in out)  :: this
        
        ! local variables
        real(KREAL)  :: tmp_sigma_s(ns%state%ng, ns%state%ng, ns%deduce%scat_xs)
        integer  :: ig, iig, iz, ia
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                do ig = 1, ns%state%ng
                    do iig = 1, ns%state%ng
                        tmp_sigma_s(ig, iig, :) = this%matrixs(iz,ia)%sigma_s(ig, iig, :)
                    end do
                end do
                
                do ig = 1, ns%state%ng
                    do iig = 1, ns%state%ng
                        this%matrixs(iz,ia)%sigma_s(ig, iig, :) = tmp_sigma_s(iig, ig, :)
                    end do
                end do
            end do
        end do
    
    end subroutine Transpose_scatter_xsec
    
    !$
    !===============================================================================================
    ! wether contains fission material
    !===============================================================================================
    function Is_contain_fission_material (this)  result (is_true)
        
        class(CrossSection), intent(in)  :: this
        logical  :: is_true
        
        real(KREAL)  :: fission
        integer  :: ig, ia, iz
        
        is_true = .TRUE.                                                        ! default is true
        fission = 0.0
        
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    fission = fission + this%matrixs(iz, ia)%sigma_f_nu(ig)
                end do
            end do
        end do
        
        if (ABS(fission-0.0) < EPS_EQUAL)  then
            is_true = .FALSE.
        end if
    
    end function Is_contain_fission_material
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_CrossSection (this, unit_, zone, layer)
    
        class(CrossSection), intent(in)  :: this
        integer, intent(in)            :: unit_
        integer, intent(in), optional  :: zone
        integer, intent(in), optional  :: layer
        
        integer  :: zone_start, zone_end
        integer  :: layer_start, layer_end
        integer  :: ia, iz, ig, ic
        
        if (present(zone) .and. present(layer))  then
            zone_start = zone
            zone_end = zone
            layer_start = layer
            layer_end = layer
        else
            zone_start = 1
            zone_end = ns%state%zone
            layer_start = 1
            layer_end = ns%state%layer
        end if
        
        do ia = layer_start, layer_end
            do iz = zone_start, zone_end
                write(unit=unit_, fmt=*) TRIM(CHAR_SSUBMARK)
                write(unit=unit_, fmt=*) 'layer is:', ia, 'zone is :', iz
                
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%chi_steady(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%sigma_t(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%sigma_f_nu(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%sigma_f_kappa(:)
                
                do ic = 1, ns%deduce%scat_xs
                    do ig = 1, ns%state%ng
                        write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%sigma_s(ig, :, ic)
                    end do
                end do
            end do
        end do
    
    end subroutine Print_CrossSection

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_KineticsParameterInput (this)
        
        class(KineticsParameterInput), intent(in out)  :: this
        integer  :: i_allocate
        integer  :: i
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%mats(ns%state%mat), stat=i_allocate)
        do i = 1, SIZE(this%mats, dim=1)
            call this%mats(i)%alloc ()
        end do
        
        this%memory = REAL_ZERO
        this%memory = this%mats(1)%memory * SIZE(this%mats)
    
    end subroutine Allocate_KineticsParameterInput
    
    !$
    !===============================================================================================
    ! finalizer for class of KineticsParameterInput
    !===============================================================================================
    subroutine Free_KineticsParameterInput (this)
        
        class(KineticsParameterInput), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%mats))  then
            do i = 1, SIZE(this%mats, dim=1)
                call this%mats(i)%clean ()
            end do
            deallocate(this%mats)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_KineticsParameterInput
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_KineticsParameter_zero (this)
    
        class(KineticsParameter), intent(in out)  :: this
        integer  :: i, j
        
        do i = 1, SIZE(this%matrixs, dim=1)
            do j = 1, SIZE(this%matrixs, dim=2)
                call this%matrixs(i, j)%zero ()
            end do
        end do
        
    end subroutine Set_KineticsParameter_zero
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_KineticsParameter (this)
    
        class(KineticsParameter), intent(in out)  :: this
        integer  :: i_allocate
        integer  :: i, j
        
        ! check allocated status first
        call this%clean ()
        
        ! allocate itself first
        allocate(this%matrixs(ns%state%zone, ns%state%layer), stat=i_allocate)
        
        ! allocate nested type object
        do i = 1, SIZE(this%matrixs, dim=1) 
            do j = 1, SIZE(this%matrixs, dim=2)
                call this%matrixs(i, j)%alloc ()
            end do
        end do
        
        this%memory = REAL_ZERO
        this%memory = this%matrixs(1, 1)%memory * SIZE(this%matrixs)
    
    end subroutine Allocate_KineticsParameter
    
    !$
    !===============================================================================================
    ! finalizer for class of KineticsParameter
    !===============================================================================================
    subroutine Free_KineticsParameter (this)
    
        class(KineticsParameter), intent(in out)  :: this
        integer  :: i, j
        
        if (allocated(this%matrixs))  then
            ! free the nested type objects
            do i = 1, SIZE(this%matrixs, dim=1)
                do j = 1, SIZE(this%matrixs, dim=2)
                    call this%matrixs(i, j)%clean ()
                end do
            end do
            ! free it self
            deallocate(this%matrixs)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_KineticsParameter
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_KineticsParameter (this, unit_, zone, layer)
    
        class(KineticsParameter), intent(in)  :: this
        integer, intent(in)            :: unit_
        integer, intent(in), optional  :: zone
        integer, intent(in), optional  :: layer
        
        integer  :: zone_start, zone_end
        integer  :: layer_start, layer_end
        integer  :: ia, iz, ip
        
        if (present(zone) .and. present(layer))  then
            zone_start = zone
            zone_end = zone
            layer_start = layer
            layer_end = layer
        else
            zone_start = 1
            zone_end = ns%state%zone
            layer_start = 1
            layer_end = ns%state%layer
        end if
        
        do ia = layer_start, layer_end
            do iz = zone_start, zone_end
                write(unit=unit_, fmt=*) '______________________________________________'
                write(unit=unit_, fmt=*) 'layer is:', ia, 'zone is :', iz
                
                do ip = 1, nt%state%dg
                    write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%chi_delay(ip, :)
                end do
                
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%velocity(:)
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%beta
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%matrixs(iz,ia)%lambda
            end do
        end do
    
    end subroutine Print_KineticsParameter
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
        
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_Material (this)
        
        class(Material), intent(in out)  :: this
        
        integer  :: i_allocate
        integer  :: i
        
        ! check allocated status first
        call this%clean ()

        allocate(this%loading(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%libs(ns%state%mat), stat=i_allocate)
        
        allocate(this%mask_core(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%mask_axial(ns%state%layer), stat=i_allocate)
        allocate(this%mask_radial(ns%state%zone), stat=i_allocate)
        
        ! initalize value
        this%loading  = INT_ZERO
        do i = 1, SIZE(this%libs)
            this%libs(i)%ID      = INT_ZERO
            this%libs(i)%name    = CHAR_NULL
            this%libs(i)%is_CR   = .FALSE.
            this%libs(i)%is_fuel = .TRUE.
        end do
        
        this%mask_core = .TRUE.
        this%mask_axial = .TRUE.
        this%mask_radial = .TRUE.
    
    end subroutine Allocate_Material
    
    !$
    !===============================================================================================
    ! set material loading information
    !===============================================================================================
    subroutine Set_Material_by_row (this, assign, config)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: assign(:)
        integer, intent(in)  :: config(:, :)
        
        integer  :: ia, i

        do ia = 1, ns%state%layer
            this%loading(:, ia) = config(:, assign(ia))
        end do
    
    end subroutine Set_Material_by_row
    
    !$
    !===============================================================================================
    ! set material loading information
    !===============================================================================================
    subroutine Set_Material_by_column (this, assign, config)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: assign(:)
        integer, intent(in)  :: config(:, :)
        
        integer  :: iz, i
        
        do iz = 1, ns%state%zone
            this%loading(iz, :) = config(assign(iz), :)
        end do
    
    end subroutine Set_Material_by_column
    
    !$
    !===============================================================================================
    ! wether contains fission material, "TRUE" means active fuel zone
    !===============================================================================================
    subroutine Set_Material_mask (this, xsec)
        
        class(Material), intent(in out)  :: this
        type(CrossSection), intent(in)   :: xsec
        real(KREAL)  :: fission
        integer      :: ia, iz
        
        this%mask_core = .TRUE.
        do ia = 1, SIZE(xsec%matrixs, dim=2)
            do iz = 1, SIZE(xsec%matrixs, dim=1)
                fission = SUM(xsec%matrixs(iz, ia)%sigma_f_nu(:))
                
                if (ABS(fission) < EPS_ZERO)  then
                    this%mask_core(iz, ia) = .FALSE.
                end if
            end do
        end do
        
        this%mask_axial = .FALSE.
        do ia = 1, SIZE(this%mask_core, dim=2)
            if (ANY(this%mask_core(:, ia)))  then
                this%mask_axial(ia) = .TRUE.
            end if
        end do
        
        this%mask_radial = .FALSE.
        do iz = 1, SIZE(this%mask_core, dim=1)
            if (ANY(this%mask_core(iz, :)))  then
                this%mask_radial(iz) = .TRUE.
            end if
        end do
    
    end subroutine Set_Material_mask
    
    !$
    !===============================================================================================
    ! rerange material configuration for xsec import
    !===============================================================================================
    subroutine Fix_Material_import (this)
    
        class(Material), intent(in out)  :: this
        
        ! local
        integer  :: i_allocate
        integer  :: ia, iz, im
        
        if (allocated(this%libs))   then
            deallocate(this%libs)
            allocate(this%libs(ns%state%mat), stat=i_allocate)
        end if
        
        im = 0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = im + 1
                this%loading(iz, ia) = im
                this%libs(im)%ID = im
                this%libs(im)%name = ' '
                this%libs(im)%file = ' '
                this%libs(im)%is_CR = .FALSE.
            end do
        end do
    
    end subroutine Fix_Material_import
    
    !$
    !===============================================================================================
    ! finalizer for class of Material
    !===============================================================================================
    subroutine Free_Material (this)
        
        class(Material), intent(in out)  :: this

        if (allocated(this%loading))        deallocate(this%loading)
        if (allocated(this%libs))           deallocate(this%libs)
        
        if (allocated(this%mask_core))      deallocate(this%mask_core)
        if (allocated(this%mask_axial))     deallocate(this%mask_axial)
        if (allocated(this%mask_radial))    deallocate(this%mask_radial)
    
    end subroutine Free_Material
    
    !$
    !===============================================================================================
    ! get material name by ID
    !===============================================================================================
    function Get_material_name (this, ID)  result(name)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: ID
        character(len=MAX_WORD_LEN)  :: name
        
        integer  :: i
        
        do i = 1, SIZE(this%libs)
            if (this%libs(i)%ID == ID)  then
                name = this%libs(i)%name
                return
            end if
        end do
        
        ! when can not find result
        call a_warning%set (INFO_LIST_XSEC, 'can not find the target material name')
        call a_warning%print (OUTPUT_UNIT)
        
    end function Get_material_name
    
    !$
    !===============================================================================================
    ! get material file by ID
    !===============================================================================================
    function Get_material_file (this, ID)  result(file)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: ID
        character(len=MAX_WORD_LEN)  :: file
        
        integer  :: i
        
        do i = 1, SIZE(this%libs)
            if (this%libs(i)%ID == ID)  then
                file = this%libs(i)%file
                return
            end if
        end do
        
        ! when can not find result
        call a_warning%set (INFO_LIST_XSEC, 'can not find the target material file')
        call a_warning%print (OUTPUT_UNIT)
        
    end function Get_material_file
    
    !$
    !===============================================================================================
    ! get material ID by name
    !===============================================================================================
    function Get_material_ID (this, name)  result(ID)
        
        class(Material), intent(in out)  :: this
        character(len=*), intent(in)  :: name
        integer  :: ID
        
        integer  :: i
        
        do i = 1, SIZE(this%libs)
            if (TRIM(ADJUSTL(this%libs(i)%name)) == TRIM(ADJUSTL(name)))  then
                ID = this%libs(i)%ID
                return
            end if
        end do
        
        ! when can not find result
        call a_warning%set (INFO_LIST_XSEC, 'can not find the target material ID')
        call a_warning%print (OUTPUT_UNIT)
        
    end function Get_material_ID
    
    !$
    !===============================================================================================
    ! material is fuel type ?
    !===============================================================================================
    function Is_materil_fuel (this, ID)  result(is_true)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: ID
        logical  :: is_true
        
        is_true = this%libs(ID)%is_fuel
    
    end function Is_materil_fuel 
    
    !$
    !===============================================================================================
    ! material is Control Rod type ?
    !===============================================================================================
    function Is_material_CR (this, ID)  result(is_true)
        
        class(Material), intent(in out)  :: this
        integer, intent(in)  :: ID
        logical  :: is_true
        
        is_true = this%libs(ID)%is_CR
    
    end function Is_material_CR
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
        
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_external_source_info_tp (this)
        
        class(external_source_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%intensity(ns%state%ng), stat=i_allocate)
        
        this%intensity = REAL_ZERO

        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%intensity)
        
    end subroutine Allocate_external_source_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of external_source_info_tp
    !===============================================================================================
    subroutine Free_external_source_info_tp (this)
        
        class(external_source_info_tp), intent(in out)  :: this
        
        if (allocated(this%intensity))          deallocate(this%intensity)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_external_source_info_tp
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_ExternalSource (this)
        
        class(ExternalSource), intent(in out)  :: this
        
        integer  :: i_allocate
        integer  :: i, j, k
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%adding(ns%state%zone, ns%state%layer), stat=i_allocate)
        
        allocate(this%kinds(ns%state%source), stat=i_allocate)
        do i = 1, SIZE(this%kinds)
            call this%kinds(i)%alloc ()
        end do
        
        allocate(this%matrixs(ns%state%zone, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%matrixs, dim=1)
            do j = 1, SIZE(this%matrixs, dim=2)
                call this%matrixs(i, j)%alloc () 
            end do
        end do
        
        ! NOTE, per nodal, for the flux is per nodal
        allocate(this%iter_scalar(ns%state%nodal, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%iter_scalar, dim=1)
            do j = 1, SIZE(this%iter_scalar, dim=2)
                call this%iter_scalar(i, j)%alloc ()
            end do
        end do
    
        this%adding = INT_ZERO
        this%is_adjust = .FALSE.
        this%iter_normal = 1.0
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%adding)
        if (SIZE(this%kinds) > 0)  then
            this%memory = this%memory + SIZE(this%kinds) * this%kinds(1)%memory
        end if
        if (SIZE(this%kinds) > 0)  then
            this%memory = this%memory + SIZE(this%matrixs) * this%matrixs(1, 1)%memory
        end if
        if (SIZE(this%kinds) > 0)  then
            this%memory = this%memory + SIZE(this%iter_scalar) * this%iter_scalar(1, 1)%memory
        end if
    
    end subroutine Allocate_ExternalSource
    
    !$
    !===============================================================================================
    ! set  
    !===============================================================================================
    subroutine Set_ExternalSource_by_row (this, assign, config)
        
        class(ExternalSource), intent(in out)  :: this
        integer, intent(in)  :: assign(:)
        integer, intent(in)  :: config(:, :)
        
        integer  :: ia, i

        do ia = 1, ns%state%layer
            this%adding(:, ia) = config(:, assign(ia))
        end do
    
    end subroutine Set_ExternalSource_by_row
    
    !$
    !===============================================================================================
    ! set  
    !===============================================================================================
    subroutine Set_ExternalSource_by_column (this, assign, config)
        
        class(ExternalSource), intent(in out)  :: this
        integer, intent(in)  :: assign(:)
        integer, intent(in)  :: config(:, :)
        
        integer  :: iz, i
        
        do iz = 1, ns%state%zone
            this%adding(iz, :) = config(assign(iz), :)
        end do
    
    end subroutine Set_ExternalSource_by_column
    
    !$
    !===============================================================================================
    ! finalizer for class of ExternalSource
    !===============================================================================================
    subroutine Free_ExternalSource (this)
        
        class(ExternalSource), intent(in out)  :: this
        integer  :: i, j, k
        
        if (allocated(this%adding))             deallocate(this%adding)
        
        if (allocated(this%kinds))  then
            do i = 1, SIZE(this%kinds)
                call this%kinds(i)%clean ()
            end do
            ! free itself
            deallocate(this%kinds)
        end if
        
        if (allocated(this%matrixs))  then
            do i = 1, SIZE(this%matrixs, dim=1)
                do j = 1, SIZE(this%matrixs, dim=2)
                    call this%matrixs(i, j)%clean 
                end do
            end do
            ! free itself
            deallocate(this%matrixs)
        end if
    
        if (allocated(this%iter_scalar))  then
            do i = 1, SIZE(this%iter_scalar, dim=1)
                do j = 1, SIZE(this%iter_scalar, dim=2)
                    call this%iter_scalar(i, j)%clean 
                end do
            end do
            ! free itself
            deallocate(this%iter_scalar)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_ExternalSource
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! adjuslt all the source parameter by specific level
    !===============================================================================================
    subroutine Adjust_ExternalSource (this, factor)
    
        class(ExternalSource), intent(in out)  :: this
        real(KREAL), intent(in)  :: factor
        
        integer  :: i, j, k
        
        ! adjust input parameter
        do i = 1, SIZE(this%kinds, dim=1)
            this%kinds(i)%intensity = this%kinds(i)%intensity * factor
        end do
        
        ! adjust real used parameter
        do i = 1, SIZE(this%matrixs, dim=1)
            do j = 1, SIZE(this%matrixs, dim=2)
                this%matrixs(i,j)%intensity = this%matrixs(i,j)%intensity * factor
            end do
        end do
        
        this%is_adjust = .TRUE.
    
    end subroutine Adjust_ExternalSource
    
    !$
    !===============================================================================================
    ! normalize external source to '1'
    !===============================================================================================
    subroutine Normal_ExternalSource (this, geom, is_normal)
        
        class(ExternalSource), intent(in out)  :: this
        type(Geometry), intent(in)  :: geom
        logical, intent(in)         :: is_normal
        
        real(KREAL)  :: total, tmp
        integer  :: ir, ia, iz, ig
        
        ! get normalize factor
        total = 0.0
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ig = 1, ns%state%ng
                    total = total + this%iter_scalar(ir, ia)%intensity(ig)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do
        
        ! @iter_normal
        if (is_normal )  then
            this%iter_normal = total
            tmp = total
        else 
            this%iter_normal = 1.0
            tmp = 1.0
        end if
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                this%iter_scalar(ir, ia)%intensity = this%iter_scalar(ir, ia)%intensity / tmp
            end do
        end do
    
    end subroutine Normal_ExternalSource
    
    !$
    !===============================================================================================
    ! get strength for source really use in iteration
    !===============================================================================================
    subroutine Get_ExternalSource_strength (this, geom, strength)
        
        class(ExternalSource), intent(in out)  :: this
        type(Geometry), intent(in)  :: geom
        real(KREAL), intent(in out)  :: strength
    
        integer  :: ir, ia, ig
        
        strength = 0.0
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                do ig = 1, ns%state%ng
                    strength = strength + this%iter_scalar(ir,ia)%intensity(ig) * geom%area(ir) * geom%height(ia)
                end do
            end do
        end do
        
    end subroutine Get_ExternalSource_strength
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_ExternalSource (this, unit_)
    
        class(ExternalSource), intent(in)  :: this
        integer, intent(in)  :: unit_
        integer  :: ia, iz, ig, ic, ir, is
        
        write(unit=unit_, fmt=*)  ' '
        write(unit=unit_, fmt=*)  '==========================================================='
        write(unit=unit_, fmt=*)  'input status:'
        do ia = 1, SIZE(this%kinds)
            write(unit=unit_, fmt="(1x, *(ES20.13, TR3))") this%kinds(ia)%intensity(:)
        end do
        
        write(unit=unit_, fmt=*)  ' '
        write(unit=unit_, fmt=*)  '==========================================================='
        write(unit=unit_, fmt=*)  'current status:'
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                write(unit=unit_, fmt="(1x, *(ES20.13, TR3))") this%matrixs(iz,ia)%intensity(:)
            end do
        end do
        
        write(unit=unit_, fmt=*)  ' '
        write(unit=unit_, fmt=*)  '==========================================================='
        write(unit=unit_, fmt=*)  'iteration status:'
        write(unit=unit_, fmt=*)  'is adjusted ?', this%is_adjust
        write(unit=unit_, fmt=*)  'normal value:', this%iter_normal
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt="(1x, *(ES13.6, TR3))") this%iter_scalar(ir,ia)%intensity(:)
            end do
        end do
    
    end subroutine Print_ExternalSource
    
end module material_header
