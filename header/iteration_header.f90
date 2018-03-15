!$
!===================================================================================================
!
!   class for iteration paramters
!   ----------------------------------------------------------------------------
!   | (5)             5/4/3/2/1--top/down/v/u/x surface
!   | (2,2,3)         left and right, order 1 & 2, direction for XUV
!   | (8)             j1*2, j2*2, j3*2, upper and lower
!   ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          IterationSource
!                               IterationFlux
!                               IterationCounter
!                               AnisotropicScatterFlux
!                               IterationCriterion
!
!===================================================================================================
module iteration_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV

    use geometry_header,        only : Meshing, Boundary
    use material_header,        only : CrossSection
    use vector_operation
    
    implicit none 
    private
    public  :: IterationSource, IterationFlux, IterationCounter
    public  :: AnisotropicScatterFlux, IterationCriterion
    
    ! --------------------------------------------------------------------------
    ! type for source during iteration
    type, private  :: source_moments_info_tp
        real(KREAL), public, allocatable  :: total_moment(:, :, :)
        real(KREAL), public, allocatable  :: out_group_moment(:, :, :)
    contains
        procedure, public  :: alloc => Alloc_source_moments_info_tp
        procedure, public  :: clean => Free_source_moments_info_tp
    end type source_moments_info_tp
    
    ! type for fission source infor during iteration
    type, private  :: fission_source_info_tp
        real(KREAL), public, allocatable  :: moment(:, :)
        real(KREAL), public, allocatable  :: old(:)
    contains
        procedure, public  :: alloc => Alloc_fission_source_info_tp
        procedure, public  :: clean => Free_fission_source_info_tp
    end type fission_source_info_tp
    
    ! type for flux information during iteration
    type, private  :: flux_info_tp
        real(KREAL), public, allocatable  :: moment(:, :, :)
        real(KREAL), public, allocatable  :: moment_omp(:, :, :)
        real(KREAL), public, allocatable  :: old(:, :)
    contains
        procedure, public  :: alloc => Allocate_flux_info_tp
        procedure, public  :: clean => Free_flux_info_tp
    end type flux_info_tp
    
    ! type for flux distribution during iteration
    type, private  :: flux_distribution_info_tp
        real(KREAL), public, allocatable  :: nodal(:, :, :, :)
        real(KREAL), public, allocatable  :: surface(:, :, :, :)
        real(KREAL), public, allocatable  :: point(:, :, :, :)
        
        integer, public                       :: count
        real(KREAL), public, allocatable  :: rad_surf(:, :, :, :)           ! outer surface
        real(KREAL), public, allocatable  :: axi_surf(:, :, :, :)           ! top & bottom surface
    contains
        procedure, public  :: alloc => Alloc_flux_distribution_info_tp
        procedure, public  :: alloc_bound => Alloc_flux_distribution_info_tp_bound
        procedure, public  :: clean => Free_flux_distribution_info_tp
    end type flux_distribution_info_tp
        
    ! --------------------------------------------------------------------------
    ! type for IterationSource used to hold source information in iteration
    type IterationSource
        type(source_moments_info_tp), public  :: info                          ! source information, including moment
        type(fission_source_info_tp), public  :: fission                       ! fission source information
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Alloc_IterationSource
        procedure, public  :: transit => Transit_old_source_value
        procedure, public  :: clean =>  Free_IterationSource
    end type IterationSource
            
    ! type for IterationFlux used to hold flux information in iteration
    type IterationFlux
        type(flux_info_tp), public                :: info                       ! flux information, including moment
        type(flux_distribution_info_tp), public   :: dist                       ! iteration flux distribution per nodal
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Alloc_IterationFlux
        procedure, public  :: alloc_bound => Alloc_IterationFlux_bound
        procedure, public  :: transit => Transit_old_flux_value
        procedure, public  :: cycle_in => Cycle_in_iteration
        procedure, public  :: clean =>  Free_IterationFlux
    end type IterationFlux
    
    ! type for iteration counter
    type  IterationCounter
        integer, public          :: in                                          ! counter for inner iteration
        integer, public          :: out                                         ! counter for outer iteration 
        real(KREAL), public      :: eigenvalue                                  ! eigenvalue of power iteration
        real(KREAL), public      :: ks                                          ! ks for subcritical iteration
        real(KREAL), public      :: ksub                                        ! ksub for subcritical iteration
        real(KREAL), public      :: coeff_FSP     = 1.0D0                       ! use for Ks method
        real(KREAL), public      :: kcritical     = 1.0D0                       ! keff value after critical adjusted 
        
        integer, allocatable     :: ng_start(:)
        integer, allocatable     :: ng_end(:)
        integer, allocatable     :: ng_step(:)
    contains
        procedure, public  :: alloc => Allocate_IterationCounter
        procedure, public  :: clean => Free_IterationCounter
        procedure, public  :: set_upscatter => Set_upscatter_sequence
        procedure, public  :: set_downscatter => Set_upscatter_reverse
        procedure, private :: Equal_IterationCounter
        generic,   public  :: assignment (=) => Equal_IterationCounter
    end type IterationCounter
    
    ! --------------------------------------------------------------------------
    ! type for anisotropic scatter flux per group
    type, private  :: anisotropic_scatter_tp
        real(KREAL), public, allocatable  :: aniso_zero(:, :, :)
        real(KREAL), public, allocatable  :: aniso_cos(:, :, :, :)
        real(KREAL), public, allocatable  :: aniso_sin(:, :, :, :)
    contains
        procedure, public  :: alloc => Allocate_anisotropic_scatter_tp
        procedure, public  :: clean => Free_anisotropic_scatter_tp
    end type anisotropic_scatter_tp

    ! type for AnisotropicScatterFlux flux
    type  AnisotropicScatterFlux
        type(anisotropic_scatter_tp), public, allocatable  :: ngs(:)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_AnisotropicScatterFlux
        procedure, public  :: value_zero => Set_AnisotropicScatterFlux_value_to_zero
        procedure, public  :: moment_zero => Set_AnisotropicScatterFlux_moment_to_zero
        procedure, public  :: update => Update_AnisotropicScatterFlux
        procedure, public  :: weight => Weight_AnisotropicScatterFlux
        procedure, public  :: clean =>  Free_AnisotropicScatterFlux
    end type AnisotropicScatterFlux
    
    ! type for criteria of the iteration convergence
    type  IterationCriterion
        integer, public             :: error_type           = 4                 ! norm type selection of vector
        integer, public             :: max_inner            = 10                ! max iteration number of inner cycle
        integer, public             :: max_outer            = 500               ! max iteration number of outer cycle
        real(KREAL), public         :: error_inner_flux     = 5.0E-6            ! error for flux in inner iteration
        real(KREAL), public         :: error_outer_flux     = 1.0E-5            ! error for flux in outer iteration
        real(KREAL), public         :: error_fission_rate   = 1.0E-5            ! errir for fisssion rate in outer iteration
        real(KREAL), public         :: error_eigen          = 1.0E-5            ! error for eigenvalue
    contains
        procedure, public  :: check_inner => Check_inner_iteration
        procedure, public  :: check_outer => Check_outer_iteration
        procedure, private :: Equal_IterationCriterion
        generic,   public  :: assignment (=) => Equal_IterationCriterion
    end type IterationCriterion
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_source_moments_info_tp (this)
        
        class(source_moments_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%total_moment(0:8, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%out_group_moment(0:8, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        
        this%total_moment      = REAL_ZERO
        this%out_group_moment  = REAL_ZERO
    
    end subroutine Alloc_source_moments_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of source_moments_info_tp
    !===============================================================================================
    subroutine Free_source_moments_info_tp (this)
        
        class(source_moments_info_tp), intent(in out)  :: this
        
        if (allocated(this%total_moment))              deallocate(this%total_moment)
        if (allocated(this%out_group_moment))          deallocate(this%out_group_moment)
    
    end subroutine Free_source_moments_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_fission_source_info_tp (this)
        
        class(fission_source_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%moment(0:8, ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%old(ns%deduce%nodal_total), stat=i_allocate)
        
        this%moment  = REAL_ZERO
        this%old     = REAL_ZERO
    
    end subroutine Alloc_fission_source_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of fission_source_info_tp
    !===============================================================================================
    subroutine Free_fission_source_info_tp (this)
        
        class(fission_source_info_tp), intent(in out)  :: this
        
        if (allocated(this%moment))             deallocate(this%moment)
        if (allocated(this%old))                deallocate(this%old)
    
    end subroutine Free_fission_source_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_flux_info_tp(this)
        
        class(flux_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%moment(0:8, ns%deduce%nodal_total, ns%state%ng), stat=i_allocate)
        allocate(this%moment_omp(1:8, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%old(ns%deduce%nodal_total, ns%state%ng), stat=i_allocate)
        
        this%moment      = REAL_ZERO
        this%moment_omp  = REAL_ZERO
        this%old         = REAL_ZERO
    
    end subroutine Allocate_flux_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of flux_info_tp
    !===============================================================================================
    subroutine Free_flux_info_tp(this)
    
        class(flux_info_tp), intent(in out)  :: this
        
        if (allocated(this%moment))             deallocate(this%moment)
        if (allocated(this%moment_omp))         deallocate(this%moment_omp)
        if (allocated(this%old))                deallocate(this%old)
    
    end subroutine Free_flux_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_flux_distribution_info_tp (this)
        
        class(flux_distribution_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allcoated status first
        call this%clean ()
        
        allocate(this%nodal(ns%state%nodal, ns%state%layer, ns%deduce%direction, ns%state%ng), stat=i_allocate)
        allocate(this%surface(5, ns%state%nodal, ns%state%layer, ns%deduce%direction), stat=i_allocate)
        allocate(this%point(ns%state%point, ns%state%layer, ns%deduce%direction, ns%state%ng), stat=i_allocate)
        
        this%nodal    = REAL_ZERO
        this%surface  = REAL_ZERO
        this%point    = REAL_ZERO
    
    end subroutine Alloc_flux_distribution_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_flux_distribution_info_tp_bound (this, bound)
    
        class(flux_distribution_info_tp), intent(in out)  :: this
        type(Boundary), intent(in)  :: bound
        
        integer  :: i_allocate
        integer  :: ir, iy
        
        ! get count
        this%count = 0
        do ir = 1, ns%state%nodal
            do iy = 1, 3
                if (ABS(bound%nodal(iy,ir)-bound%INNER) >= EPS_ZERO)  then
                    this%count = this%count + 1
                end if
            end do
        end do
        
        ! check allocated status first
        if (allocated(this%rad_surf))       deallocate(this%rad_surf)
        if (allocated(this%axi_surf))       deallocate(this%axi_surf)
        
        allocate(this%rad_surf(this%count, ns%state%layer, ns%deduce%direction, ns%state%ng), stat=i_allocate)
        allocate(this%axi_surf(ns%state%nodal, 2, ns%deduce%direction, ns%state%ng), stat=i_allocate)
        
        this%rad_surf = REAL_ZERO
        this%axi_surf = REAL_ZERO
        
    end subroutine Alloc_flux_distribution_info_tp_bound

    !$
    !===============================================================================================
    ! finalizer for class of flux_distribution_info_tp
    !===============================================================================================
    subroutine Free_flux_distribution_info_tp (this)
        
        class(flux_distribution_info_tp), intent(in out)  :: this
        
        if (allocated(this%nodal))          deallocate(this%nodal)
        if (allocated(this%surface))        deallocate(this%surface)
        if (allocated(this%point))          deallocate(this%point)
        
        if (allocated(this%rad_surf))       deallocate(this%rad_surf)
        if (allocated(this%axi_surf))       deallocate(this%axi_surf)
    
    end subroutine Free_flux_distribution_info_tp
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_IterationSource (this)
        
        class(IterationSource), intent(in out)  :: this
        
        ! check allocated status first
        call this%clean ()
        
        call this%info%alloc ()
        call this%fission%alloc ()
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%info%total_moment)
        this%memory = this%memory + REAL_BYTE * SIZE(this%info%out_group_moment)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%fission%moment)
        this%memory = this%memory + REAL_BYTE * SIZE(this%fission%old)
    
    end subroutine Alloc_IterationSource
    
    !$
    !===============================================================================================
    ! finalizer for class of IterationSource
    !===============================================================================================
    subroutine Free_IterationSource (this)
    
        class(IterationSource), intent(in out)  :: this
        
        call this%info%clean ()
        call this%fission%clean ()
        
        this%memory = REAL_ZERO
        
    end subroutine Free_IterationSource
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_IterationFlux (this)
    
        class(IterationFlux), intent(in out)  ::  this
        
        ! check allocated status first
        call this%clean ()
        
        call this%info%alloc ()
        call this%dist%alloc ()
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%info%moment)
        this%memory = this%memory + REAL_BYTE * SIZE(this%info%moment_omp)
        this%memory = this%memory + REAL_BYTE * SIZE(this%info%old)

        this%memory = this%memory + REAL_BYTE * SIZE(this%dist%nodal)
        this%memory = this%memory + REAL_BYTE * SIZE(this%dist%surface)
        this%memory = this%memory + REAL_BYTE * SIZE(this%dist%point)
        
    end subroutine Alloc_IterationFlux
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_IterationFlux_bound (this, bound)
        
        class(IterationFlux), intent(in out)  ::  this
        type(Boundary), intent(in)  :: bound
        
        call this%dist%alloc_bound (bound)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%dist%rad_surf)
        this%memory = this%memory + REAL_BYTE * SIZE(this%dist%axi_surf)
    
    end subroutine Alloc_IterationFlux_bound
    
    !$
    !===============================================================================================
    ! finalizer for class of IterationFlux
    !===============================================================================================
    subroutine Free_IterationFlux (this)
        
        class(IterationFlux), intent(in out)  :: this
        
        call this%info%clean ()
        call this%dist%clean ()
        
        this%memory = REAL_ZERO
    
    end subroutine Free_IterationFlux

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Transit_old_source_value (this)
        
        class(IterationSource), intent(in out)  :: this
        
        integer  :: il
        
        do il = 1, ns%deduce%nodal_total
            this%fission%old(il) = this%fission%moment(0, il)
            this%fission%moment(0, il) = REAL_ZERO
        end do
    
    end subroutine Transit_old_source_value
    
    !$
    !===============================================================================================
    ! transit to old value, and reset new value to zero
    !===============================================================================================
    subroutine Transit_old_flux_value (this, ig)
    
        class(IterationFlux), intent(in out)  :: this
        integer, intent(in)  :: ig
        
        integer  :: il
        
        do il = 1, ns%deduce%nodal_total
            this%info%old(il,ig) = this%info%moment(0,il,ig)
            this%info%moment(0,il,ig) = REAL_ZERO
        end do
        
    end subroutine Transit_old_flux_value
    
    !$
    !===============================================================================================
    ! generate error of inner iteration
    !===============================================================================================
    subroutine Cycle_in_iteration (this, criteria, ig, error_inner)
    
        class(IterationFlux), intent(in out)  :: this
        type(IterationCriterion), intent(in)  :: criteria
        integer, intent(in)  :: ig
        real(KREAL), intent(in out)  :: error_inner                         ! error of flux for inner iteration
        
        real(KREAL)  :: old_flux(ns%deduce%nodal_total)
        real(KREAL)  :: new_flux(ns%deduce%nodal_total)
        
        error_inner = 0.0       
        old_flux = this%info%moment(0, :, ig)
        new_flux = this%info%old(:, ig)
        
        error_inner = get_vector_error (old_flux, new_flux, criteria%error_type)
    
    end subroutine Cycle_in_iteration
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_IterationCounter (this)
        
        class(IterationCounter), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%ng_start(1), stat=i_allocate)
        allocate(this%ng_end(1), stat=i_allocate)
        allocate(this%ng_step(1), stat=i_allocate)
        
        this%ng_start = 1
        this%ng_end   = ns%state%ng
        this%ng_step  = 1
    
    end subroutine Allocate_IterationCounter
    
    !$
    !===============================================================================================
    ! finalizer for class of IterationCounter
    !===============================================================================================
    subroutine Free_IterationCounter (this)
        
        class(IterationCounter), intent(in out)  :: this
    
        if (allocated(this%ng_start))               deallocate(this%ng_start)
        if (allocated(this%ng_end))                 deallocate(this%ng_end)
        if (allocated(this%ng_step))                deallocate(this%ng_step)
    
    end subroutine Free_IterationCounter
    
    !$
    !===============================================================================================
    ! set energy group iteration sequence for forward calculation
    !===============================================================================================
    subroutine Set_upscatter_sequence (this, xsec)
        
        class(IterationCounter), intent(in out)  :: this
        type(CrossSection), intent(in)  :: xsec
        
        ! local variables
        logical  :: is_upscatter
        integer  :: group_index, tmp_index
        integer  :: iz, ia, i_scat
        integer  :: i, j
        integer  :: i_allocate
        
        ! single group, do nothing
        if (ns%state%ng == 1)  then
            call this%alloc ()
            this%ng_start = 1
            this%ng_end   = ns%state%ng
            this%ng_step  = 1
            return
        end if
        
        is_upscatter = .FALSE.
        group_index = ns%state%ng
        tmp_index = ns%state%ng
        
        do i_scat = 1, ns%deduce%scat_xs
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    loop: do i = 2, ns%state%ng
                        do j = 1, i-1
                            if (xsec%matrixs(iz, ia)%sigma_s(i,j,i_scat) >= EPS_ZERO)  then
                                tmp_index = MAX(i, j)
                                is_upscatter = .TRUE.
                                exit loop
                            end if
                        end do
                    end do loop
                    
                    if (tmp_index <= group_index)  then
                        group_index = tmp_index
                    end if
                end do
            end do
        end do
        
        ! with no upscatter, do nothing
        if (.NOT. is_upscatter)  then
            call this%alloc ()
            this%ng_start = 1
            this%ng_end   = ns%state%ng
            this%ng_step  = 1
            return
        end if

        ! upscatter cycle is false, do nothing
        if (.NOT. ns%method%is_upscatter_cycle)  then
            call this%alloc ()
            this%ng_start = 1
            this%ng_end   = ns%state%ng
            this%ng_step  = 1
            return
        end if
        
        ! number of cycle is zero, do nothing
        if (ns%method%is_upscatter_cycle  .and. ns%method%n_upscatter_cycle == 0)  then
            call this%alloc ()
            this%ng_start = 1
            this%ng_end   = ns%state%ng
            this%ng_step  = 1
            return
        end if
        
        ! set group index
        call this%clean ()
        allocate(this%ng_start(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        allocate(this%ng_end(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        allocate(this%ng_step(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        
        do i = 1, SIZE(this%ng_start)
            if (i == 1)  then
                this%ng_start(i) = 1
                this%ng_end(i)   = ns%state%ng
                this%ng_step(i)  = 1
            else
                this%ng_start(i) = group_index
                this%ng_end(i)   = ns%state%ng
                this%ng_step(i)  = 1
            end if
        end do
        
    end subroutine Set_upscatter_sequence
    
    !$
    !===============================================================================================
    ! set energy group iteration sequence for adjoint calculation
    !===============================================================================================
    subroutine Set_upscatter_reverse (this, xsec)
        
        class(IterationCounter), intent(in out)  :: this
        type(CrossSection), intent(in)  :: xsec
        
        ! local variables
        logical  :: is_upscatter
        integer  :: group_index, tmp_index
        integer  :: iz, ia, i_scat
        integer  :: i, j
        integer  :: i_allocate
        
        ! single group, do nothing
        if (ns%state%ng == 1)  then
            call this%alloc ()
            this%ng_start = ns%state%ng
            this%ng_end   = 1
            this%ng_step  = -1
            return
        end if
        
        is_upscatter = .FALSE.
        group_index = 1
        tmp_index = 1
        
        do i_scat = 1, ns%deduce%scat_xs
            do ia = 1, ns%state%layer
                do iz = 1, ns%state%zone
                    loop: do i = ns%state%ng-1, 1, -1
                        do j = ns%state%ng, i+1, -1
                            if (xsec%matrixs(iz, ia)%sigma_s(i,j,i_scat) >= EPS_ZERO)  then
                                tmp_index = MAX(i, j)
                                is_upscatter = .TRUE.
                                exit loop
                            end if
                        end do
                    end do loop
                    
                    if (tmp_index >= group_index)  then
                        group_index = tmp_index
                    end if
                end do
            end do
        end do
        
        ! with no upscatter, do nothing
        if (.NOT. is_upscatter)  then
            call this%alloc ()
            this%ng_start = ns%state%ng
            this%ng_end   = 1
            this%ng_step  = -1
            return
        end if

        ! upscatter cycle is false, do nothing
        if (.NOT. ns%method%is_upscatter_cycle)  then
            call this%alloc ()
            this%ng_start = ns%state%ng
            this%ng_end   = 1
            this%ng_step  = -1
            return
        end if
        
        ! number of cycle is zero, do nothing
        if (ns%method%is_upscatter_cycle  .and. ns%method%n_upscatter_cycle == 0)  then
            call this%alloc ()
            this%ng_start = ns%state%ng
            this%ng_end   = 1
            this%ng_step  = -1
            return
        end if
        
        ! set group index
        call this%clean ()
        allocate(this%ng_start(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        allocate(this%ng_end(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        allocate(this%ng_step(ns%method%n_upscatter_cycle+1), stat=i_allocate)
        
        do i = 1, SIZE(this%ng_start)
            if (i == 1)  then
                this%ng_start(i) = ns%state%ng
                this%ng_end(i)   = 1
                this%ng_step(i)  = -1
            else
                this%ng_start(i) = group_index
                this%ng_end(i)   = 1
                this%ng_step(i)  = -1
            end if
        end do
    
    end subroutine Set_upscatter_reverse
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Equal_IterationCounter (left, right)
            
        class(IterationCounter), intent(in out)  :: left
        type(IterationCounter), intent(in)  :: right
        
        integer  :: i_allocate
        
        left%in          = right%in
        left%out         = right%out
        left%eigenvalue  = right%eigenvalue
        left%ks          = right%ks
        left%ksub        = right%ksub
        left%coeff_FSP   = right%coeff_FSP
        
        call left%clean ()
        
        allocate(left%ng_start(SIZE(right%ng_start)), stat=i_allocate)
        allocate(left%ng_end(SIZE(right%ng_end)), stat=i_allocate)
        allocate(left%ng_step(SIZE(right%ng_step)), stat=i_allocate)
        
        left%ng_start  = right%ng_start
        left%ng_end    = right%ng_end
        left%ng_step   = right%ng_step 
    
    end subroutine Equal_IterationCounter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_anisotropic_scatter_tp (this)
        
        class(anisotropic_scatter_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allcoated status first
        call this%clean ()
        
        allocate(this%aniso_zero(ns%state%scat_order, 0:8, ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%aniso_cos(ns%state%scat_order, ns%state%scat_order, 0:8, ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%aniso_sin(ns%state%scat_order, ns%state%scat_order, 0:8, ns%deduce%nodal_total), stat=i_allocate)
        
        this%aniso_zero  = REAL_ZERO
        this%aniso_cos   = REAL_ZERO
        this%aniso_sin   = REAL_ZERO
    
    end subroutine Allocate_anisotropic_scatter_tp

    !$
    !===============================================================================================
    ! finalizer for class of anisotropic_scatter_tp
    !===============================================================================================
    subroutine Free_anisotropic_scatter_tp (this)
        
        class(anisotropic_scatter_tp), intent(in out)  :: this
        
        if (allocated(this%aniso_zero))          deallocate(this%aniso_zero)
        if (allocated(this%aniso_cos))           deallocate(this%aniso_cos)
        if (allocated(this%aniso_sin))           deallocate(this%aniso_sin)
    
    end subroutine Free_anisotropic_scatter_tp

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_AnisotropicScatterFlux (this)
            
        class(AnisotropicScatterFlux), intent(in out)  :: this
        
        integer  :: i_allocate
        integer  :: i
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%ngs(ns%state%ng), stat=i_allocate)
        
        do i = 1, SIZE(this%ngs, dim=1)
            call this%ngs(i)%alloc ()
        end do
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%ngs(1)%aniso_zero)
        this%memory = this%memory + REAL_BYTE * SIZE(this%ngs(1)%aniso_cos)
        this%memory = this%memory + REAL_BYTE * SIZE(this%ngs(1)%aniso_sin)

        this%memory = this%memory * SIZE(this%ngs)
        
    end subroutine Allocate_AnisotropicScatterFlux
    
    !$
    !===============================================================================================
    ! finalizer for class of AnisotropicScatterFlux
    !===============================================================================================
    subroutine Free_AnisotropicScatterFlux (this)
        
        class(AnisotropicScatterFlux), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%ngs))  then
            do i = 1, SIZE(this%ngs, dim=1)
                call this%ngs(i)%clean ()
            end do
            ! free it self
            deallocate(this%ngs)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_AnisotropicScatterFlux
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Update_AnisotropicScatterFlux (this)
        
        class(AnisotropicScatterFlux), intent(in out)  :: this
    
    end subroutine Update_AnisotropicScatterFlux
    
    !$
    !===============================================================================================
    ! set value(zero order moment) to zero
    !===============================================================================================
    subroutine Set_AnisotropicScatterFlux_value_to_zero (this, ig)
        
        class(AnisotropicScatterFlux), intent(in out)  :: this
        integer, intent(in)  :: ig
        integer  :: i
        
        ! zero order means real value
        this%ngs(ig)%aniso_zero(:, 0, :)     = REAL_ZERO
        this%ngs(ig)%aniso_cos(:, :, 0, :)   = REAL_ZERO
        this%ngs(ig)%aniso_sin(:, :, 0, :)   = REAL_ZERO
    
    end subroutine Set_AnisotropicScatterFlux_value_to_zero
    
    !$
    !===============================================================================================
    ! set moment per order to zero
    !===============================================================================================
    subroutine Set_AnisotropicScatterFlux_moment_to_zero (this, ig)
        
        class(AnisotropicScatterFlux), intent(in out)  :: this
        integer, intent(in)  :: ig
        integer  :: i
        
        ! order 1 to 8 means moments
        this%ngs(ig)%aniso_zero(:, 1:8, :)     = REAL_ZERO
        this%ngs(ig)%aniso_cos(:, :, 1:8, :)   = REAL_ZERO
        this%ngs(ig)%aniso_sin(:, :, 1:8, :)   = REAL_ZERO
    
    end subroutine Set_AnisotropicScatterFlux_moment_to_zero
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Weight_AnisotropicScatterFlux (this)
    
        class(AnisotropicScatterFlux), intent(in out)  :: this
    
    end subroutine Weight_AnisotropicScatterFlux
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_inner_iteration (this, count_in, error_inner, is_passing)
        
        class(IterationCriterion), intent(in out)  :: this
        integer, intent(in)          :: count_in
        real(KREAL), intent(in)  :: error_inner
        logical, intent(out)         :: is_passing
        
        is_passing = .TRUE.
        
        if (error_inner>this%error_inner_flux .AND. count_in<this%max_inner)  then
            is_passing = .FALSE.
        end if
    
    end subroutine Check_inner_iteration
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_outer_iteration (this, count_out, error_eigen, error_flux, error_fission_rate, is_passing)
        
        class(IterationCriterion), intent(in out)  :: this
        integer, intent(in)          :: count_out
        real(KREAL), intent(in)  :: error_eigen
        real(KREAL), intent(in)  :: error_flux
        real(KREAL), intent(in)  :: error_fission_rate
        logical, intent(out)  :: is_passing
            
        is_passing = .TRUE.
        
        if ((error_eigen>this%error_eigen.or.error_flux>this%error_outer_flux) .and. count_out<this%max_outer) then
            is_passing = .FALSE.
        end if
        
    end subroutine Check_outer_iteration
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Equal_IterationCriterion (left, right)
        
        class(IterationCriterion), intent(in out)  :: left
        type(IterationCriterion), intent(in)       :: right
        
        left%error_type          = right%error_type          
        left%max_inner           = right%max_inner           
        left%max_outer           = right%max_outer           
                                   
        left%error_inner_flux    = right%error_inner_flux    
        left%error_outer_flux    = right%error_outer_flux    
        left%error_fission_rate  = right%error_fission_rate  
        left%error_eigen         = right%error_eigen         
    
    end subroutine Equal_IterationCriterion
    
end module iteration_header
