!$
!===================================================================================================
!
!   class for parameter distribution and transient time state container
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          TimeListParameter
!                               DistributionParameter
!                               GroupsFlux
!
!===================================================================================================
module contain_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use stastics,                   only : stastics_max_value, stastics_average_value, stastics_max_id
    use geometry_header,            only : Geometry, Meshing
    use quadrature_header,          only : QuadratureSet
    use feedback_header,            only : FeedbackParameter
    
    implicit none 
    private
    public  :: TimeListParameter, DistributionParameter, GroupsFlux
    
    ! --------------------------------------------------------------------------
    ! type for containing integral parameters per time step
    type  TimeListParameter
        real(KREAL), public               :: power                          ! power level of current time point (Watt)
        real(KREAL), public, allocatable  :: precursor(:)                   ! integrated precursor 

        real(KREAL), public               :: normal_factor                  ! normal factor for power level
        real(KREAL), public               :: initial_power                  ! save for initial power 
        
        real(KREAL), public               :: reactivity = REAL_ONE          ! dynamic reactivity 
        real(KREAL), public               :: beta       = REAL_ONE          ! effect delayed neutron fraction
        real(KREAL), public               :: fxy_nodal  = REAL_ONE
        integer, public                   :: idx_nodal  = INT_ONE
        real(KREAL), public               :: fxy_FA     = REAL_ONE
        integer, public                   :: idx_FA     = INT_ONE
        real(KREAL), public               :: fz         = REAL_ONE
        integer, public                   :: idx_z      = INT_ONE
        real(KREAL), public               :: fq         = REAL_ONE
        real(KREAL), public               :: Tm_max     = CKELVIN
        real(KREAL), public               :: Tm_outlet  = CKELVIN
        real(KREAL), public               :: Tm_avg     = CKELVIN
        real(KREAL), public               :: Tf_max     = CKELVIN
        real(KREAL), public               :: Tf_avg     = CKELVIN
    contains
        procedure, public  :: alloc => Allocate_TimeList
        procedure, public  :: clean => Free_TimeList
        procedure, public  :: set => Set_Timelist
    end type TimeListParameter
    
    ! type for containing distribution parameter per time step
    type  DistributionParameter
        real(KREAL), public, allocatable  :: matrix(:, :)
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public :: alloc => Allocate_Distribution
        procedure, public :: clean => Free_Distribution
        procedure, public :: set => Set_Distribution
        procedure, public :: nodal => Convert_to_nodal
        procedure, public :: zone_layer => Convert_to_zone_layer
        procedure, public :: zone => Convert_to_zone
        procedure, public :: layer => Convert_to_layer
        procedure, public :: fq_zone => Get_zone_fq_factor
        procedure, public :: fq_zone_layer => Get_zone_layer_fq_factor
        procedure, public :: factor => Get_distribution_factor
    end type DistributionParameter
    
    ! --------------------------------------------------------------------------
    ! type for flux distribution per group 
    type, private  :: flux_distribution_tp
        real(KREAL), public, allocatable  :: scalar(:, :)
        real(KREAL), public, allocatable  :: angular(:, :, :)
    contains
        procedure, public  :: alloc => Allocate_flux_distribution_tp
        procedure, public  :: clean => Free_flux_distribution_tp
    end type flux_distribution_tp
    
    ! type for flux distribution
    type  GroupsFlux 
        type(flux_distribution_tp), public, allocatable  :: ngs(:)
        logical, public          :: is_angular = .FALSE.
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_GroupsFlux
        procedure, public  :: set_scalar => Set_GroupsFlux_scalar
        procedure, public  :: sca2ang => Override_scalar2angular
        procedure, public  :: clean =>  Free_GroupsFlux
        procedure, private :: Equal_GroupsFlux
        generic,  public   :: assignment (=) => Equal_GroupsFlux
        procedure, public  :: fix => Fix_GroupsFlux_angular
    end type GroupsFlux
    
    ! private the real function name
    private  :: Allocate_TimeList, Free_TimeList, Set_Timelist
    private  :: Allocate_Distribution, Set_Distribution, Free_Distribution
    private  :: Convert_to_nodal, Convert_to_zone_layer, Convert_to_zone, Convert_to_layer
    private  :: Get_zone_fq_factor, Get_zone_layer_fq_factor
    
    private  :: Allocate_flux_distribution_tp, Free_flux_distribution_tp
    private  :: Allocate_GroupsFlux, Free_GroupsFlux, Set_GroupsFlux_scalar, Equal_GroupsFlux, Fix_GroupsFlux_angular
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_TimeList (this)
    
        class(TimeListParameter), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%precursor(nt%state%dg), stat=i_allocate)
        
        this%precursor  = REAL_ZERO
    
    end subroutine Allocate_TimeList
    
    !$
    !===============================================================================================
    ! finalizer for class of TimeListParameter
    !===============================================================================================
    subroutine Free_TimeList (this)
        
        class(TimeListParameter), intent(in out)  :: this
        
        if (allocated(this%precursor))       deallocate(this%precursor)
        
    end subroutine Free_TimeList
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_Timelist (this, geom, mesh, dist_power, self_fdbk, reactivity, beta)
        
        class(TimeListParameter), intent(in out)  :: this
        type(Geometry), intent(in)                :: geom
        type(Meshing), intent(in)                 :: mesh
        type(DistributionParameter), intent(in)   :: dist_power
        type(FeedbackParameter), intent(in)       :: self_fdbk
        real(KREAL), intent(in)  :: reactivity
        real(KREAL), intent(in)  :: beta
        
        real(KREAL)  :: fxy_nodal, fxy_FA, fz
        integer      :: idx_nodal, idx_FA, idx_z
        
        this%reactivity = reactivity
        this%beta = beta
        
        call dist_power%factor (mesh, geom, fxy_nodal, fxy_FA, fz, idx_nodal, idx_FA, idx_z)
        this%fxy_nodal = fxy_nodal
        this%fxy_FA = fxy_FA
        this%fz = fz
        this%idx_nodal = idx_nodal
        this%idx_FA = idx_FA
        this%idx_z = idx_z
        
        if (ns%feedback%is_feedback)  then
            this%Tm_outlet = self_fdbk%Tm%out
            this%Tm_max = self_fdbk%Tm%max
            this%Tm_avg = self_fdbk%Tm%avg
            this%Tf_max = self_fdbk%Tf%max
            this%Tf_avg = self_fdbk%Tf%avg
        end if
    
    end subroutine Set_Timelist
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_Distribution (this)
        
        class(DistributionParameter), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%matrix(ns%state%nodal, ns%state%layer), stat=i_allocate)
        
        this%matrix  = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%matrix)
    
    end subroutine Allocate_Distribution
    
    !$
    !===============================================================================================
    ! finalizer for class of DistributionParameter
    !===============================================================================================
    subroutine Free_Distribution (this)
        
        class(DistributionParameter), intent(in out)  :: this
        
        if (allocated(this%matrix))             deallocate(this%matrix)
        
        this%memory = REAL_ZERO
        
    end subroutine Free_Distribution
    
    !$
    !===============================================================================================
    ! assignment matrix to this class
    !===============================================================================================
    subroutine Set_Distribution (this, matrix)
        
        class(DistributionParameter), intent(in out)  :: this
        real(KREAL), intent(in)  :: matrix(:, :)
        
        integer  :: i, j
        integer  :: ii, jj
        do i = lbound(matrix, dim=1), ubound(matrix, dim=1) 
            do j = lbound(matrix, dim=2), ubound(matrix, dim=2)
                ii = i - lbound(matrix, dim=1) + 1                              ! low bound always begin with 1 
                jj = j - lbound(matrix, dim=2) + 1
                this%matrix(ii, jj) = matrix(i, j)
            end do
        end do
        
    end subroutine Set_Distribution

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! convert(nodal, layer) to (zone, layer) by is_function
    ! is_function = .TRUE.  : two sub zone will summed to generate a new zone
    ! is_function = .FALSE. : two sub zone will [volume-averaged] to generate a new zone
    !                         this is suit for the following subroutine names  Convert_to_*
    !===============================================================================================
    subroutine Convert_to_zone_layer (this, mesh, geom, is_function, zone_layer)
        
        class(DistributionParameter), intent(in) :: this  
        type(Meshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        logical, intent(in)         :: is_function
        real(KREAL), intent(in out)   :: zone_layer(:, :)                   ! the target holder
        
        integer  :: ia, ir, iz, im, ig
        integer  :: i_allocate
        
        ! this is a distribution, here take area into consideration only, not volume
        if (.NOT. is_function)  then
            zone_layer = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    zone_layer(iz, ia) = zone_layer(iz, ia) + this%matrix(ir, ia)*geom%area(ir)
                end do
            end do
            do iz = 1, ns%state%zone
                zone_layer(iz, :) = zone_layer(iz, :) / geom%zone_area(iz)
            end do
            
        ! this is a function, add directely
        else
            zone_layer = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    zone_layer(iz, ia) = zone_layer(iz, ia) + this%matrix(ir, ia)
                end do
            end do
        end if
    
    end subroutine Convert_to_zone_layer
    
    !$
    !===============================================================================================
    ! convert(nodal, layer) to (nodal) by is_function
    !===============================================================================================
    subroutine Convert_to_nodal (this, mesh, geom, is_function, nodal)
    
        class(DistributionParameter), intent(in) :: this  
        type(Meshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        logical, intent(in)         :: is_function
        real(KREAL), intent(in out)   :: nodal(:)                           ! the target holder

        integer  :: ia, ir, iz, im, ig
        integer  :: i_allocate
        
        ! this is a distribution
        if (.NOT. is_function)  then
            nodal = 0.0
            do ia = 1, ns%state%layer
                nodal(:) = nodal(:) + this%matrix(:, ia) *geom%height(ia)
            end do
            nodal(:) = nodal(:) / SUM(geom%height)
        
        ! this is a function
        else
            nodal = 0.0
            do ia = 1, ns%state%layer
                nodal(:) = nodal(:) + this%matrix(:, ia)
            end do
        end if
        
    end subroutine Convert_to_nodal
    
    !$
    !===============================================================================================
    ! convert(nodal, layer) to (zone) by is_function
    !===============================================================================================
    subroutine Convert_to_zone (this, mesh, geom, is_function, zone)
        
        class(DistributionParameter), intent(in) :: this
        type(Meshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        logical, intent(in)  :: is_function
        real(KREAL), intent(in out)   :: zone(:)                            ! the target holder
        
        real(KREAL), allocatable  :: zone_layer(:, :)                       ! tmp holder
        integer  :: ia, ir, iz, im, ig
        integer  :: i_allocate
        
        allocate(zone_layer(ns%state%zone, ns%state%layer), stat=i_allocate)
        
        ! this is a distribution
        if (.NOT. is_function)  then
            zone_layer = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    zone_layer(iz, ia) = zone_layer(iz, ia) + this%matrix(ir, ia)*geom%area(ir)
                end do
            end do
            do iz = 1, ns%state%zone
                zone_layer(iz, :) = zone_layer(iz, :) / geom%zone_area(iz)
            end do
            
            zone = 0.0
            do ia = 1, ns%state%layer
                zone(:) = zone(:) + zone_layer(:, ia) *geom%height(ia)
            end do
            zone(:) = zone(:) / SUM(geom%height)
            
        ! this is a function, add directely
        else
            zone_layer = 0.0
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    zone_layer(iz, ia) = zone_layer(iz, ia) + this%matrix(ir, ia)
                end do
            end do
            
            zone = 0.0
            do ia = 1, ns%state%layer
                zone(:) = zone(:) + zone_layer(:, ia)
            end do
        end if
        
        deallocate(zone_layer)
        
    end subroutine Convert_to_zone
    
    !$
    !===============================================================================================
    ! convert(nodal, layer) to (layer) by is_function
    !===============================================================================================
    subroutine Convert_to_layer (this, mesh, geom, is_function, layer)
        
        class(DistributionParameter), intent(in) :: this
        type(Meshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        logical, intent(in)  :: is_function
        real(KREAL), intent(in out)   :: layer(:)                           ! the target holder
        
        integer  :: ia, ir, iz, im, ig 
        integer  :: i_allocate
        
        ! this is a distribution
        if (.NOT. is_function)  then
            layer = 0.0
            do ir = 1, ns%state%nodal
                layer(:) = layer(:) + this%matrix(ir, :) * geom%area(ir)
            end do
            layer(:) = layer(:) / SUM(geom%zone_area)
            
        ! this is a function
        else
            layer = 0.0
            do ir = 1, ns%state%layer
                layer(:) = layer(:) + this%matrix(ir, :)
            end do 
        end if
        
    end subroutine Convert_to_layer
    
    !$
    !===============================================================================================
    ! get peak factor per zone
    !===============================================================================================
    subroutine Get_zone_fq_factor (this, mesh, geom, fq_zone)
        
        class(DistributionParameter), intent(in) :: this
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        real(KREAL), intent(in out) :: fq_zone(:)
        
        ! local varibles
        type  :: zone_count_tp
            integer          :: count  = 0
            real(KREAL)  :: maxima = 0.0
            real(KREAL)  :: total  = 0.0
        end type zone_count_tp
        
        type(zone_count_tp)  :: zone(ns%state%zone)
        real(KREAL)      :: nodal(ns%state%nodal)
        
        integer  :: iz, ir
        
        call this%nodal (mesh, geom, .FALSE., nodal)
        do ir = 1, ns%state%nodal
            iz = mesh%zone(ir)
            
            zone(iz)%count = zone(iz)%count + 1
            zone(iz)%total = zone(iz)%total + nodal(ir)*geom%area(ir)
            
            if (nodal(ir) >= zone(iz)%maxima)  then
                zone(iz)%maxima = nodal(ir)
            end if
        end do
        
        do iz = 1, ns%state%zone
            if (zone(iz)%maxima < EPS_ZERO)  then
                fq_zone(iz) = 1.0
            else
!                fq_zone(iz) = zone(iz)%maxima / (zone(iz)%total / zone(iz)%count)
                fq_zone(iz) = zone(iz)%maxima / (zone(iz)%total / geom%zone_area(iz))
            end if
        end do
    
    end subroutine Get_zone_fq_factor
    
    !$
    !===============================================================================================
    ! get peak factor per zone per layer
    !===============================================================================================
    subroutine Get_zone_layer_fq_factor (this, mesh, geom, fq_zone_layer)
        
        class(DistributionParameter), intent(in)  :: this
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        real(KREAL), intent(in out) :: fq_zone_layer(:, :)
        
        ! local variables
        type  :: zone_count_tp
            integer          :: count  = 0
            real(KREAL)  :: maxima = 0.0
            real(KREAL)  :: total  = 0.0
        end type zone_count_tp
        
        type(zone_count_tp)  :: zone_layer(ns%state%zone, ns%state%layer)
        integer  :: ia, iz, ir
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                
                zone_layer(iz, ia)%count = zone_layer(iz, ia)%count + 1
                zone_layer(iz, ia)%total = zone_layer(iz, ia)%total + this%matrix(ir, ia)*geom%area(ir)
                
                if (this%matrix(ir, ia) >= zone_layer(iz, ia)%maxima)  then
                    zone_layer(iz, ia)%maxima = this%matrix(ir, ia)
                end if
            end do
        end do
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                if (zone_layer(iz, ia)%maxima < EPS_ZERO)  then
                    fq_zone_layer(iz, ia) = 1.0
                else
!                    fq_zone_layer(iz, ia) = zone_layer(iz, ia)%maxima / (zone_layer(iz, ia)%total / zone_layer(iz, ia)%count)
                    fq_zone_layer(iz, ia) = zone_layer(iz, ia)%maxima / (zone_layer(iz, ia)%total / geom%zone_area(iz))
                end if
            end do
        end do
    
    end subroutine Get_zone_layer_fq_factor
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_distribution_factor (this, mesh, geom, fxy_nodal, fxy_FA, fz, idx_nodal, idx_FA, idx_z)
        
        class(DistributionParameter), intent(in)  :: this
        type(Meshing), intent(in)                :: mesh
        type(Geometry), intent(in)               :: geom
        real(KREAL), intent(in out)              :: fxy_nodal
        real(KREAL), intent(in out)              :: fxy_FA
        real(KREAL), intent(in out)              :: fz
        integer, intent(in out)                  :: idx_nodal
        integer, intent(in out)                  :: idx_FA
        integer, intent(in out)                  :: idx_z
        
        real(KREAL), allocatable  :: nodal(:)
        real(KREAL), allocatable  :: zone(:)
        real(KREAL), allocatable  :: layer(:)
        real(KREAL)  :: Rmax, Ravg
        integer      :: na, nr, nz
        integer      :: ixy, iz, ia
        integer      :: i_allocate
        
        nr = SIZE(this%matrix, dim=1)
        na = SIZE(this%matrix, dim=2)
        nz = SIZE(geom%zone_area)
        
        ! fxy_nodal
        allocate(nodal(nr), stat=i_allocate)
        nodal = REAL_ZERO
        do ia = 1, na
            nodal(:) = nodal(:) + geom%height(ia) * this%matrix(:, ia)
        end do
        nodal = nodal / SUM(geom%height)
        
        call stastics_max_id (nodal, ixy)
        Rmax = stastics_max_value (nodal)
        Ravg = stastics_average_value (nodal, geom%area, is_zero=.FALSE.)
        fxy_nodal = Rmax / Ravg
        idx_nodal = ixy
        if (allocated(nodal))           deallocate(nodal)
        
        ! fxy_FA
        allocate(zone(nz), stat=i_allocate)
        zone = REAL_ZERO
        call this%zone(mesh, geom, .FALSE., zone)
        
        call stastics_max_id (zone, ixy)
        Rmax = stastics_max_value (zone)
        Ravg = stastics_average_value (zone, geom%zone_area, is_zero=.FALSE.)
        fxy_FA = Rmax / Ravg
        idx_FA = ixy
        if (allocated(zone))            deallocate(zone)
        
        ! fz
        allocate(layer(na), stat=i_allocate)
        layer = REAL_ZERO
        call this%layer(mesh, geom, .FALSE., layer)
        
        call stastics_max_id (layer, iz)
        Rmax = stastics_max_value (layer)
        Ravg = stastics_average_value (layer, geom%height, is_zero=.FALSE.)
        fz = Rmax / Ravg
        idx_z = iz
        if (allocated(layer))           deallocate(layer)
    
    end subroutine Get_distribution_factor
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_flux_distribution_tp (this, is_angular)
        
        class(flux_distribution_tp), intent(in out)  :: this
        logical, intent(in)                          :: is_angular
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%scalar(ns%state%nodal, ns%state%layer), stat=i_allocate)
        this%scalar = REAL_ZERO
        
        if (is_angular )  then
            allocate(this%angular(ns%state%nodal, ns%state%layer, ns%deduce%direction), stat=i_allocate)
            this%angular = REAL_ZERO
        end if
    
    end subroutine Allocate_flux_distribution_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of flux_distrubution_tp
    !===============================================================================================
    subroutine Free_flux_distribution_tp (this)
        
        class(flux_distribution_tp), intent(in out)  :: this
        
        if (allocated(this%scalar))             deallocate(this%scalar)
        if (allocated(this%angular))            deallocate(this%angular)
    
    end subroutine Free_flux_distribution_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_GroupsFlux (this, is_angular)
        
        class(GroupsFlux), intent(in out)  :: this
        logical, intent(in)                :: is_angular
        
        integer  :: i_allocate
        integer  :: i
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%ngs(ns%state%ng), stat=i_allocate)
        
        do i = 1, SIZE(this%ngs, dim=1)
            call this%ngs(i)%alloc (is_angular)
        end do
                
        this%is_angular = is_angular
        
        ! get memory size
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%ngs(1)%scalar)
        if (is_angular )  then
            this%memory = this%memory + REAL_BYTE * SIZE(this%ngs(1)%angular)
        end if
        
        this%memory = this%memory * SIZE(this%ngs)
    
    end subroutine Allocate_GroupsFlux
    
    !$
    !===============================================================================================
    ! finalizer for classs of GroupsFlux
    !===============================================================================================
    subroutine Free_GroupsFlux (this)
    
        class(GroupsFlux), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%ngs))  then
            do i = 1, SIZE(this%ngs, dim=1)
                call this%ngs(i)%clean ()
            end do
            ! free itself
            deallocate(this%ngs)
        end if
        
        this%memory = REAL_ZERO
    
    end subroutine Free_GroupsFlux
    
    !$
    !===============================================================================================
    ! set scalar flux by angular flux
    !===============================================================================================
    subroutine Set_GroupsFlux_scalar (this, quad)
        
        class(GroupsFlux), intent(in out)  :: this
        type(QuadratureSet), intent(in)    :: quad
        
        integer  :: ig, is, ia, ir
        
        if (.NOT. this%is_angular)  then
            return
        end if
        
        do ig = 1, ns%state%ng
            this%ngs(ig)%scalar = 0.0
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    do is = 1, ns%deduce%direction
                        this%ngs(ig)%scalar(ir, ia) = this%ngs(ig)%scalar(ir, ia) + this%ngs(ig)%angular(ir, ia, is) * quad%directions(is)%wmu
                    end do
                end do
            end do
        end do
        
    end subroutine Set_GroupsFlux_scalar
    
    !$
    !===============================================================================================
    ! override angular flux by scalar flux
    !===============================================================================================
    subroutine Override_scalar2angular (this)
        
        class(GroupsFlux), intent(in out)  :: this
        
        integer  :: ig, is, ia, ir 
        
        do ig = 1, ns%state%ng
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    this%ngs(ig)%angular(ir, ia, :) = this%ngs(ig)%scalar(ir, ia)
                end do 
            end do 
        end do 
    
    end subroutine Override_scalar2angular
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Equal_GroupsFlux (left, right)
        
        class(GroupsFlux), intent(in out)  :: left
        type(GroupsFlux), intent(in)       :: right
        
        integer  :: ig
        
        do ig = 1, SIZE(left%ngs, dim=1)
            left%ngs(ig)%scalar = right%ngs(ig)%scalar
        end do

        if (left%is_angular )  then
            do ig = 1, SIZE(left%ngs, dim=1)
                left%ngs(ig)%angular = right%ngs(ig)%angular
            end do
        end if
        
    end subroutine Equal_GroupsFlux

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Fix_GroupsFlux_angular (this)
        
        class(GroupsFlux), intent(in out)  :: this
        
        real(KREAL), allocatable  :: tmp(:, :, :)
        integer  :: up_start, up_end, down_start, down_end
        integer  :: octant
        integer  :: ig
        integer  :: i_allocate
        
        octant = SIZE(this%ngs(1)%angular, dim=3) / 8
        allocate(tmp(ns%state%nodal, ns%state%layer, octant), stat=i_allocate)
        
        if (.TRUE.)  then 
            do ig = 1, SIZE(this%ngs)
                up_start = 0*octant + 1
                up_end = 1*octant
                down_start = 4*octant + 1
                down_end = 5*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
                
                up_start = 1*octant + 1
                up_end = 2*octant
                down_start = 5*octant + 1
                down_end = 6*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
                
                up_start = 2*octant + 1
                up_end = 3*octant
                down_start = 6*octant + 1
                down_end = 7*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
            
                up_start = 3*octant + 1
                up_end = 4*octant
                down_start = 7*octant + 1
                down_end = 8*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
            end do
        end if  
        
        ! 1<->8, 2<->7, 3<->6, 4<->5 (re-index direction, to comparison)
        if (.FALSE.)  then 
            do ig = 1, SIZE(this%ngs)
                up_start = 0*octant + 1
                up_end = 1*octant
                down_start = 7*octant + 1
                down_end = 8*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
                
                up_start = 1*octant + 1
                up_end = 2*octant
                down_start = 6*octant + 1
                down_end = 7*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
                
                up_start = 2*octant + 1
                up_end = 3*octant
                down_start = 5*octant + 1
                down_end = 6*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
            
                up_start = 3*octant + 1
                up_end = 4*octant
                down_start = 4*octant + 1
                down_end = 5*octant
                tmp(:,:,1:octant) = this%ngs(ig)%angular(:,:,up_start:up_end)
                this%ngs(ig)%angular(:,:,up_start:up_end) = this%ngs(ig)%angular(:,:,down_start:down_end)
                this%ngs(ig)%angular(:,:,down_start:down_end) = tmp(:,:,1:octant)
            end do
        end if 
        
        if (allocated(tmp))         deallocate(tmp)
        
    end subroutine Fix_GroupsFlux_angular
    
end module contain_header
