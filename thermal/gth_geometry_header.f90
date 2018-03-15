!$
!===================================================================================================
!
!   module for geometry information for thermal calculation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          ThermalScale
!                               ThermalGeometry
!                               ThermalAssemblyGeometry
!
!===================================================================================================
module gth_geometry_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private 
    public  :: ThermalScale, ThermalGeometry, ThermalAssemblyGeometry
    
    ! --------------------------------------------------------------------------
    ! type for TH scale parameter
    type  ThermalScale 
        integer, public         :: nf                                           ! number of pellet sections in fuel conduction calculaton
        integer, public         :: nc                                           ! number of cladding sections in fuel conduction calculaton
        integer, public         :: n_mesh                                       ! total number of radial mesh point

        integer, public         :: nr                                           ! number of radial regions
        integer, public         :: na                                           ! number of axial meshs
        integer, public         :: na_start                                     ! active core start index
        integer, public         :: na_end                                       ! active core end index
        integer, public         :: n_assm_geom                                  ! number of assembly type by geometry
    contains
        procedure, public  :: mesh => Set_conduction_mesh
        procedure, public  :: set => Set_ThermalScale
    end type ThermalScale
    
    ! type for geometry information use in TH 
    type  ThermalGeometry
        real(KREAL), public, allocatable     :: height(:)                       ! (m) --axial mesh length of each mesh
        integer, public, allocatable         :: geom_type(:)                    ! geometry type per channel
        
        integer, public                      :: coolant_type       = 1          ! coolant type
        integer, public, allocatable         :: cladding_type(:)                ! cladding per geometry type
        integer, public, allocatable         :: gap_type(:)                     ! gap per geometry type
        integer, public, allocatable         :: fuel_type(:)                    ! fuel per geometry type
        real(KREAL), public, allocatable     :: fuel_TRU(:)                     ! weight fraction of TRU per geometry type
    contains
        procedure, public  :: alloc => Alloc_ThermalGeometry
        procedure, public  :: clean => Free_ThermalGeometry
        procedure, public  :: set_height => Set_ThermalGeometry_height
    end type ThermalGeometry

    ! type for assembly type by geometry
    type  ThermalAssemblyGeometry
        integer, public          :: n_pin                                       ! number of pin per assembly
        integer, public          :: n_fuelpin                                   ! number of fuel pin per assembly
        real(KREAL), public      :: pitch                                       ! (m)--pin pitch
        real(KREAL), public      :: rod                                         ! (m)--rod radius (contain clad)
        real(KREAL), public      :: cladth                                      ! (m)--clad thickness
        real(KREAL), public      :: bond                                        ! (m)--bond thickness
        real(KREAL), public      :: pellet                                      ! (m)--fuel pellet radius
        real(KREAL), public      :: hole                                        ! (m)--inner hole radius
                                 
        real(KREAL), public      :: area                                        ! (m^2)--fuel area
        real(KREAL), public      :: flowarea                                    ! (m^2)--flow area per pin
        logical, public          :: is_hexagonal       = .TRUE.                 ! is hexagonal lattice
        
        real(KREAL), public      :: P2D                                         ! ratio of pitch to diameter
        real(KREAL), public      :: Dh                                          ! equivalant diameter of coolant channel
        
        real(KREAL), public, allocatable     :: x_point(:)                      ! radial mesh for conductive calculation
        real(KREAL), public, allocatable     :: x_surface(:)                    ! radial surface for conductive calculation
        real(KREAL), public      :: df                                          ! length per fuel calculation section
        real(KREAL), public      :: dg                                          ! length of gap
        real(KREAL), public      :: dc                                          ! length of coolant calculation section
    contains
        procedure, public  :: alloc => Alloc_TypeGeometryAssembly
        procedure, public  :: clean => Free_TypeGeometryAssembly
        procedure, public  :: set => Set_TypeGeometryAssembly
    end type ThermalAssemblyGeometry
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Set_conduction_mesh, Set_ThermalScale
    private  :: Alloc_ThermalGeometry, Free_ThermalGeometry, Set_ThermalGeometry_height
    private  :: Alloc_TypeGeometryAssembly, Free_TypeGeometryAssembly, Set_TypeGeometryAssembly
    
    real(KREAL), parameter      :: CM2METER_ = 0.01D0                           ! cent-meter to meter
    
contains
    !$
    !===============================================================================================
    ! get radial mesh point
    !===============================================================================================
    subroutine Set_conduction_mesh (this, nf, nc)
        
        class(ThermalScale), intent(in out)  :: this
        integer, intent(in)  :: nf
        integer, intent(in)  :: nc
        
        this%nf = nf
        this%nc = nc
        
        this%n_mesh = this%nf + this%nc + 2
    
    end subroutine Set_conduction_mesh

    !$
    !===============================================================================================
    ! get axial calculation mesh point
    !===============================================================================================
    subroutine Set_ThermalScale (this, zone, layer, top_layer, bottom_layer)
        
        class(ThermalScale), intent(in out)  :: this
        integer, intent(in)  :: zone
        integer, intent(in)  :: layer
        integer, intent(in)  :: top_layer
        integer, intent(in)  :: bottom_layer
        
        this%nr = zone
        this%na = layer
        
        this%na_start = bottom_layer + 1
        this%na_end = this%na - top_layer
        
        this%n_mesh = this%nf + this%nc + 2
        
    end subroutine Set_ThermalScale

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_ThermalGeometry (this, nth)
        
        class(ThermalGeometry), intent(in out)  :: this
        type(ThermalScale), intent(in)          :: nth
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%height(nth%na), stat=i_allocate)
        allocate(this%geom_type(nth%nr), stat=i_allocate)

        allocate(this%cladding_type(nth%n_assm_geom), stat=i_allocate)
        allocate(this%gap_type(nth%n_assm_geom), stat=i_allocate)
        allocate(this%fuel_type(nth%n_assm_geom), stat=i_allocate)
        allocate(this%fuel_TRU(nth%n_assm_geom), stat=i_allocate)
        
        this%height     = 0.1D0
        this%geom_type  = 1
        
        this%coolant_type  = 1                                                  ! default is parcs-water
        this%cladding_type = 1
        this%gap_type      = 1
        this%fuel_type     = 1
        this%fuel_TRU      = 0.382D0
    
    end subroutine Alloc_ThermalGeometry
    
    !$
    !===============================================================================================
    ! finalizer for class of ThermalGeometry
    !===============================================================================================
    subroutine Free_ThermalGeometry (this)
        
        class(ThermalGeometry), intent(in out)  :: this
        
        if (allocated(this%height))             deallocate(this%height)
        if (allocated(this%geom_type))          deallocate(this%geom_type)
        
        if (allocated(this%cladding_type))      deallocate(this%cladding_type)
        if (allocated(this%gap_type))           deallocate(this%gap_type)
        if (allocated(this%fuel_type))          deallocate(this%fuel_type)
        if (allocated(this%fuel_TRU))           deallocate(this%fuel_TRU)
        
    end subroutine Free_ThermalGeometry
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_ThermalGeometry_height (this, height)
        
        class(ThermalGeometry), intent(in out)  :: this
        real(KREAL), intent(in)                 :: height(:)
        
        integer  :: ia
        
        this%height = height
        
        ! change unit from cent-meter(input) to meter(real used)
        this%height = this%height * CM2METER_
    
    end subroutine Set_ThermalGeometry_height
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_TypeGeometryAssembly (this, nth)
        
        class(ThermalAssemblyGeometry), intent(in out)  :: this
        type(ThermalScale), intent(in)                  :: nth
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%x_point(1: nth%n_mesh), stat=i_allocate)
        allocate(this%x_surface(0: nth%n_mesh), stat=i_allocate)
        
        this%x_point   = 0.0D0
        this%x_surface = 0.0D0
        
    end subroutine Alloc_TypeGeometryAssembly
    
    !$
    !===============================================================================================
    ! finalizer for class of ThermalAssemblyGeometry
    !===============================================================================================
    subroutine Free_TypeGeometryAssembly (this)
        
        class(ThermalAssemblyGeometry), intent(in out)  :: this
        
        if (allocated(this%x_point))             deallocate(this%x_point)
        if (allocated(this%x_surface))           deallocate(this%x_surface)
    
    end subroutine Free_TypeGeometryAssembly
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_TypeGeometryAssembly (this, nth)
        
        class(ThermalAssemblyGeometry), intent(in out)  :: this
        type(ThermalScale), intent(in)          :: nth
        
        integer          :: i
        
        ! flowarea
        if (this%is_hexagonal)  then
            this%flowarea = SQRT(3.0D0)*this%pitch**2/2.0D0 - PI*this%rod**2    ! flow area per rod
        else
            this%flowarea = this%pitch**2 - PI*this%rod**2
        end if
        
        ! fuel area
        this%area = PI*(this%rod - this%cladth - this%bond)**2 - PI*(this%hole)**2
        
        ! pitch to diameter
        this%P2D = this%pitch / (2.0D0*this%rod)
        
        ! equivalant diameter of coolant channel
        this%Dh = 4.0D0*this%flowarea / (2.0D0*PI*this%rod)
        
        ! fuel pellet
        this%pellet = this%rod - this%cladth - this%bond
        
        ! radial mesh
        this%df = (this%rod - this%cladth - this%bond - this%hole) / nth%nf
        this%dg = this%bond
        this%dc = this%cladth / nth%nc
        
        do i = 1, nth%n_mesh
            if (i <= nth%nf + 1)  then
                this%x_point(i) = this%hole+ (i-1) * this%df
            else if (i == nth%nf + 2)  then
                this%x_point(i) = this%x_point(i-1) + this%dg
            else
                this%x_point(i) = this%x_point(i-1) + this%dc
            end if
        end do
        
        ! radial surface
        do i = 0, nth%n_mesh
            if (i == 0)  then
                this%x_surface(i) = this%x_point(1)
            else if (i <= nth%nf)  then
                this%x_surface(i) = this%x_point(i) + 0.5D0 * this%df
            else if (i == nth%nf + 1)  then
                this%x_surface(i) = this%x_point(i) + 0.5D0 * this%dg
            else if (i <= nth%nf + nth%nc + 1)  then
                this%x_surface(i) = this%x_point(i) + 0.5D0 * this%dc
            else if (i == nth%nf + nth%nc + 2)  then
                this%x_surface(i) = this%x_point(i)
            end if
        end do
        
        ! change unit from cent-meter(input) to meter(real used)
        this%pitch = this%pitch * CM2METER_
        this%rod = this%rod * CM2METER_
        this%cladth = this%cladth * CM2METER_
        this%bond = this%bond * CM2METER_
        this%pellet = this%pellet * CM2METER_
        this%hole = this%hole * CM2METER_
        
        this%area = this%area * CM2METER_**2
        this%flowarea = this%flowarea * CM2METER_**2
        
        this%x_point = this%x_point * CM2METER_
        this%x_surface = this%x_surface * CM2METER_
        this%df = this%df * CM2METER_
        this%dg = this%dg * CM2METER_
        this%dc = this%dc * CM2METER_
        
        ! NOTE:
        this%Dh = this%Dh * CM2METER_
        
        ! ----------------------------------------------------------------------
        ! NOTE: pin axial mesh point
        ! 
        ! !<-----------pellet------------------>|<--gap-->|<-----clad---------->|
        ! .___________._______.______._________.__________.________.____________.
        ! 1           2       ...            nf+1       nf+2               nf+nc+2
        ! |<---df---->|                         |<--dg--->|        |<---dc---- >|
        !
        ! ----------------------------------------------------------------------
        
    end subroutine Set_TypeGeometryAssembly
    
end module gth_geometry_header
