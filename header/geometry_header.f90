!$
!===================================================================================================
!
!   class for geometry, meshing and boundary
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          Meshing
!                               Geometry
!                               Boundary
!                               VTKMeshing
!
!===================================================================================================
module geometry_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
        
    implicit none 
    private
    public  :: Meshing, Geometry, Boundary, VTKMeshing, coordinate_info_tp
    
    ! --------------------------------------------------------------------------
    ! type for meshing of geometry
    type  Meshing
        integer, public               :: meshtype = 0                           ! (0-ansys/ 1-rec/ 2-hex)
        character(len=MAX_WORD_LEN)   :: xyfile   = './input.coordinate'        ! 
        character(len=MAX_WORD_LEN)   :: iifile   = './input.mesh'              ! 
        integer, public               :: MESH_ANSYS = 0
        integer, public               :: MESH_REC   = 1
        integer, public               :: MESH_HEX   = 2
        
        integer, public, allocatable  :: zone(:)                                ! material zone
        integer, public, allocatable  :: point(:, :)                            ! global point ID per point

        integer, public, allocatable  :: nearby_nodal(:, :)                     ! adjacent nodal ID per edge
        integer, public, allocatable  :: localID(:, :)                          ! local ID per dege
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_Meshing
        procedure, public  :: set => Set_Meshing
        procedure, public  :: clean => Free_Meshing
    end type Meshing
        
    ! type for geometry scale information
    type  Geometry
        real(KREAL), public               :: dx0                                ! base dimension
        real(KREAL), public               :: dx1                                ! expansion dimension
        real(KREAL), public, allocatable  :: coordinate(:, :)                   ! 2D coordinate per nodal point
        real(KREAL), public, allocatable  :: height(:)                          ! height per axial layer
                                                                                    
        real(KREAL), public, allocatable  :: area(:)                            ! area per nodal
        real(KREAL), public, allocatable  :: zone_area(:)                       ! area per zone
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_Geometry
        procedure, public  :: get_area => Get_area_per_nodal
        procedure, public  :: get_zone_area => Get_area_per_zone
        procedure, public  :: clean =>  Free_Geometry
    end type Geometry
    
    ! type for geometry boundary condition
    type  Boundary
        real(KREAL), public, allocatable   :: radial(:)                         ! radial boundary condition per segment
        real(KREAL), public, allocatable   :: axial(:)                          ! axial boundary conditon: lower(1), uper(2)
        
        real(KREAL), public, allocatable   :: nodal(:, :)                       ! boundary condition of three side per nodal
        integer, public, allocatable   :: segmentID(:, :)                       ! boundary segment ID of three side
        
        real(KREAL), public          :: INNER   = -1.0D0
        real(KREAL), public          :: VACCUM  =  0.0D0
        real(KREAL), public          :: REFLECT =  1.0D0
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_Boundary
        procedure, public  :: set => Set_nodal_boundary
        procedure, public  :: inbc => Set_inner_boundary 
        procedure, public  :: clean => Free_Boundary
    end type Boundary
    
    ! --------------------------------------------------------------------------
    ! type for point coordinate
    type  :: coordinate_info_tp
        real(KREAL), public                             :: x                    ! first coordinate
        real(KREAL), public                             :: y                    ! second coordinate
        integer, public                                 :: index                ! index among all point
    contains
        procedure, public  :: is_same => Is_coordinate_same
        procedure, public  :: is_line => Is_coordinate_line
        generic,   public  :: assignment (=) => Equal_coordinate_info_tp
        procedure          :: Equal_coordinate_info_tp
    end type coordinate_info_tp
    
    ! type for mesh information of zone
    type, private  :: vtk_zone_info_tp
        integer, public                                 :: mesh_type       = 7  ! mesh type in vtk define (vtk-polygon)
        integer, public                                 :: nodal                ! number of nodal per zone
        integer, public                                 :: edge                 ! number of boundary edge
        integer, public                                 :: point                ! number of boundary point
        type(coordinate_info_tp), allocatable, public   :: coordinates(:)
    contains
        procedure, public  :: add => Add_point_info
        procedure, public  :: delete => Delete_point_info
        procedure, public  :: sort => Sort_point_info
        procedure, public  :: line => Check_line_info
    end type vtk_zone_info_tp
    
    ! type for vtk zone
    type  VTKMeshing
        integer, public                                 :: total_point
        integer, public                                 :: total_zone
        type(vtk_zone_info_tp), allocatable, public     :: zones(:)
        integer, allocatable, public                    :: mapping(:)           ! print mask, '0' means none
    contains
        procedure, public  :: alloc => Allocate_VTKMeshing
        procedure, public  :: clean => Free_VTKMeshing
        procedure, public  :: set => Set_VTKMeshing
    end type VTKMeshing
    
    ! --------------------------------------------------------------------------
    ! private the real funciton name 
    private  :: Allocate_Meshing, Set_Meshing, Free_Meshing
    private  :: Allocate_Geometry, Get_area_per_nodal, Get_area_per_zone, Free_Geometry
    private  :: Allocate_Boundary, Set_nodal_boundary, Free_Boundary
    
    private  :: Is_coordinate_same, Is_coordinate_line, Equal_coordinate_info_tp
    private  :: Add_point_info, Delete_point_info, Sort_point_info, Check_line_info
    private  :: Allocate_VTKMeshing, Free_VTKMeshing, Set_VTKMeshing
    
    real(KREAL), parameter  :: SAME_LIMIT = 1.0D-3                          ! judage coordinate relationship
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_Meshing (this)
        
        class(Meshing), intent(in out)  :: this
        
        integer  :: i_allocate
        integer, parameter    :: N_EDGE  = 3
        integer, parameter    :: N_POINT = 3
        
        ! check for allocated status first
        call this%clean ()
        
        ! allocate memory by global parameter
        allocate(this%zone(ns%state%nodal), stat=i_allocate)
        allocate(this%point(N_POINT, ns%state%nodal), stat=i_allocate)
        allocate(this%nearby_nodal(N_EDGE, ns%state%nodal), stat=i_allocate)
        allocate(this%localID(N_EDGE, ns%state%nodal), stat=i_allocate)
    
        ! initialize value
        this%zone         = INT_ZERO
        this%point        = INT_ZERO
        this%nearby_nodal = INT_ZERO
        this%localID      = INT_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + INT_BYTE * SIZE(this%zone)
        this%memory = this%memory + INT_BYTE * SIZE(this%point)
        this%memory = this%memory + INT_BYTE * SIZE(this%nearby_nodal)
        this%memory = this%memory + INT_BYTE * SIZE(this%localID)
        
    end subroutine Allocate_Meshing
    
    !$
    !===============================================================================================
    ! finalizer for class of Meshing
    !===============================================================================================
    subroutine Free_Meshing (this)
    
        class(Meshing), intent(in out)  :: this
        
        if (allocated(this%zone))               deallocate(this%zone)
        if (allocated(this%point))              deallocate(this%point)
        if (allocated(this%nearby_nodal))       deallocate(this%nearby_nodal)
        if (allocated(this%localID))            deallocate(this%localID)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_Meshing
    
    !$
    !===============================================================================================
    ! define nearby_nodal and nodal point's local ID
    !===============================================================================================
    subroutine Set_Meshing (this)
    
        class(Meshing), intent(in out)  :: this
        
        ! local variables
        integer  :: j, i, j4, j5, j6, i4
        integer  :: i1, i2, i3, j1, j2, j3
        
        ! nodal illustration:
        ! ----------------------------------------------------------------------
        !               i1
        !              /  \
        !            j3    j2
        !            /      \
        !          i2--j1---i3  
        !
        
        first_nodal: do i = 1, ns%state%nodal
            ! initialize encoding
            ! if i is the boundary nodal, i4 is the boundary edge
            do i4 = 1, 3
                this%localID(i4,i) = i4
                this%nearby_nodal(i4,i) = i
            end do
            i1 = this%point(1, i)
            i2 = this%point(2, i)
            i3 = this%point(3, i)
            
            ! search for the whole nodals to find the adjacent three nodal
            second_nodal: do j = 1 , ns%state%nodal
                if (j /= i) then
                    do j4 = 1, 3
                    
                        ! have one point is the same (i2)
                        if(this%point(j4,j) == i2) then
                            j2 = j4
                            do j5 = 1, 3
                                ! have another point is the same (i3)
                                if (this%point(j5,j) == i3) then
                                    j3 = j5
                                    ! point i2 and point i3 is the same, then the nodal is adjacent to line 1
                                    this%nearby_nodal(1, i) = j
                                    do j6 = 1, 3
                                        ! local encode j2 and j3 is used, only left j6
                                        if (j6/=j2 .and. j6/=j3)  this%localID(1,i) = j6
                                    end do
                                end if
                            end do
                        end if
                        
                        if (this%point(j4,j) == i3) then
                            j3 = j4
                            do j5 = 1, 3
                                if(this%point(j5,j) == i1) then
                                    j1 = j5
                                    this%nearby_nodal(2,i) = j
                                    do j6 = 1, 3
                                        if (j6/=j1 .and. j6/=j3)  this%localID(2,i) = j6
                                    end do
                                end if
                            end do
                        end if
                        
                        if (this%point(j4,j) == i1) then
                            j1 = j4
                            do j5 = 1, 3
                                if (this%point(j5,j) == i2) then
                                    j2 = j5
                                    this%nearby_nodal(3,i) = j
                                    do j6 = 1, 3
                                        if (j6/=j1 .and. j6/=j2)  this%localID(3,i) = j6
                                    end do
                                end if
                            end do
                        end if
                        
                    end do
                end if
            end do second_nodal
        end do first_nodal
                
    end subroutine Set_Meshing
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_Geometry (this)
        
        class(Geometry), intent(in out)  :: this
        
        integer  :: i_allocate
        integer, parameter   :: N_DIMENSION = 2
        
        ! check allcoated status first
        call this%clean ()
    
        ! allcoate memory for this type
        allocate(this%coordinate(N_DIMENSION, ns%state%point), stat=i_allocate)
        allocate(this%height(ns%state%layer), stat=i_allocate)
        allocate(this%area(ns%state%nodal), stat=i_allocate)
        allocate(this%zone_area(ns%state%zone), stat=i_allocate)
        
        ! initalize value
        this%dx0        = REAL_ONE
        this%dx1        = REAL_ONE
        
        this%coordinate = REAL_ZERO
        this%height     = REAL_ZERO
        this%area       = REAL_ZERO
        this%zone_area  = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%coordinate)
        this%memory = this%memory + REAL_BYTE * SIZE(this%height)
        this%memory = this%memory + REAL_BYTE * SIZE(this%area)
        this%memory = this%memory + REAL_BYTE * SIZE(this%zone_area)
    
    end subroutine Allocate_Geometry
    
    !$
    !===============================================================================================
    ! finalizer for class of Geometry
    !===============================================================================================
    subroutine Free_Geometry (this)
        
        class(Geometry), intent(in out)  :: this
        
        if (allocated(this%coordinate))    deallocate(this%coordinate)
        if (allocated(this%height))        deallocate(this%height)
        if (allocated(this%area))          deallocate(this%area)
        if (allocated(this%zone_area))     deallocate(this%zone_area)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_Geometry

    !$
    !===============================================================================================
    ! NOTE: mesh is also been changed if area if negative
    !===============================================================================================
    subroutine Get_area_per_nodal (this, mesh)
        
        class(Geometry), intent(in out)  :: this
        type(Meshing), intent(in out)    :: mesh
        
        integer  :: i
        integer  :: i_tmp

        ! re-arrange according to counter-clockwise, and get areas
        do i = 1, ns%state%nodal
            associate (x1 => this%coordinate(1, mesh%point(1, i)), y1 => this%coordinate(2, mesh%point(1, i)),     &
                       x2 => this%coordinate(1, mesh%point(2, i)), y2 => this%coordinate(2, mesh%point(2, i)),     &
                       x3 => this%coordinate(1, mesh%point(3, i)), y3 => this%coordinate(2, mesh%point(3, i)))
                this%area(i) = (x2*y3-x3*y2+x3*y1-x1*y3+x1*y2-x2*y1) / 2.0
            end associate
            
            if (this%area(i) < 0.0)  then
                i_tmp = mesh%point(2, i)
                mesh%point(2, i) = mesh%point(3, i)
                mesh%point(3, i) = i_tmp
                this%area(i) = -this%area(i)
            end if
        end do
    
    end subroutine Get_area_per_nodal
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Get_area_per_zone (this, mesh)
        
        class(Geometry), intent(in out)  :: this
        type(Meshing), intent(in out)    :: mesh
        integer  :: i, iz
        
        this%zone_area = 0.0
        do i = 1, ns%state%nodal
            iz = mesh%zone(i)
            this%zone_area(iz) = this%zone_area(iz) + this%area(i)
        end do
    
    end subroutine Get_area_per_zone
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_Boundary (this)
        
        class(Boundary), intent(in out)  :: this
        
        integer  :: i_allocate
        integer, parameter  :: UP_AND_DOWN = 2
        integer, parameter  :: N_EDGE      = 3
        
        ! check allocated status first
        call this%clean ()
        
        ! allocate memory
        allocate(this%radial(ns%state%segment), stat=i_allocate)
        allocate(this%axial(UP_AND_DOWN), stat=i_allocate)

        allocate(this%nodal(N_EDGE, ns%state%nodal), stat=i_allocate)
        allocate(this%segmentID(N_EDGE, ns%state%nodal), stat=i_allocate)
        
        ! initialize value
        this%radial     = this%VACCUM
        this%axial      = this%VACCUM
        this%nodal      = this%INNER
        this%segmentID  = INT_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + INT_BYTE * SIZE(this%radial)
        this%memory = this%memory + INT_BYTE * SIZE(this%axial)
        this%memory = this%memory + INT_BYTE * SIZE(this%nodal)
        this%memory = this%memory + INT_BYTE * SIZE(this%segmentID)
    
    end subroutine Allocate_Boundary
    
    !$
    !===============================================================================================
    ! finalizer for class of Boundary
    !===============================================================================================
    subroutine Free_Boundary (this)
        
        class(Boundary), intent(in out)  :: this
        
        if (allocated(this%radial))             deallocate(this%radial)
        if (allocated(this%axial))              deallocate(this%axial)
        if (allocated(this%nodal))              deallocate(this%nodal)
        if (allocated(this%segmentID))          deallocate(this%segmentID)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_Boundary
    
    !$
    !===============================================================================================
    ! define boundary condition
    !===============================================================================================
    subroutine Set_nodal_boundary (this, mesh, geom, nd, i)
        
        class(Boundary), intent(in out)  :: this
        type(Meshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        integer, intent(in)  :: nd(3)                                           ! point information per segment
        integer, intent(in)  :: i                                               ! segment index
        
        ! local variables
        integer  :: j, k
        integer  :: k2, k3
        real(KREAL) :: x_end, x_start, y_end, y_start
        real(KREAL) :: xl, yl, xr, yr, xm, ym
        logical     :: lefe_in, right_in, center_in
            
        if (nd(1) /= nd(2)) then
            do j = 1, ns%state%nodal
                do k = 1, 3
                    ! this nodal
                    if (mesh%nearby_nodal(k,j) == j)  then
                        k2 = (-2)**(k/3) + k
                        k3 = (-2)**(k2/3) + k2
                        
                        ! the boundary is made of nd(1) and nd(2)
                        if (mesh%point(k2,j)==nd(1) .AND. mesh%point(k3,j)==nd(2)) then
                            this%nodal(k,j) = this%radial(i)
                            this%segmentID(k,j)  = i
                        end if
                    end if
                end do
            end do
        end if
        
        ! NOTE: use self-made function to ruturn with the same kind
        x_end   = Self_max_real (geom%coordinate(1,nd(3)), geom%coordinate(1,nd(2)))
        x_start = Self_min_real (geom%coordinate(1,nd(3)), geom%coordinate(1,nd(2)))
        y_end   = Self_max_real (geom%coordinate(2,nd(3)), geom%coordinate(2,nd(2)))
        y_start = Self_min_real (geom%coordinate(2,nd(3)), geom%coordinate(2,nd(2)))
        
        do j = 1, ns%state%nodal
            do k = 1, 3
                ! boundary nodal
                if (mesh%nearby_nodal(k,j) == j) then
                    k2 = (-2)**(k/3) + k
                    k3 = (-2)**(k2/3) + k2
                    xl = geom%coordinate(1, mesh%point(k2,j))
                    yl = geom%coordinate(2, mesh%point(k2,j))
                    xr = geom%coordinate(1, mesh%point(k3,j))
                    yr = geom%coordinate(2, mesh%point(k3,j))
                    xm = (xl + xr) / 2.0
                    ym = (yl + yr) / 2.0
                    
                    lefe_in = Is_in_scope (xl, yl, x_start, x_end, y_start, y_end)
                    right_in = Is_in_scope (xr, yr, x_start, x_end, y_start, y_end)
                    center_in = Is_in_scope (xm, ym, x_start, x_end, y_start, y_end)
                    
                    if (lefe_in .AND. right_in .AND. center_in)  then 
                        this%nodal(k, j) = this%radial(i)
                        this%segmentID(k, j) = i
                    end if
                end if
            end do
        end do
        
    end subroutine Set_nodal_boundary
    
    !$
    !===============================================================================================
    ! define boundary condition
    !===============================================================================================
    subroutine Set_inner_boundary (this, mesh)
        
        class(Boundary), intent(in out)  :: this
        type(Meshing), intent(in)   :: mesh
        
        integer  :: i, j
        
        do i = 1, ns%state%nodal
            do j = 1, 3
                if (mesh%nearby_nodal(j,i) /= i)  then
                    this%nodal(j,i) = this%INNER
                    this%segmentID(j,i) = 0
                end if 
            end do 
        end do 
    
    end subroutine Set_inner_boundary
        
    !$
    !===============================================================================================
    ! return max one of two input with the same type
    !   the intrinsic max return default type always
    !===============================================================================================
    function Self_max_real (one, two)  result(output)
        
        real(KREAL), intent(in)  :: one
        real(KREAL), intent(in)  :: two
        real(KREAL)  :: output
        
        if (one >= two)  then
            output = one 
        else 
            output = two
        end if
    
    end function Self_max_real
        
    !$
    !===============================================================================================
    ! return min one of two input with the same type
    !===============================================================================================
    function Self_min_real (one, two)  result(output)
        
        real(KREAL), intent(in)  :: one
        real(KREAL), intent(in)  :: two
        real(KREAL)  :: output
        
        if (one <= two)  then
            output = one 
        else 
            output = two
        end if
    
    end function Self_min_real
    
    !$
    !===============================================================================================
    ! is in scope
    !===============================================================================================
    function Is_in_scope (xin, yin, x_start, x_end, y_start, y_end)  result(output)
        
        real(KREAL), intent(in)  :: xin
        real(KREAL), intent(in)  :: yin
        real(KREAL), intent(in)  :: x_start
        real(KREAL), intent(in)  :: x_end
        real(KREAL), intent(in)  :: y_start
        real(KREAL), intent(in)  :: y_end
        logical  :: output
        
        real(KREAL), parameter  :: EPS_ADDED = 1.0D-3                           ! NOTE: this is used for reading data error
        
        output = .FALSE.
        if ((x_start-EPS_ADDED<xin) .AND. (xin<x_end+EPS_ADDED) .AND. (y_start-EPS_ADDED<yin) .AND. (yin<y_end+EPS_ADDED))  then
            output = .TRUE.
        end if 
    
    end function Is_in_scope    
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! is the coordinate the same one ? by coordinate value
    !===============================================================================================
    function Is_coordinate_same (this, second)  result (is_true)
        
        class(coordinate_info_tp), intent(in)  :: this
        type(coordinate_info_tp), intent(in)   :: second
        logical  :: is_true
        
        is_true = .FALSE.
        if (ABS(this%x-second%x)<=SAME_LIMIT .and. ABS(this%y-second%y)<=SAME_LIMIT)  then
            is_true = .TRUE.
            return
        end if
    
    end function Is_coordinate_same
    
    !$
    !===============================================================================================
    ! is the coordinate along a line ? by coordinate value
    !===============================================================================================
    function Is_coordinate_line (this, second, third)  result (is_true)
        
        class(coordinate_info_tp), intent(in out)  :: this
        type(coordinate_info_tp), intent(in)       :: second
        type(coordinate_info_tp), intent(in)       :: third
        logical  :: is_true
        
        ! local
        real(KREAL)  :: tmp
        
        is_true = .FALSE.
        if (ABS(this%x-second%x)<=SAME_LIMIT .and. ABS(this%x-third%x)<=SAME_LIMIT)  then
            is_true = .TRUE.
            return
        end if
        
        if (ABS(this%y-second%y)<=SAME_LIMIT .and. ABS(this%y-third%y)<=SAME_LIMIT)  then
            is_true = .TRUE.
            return
        end if
        
        tmp = (third%y-this%y)/(third%x-this%x) - (second%y-this%y)/(second%x-second%x)
        if (ABS(tmp) <= SAME_LIMIT)  then
            is_true = .TRUE.
            return
        end if
    
    end function Is_coordinate_line 
    
    !$
    !===============================================================================================
    ! assignment implement for coordinate_info_tp
    !===============================================================================================
    subroutine Equal_coordinate_info_tp (left, right)
        
        class(coordinate_info_tp), intent(in out)  :: left
        type(coordinate_info_tp), intent(in)       :: right
        
        left%x = right%x
        left%y = right%y
        
        left%index = right%index
    
    end subroutine Equal_coordinate_info_tp

    !$
    !===============================================================================================
    ! add a point coordinate, if it is not added
    !===============================================================================================
    subroutine Add_point_info (this, x, y, index)
    
        class(vtk_zone_info_tp), intent(in out)  :: this
        real(KREAL), intent(in)  :: x
        real(KREAL), intent(in)  :: y
        integer, intent(in)          :: index
        
        ! local 
        type(coordinate_info_tp), allocatable  :: tmp(:)
        type(coordinate_info_tp)               :: point
        integer  :: count
        integer  :: i, i_allocate
        logical  :: is_added = .FALSE.
        
        point%x = x
        point%y = y
        point%index = index
        
        if (allocated(this%coordinates))  then
            ! is existed ?
            do i = 1, SIZE(this%coordinates)
                is_added = this%coordinates(i)%is_same (point)
                if (is_added)  then
                    return
                end if
            end do
        
            this%point = this%point + 1
            count = SIZE(this%coordinates)
            allocate(tmp(count), stat=i_allocate)
            tmp = this%coordinates
            
            deallocate(this%coordinates)
            allocate(this%coordinates(count+1), stat=i_allocate)
            this%coordinates(1:count) = tmp
            this%coordinates(count+1) = point
            deallocate(tmp)
            
        ! this is the first point
        else
            this%point = 1
            allocate(this%coordinates(1), stat=i_allocate)
            this%coordinates(1) = point
        end if
    
    end subroutine Add_point_info
    
    !$
    !===============================================================================================
    ! delete a point by index, left more than one
    !===============================================================================================
    subroutine Delete_point_info (this, index)
    
        class(vtk_zone_info_tp), intent(in out)  :: this
        integer, intent(in)                      :: index
    
        ! local
        type(coordinate_info_tp), allocatable    :: tmp(:)
        integer  :: count
        integer  :: i, i_allocate
        
        count = SIZE(this%coordinates)
        allocate(tmp(count-1), stat=i_allocate)
        
        if (index == 1)  then
            tmp(1:) = this%coordinates(2:count)
        else if (index == count)  then
            tmp(1:) = this%coordinates(1:count-1)
        else
            tmp(1:index-1) = this%coordinates(1:index-1)
            tmp(index:)    = this%coordinates(index+1:count)
        end if
        
        deallocate(this%coordinates)
        allocate(this%coordinates(count-1), stat=i_allocate)
        this%coordinates = tmp
        deallocate(tmp)
        
    end subroutine Delete_point_info

    !$
    !===============================================================================================
    ! sort the point anti-clockwise, starting at [-0.5*pi, 1.5*pi)
    !===============================================================================================
    subroutine Sort_point_info (this)
        
        class(vtk_zone_info_tp), intent(in out)  :: this
        
        ! local 
        type(coordinate_info_tp), allocatable  :: copy(:)
        real(KREAL), allocatable               :: angular(:)
        type(coordinate_info_tp)               :: tmp
        
        real(KREAL)  :: value
        real(KREAL)  :: center_x, center_y, x, y
        integer  :: i, j, i_allocate
        
        ! get the shape center, make sure that the angular within [-0.5*pi, 1.5*pi]
        center_x = 0.0; center_y = 0.0;
        do i = 1, this%point
            center_x = center_x + this%coordinates(i)%x
            center_y = center_y + this%coordinates(i)%y
        end do
        center_x = center_x / this%point
        center_y = center_y / this%point
        
        allocate(copy(this%point), stat=i_allocate)
        allocate(angular(this%point), stat=i_allocate)
        
        ! get angular info
        do i = 1, this%point
            x = this%coordinates(i)%x - center_x
            y = this%coordinates(i)%y - center_y
            value = y / SQRT(x*x + y*y)
            if (value >=  1.0D0)  value = 1.0D0
            if (value <= -1.0D0)  value = -1.0D0
            
            angular(i) = ASIN(value)
            if (x < 0.0)  then
                angular(i) = PI - ASIN(value)
            end if
            angular(i) = angular(i) / PI
        end do
        
        ! sort coordinate and angular at the same time
        copy = this%coordinates
        do i = 1, this%point-1
            do j = i+1, this%point
                if (angular(j) <= angular(i))  then
                    value = angular(i)
                    angular(i) = angular(j)
                    angular(j) = value
                    
                    tmp = copy(i)
                    copy(i) = copy(j)
                    copy(j) = tmp
                end if
            end do
        end do
        
        this%coordinates = copy
        
        if (allocated(copy))        deallocate(copy)
        if (allocated(angular))     deallocate(angular)
        
    end subroutine Sort_point_info
    
    !$
    !===============================================================================================
    ! delete the second point if three countinuous points along a line
    !===============================================================================================
    subroutine Check_line_info (this) 
        
        class(vtk_zone_info_tp), intent(in out)  :: this
        
        ! local
        type(coordinate_info_tp)  :: second
        type(coordinate_info_tp)  :: third
        integer  :: i1, i2, i3, i, j, point
        logical  :: is_line 
        
        point = this%point
        do i = 1, point
            if (i == 1)  then
                i1 = 0
            end if
            i1 = i1 + 1                                                         ! i1 increase
            i2 = i1 + 1                                                         ! i2 is the next point to i1
            i3 = i1 + 2
            
            if (i2 > this%point)  i2 = i2 - this%point                          ! back to starting point
            if (i3 > this%point)  i3 = i3 - this%point
            
            second = this%coordinates(i2)
            third = this%coordinates(i3)
            is_line = this%coordinates(i1)%is_line (second, third)
            
            if (is_line)  then
                call this%delete (i2)                                           ! delete the center point
                i1 = i1 - 1
                this%point = this%point - 1
                this%edge = this%edge - 1
            end if
        end do
    
    end subroutine Check_line_info

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_VTKMeshing (this)
        
        class(VTKMeshing), intent(in out)  :: this
        integer  :: i, i_allocate
        
        ! check for allocated first
        call this%clean ()
        
        allocate(this%zones(ns%state%zone), stat=i_allocate)
        allocate(this%mapping(ns%state%zone), stat=i_allocate)
        
        do i = 1, SIZE(this%zones)
            this%zones(i)%nodal = 0
            this%zones(i)%edge  = 0
            this%zones(i)%point = 0
        end do
        
        this%mapping = 1
        if (SIZE(this%mapping)>=2)  then
            this%mapping(1:2)= 1
        end if 
        this%total_zone  = 0
        this%total_point = 0
        
    end subroutine Allocate_VTKMeshing
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_VTKMeshing (this)
        
        class(VTKMeshing), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%mapping))    deallocate(this%mapping)
        
        if (allocated(this%zones))  then
            do i = 1, SIZE(this%zones)
                if (allocated(this%zones(i)%coordinates))  then
                    deallocate(this%zones(i)%coordinates)
                end if
            end do
            
            ! free itself
            deallocate(this%zones)
        end if
        
    end subroutine Free_VTKMeshing
    
    !$
    !===============================================================================================
    ! set zone information for vtk plotting without fine mesh
    !===============================================================================================
    subroutine Set_VTKMeshing (this, mesh, geom)
    
        class(VTKMeshing), intent(in out)   :: this
        type(Meshing), intent(in)           :: mesh
        type(Geometry), intent(in)          :: geom
    
        ! local
        real(KREAL)  :: x, y
        integer  :: index, j
        integer  :: ir, iz, i1, i2, i3, id
        
        do ir = 1, ns%state%nodal
            iz = mesh%zone(ir)
            this%zones(iz)%nodal = this%zones(iz)%nodal + 1
            
            do i1 = 1, 3
                i2 = i1 + (-2)**(i1/3)
                i3 = i2 + (-2)**(i2/3)
                id = mesh%nearby_nodal(i1, ir)
                
                ! 'i1' is edge of zone
                if ((id == ir) .or. (mesh%zone(id) /= iz)) then
                    this%zones(iz)%edge = this%zones(iz)%edge + 1
                    
                    index = mesh%point(i2, ir)
                    x = geom%coordinate(1, index)
                    y = geom%coordinate(2, index)
                    call this%zones(iz)%add (x, y, index)
                    
                    index = mesh%point(i3, ir)
                    x = geom%coordinate(1, index)
                    y = geom%coordinate(2, index)
                    call this%zones(iz)%add (x, y, index)
                end if
            end do
        end do
        
        ! sort
        this%total_zone = 0
        this%total_point = 0
        do iz = 1, ns%state%zone
            call this%zones(iz)%sort ()
            call this%zones(iz)%line ()
            
            if (this%mapping(iz) /= 0)  then
                this%total_zone = this%total_zone + 1
                this%total_point = this%total_point + this%zones(iz)%point
            end if
        end do
        
    end subroutine Set_VTKMeshing
        
end module geometry_header
