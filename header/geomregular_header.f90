!$
!===================================================================================================
!
!   class for regular geometry
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          GeomRec
!                               GeomHex
!
!===================================================================================================
module geomregular_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use string 
    
    use state_header,           only : SteadyState
    use geometry_header,        only : Meshing, Geometry, Boundary, VTKMeshing, coordinate_info_tp
    
    implicit none
    private
    public  :: GeomRec, GeomHex
    
    ! --------------------------------------------------------------------------
    ! type for rectangular geometry
    type  GeomRec
        logical, public                           :: is_mesh_define = .FALSE.
        logical, public                           :: is_bc_define = .FALSE.
        logical, public                           :: is_xdefine = .FALSE.
        logical, public                           :: is_ydefine = .FALSE.
        character(len=MAX_LINE_LEN)               :: xchar = CHAR_NULL
        character(len=MAX_LINE_LEN)               :: ychar = CHAR_NULL
        
        integer, public                           :: nx 
        integer, public                           :: ny 
        real(KREAL), public, allocatable          :: xdim(:)
        real(KREAL), public, allocatable          :: ydim(:)
        character(len=MAX_WORD_LEN), allocatable  :: str0(:)
        integer, public, allocatable              :: conf(:, :)
        integer, public, allocatable              :: rowFA(:)                   ! number of sub-FA per row
        real(KREAL), public, allocatable          :: cx(:, :)                   ! center x 
        real(KREAL), public, allocatable          :: cy(:, :)                   ! center y 
        real(KREAL), public                       :: bc_xmin
        real(KREAL), public                       :: bc_xmax
        real(KREAL), public                       :: bc_ymin
        real(KREAL), public                       :: bc_ymax
        integer, public                           :: meshsize 
    contains
        procedure, public  :: alloc => Allocate_GeomRec
        procedure, public  :: clean => Free_GeomRec
        procedure, public  :: set => Set_GeomRec_mesh 
    end type GeomRec

    ! type for hexagonal geometry
    type  GeomHex
        logical, public                           :: is_mesh_define = .FALSE.
        logical, public                           :: is_bc_define = .FALSE.
        
        integer, public                           :: degree                     ! degree region to calculated
        integer, public                           :: ring                       ! number of rings, including the center one
        real(KREAL), public                       :: pitch                      ! sub-FA pitch, in cm
        character(len=MAX_WORD_LEN), allocatable  :: str0(:)                    ! sub-FA configure in 1D string
        integer, public, allocatable              :: conf(:, :)                 ! sub-FA configure in 2D integer
        integer, public, allocatable              :: rowFA(:)                   ! number of sub-FA per row
        real(KREAL), public, allocatable          :: cx(:, :)                   ! center x 
        real(KREAL), public, allocatable          :: cy(:, :)                   ! center y 
        real(KREAL), public                       :: bc_symline
        real(KREAL), public                       :: bc_outer
        integer, public                           :: meshsize 
    contains
        procedure, public  :: alloc => Allocate_GeomHex
        procedure, public  :: clean => Free_GeomHex
        procedure, public  :: set => Set_GeomHex_mesh
    end type GeomHex
    
    ! --------------------------------------------------------------------------
    type, private  :: global_xy_tp
        integer, public                        :: point = 0
        type(coordinate_info_tp), allocatable  :: xy(:)
    contains
        procedure, public  :: add => Add_point_info
        procedure, public  :: clean => Free_point_info 
    end type global_xy_tp 
    
    type(global_xy_tp)  :: rec_xy
    type(global_xy_tp)  :: hex_xy
    
    integer, parameter  :: BY_SCAN = 0
    integer, parameter  :: BY_READ = 1 
    
    integer, parameter  :: SECTOR60   = 60
    integer, parameter  :: SECTOR90   = 90
    integer, parameter  :: SECTOR120  = 120
    integer, parameter  :: SECTOR180  = 180
    integer, parameter  :: SECTOR360  = 360
    
    integer, parameter  :: SEGMENT_LAST  = 4
    integer, parameter  :: SEGMENT_LL    = 3
    integer, parameter  :: SEGMENT_OTHER = 1
    
    integer, parameter  :: SEGMENT_XMAX  = 1  
    integer, parameter  :: SEGMENT_YMAX  = 2 
    integer, parameter  :: SEGMENT_XMIN  = 3 
    integer, parameter  :: SEGMENT_YMIN  = 4 
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_GeomRec (this)
        
        class(GeomRec), intent(in out)  :: this
        
        integer  :: ny, nx 
        integer  :: i_allocate
        
        ! check for allocated status first
        call this%clean ()
        
        ny = this%ny 
        nx = this%nx 
        
        allocate(this%xdim(nx), stat=i_allocate)
        allocate(this%ydim(ny), stat=i_allocate)
        allocate(this%str0(ny), stat=i_allocate)
        allocate(this%conf(ny, nx), stat=i_allocate)
        allocate(this%rowFA(ny), stat=i_allocate)
        allocate(this%cx(ny, nx), stat=i_allocate)
        allocate(this%cy(ny, nx), stat=i_allocate)
    
        this%xdim   = REAL_ZERO
        this%ydim   = REAL_ZERO
        this%str0 = CHAR_NULL
        this%conf = INT_ZERO
        this%rowFA = INT_ZERO
        this%cx = REAL_ZERO
        this%cy = REAL_ZERO
        
        this%bc_xmin = REAL_ZERO
        this%bc_xmax = REAL_ZERO
        this%bc_ymin = REAL_ZERO
        this%bc_ymax = REAL_ZERO
        this%meshsize = INT_ONE
        
    end subroutine Allocate_GeomRec
    
    !$
    !===============================================================================================
    ! finalizer for class of GeomRec
    !===============================================================================================
    subroutine Free_GeomRec (this)
    
        class(GeomRec), intent(in out)  :: this
        
        if (allocated(this%xdim))           deallocate(this%xdim)
        if (allocated(this%ydim))           deallocate(this%ydim)
        if (allocated(this%str0))           deallocate(this%str0)
        if (allocated(this%conf))           deallocate(this%conf)
        if (allocated(this%rowFA))          deallocate(this%rowFA)
        if (allocated(this%cx))             deallocate(this%cx)
        if (allocated(this%cy))             deallocate(this%cy)
    
    end subroutine Free_GeomRec
    
    !$
    !===============================================================================================
    ! set rec-mesh by meshsize 
    !===============================================================================================
    subroutine Set_GeomRec_mesh (this, ns, geom, mesh, bound, mesh_vtk, model)
        
        class(GeomRec), intent(in out)     :: this
        type(SteadyState), intent(in out)  :: ns
        type(Geometry), intent(in out)     :: geom
        type(Meshing), intent(in out)      :: mesh
        type(Boundary), intent(in out)     :: bound
        type(VTKMeshing), intent(in out)   :: mesh_vtk
        integer, intent(in)                :: model 
        
        logical  :: is_4side = .TRUE.
        integer  :: i, j, k 
        integer  :: ibeg, iend 
        integer  :: jbeg, jend 
        integer  :: zone_idx 
        integer  :: npt
        integer, allocatable  :: pidx(:), eidx(:)
        
        real(KREAL)  :: DX, DY
        real(KREAL)  :: ex(4), ey(4)
        real(KREAL), allocatable  :: px(:)
        real(KREAL), allocatable  :: py(:)
        
        character(len=MAX_WORD_LEN)  :: words(MAX_WORDS)                        ! hold for words after a line splited
        integer  :: n_word                                                      ! how many words of a line
        
        if (model == BY_SCAN)  then
            ns%flag%is_60degree = .FALSE.
            ns%flag%is_bevel_edge = .TRUE.
            ns%state%point = 0
            ns%state%nodal = 0
            ns%state%zone = 0
            ns%state%segment = 0
        end if 
        
        ! conf from string to integer
        this%conf = 0
        zone_idx = 0 
        do i = 1, SIZE(this%str0)
            call Split_string (this%str0(i), words, n_word)
            jbeg = String_to_int(words(1))
            jend = jbeg + n_word - 2 
            
            this%rowFA(i) = jend - jbeg + 1  
            do j = jbeg, jend 
                this%conf(i, j) = String_to_int(words(j-jbeg+2))
                if (this%conf(i, j) == 0)  then
                    cycle 
                end if 
                zone_idx = zone_idx + 1 
            end do 
        end do 
        
        if (zone_idx == this%nx*this%ny)  then
            is_4side = .TRUE.
        else 
            is_4side = .FALSE.
        end if 
        
        if (model == BY_READ)  then
            bound%nodal = this%bc_xmax 
            bound%segmentID = SEGMENT_XMAX
        end if 
        
        zone_idx = 0
        ibeg = 0; iend = 0;
        ! coarse-mesh
        if (this%meshsize == 1)  then
            if (is_4side)  then
                do i = 1, this%ny
                    do j = 1, this%nx
                        zone_idx = zone_idx + 1
                        DX = this%xdim(j)
                        DY = this%ydim(i)
                        this%cx(i, j) = SUM(this%xdim(1:j)) - 0.5*DX
                        this%cy(i, j) = SUM(this%ydim(1:i)) - 0.5*DY
                        call get_rec_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), DX, DY, px, py, npt, ex, ey) 
                        
                        call rec_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                        ibeg = iend + 1
                        iend = iend + 4 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)] ! x+
                            mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)] ! y+
                            mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)] ! x-
                            mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(2)] ! y-
                        end if 
                        
                        if (model == BY_READ)  then 
                            ! corner 
                            if ((i == 1) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+3) = this%bc_ymin
                                bound%segmentID(1, ibeg+3) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+2) = this%bc_xmin
                                bound%segmentID(1, ibeg+2) = SEGMENT_XMIN
                            else if ((i == this%ny) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+1) = this%bc_ymax
                                bound%segmentID(1, ibeg+1) = SEGMENT_YMAX
                                bound%nodal(1, ibeg+2) = this%bc_xmin
                                bound%segmentID(1, ibeg+2) = SEGMENT_XMIN
                            else if ((i == 1) .AND. (j == this%nx))  then
                                bound%nodal(1, ibeg+3) = this%bc_ymin
                                bound%segmentID(1, ibeg+3) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+0) = this%bc_xmax
                                bound%segmentID(1, ibeg+0) = SEGMENT_XMAX
                            else if ((i == this%ny) .AND. (j == this%nx))  then
                                bound%nodal(1, ibeg+1) = this%bc_ymax
                                bound%segmentID(1, ibeg+1) = SEGMENT_YMAX
                                bound%nodal(1, ibeg+0) = this%bc_xmax
                                bound%segmentID(1, ibeg+0) = SEGMENT_XMAX
                                
                            ! edge 
                            else if (i == 1)  then
                                bound%nodal(1, ibeg+3) = this%bc_ymin
                                bound%segmentID(1, ibeg+3) = SEGMENT_YMIN
                            else if (i == this%ny)  then
                                bound%nodal(1, ibeg+1) = this%bc_ymax
                                bound%segmentID(1, ibeg+1) = SEGMENT_YMAX
                            else if (j == 1)  then
                                bound%nodal(1, ibeg+2) = this%bc_xmin
                                bound%segmentID(1, ibeg+2) = SEGMENT_XMIN
                            else if (j == this%nx)  then
                                bound%nodal(1, ibeg+0) = this%bc_xmax
                                bound%segmentID(1, ibeg+0) = SEGMENT_XMAX
                            end if 
                        end if 
                        
                    end do 
                end do 
                
            
            else
                do i = 1, this%ny
                    do j = 1, this%nx
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1
                        DX = this%xdim(j)
                        DY = this%ydim(i)
                        this%cx(i, j) = SUM(this%xdim(1:j)) - 0.5*DX
                        this%cy(i, j) = SUM(this%ydim(1:i)) - 0.5*DY
                        call get_rec_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), DX, DY, px, py, npt, ex, ey) 
                        
                        call rec_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                        ibeg = iend + 1
                        iend = iend + 4 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)] ! x+
                            mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)] ! y+
                            mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)] ! x-
                            mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(2)] ! y-
                        end if 
                        
                        if (model == BY_READ)  then 
                            if ((i == 1) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+3) = this%bc_ymin
                                bound%segmentID(1, ibeg+3) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+2) = this%bc_xmin
                                bound%segmentID(1, ibeg+2) = SEGMENT_XMIN
                            else if (i == 1)  then
                                bound%nodal(1, ibeg+3) = this%bc_ymin
                                bound%segmentID(1, ibeg+3) = SEGMENT_YMIN
                            else if (j == 1)  then
                                bound%nodal(1, ibeg+2) = this%bc_xmin
                                bound%segmentID(1, ibeg+2) = SEGMENT_XMIN
                            end if 
                        end if 
                        
                    end do 
                    
                end do 
            end if 
        
        ! fine-mesh
        else if (this%meshsize == 2)  then
            if (is_4side)  then
                do i = 1, this%ny
                    do j = 1, this%nx
                        zone_idx = zone_idx + 1
                        DX = this%xdim(j)
                        DY = this%ydim(i)
                        this%cx(i, j) = SUM(this%xdim(1:j)) - 0.5*DX
                        this%cy(i, j) = SUM(this%ydim(1:i)) - 0.5*DY
                        call get_rec_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), DX, DY, px, py, npt, ex, ey) 
                        
                        call rec_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                        call rec_xy%add (ex, ey, [1, 2, 3, 4], eidx)
                        ibeg = iend + 1
                        iend = iend + 8 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(1), pidx(2), eidx(1)] ! x+
                            mesh%point(:, ibeg+1) = [pidx(1), eidx(1), pidx(3)] ! x+
                            mesh%point(:, ibeg+2) = [pidx(1), pidx(3), eidx(2)] ! y+
                            mesh%point(:, ibeg+3) = [pidx(1), eidx(2), pidx(4)] ! y+
                            mesh%point(:, ibeg+4) = [pidx(1), pidx(4), eidx(3)] ! x-
                            mesh%point(:, ibeg+5) = [pidx(1), eidx(3), pidx(5)] ! x-
                            mesh%point(:, ibeg+6) = [pidx(1), pidx(2), eidx(4)] ! y-
                            mesh%point(:, ibeg+7) = [pidx(1), eidx(4), pidx(5)] ! y-
                        end if 
                        
                        if (model == BY_READ)  then 
                            ! corner 
                            if ((i == 1) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+6: ibeg+7) = this%bc_ymin
                                bound%segmentID(1, ibeg+6: ibeg+7) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+4: ibeg+5) = this%bc_xmin
                                bound%segmentID(1, ibeg+4: ibeg+5) = SEGMENT_XMIN
                            else if ((i == this%ny) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+2: ibeg+3) = this%bc_ymax
                                bound%segmentID(1, ibeg+2: ibeg+3) = SEGMENT_YMAX
                                bound%nodal(1, ibeg+4: ibeg+5) = this%bc_xmin
                                bound%segmentID(1, ibeg+4: ibeg+5) = SEGMENT_XMIN
                            else if ((i == 1) .AND. (j == this%nx))  then
                                bound%nodal(1, ibeg+6: ibeg+7) = this%bc_ymin
                                bound%segmentID(1, ibeg+6: ibeg+7) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+0: ibeg+1) = this%bc_xmax
                                bound%segmentID(1, ibeg+0: ibeg+1) = SEGMENT_XMAX
                            else if ((i == this%ny) .AND. (j == this%nx))  then
                                bound%nodal(1, ibeg+2: ibeg+3) = this%bc_ymax
                                bound%segmentID(1, ibeg+2: ibeg+3) = SEGMENT_YMAX
                                bound%nodal(1, ibeg+0: ibeg+1) = this%bc_xmax
                                bound%segmentID(1, ibeg+0: ibeg+1) = SEGMENT_XMAX
                                
                            ! edge 
                            else if (i == 1)  then
                                bound%nodal(1, ibeg+6: ibeg+7) = this%bc_ymin
                                bound%segmentID(1, ibeg+6: ibeg+7) = SEGMENT_YMIN
                            else if (i == this%ny)  then
                                bound%nodal(1, ibeg+2: ibeg+3) = this%bc_ymax
                                bound%segmentID(1, ibeg+2: ibeg+3) = SEGMENT_YMAX
                            else if (j == 1)  then
                                bound%nodal(1, ibeg+4: ibeg+5) = this%bc_xmin
                                bound%segmentID(1, ibeg+4: ibeg+5) = SEGMENT_XMIN
                            else if (j == this%nx)  then
                                bound%nodal(1, ibeg+0: ibeg+1) = this%bc_xmax
                                bound%segmentID(1, ibeg+0: ibeg+1) = SEGMENT_XMAX
                            end if 
                        end if 
                        
                    end do 
                end do 
            
            else
                do i = 1, this%ny
                    do j = 1, this%nx
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1
                        DX = this%xdim(j)
                        DY = this%ydim(i)
                        this%cx(i, j) = SUM(this%xdim(1:j)) - 0.5*DX
                        this%cy(i, j) = SUM(this%ydim(1:i)) - 0.5*DY
                        call get_rec_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), DX, DY, px, py, npt, ex, ey) 
                        
                        call rec_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                        call rec_xy%add (ex, ey, [1, 2, 3, 4], eidx)
                        ibeg = iend + 1
                        iend = iend + 8 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(1), pidx(2), eidx(1)] ! x+
                            mesh%point(:, ibeg+1) = [pidx(1), eidx(1), pidx(3)] ! x+
                            mesh%point(:, ibeg+2) = [pidx(1), pidx(3), eidx(2)] ! y+
                            mesh%point(:, ibeg+3) = [pidx(1), eidx(2), pidx(4)] ! y+
                            mesh%point(:, ibeg+4) = [pidx(1), pidx(4), eidx(3)] ! x-
                            mesh%point(:, ibeg+5) = [pidx(1), eidx(3), pidx(5)] ! x-
                            mesh%point(:, ibeg+6) = [pidx(1), pidx(2), eidx(4)] ! y-
                            mesh%point(:, ibeg+7) = [pidx(1), eidx(4), pidx(5)] ! y-
                        end if 
                        
                        if (model == BY_READ)  then 
                            if ((i == 1) .AND. (j == 1))  then
                                bound%nodal(1, ibeg+6: ibeg+7) = this%bc_ymin
                                bound%segmentID(1, ibeg+6: ibeg+7) = SEGMENT_YMIN
                                bound%nodal(1, ibeg+4: ibeg+5) = this%bc_xmin
                                bound%segmentID(1, ibeg+4: ibeg+5) = SEGMENT_XMIN
                            else if (i == 1)  then
                                bound%nodal(1, ibeg+6: ibeg+7) = this%bc_ymin
                                bound%segmentID(1, ibeg+6: ibeg+7) = SEGMENT_YMIN
                            else if (j == 1)  then
                                bound%nodal(1, ibeg+4: ibeg+5) = this%bc_xmin
                                bound%segmentID(1, ibeg+4: ibeg+5) = SEGMENT_XMIN
                            end if 
                        end if 
                        
                    end do 
                end do 
            end if 
            
        end if 
        
        if (is_4side)  then
            ns%flag%is_bevel_edge = .FALSE.
            ns%state%segment = SEGMENT_YMIN 
        else
            ns%flag%is_bevel_edge = .TRUE.
            ns%state%segment = SEGMENT_LAST                                     ! for bevel edge case, only segment and segment-1 is adopted
        end if 
        
        ns%state%point = SIZE(rec_xy%xy)
        ns%state%nodal = iend 
        ns%state%zone = zone_idx
        
        if (model == BY_READ)  then
            do i = 1, ns%state%point
                geom%coordinate(1, i) = rec_xy%xy(i)%x
                geom%coordinate(2, i) = rec_xy%xy(i)%y
            end do 
        end if 
        
        call rec_xy%clean ()
        if (allocated(px))          deallocate(px)
        if (allocated(py))          deallocate(py)
        if (allocated(pidx))        deallocate(pidx)
        if (allocated(eidx))        deallocate(eidx)
    
    end subroutine Set_GeomRec_mesh 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_GeomHex (this)
        
        class(GeomHex), intent(in out)  :: this
        
        integer  :: ny, nx 
        integer  :: i_allocate
        
        ! check for allocated status first
        call this%clean ()
        
        if (this%degree == SECTOR360)  then
            ny = 2*this%ring - 1
        else
            ny = this%ring
        end if 
        
        if ((this%degree == SECTOR180) .OR. (this%degree == SECTOR360)) then
            nx = 2*this%ring - 1
        else
            nx = this%ring
        end if 
        
        allocate(this%str0(ny), stat=i_allocate)
        allocate(this%conf(ny, nx), stat=i_allocate)
        allocate(this%rowFA(ny), stat=i_allocate)
        allocate(this%cx(ny, nx), stat=i_allocate)
        allocate(this%cy(ny, nx), stat=i_allocate)
        
        this%str0 = CHAR_NULL
        this%conf = INT_ZERO
        this%rowFA = INT_ZERO
        this%cx = REAL_ZERO
        this%cy = REAL_ZERO
        
        this%bc_symline = REAL_ZERO
        this%bc_outer = REAL_ZERO
        this%meshsize = INT_ONE
        
    end subroutine Allocate_GeomHex
    
    !$
    !===============================================================================================
    ! finalizer for class of GeomHex
    !===============================================================================================
    subroutine Free_GeomHex (this)
    
        class(GeomHex), intent(in out)  :: this
        
        if (allocated(this%str0))           deallocate(this%str0)
        if (allocated(this%conf))           deallocate(this%conf)
        if (allocated(this%rowFA))          deallocate(this%rowFA)
        if (allocated(this%cx))             deallocate(this%cx)
        if (allocated(this%cy))             deallocate(this%cy)
        
        this%is_mesh_define = .FALSE.
        this%is_bc_define = .FALSE.
    
    end subroutine Free_GeomHex
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_GeomHex_mesh (this, ns, geom, mesh, bound, mesh_vtk, model)
    
        class(GeomHex), intent(in out)     :: this
        type(SteadyState), intent(in out)  :: ns
        type(Geometry), intent(in out)     :: geom
        type(Meshing), intent(in out)      :: mesh
        type(Boundary), intent(in out)     :: bound
        type(VTKMeshing), intent(in out)   :: mesh_vtk
        integer, intent(in)                :: model 
        
        logical  :: is_half = .FALSE.
        integer  :: i, j, k 
        integer  :: ibeg, iend 
        integer  :: jbeg, jend 
        integer  :: zone_idx 
        integer  :: npt
        integer, allocatable  :: pidx(:), eidx(:)
        
        real(KREAL)  :: PCH, LEN
        real(KREAL)  :: ex(6), ey(6)
        real(KREAL), allocatable  :: px(:)
        real(KREAL), allocatable  :: py(:)
        
        character(len=MAX_WORD_LEN)  :: words(MAX_WORDS)                        ! hold for words after a line splited
        integer  :: n_word                                                      ! how many words of a line
        
        if (model == BY_SCAN)  then
        select case (this%degree)
            case (SECTOR60, SECTOR120, SECTOR180)
                ns%flag%n_theta = this%degree 
                ns%flag%is_60degree = .TRUE.
                ns%flag%is_bevel_edge = .TRUE.
                ns%state%point = 0
                ns%state%nodal = 0
                ns%state%zone = 0
                ns%state%segment = 0
                
            case (SECTOR90, SECTOR360)
                ns%flag%n_theta = this%degree 
                ns%flag%is_60degree = .FALSE.
                ns%flag%is_bevel_edge = .TRUE.
                ns%state%point = 0
                ns%state%nodal = 0
                ns%state%zone = 0
                ns%state%segment = 0
            
            case default
            end select 
        end if 
        
        ! conf from string to integer
        this%conf = 0
        do i = 1, SIZE(this%str0)
            call Split_string (this%str0(i), words, n_word)
            jbeg = String_to_int(words(1))
            jend = jbeg + n_word - 2 
            
            this%rowFA(i) = jend - jbeg + 1  
            do j = jbeg, jend 
                this%conf(i, j) = String_to_int(words(j-jbeg+2))
            end do 
        end do 
        
        if (model == BY_READ)  then
            bound%nodal = this%bc_outer 
            bound%segmentID = SEGMENT_OTHER
        end if 
        
        zone_idx = 0 
        PCH = this%pitch
        LEN = this%pitch / SQRT(3.0)
        ibeg = 0; iend = 0; 
        ! coarse-mesh 
        if (this%meshsize == 1)  then
            select case (this%degree)
            case (SECTOR60)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH + (i-1)*0.5*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                        
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4], pidx)
                            call hex_xy%add (ex, ey, [2, 3], eidx)
                            ibeg = iend + 1
                            iend = iend + 2 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(2), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+1) = this%bc_symline
                                bound%segmentID(1, ibeg+1) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(4), pidx(1), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(4), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+3) = this%bc_symline
                                bound%segmentID(1, ibeg+3) = SEGMENT_LAST
                            end if 
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4], pidx)
                            call hex_xy%add (ex, ey, [3, 6], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(4), eidx(1), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(4), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(2), pidx(3), pidx(1)]
                                mesh%point(:, ibeg+3) = [pidx(2), pidx(1), eidx(2)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+3) = this%bc_symline
                                bound%segmentID(1, ibeg+3) = SEGMENT_LL
                            end if  
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                                mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                            end if 
                        end if 
                    end do 
                end do 
                
            case (SECTOR90)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        if (MOD(i,2) == 1)  then
                            is_half = .TRUE.
                        else
                            is_half = .FALSE.
                        end if 
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH 
                        if (.NOT. is_half)  then
                            this%cx(i, j) = this%cx(i, j) + 0.5*PCH
                        end if 
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4, 5], pidx)
                            call hex_xy%add (ex, ey, [2], eidx)
                            ibeg = iend + 1
                            iend = iend + 2 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(2), pidx(3), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+1) = this%bc_symline
                                bound%segmentID(1, ibeg+1) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(4), pidx(1), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(4), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+3) = this%bc_symline
                                bound%segmentID(1, ibeg+3) = SEGMENT_LAST
                            end if 
                            
                        else if ((j == 1) .AND. is_half)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                            ibeg = iend + 1
                            iend = iend + 3 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(4), pidx(5), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(4), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(3), pidx(1), pidx(2)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+2) = this%bc_symline
                                bound%segmentID(1, ibeg+2) = SEGMENT_LL
                            end if  
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                                mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                                
                                bound%nodal(1, ibeg+4) = this%bc_symline
                                bound%segmentID(1, ibeg+4) = SEGMENT_LL
                            end if  
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                                mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                            end if 
                        end if 
                    end do 
                end do 
                            
            case (SECTOR120)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH - (i-1)*0.5*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                        
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4, 5], pidx)
                            call hex_xy%add (ex, ey, [2, 4], eidx)
                            ibeg = iend + 1
                            iend = iend + 3 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(2), pidx(3), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(3), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+2) = this%bc_symline
                                bound%segmentID(1, ibeg+2) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(4), pidx(1), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(4), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+3) = this%bc_symline
                                bound%segmentID(1, ibeg+3) = SEGMENT_LAST
                            end if 
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 3, 4, 5], pidx)
                            call hex_xy%add (ex, ey, [1, 4], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(4), eidx(2), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(4), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(2), pidx(3), pidx(1)]
                                mesh%point(:, ibeg+3) = [pidx(2), pidx(1), eidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+3) = this%bc_symline
                                bound%segmentID(1, ibeg+3) = SEGMENT_LL
                            end if  
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                                mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                            end if 
                        end if 
                    end do 
                end do 
            
            case (SECTOR180)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH + (i-1)*0.5*PCH - (this%ring-1)*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 4 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(4), pidx(1), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(4), eidx(2), pidx(1)]
                                
                                if (j < this%ring)  then
                                    bound%nodal(1, ibeg+3) = this%bc_symline
                                    bound%segmentID(1, ibeg+3) = SEGMENT_LL
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                else if (j == this%ring)  then
                                    bound%nodal(1, ibeg+3) = this%bc_symline
                                    bound%segmentID(1, ibeg+3) = SEGMENT_LL
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                else
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                    bound%nodal(1, ibeg+3) = this%bc_symline
                                    bound%segmentID(1, ibeg+3) = SEGMENT_LAST
                                end if 
                            end if 
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                                mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                            end if 
                        end if 
                    end do 
                end do 
            
            case (SECTOR360)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-this%ring)*PCH + ABS((i-this%ring)*0.5*PCH)  
                        
                        this%cy(i, j) = (i-1)*1.5*LEN - (this%ring-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7], pidx)
                        ibeg = iend + 1
                        iend = iend + 6 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(1), pidx(2), pidx(3)]
                            mesh%point(:, ibeg+1) = [pidx(1), pidx(3), pidx(4)]
                            mesh%point(:, ibeg+2) = [pidx(1), pidx(4), pidx(5)]
                            mesh%point(:, ibeg+3) = [pidx(1), pidx(5), pidx(6)]
                            mesh%point(:, ibeg+4) = [pidx(1), pidx(6), pidx(7)]
                            mesh%point(:, ibeg+5) = [pidx(1), pidx(7), pidx(2)]
                        end if 
                    end do 
                end do 
                
            case default
            
            end select 
        
        ! fine-mesh 
        else if (this%meshsize == 2)  then
            select case (this%degree)
            case (SECTOR60)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH + (i-1)*0.5*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                        
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4], pidx)
                            call hex_xy%add (ex, ey, [2, 3], eidx)
                            ibeg = iend + 1
                            iend = iend + 2 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(2), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+1) = this%bc_symline
                                bound%segmentID(1, ibeg+1) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6, 9], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+4) = [pidx(5), pidx(4), eidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+5) = this%bc_symline
                                bound%segmentID(1, ibeg+5) = SEGMENT_LAST
                            end if 
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4, 8], pidx)
                            call hex_xy%add (ex, ey, [3, 6], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), eidx(1), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), pidx(4), eidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+4) = [pidx(5), eidx(2), pidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), pidx(1), eidx(2)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+5) = this%bc_symline
                                bound%segmentID(1, ibeg+5) = SEGMENT_LL
                            end if   
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                            ibeg = iend + 1
                            iend = iend + 10 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                                mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                                mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                                mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                                mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                                mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                            end if 
                        end if
                        
                    end do 
                end do 
                
            case (SECTOR90)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        if (MOD(i,2) == 1)  then
                            is_half = .TRUE.
                        else
                            is_half = .FALSE.
                        end if 
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH 
                        if (.NOT. is_half)  then
                            this%cx(i, j) = this%cx(i, j) + 0.5*PCH
                        end if 
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4, 5, 9], pidx)
                            call hex_xy%add (ex, ey, [2], eidx)
                            ibeg = iend + 1
                            iend = iend + 3 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(2), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(2), pidx(4), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(2), pidx(3), pidx(4)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+1) = this%bc_symline
                                bound%segmentID(1, ibeg+1) = SEGMENT_LL
                                bound%nodal(1, ibeg+2) = this%bc_symline
                                bound%segmentID(1, ibeg+2) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6, 9], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+4) = [pidx(5), pidx(4), eidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+5) = this%bc_symline
                                bound%segmentID(1, ibeg+5) = SEGMENT_LAST
                            end if  
                            
                        else if ((j == 1) .AND. is_half)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5], pidx)
                            ibeg = iend + 1
                            iend = iend + 3 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(4), pidx(5), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(3), pidx(4), pidx(1)]
                                mesh%point(:, ibeg+2) = [pidx(3), pidx(1), pidx(2)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+2) = this%bc_symline
                                bound%segmentID(1, ibeg+2) = SEGMENT_LL
                            end if  
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                            ibeg = iend + 1
                            iend = iend + 10 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                                mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                                mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                                mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                                mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                                mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                                
                                bound%nodal(1, ibeg+4) = this%bc_symline
                                bound%segmentID(1, ibeg+4) = SEGMENT_LL
                            end if 
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                            ibeg = iend + 1
                            iend = iend + 10 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                                mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                                mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                                mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                                mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                                mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                            end if 
                        end if
                    end do 
                end do 
                            
            case (SECTOR120)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH - (i-1)*0.5*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                        
                        if ((i == 1) .AND. (j == 1))  then
                            call hex_xy%add (px, py, [1, 4, 5, 11], pidx)
                            call hex_xy%add (ex, ey, [2, 4], eidx)
                            ibeg = iend + 1
                            iend = iend + 5 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(4), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(4), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(4), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(4), pidx(3), eidx(2)]
                                mesh%point(:, ibeg+4) = [pidx(4), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+4) = this%bc_symline
                                bound%segmentID(1, ibeg+4) = SEGMENT_LL
                            end if 
                            
                        else if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6, 9], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+4) = [pidx(5), pidx(4), eidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), eidx(2), pidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                bound%nodal(1, ibeg+5) = this%bc_symline
                                bound%segmentID(1, ibeg+5) = SEGMENT_LAST
                            end if 
                            
                        else if (j == 1)  then
                            call hex_xy%add (px, py, [1, 3, 4, 5, 11], pidx)
                            call hex_xy%add (ex, ey, [1, 4], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), eidx(2), pidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), pidx(4), eidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+4) = [pidx(5), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), pidx(1), eidx(1)]
                                
                                bound%nodal(1, ibeg+0) = this%bc_symline
                                bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                bound%nodal(1, ibeg+5) = this%bc_symline
                                bound%segmentID(1, ibeg+5) = SEGMENT_LL
                            end if   
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                            ibeg = iend + 1
                            iend = iend + 10 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                                mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                                mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                                mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                                mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                                mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                            end if 
                        end if 
                    end do 
                end do 
            
            case (SECTOR180)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-1)*PCH + (i-1)*0.5*PCH - (this%ring-1)*PCH
                        this%cy(i, j) = (i-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        if (i == 1)  then
                            call hex_xy%add (px, py, [1, 4, 5, 6, 9], pidx)
                            call hex_xy%add (ex, ey, [2, 5], eidx)
                            ibeg = iend + 1
                            iend = iend + 6 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(5), pidx(1), eidx(1)]
                                mesh%point(:, ibeg+1) = [pidx(5), eidx(1), pidx(2)]
                                mesh%point(:, ibeg+2) = [pidx(5), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+3) = [pidx(5), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+4) = [pidx(5), pidx(4), eidx(2)]
                                mesh%point(:, ibeg+5) = [pidx(5), eidx(2), pidx(1)]
                                
                                if (j < this%ring)  then
                                    bound%nodal(1, ibeg+5) = this%bc_symline
                                    bound%segmentID(1, ibeg+5) = SEGMENT_LL
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LL
                                else if (j == this%ring)  then
                                    bound%nodal(1, ibeg+5) = this%bc_symline
                                    bound%segmentID(1, ibeg+5) = SEGMENT_LL
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                else
                                    bound%nodal(1, ibeg+0) = this%bc_symline
                                    bound%segmentID(1, ibeg+0) = SEGMENT_LAST
                                    bound%nodal(1, ibeg+5) = this%bc_symline
                                    bound%segmentID(1, ibeg+5) = SEGMENT_LAST
                                end if 
                            end if 
                            
                        else 
                            call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                            ibeg = iend + 1
                            iend = iend + 10 
                            
                            if (model == BY_READ)  then 
                                mesh%zone(ibeg:iend) = zone_idx
                                mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                                mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                                mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                                mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                                mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                                mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                                mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                                mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                                mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                                mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                            end if 
                        end if 
                    end do 
                end do 
            
            case (SECTOR360)
                do i = 1, SIZE(this%conf, dim=1)
                    do j = 1, SIZE(this%conf, dim=2)
                        if (this%conf(i,j) == 0)  cycle
                        
                        zone_idx = zone_idx + 1 
                        this%cx(i, j) = (j-this%ring)*PCH + ABS((i-this%ring)*0.5*PCH)  
                        
                        this%cy(i, j) = (i-1)*1.5*LEN - (this%ring-1)*1.5*LEN
                        call get_hex_meshsize (this%meshsize, this%cx(i, j), this%cy(i, j), PCH, LEN, px, py, npt, ex, ey)
                
                        call hex_xy%add (px, py, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], pidx)
                        ibeg = iend + 1
                        iend = iend + 10 
                        
                        if (model == BY_READ)  then 
                            mesh%zone(ibeg:iend) = zone_idx
                            mesh%point(:, ibeg+0) = [pidx(8), pidx(2), pidx(3)]
                            mesh%point(:, ibeg+1) = [pidx(8), pidx(3), pidx(4)]
                            mesh%point(:, ibeg+2) = [pidx(9), pidx(4), pidx(5)]
                            mesh%point(:, ibeg+3) = [pidx(9), pidx(5), pidx(6)]
                            mesh%point(:, ibeg+4) = [pidx(10), pidx( 6), pidx( 7)]
                            mesh%point(:, ibeg+5) = [pidx(10), pidx( 7), pidx( 2)]
                            mesh%point(:, ibeg+6) = [pidx( 8), pidx(10), pidx( 2)]
                            mesh%point(:, ibeg+7) = [pidx( 9), pidx( 8), pidx( 4)]
                            mesh%point(:, ibeg+8) = [pidx(10), pidx( 9), pidx( 6)]
                            mesh%point(:, ibeg+9) = [pidx( 8), pidx( 9), pidx(10)]
                        end if
                    end do 
                end do 
                
            case default
            
            end select 
                
        end if 
        
        ns%state%point = SIZE(hex_xy%xy)
        ns%state%nodal = iend 
        ns%state%zone = zone_idx
        ns%state%segment = SEGMENT_LAST                                         ! for bevel edge case, only segment and segment-1 is adopted
        
        if (model == BY_READ)  then
            do i = 1, ns%state%point
                geom%coordinate(1, i) = hex_xy%xy(i)%x
                geom%coordinate(2, i) = hex_xy%xy(i)%y
            end do 
        end if 
        
        call hex_xy%clean ()
        if (allocated(px))          deallocate(px)
        if (allocated(py))          deallocate(py)
        if (allocated(pidx))        deallocate(pidx)
        if (allocated(eidx))        deallocate(eidx)
        
    end subroutine Set_GeomHex_mesh
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! meshsize = 1/2
    !===============================================================================================
    subroutine get_rec_meshsize (meshsize, cx, cy, DX, DY, px, py, npt, ex, ey)
        
        integer, intent(in)      :: meshsize 
        real(KREAL), intent(in)  :: cx
        real(KREAL), intent(in)  :: cy
        real(KREAL), intent(in)  :: DX
        real(KREAL), intent(in)  :: DY
        real(KREAL), intent(in out), allocatable  :: px(:)
        real(KREAL), intent(in out), allocatable  :: py(:)
        integer, intent(out)         :: npt
        real(KREAL), intent(in out)  :: ex(:)
        real(KREAL), intent(in out)  :: ey(:)
        
        integer  :: nzone
        integer  :: i_allocate
        
        ex(1) = cx + 0.5*DX;  ey(1) = cy + 0.0*DY;
        ex(2) = cx + 0.0*DX;  ey(2) = cy + 0.5*DY;
        ex(3) = cx - 0.5*DX;  ey(3) = cy + 0.0*DY;
        ex(4) = cx + 0.0*DX;  ey(4) = cy - 0.5*DY;
        
        if (allocated(px))      deallocate(px)
        if (allocated(py))      deallocate(py)
        
        select case (meshsize)
        case (1)
            npt = 5
            nzone = 4
            allocate(px(npt), stat=i_allocate)
            allocate(py(npt), stat=i_allocate)
            px(1) = cx + 0.0*DX;  py(1) = cy + 0.0*DY;
            
            px(2) = cx + 0.5*DX;  py(2) = cy - 0.5*DY;
            px(3) = cx + 0.5*DX;  py(3) = cy + 0.5*DY;
            px(4) = cx - 0.5*DX;  py(4) = cy + 0.5*DY;
            px(5) = cx - 0.5*DX;  py(5) = cy - 0.5*DY;
            
        case (2)
            npt = 5
            nzone = 8
            allocate(px(npt), stat=i_allocate)
            allocate(py(npt), stat=i_allocate)
            px(1) = cx + 0.0*DX;  py(1) = cy + 0.0*DY;
            
            px(2) = cx + 0.5*DX;  py(2) = cy - 0.5*DY;
            px(3) = cx + 0.5*DX;  py(3) = cy + 0.5*DY;
            px(4) = cx - 0.5*DX;  py(4) = cy + 0.5*DY;
            px(5) = cx - 0.5*DX;  py(5) = cy - 0.5*DY;
        
        case default
        
        end select 
        
    end subroutine get_rec_meshsize 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine get_hex_meshsize (meshsize, cx, cy, PCH, LEN, px, py, npt, ex, ey)
        
        integer, intent(in)      :: meshsize 
        real(KREAL), intent(in)  :: cx
        real(KREAL), intent(in)  :: cy
        real(KREAL), intent(in)  :: PCH
        real(KREAL), intent(in)  :: LEN
        real(KREAL), intent(in out), allocatable  :: px(:)
        real(KREAL), intent(in out), allocatable  :: py(:)
        integer, intent(out)         :: npt
        real(KREAL), intent(in out)  :: ex(:)
        real(KREAL), intent(in out)  :: ey(:)
        
        integer  :: nzone
        integer  :: i_allocate
        
        ex(1) = cx + 0.5*PCH*COS60;  ey(1) = cy - 0.5*PCH*SIN60;
        ex(2) = cx + 0.5*PCH;        ey(2) = cy + 0.0*PCH;
        ex(3) = cx + 0.5*PCH*COS60;  ey(3) = cy + 0.5*PCH*SIN60;
        ex(4) = cx - 0.5*PCH*COS60;  ey(4) = cy + 0.5*PCH*SIN60;
        ex(5) = cx - 0.5*PCH;        ey(5) = cy + 0.0*PCH;
        ex(6) = cx - 0.5*PCH*COS60;  ey(6) = cy - 0.5*PCH*SIN60;
        
        if (allocated(px))      deallocate(px)
        if (allocated(py))      deallocate(py)
        
        select case (meshsize)
        case (1)
            npt = 7
            nzone = 6
            allocate(px(npt), stat=i_allocate)
            allocate(py(npt), stat=i_allocate)
            px(1) = cx + 0.0*PCH;  py(1) = cy + 0.0*LEN;
            
            px(2) = cx + 0.0*PCH;  py(2) = cy - 1.0*LEN;
            px(3) = cx + 0.5*PCH;  py(3) = cy - 0.5*LEN;
            px(4) = cx + 0.5*PCH;  py(4) = cy + 0.5*LEN;
            px(5) = cx + 0.0*PCH;  py(5) = cy + 1.0*LEN;
            px(6) = cx - 0.5*PCH;  py(6) = cy + 0.5*LEN;
            px(7) = cx - 0.5*PCH;  py(7) = cy - 0.5*LEN;
            
        case (2)
            npt = 13
            nzone = 10
            allocate(px(npt), stat=i_allocate)
            allocate(py(npt), stat=i_allocate)
            px(1) = cx + 0.0*PCH;  py(1) = cy + 0.0*LEN;
            
            px(2) = cx + 0.0*PCH;  py(2) = cy - 1.0*LEN;
            px(3) = cx + 0.5*PCH;  py(3) = cy - 0.5*LEN;
            px(4) = cx + 0.5*PCH;  py(4) = cy + 0.5*LEN;
            px(5) = cx + 0.0*PCH;  py(5) = cy + 1.0*LEN;
            px(6) = cx - 0.5*PCH;  py(6) = cy + 0.5*LEN;
            px(7) = cx - 0.5*PCH;  py(7) = cy - 0.5*LEN;
            
            px(8) = cx + (1.0/6.0)*PCH;  py(8) = cy - (1.0/6.0)*LEN;
            px(9) = cx + 0.0*PCH;        py(9) = cy + (1.0/3.0)*LEN;
            px(10)= cx - (1.0/6.0)*PCH;  py(10)= cy - (1.0/6.0)*LEN;

            px(11)= cx + (1.0/6.0)*PCH;  py(11)= cy + (1.0/6.0)*LEN;
            px(12)= cx + 0.0*PCH;        py(12)= cy - (1.0/3.0)*LEN;
            px(13)= cx - (1.0/6.0)*PCH;  py(13)= cy + (1.0/6.0)*LEN;
                                      
        case default
        
        end select
        
    end subroutine get_hex_meshsize 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Add_point_info (this, px, py, idx, newidx)
        
        class(global_xy_tp), intent(in out)  :: this
        real(KREAL), intent(in)              :: px(:)
        real(KREAL), intent(in)              :: py(:)
        integer, intent(in)                  :: idx(:)
        integer, intent(in out), allocatable :: newidx(:) 
        
        type(coordinate_info_tp), allocatable  :: tmp(:)
        type(coordinate_info_tp)               :: point
        integer  :: cnt
        integer  :: i, j, i_allocate
        logical  :: is_added = .FALSE.
        
        outer1: do i = 1, SIZE(idx)
            point%x = px(idx(i))
            point%y = py(idx(i))
            
            if (allocated(this%xy))  then
                ! is existed ?
                inner1: do j = 1, SIZE(this%xy)
                    is_added = this%xy(j)%is_same (point)
                    if (is_added)  then
                        cycle outer1
                    end if
                end do inner1 
            
                this%point = this%point + 1
                cnt = SIZE(this%xy)
                allocate(tmp(cnt), stat=i_allocate)
                tmp = this%xy
                
                deallocate(this%xy)
                allocate(this%xy(cnt+1), stat=i_allocate)
                this%xy(1:cnt) = tmp
                this%xy(cnt+1) = point
                deallocate(tmp)
                
            ! this is the first point
            else
                this%point = 1
                allocate(this%xy(1), stat=i_allocate)
                this%xy(1) = point
            end if
        end do outer1 
        
        ! get-index 
        if (allocated(newidx))      deallocate(newidx)
        allocate(newidx(SIZE(idx)), stat=i_allocate)
        newidx = 0 
        
        outer2: do i = 1, SIZE(idx)
            point%x = px(idx(i))
            point%y = py(idx(i))
            
            if (allocated(this%xy))  then
                inner2: do j = 1, SIZE(this%xy)
                    is_added = this%xy(j)%is_same (point)
                    if (is_added)  then
                        newidx(i) = j
                        cycle outer2
                    end if
                end do inner2
            end if
        end do outer2
    
    end subroutine Add_point_info
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_point_info (this)
        
        class(global_xy_tp), intent(in out)  :: this
        
        if (allocated(this%xy))         deallocate(this%xy)
    
    end subroutine Free_point_info
    
end module geomregular_header
