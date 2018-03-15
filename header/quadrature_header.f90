!$
!===================================================================================================
!
!   class for quadrature sets
!   NOTE: generially, the total weight for all directions can be normalled to 1, 4, 8 or 4*pi,
!         here we choose 1, that means the 4*pi is integratd into angular flux
!         see page 191 in J.J Stamm'ler
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          QuadratureSet
!
!===================================================================================================
module quadrature_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use state_header,               only : SteadyState
    use geometry_header,            only : Meshing, Geometry
    
    implicit none 
    private
    public  :: QuadratureSet
    
    ! --------------------------------------------------------------------------
    ! type for quadrature set
    type, private  :: quadrature_set_info_tp
        real(KREAL), public               :: wmu                                ! weight per ordinates
        real(KREAL), public               :: xmu(3)                             ! project to global coordinates
        real(KREAL), public, allocatable  :: projection(:, :)                   ! project for all nodal
    contains
        procedure, private  :: alloc => Allocate_quadrature_set_info_tp
        procedure, private  :: clean => Free_quadrature_set_info_tp
    end type quadrature_set_info_tp
    
    !---------------------------------------------------------------------------
    ! type for the really quadrature
    type  QuadratureSet 
        type(quadrature_set_info_tp), public, allocatable  :: directions(:)     ! total direction
        real(KREAL), public  :: memory = REAL_ZERO
        integer, public, allocatable                       :: is_z(:)           ! symmetry with Z line
        integer, public, allocatable                       :: is_x(:)           ! symmetry with X line
        integer, public, allocatable                       :: is_y(:)           ! symmetry with Y line
        integer, public, allocatable                       :: is_symmetry(:)    ! symmetry with 60's line 
        integer, public, allocatable                       :: is_order(:)       ! sweep order 
    contains
        procedure, public  :: alloc => Allocate_QuadratureSet
        procedure, public  :: set_weight => Set_QuadratureSet_weight
        procedure, public  :: set_project => Set_QuadratureSet_project
        procedure, public  :: negate => Negate_QuadratureSet
        procedure, public  :: clean =>  Free_QuadratureSet
        procedure, public  :: set_sym => Set_QuadratureSet_symmetry
        procedure, public  :: set_order => Set_QuadratureSet_order
    end type QuadratureSet
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Allocate_quadrature_set_info_tp, Free_quadrature_set_info_tp
    private  :: Allocate_QuadratureSet, Free_QuadratureSet
    private  :: Set_QuadratureSet_weight, Set_QuadratureSet_project, Negate_QuadratureSet
    
    integer, parameter  :: SYMMETRY_DNTR            = 1
    integer, parameter  :: SYMMETRY_UPDATE          = 2
    integer, parameter  :: FULL_LEVEL_DNTR          = 1
    integer, parameter  :: FULL_LEVEL_HEBERT        = 2
    integer, parameter  :: FULL_PENTRAN_EQUAl       = 3
    integer, parameter  :: FULL_PENTRAN_CHEBYSHEV   = 4
    
    ! selection
    integer, parameter  :: SYMMETRY_SELECTION       = 2
    integer, parameter  :: FULL_SELECTION           = 4

contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_quadrature_set_info_tp (this)
        
        class(quadrature_set_info_tp), intent(in out)  :: this
       
        integer  :: i_allocate
        integer, parameter  :: N_EDGE = 3                                       ! three edge per triangular
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%projection(N_EDGE, ns%state%nodal), stat=i_allocate)
        
        this%xmu         = REAL_ZERO
        this%wmu         = REAL_ZERO
        this%projection  = REAL_ZERO
    
    end subroutine Allocate_quadrature_set_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of quadrature_set_info_tp
    !===============================================================================================
    subroutine Free_quadrature_set_info_tp (this)
        
        class(quadrature_set_info_tp), intent(in out)  :: this
        
        if (allocated(this%projection))     deallocate(this%projection)
        
    end subroutine Free_quadrature_set_info_tp

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_QuadratureSet (this)
        
        class(QuadratureSet), intent(in out)  :: this
        
        integer  :: i_allocate
        integer  :: i, nd 
        
        nd = ns%deduce%direction
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%directions(nd), stat=i_allocate)
        
        do i = 1, nd
            call this%directions(i)%alloc ()
        end do
        
        allocate(this%is_z(nd), stat=i_allocate)
        allocate(this%is_y(nd), stat=i_allocate)
        allocate(this%is_x(nd), stat=i_allocate)
        allocate(this%is_symmetry(nd), stat=i_allocate)
        allocate(this%is_order(nd), stat=i_allocate)
        
        this%is_z = INT_ONE
        this%is_y = INT_ONE
        this%is_x = INT_ONE
        this%is_symmetry = INT_ONE
        this%is_order = INT_ONE
        
        this%memory = REAL_ZERO
        if (SIZE(this%directions) >= 1)  then
            this%memory = this%memory + REAL_BYTE * SIZE(this%directions) * SIZE(this%directions(1)%xmu)
            this%memory = this%memory + REAL_BYTE * SIZE(this%directions) * SIZE(this%directions(1)%projection)
        end if 
        
    end subroutine Allocate_QuadratureSet
    
    !$
    !===============================================================================================
    ! finalizer for class of QuadratureSet
    !===============================================================================================
    subroutine Free_QuadratureSet (this)
        
        class(QuadratureSet), intent(in out)  :: this
        integer  :: i
        
        if (allocated(this%directions))  then
            do i = 1, SIZE(this%directions)
                call this%directions(i)%clean ()
            end do
            ! free itself
            deallocate(this%directions)
        end if 
        
        this%memory = REAL_ZERO
        
    end subroutine Free_QuadratureSet
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_QuadratureSet_weight (this)
        
        class(QuadratureSet), intent(in out)  :: this
    
        ! local varibles
        integer, parameter  :: n_octant = 8                                     ! number of octant
        integer             :: per_octant
        
        per_octant = ns%deduce%direction / n_octant
    
        quad_sets: select case(ns%flag%is_60degree)
        ! case=1, 60 degree
        ! SN        --S2, S4, S6, S8, S10, S12, S16
        ! direction -- 3,  6,  9, 12,  15,  18,  24
        case(.TRUE.)  quad_sets
            
            select case (SYMMETRY_SELECTION)
            case (SYMMETRY_DNTR)
                call Guass_chebyshev_symmetry_DNTR (this, per_octant)
            case (SYMMETRY_UPDATE)
                call Guass_chebyshev_symmetry (this, per_octant)
            end select
    
        ! symmetry quadrature sets
        ! SN         -- S2, S4, S6, S8, S10, S12, S14, S16
        ! direction  -- (1,  3,  6, 10,  15,  21,  28,  36)*8
        case(.FALSE.)  quad_sets
            
            select case (FULL_SELECTION)
            case (FULL_LEVEL_DNTR)
                call Level_symmetry_DNTR (this, per_octant)
            case (FULL_LEVEL_HEBERT)
                call Level_symmetry_Hebert (this, per_octant)
            case (FULL_PENTRAN_EQUAl)
                call Guass_Equal_PENTRAN (this, per_octant)
            case (FULL_PENTRAN_CHEBYSHEV)
                call Guass_ChebyShev_PENTRAN (this, per_octant)
            end select
            
        end select quad_sets
            
        ! extends the weight to all the ordinates of the octant
        call Extend_other_octant (this, per_octant, n_octant)
    
    end subroutine Set_QuadratureSet_weight
    
    !$
    !===============================================================================================
    ! negate the projection to get adjoint flux
    !===============================================================================================
    subroutine Negate_QuadratureSet (this)
        
        class(QuadratureSet), intent(in out)  :: this
        integer  :: is
        
        do is = 1, SIZE(this%directions)
            this%directions(is)%projection(:, :) = - this%directions(is)%projection(:, :)
            this%directions(is)%xmu(:) = - this%directions(is)%xmu(:)
        end do
    
    end subroutine Negate_QuadratureSet
    
    !$
    !===============================================================================================
    ! set the value of projetion per direction per nodal per edge
    !===============================================================================================
    subroutine Set_QuadratureSet_project (this, mesh, geom)
        
        class(QuadratureSet), intent(in out)  :: this
        type(MEshing), intent(in)   :: mesh
        type(Geometry), intent(in)  :: geom
        
        integer  :: i, j, k, is
        real(KREAL), parameter  :: EPS_PROJECTION = 1.0E-10                 ! seperate in and out when direction parallel with the edge

        real(KREAL)  :: tmp_var(2, 3, ns%state%nodal)                       ! temporary variable to get projection
        real(KREAL)  :: x1, y1, x2, y2, x3, y3, ux, uy
        
        do i = 1, ns%state%nodal
            x1 = geom%coordinate(1, mesh%point(1,i))
            y1 = geom%coordinate(2, mesh%point(1,i))
            x2 = geom%coordinate(1, mesh%point(2,i))
            y2 = geom%coordinate(2, mesh%point(2,i))
            x3 = geom%coordinate(1, mesh%point(3,i))
            y3 = geom%coordinate(2, mesh%point(3,i))
            tmp_var(1,1,i) = (y3-y2) / (2.0*geom%area(i))                       ! X partial value of x direction of nodal i
            tmp_var(2,1,i) = (x2-x3) / (2.0*geom%area(i))                       ! Y                  x                    i
            tmp_var(1,2,i) = (y1-y3) / (2.0*geom%area(i))                       ! X                  u                    i  
            tmp_var(2,2,i) = (x3-x1) / (2.0*geom%area(i))                       ! Y                  u                    i  
            tmp_var(1,3,i) = (y2-y1) / (2.0*geom%area(i))                       ! X                  v                    i  
            tmp_var(2,3,i) = (x1-x2) / (2.0*geom%area(i))                       ! Y                  v                    i  
        end do

        do is = 1, ns%deduce%direction
            ! cosine of direction ID is with x and axis, it means the X and Y partial value in radial
            ux = this%directions(is)%xmu(1)
            uy = this%directions(is)%xmu(2)
            do j = 1, ns%state%nodal
                ! in radial, get the inner product of specific direction and the XUV of nodal
                ! result 0 means this direction is normal to radial direction
                do k = 1, 3
                    this%directions(is)%projection(k,j) = tmp_var(1,k,j)*ux + tmp_var(2,k,j)*uy
                end do
            end do
        end do
        
        ! get the incoming and outgoing when the direction parallel with the edge
        do is = 1, ns%deduce%direction
            do j = 1, ns%state%nodal
                do k = 1, 3
                    ! NOTE--use 'EPS_PROJECTION/2.0' as limit, so it will not change twice
                    if (ABS(this%directions(is)%projection(k,j)) < EPS_PROJECTION/2.0) then
                        if ((ns%deduce%direction/4<is.AND.is<=ns%deduce%direction/2) .OR. (3*ns%deduce%direction/4<is.AND.is<=ns%deduce%direction))  then
                            ! incoming edge
                            this%directions(is)%projection(k,j) = -EPS_PROJECTION
                            ! this edge is outgoing in the neighbour
                            this%directions(is)%projection(mesh%localID(k,j), mesh%nearby_nodal(k,j)) = EPS_PROJECTION
                        else
                            this%directions(is)%projection(k,j) = EPS_PROJECTION
                            this%directions(is)%projection(mesh%localID(k,j), mesh%nearby_nodal(k,j)) = -EPS_PROJECTION
                        end if
                    end if
                end do
            end do
        end do
                
    end subroutine Set_QuadratureSet_project
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! extends the weight to all the ordinates of the octant
    !===============================================================================================
    subroutine Extend_other_octant (quad, per_octant, n_octant)
    
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        integer, intent(in)  :: n_octant
        
        integer  :: i_octant, iy
        
        do i_octant = 2, n_octant
            associate (istart => per_octant*(i_octant-1) + 1,   &
                      &  iend => per_octant*i_octant,           &
                      &  f_end => per_octant)
                quad%directions(istart: iend)%wmu = quad%directions(1: f_end)%wmu
                
                do iy = 1, 3
                    quad%directions(istart: iend)%xmu(iy) = quad%directions(1: f_end)%xmu(iy)
                end do
            end associate
        end do
        quad%directions(:)%wmu = quad%directions(:)%wmu / n_octant
    
        do i_octant = 2, n_octant
            associate (istart => per_octant*(i_octant-1) + 1,   &
                      &  iend => per_octant*i_octant,           &
                      &  f_end => per_octant)
                
                select case(i_octant)
                case(2) 
                    quad%directions(istart:iend)%xmu(1) = -quad%directions(1:f_end)%xmu(1)
                case(3)
                    quad%directions(istart:iend)%xmu(2) = -quad%directions(1:f_end)%xmu(2)
                case(4)
                    quad%directions(istart:iend)%xmu(1) = -quad%directions(1:f_end)%xmu(1)
                    quad%directions(istart:iend)%xmu(2) = -quad%directions(1:f_end)%xmu(2)
                case(5)
                    quad%directions(istart:iend)%xmu(3) = -quad%directions(1:f_end)%xmu(3)
                case(6)
                    quad%directions(istart:iend)%xmu(1) = -quad%directions(1:f_end)%xmu(1)
                    quad%directions(istart:iend)%xmu(3) = -quad%directions(1:f_end)%xmu(3)
                case(7)
                    quad%directions(istart:iend)%xmu(2) = -quad%directions(1:f_end)%xmu(2)
                    quad%directions(istart:iend)%xmu(3) = -quad%directions(1:f_end)%xmu(3)
                case(8)
                    quad%directions(istart:iend)%xmu(1) = -quad%directions(1:f_end)%xmu(1)
                    quad%directions(istart:iend)%xmu(2) = -quad%directions(1:f_end)%xmu(2)
                    quad%directions(istart:iend)%xmu(3) = -quad%directions(1:f_end)%xmu(3)
                end select
            end associate
        end do
    
    end subroutine Extend_other_octant
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Guass_chebyshev_symmetry_DNTR (quad, per_octant)
    
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
            
        if(ns%state%sn == 2) then
            quad%directions(1)%xmu(1) = 0.7886751346 
            quad%directions(2)%xmu(1) = 0.5773502692 
            quad%directions(3)%xmu(1) = 0.2113248655
            quad%directions(1)%xmu(2) = 0.2113248655
            quad%directions(2)%xmu(2) = 0.5773502692
            quad%directions(3)%xmu(2) = 0.7886751346
            quad%directions(1)%xmu(3) = 0.5773502692
            quad%directions(2)%xmu(3) = 0.5773502692
            quad%directions(3)%xmu(3) = 0.5773502692
            quad%directions(1)%wmu = 0.3333333333
            quad%directions(2)%wmu = 0.3333333333
            quad%directions(3)%wmu = 0.3333333333
        else if(ns%state%sn == 4) then
            quad%directions(1)%xmu(1) = 0.4910516985
            quad%directions(2)%xmu(1) = 0.3594747924
            quad%directions(3)%xmu(1) = 0.1315769061
            quad%directions(4)%xmu(1) = 0.9083878357
            quad%directions(5)%xmu(1) = 0.6649860488
            quad%directions(6)%xmu(1) = 0.2434017871
            quad%directions(1)%xmu(2) = 0.1315769061
            quad%directions(2)%xmu(2) = 0.3594747924
            quad%directions(3)%xmu(2) = 0.4910516985
            quad%directions(4)%xmu(2) = 0.2434017871
            quad%directions(5)%xmu(2) = 0.6649860488
            quad%directions(6)%xmu(2) = 0.9083878357
            quad%directions(1)%xmu(3) = 0.8611363116
            quad%directions(2)%xmu(3) = 0.8611363116
            quad%directions(3)%xmu(3) = 0.8611363116
            quad%directions(4)%xmu(3) = 0.3399810436
            quad%directions(5)%xmu(3) = 0.3399810436
            quad%directions(6)%xmu(3) = 0.3399810436
            quad%directions(1)%wmu = 0.1159516151
            quad%directions(2)%wmu = 0.1159516151
            quad%directions(3)%wmu = 0.1159516151
            quad%directions(4)%wmu = 0.2173817182
            quad%directions(5)%wmu = 0.2173817182
            quad%directions(6)%wmu = 0.2173817182
        else if(ns%state%sn == 6) then
            quad%directions(1)%xmu(1) = 0.3489394248
            quad%directions(2)%xmu(1) = 0.2554413876
            quad%directions(3)%xmu(1) = 0.0934980370
            quad%directions(4)%xmu(1) = 0.7246389115
            quad%directions(5)%xmu(1) = 0.5304725003
            quad%directions(6)%xmu(1) = 0.1941664112
            quad%directions(7)%xmu(1) = 0.9380233384
            quad%directions(8)%xmu(1) = 0.6866807425
            quad%directions(9)%xmu(1) = 0.2513425960
            quad%directions(1)%xmu(2) = 0.0934980370
            quad%directions(2)%xmu(2) = 0.2554413876
            quad%directions(3)%xmu(2) = 0.3489394248
            quad%directions(4)%xmu(2) = 0.1941664112
            quad%directions(5)%xmu(2) = 0.5304725003
            quad%directions(6)%xmu(2) = 0.7246389115
            quad%directions(7)%xmu(2) = 0.2513425960
            quad%directions(8)%xmu(2) = 0.6866807425
            quad%directions(9)%xmu(2) = 0.9380233384
            quad%directions(1)%xmu(3) = 0.9324695142
            quad%directions(2)%xmu(3) = 0.9324695142
            quad%directions(3)%xmu(3) = 0.9324695142
            quad%directions(4)%xmu(3) = 0.6612093865
            quad%directions(5)%xmu(3) = 0.6612093865
            quad%directions(6)%xmu(3) = 0.6612093865
            quad%directions(7)%xmu(3) = 0.2386191861
            quad%directions(8)%xmu(3) = 0.2386191861
            quad%directions(9)%xmu(3) = 0.2386191861
            quad%directions(1)%wmu = 0.0571081641
            quad%directions(2)%wmu = 0.0571081641
            quad%directions(3)%wmu = 0.0571081641
            quad%directions(4)%wmu = 0.1202538577
            quad%directions(5)%wmu = 0.1202538577
            quad%directions(6)%wmu = 0.1202538577
            quad%directions(7)%wmu = 0.1559713115
            quad%directions(8)%wmu = 0.1559713115
            quad%directions(9)%wmu = 0.1559713115 
        else if(ns%state%sn == 8) then
            quad%directions(1)%xmu(1) = 0.2694974453
            quad%directions(2)%xmu(1) = 0.1972858225
            quad%directions(3)%xmu(1) = 0.0722116229
            quad%directions(4)%xmu(1) = 0.5838240788
            quad%directions(5)%xmu(1) = 0.4273888884
            quad%directions(6)%xmu(1) = 0.1564351904
            quad%directions(7)%xmu(1) = 0.8217841742
            quad%directions(8)%xmu(1) = 0.6015877684
            quad%directions(9)%xmu(1) = 0.2201964057
            quad%directions(10)%xmu(1) = 0.9495359079
            quad%directions(11)%xmu(1) = 0.6951085282
            quad%directions(12)%xmu(1) = 0.2544273797
            quad%directions(1)%xmu(2) = 0.0722116229
            quad%directions(2)%xmu(2) = 0.1972858225
            quad%directions(3)%xmu(2) = 0.2694974453
            quad%directions(4)%xmu(2) = 0.1564351904
            quad%directions(5)%xmu(2) = 0.4273888884
            quad%directions(6)%xmu(2) = 0.5838240788
            quad%directions(7)%xmu(2) = 0.2201964057
            quad%directions(8)%xmu(2) = 0.6015877684
            quad%directions(9)%xmu(2) = 0.8217841742
            quad%directions(10)%xmu(2) = 0.2544273797
            quad%directions(11)%xmu(2) = 0.6951085282
            quad%directions(12)%xmu(2) = 0.9495359079
            quad%directions(1)%xmu(3) = 0.9602898565
            quad%directions(2)%xmu(3) = 0.9602898565
            quad%directions(3)%xmu(3) = 0.9602898565
            quad%directions(4)%xmu(3) = 0.7966664774
            quad%directions(5)%xmu(3) = 0.7966664774
            quad%directions(6)%xmu(3) = 0.7966664774
            quad%directions(7)%xmu(3) = 0.5255324099
            quad%directions(8)%xmu(3) = 0.5255324099
            quad%directions(9)%xmu(3) = 0.5255324099
            quad%directions(10)%xmu(3) = 0.1834346425
            quad%directions(11)%xmu(3) = 0.1834346425
            quad%directions(12)%xmu(3) = 0.1834346425
            quad%directions(1)%wmu = 0.0337428454
            quad%directions(2)%wmu = 0.0337428454
            quad%directions(3)%wmu = 0.0337428454
            quad%directions(4)%wmu = 0.0741270115
            quad%directions(5)%wmu = 0.0741270115
            quad%directions(6)%wmu = 0.0741270115
            quad%directions(7)%wmu = 0.1045688819
            quad%directions(8)%wmu = 0.1045688819
            quad%directions(9)%wmu = 0.1045688819
            quad%directions(10)%wmu = 0.1208945945
            quad%directions(11)%wmu = 0.1208945945
            quad%directions(12)%wmu = 0.1208945945
        else if(ns%state%sn == 10) then
            quad%directions(1)%xmu(1) = 0.2192163794
            quad%directions(2)%xmu(1) = 0.1604775276
            quad%directions(3)%xmu(1) = 0.0587388518
            quad%directions(4)%xmu(1) = 0.4845688685
            quad%directions(5)%xmu(1) = 0.3547290315
            quad%directions(6)%xmu(1) = 0.1298398370
            quad%directions(7)%xmu(1) = 0.7087570108
            quad%directions(8)%xmu(1) = 0.5188461421
            quad%directions(9)%xmu(1) = 0.1899108688
            quad%directions(10)%xmu(1) = 0.8704961020
            quad%directions(11)%xmu(1) = 0.6372473744
            quad%directions(12)%xmu(1) = 0.2332487274
            quad%directions(13)%xmu(1) = 0.9551616673
            quad%directions(14)%xmu(1) = 0.6992268699
            quad%directions(15)%xmu(1) = 0.2559347973
            quad%directions(1)%xmu(2) = 0.0587388518
            quad%directions(2)%xmu(2) = 0.1604775276
            quad%directions(3)%xmu(2) = 0.2192163794
            quad%directions(4)%xmu(2) = 0.1298398370
            quad%directions(5)%xmu(2) = 0.3547290315
            quad%directions(6)%xmu(2) = 0.4845688685
            quad%directions(7)%xmu(2) = 0.1899108688
            quad%directions(8)%xmu(2) = 0.5188461421
            quad%directions(9)%xmu(2) = 0.7087570108
            quad%directions(10)%xmu(2) = 0.2332487274
            quad%directions(11)%xmu(2) = 0.6372473744
            quad%directions(12)%xmu(2) = 0.8704961020
            quad%directions(13)%xmu(2) = 0.2559347973
            quad%directions(14)%xmu(2) = 0.6992268699
            quad%directions(15)%xmu(2) = 0.9551616673
            quad%directions(1)%xmu(3) = 0.9739065285
            quad%directions(2)%xmu(3) = 0.9739065285
            quad%directions(3)%xmu(3) = 0.9739065285
            quad%directions(4)%xmu(3) = 0.8650633667
            quad%directions(5)%xmu(3) = 0.8650633667
            quad%directions(6)%xmu(3) = 0.8650633667
            quad%directions(7)%xmu(3) = 0.6794095683
            quad%directions(8)%xmu(3) = 0.6794095683
            quad%directions(9)%xmu(3) = 0.6794095683
            quad%directions(10)%xmu(3) = 0.4333953941
            quad%directions(11)%xmu(3) = 0.4333953941
            quad%directions(12)%xmu(3) = 0.4333953941
            quad%directions(13)%xmu(3) = 0.1488743390
            quad%directions(14)%xmu(3) = 0.1488743390
            quad%directions(15)%xmu(3) = 0.1488743390
            quad%directions(1)%wmu = 0.0222237815
            quad%directions(2)%wmu = 0.0222237815
            quad%directions(3)%wmu = 0.0222237815
            quad%directions(4)%wmu = 0.0498171163
            quad%directions(5)%wmu = 0.0498171163
            quad%directions(6)%wmu = 0.0498171163
            quad%directions(7)%wmu = 0.0730287875
            quad%directions(8)%wmu = 0.0730287875
            quad%directions(9)%wmu = 0.0730287875
            quad%directions(10)%wmu = 0.0897555731
            quad%directions(11)%wmu = 0.0897555731
            quad%directions(12)%wmu = 0.0897555731
            quad%directions(13)%wmu = 0.0985080749
            quad%directions(14)%wmu = 0.0985080749
            quad%directions(15)%wmu = 0.0985080749
        else if(ns%state%sn == 12) then
            quad%directions(1)%xmu(1) = 0.1846377296
            quad%directions(2)%xmu(1) = 0.1351641990
            quad%directions(3)%xmu(1) = 0.0494735304
            quad%directions(4)%xmu(1) = 0.4127250690
            quad%directions(5)%xmu(1) = 0.3021357201
            quad%directions(6)%xmu(1) = 0.1105893489
            quad%directions(7)%xmu(1) = 0.6164165018
            quad%directions(8)%xmu(1) = 0.4512481980
            quad%directions(9)%xmu(1) = 0.1651683039
            quad%directions(10)%xmu(1) = 0.7817781530
            quad%directions(11)%xmu(1) = 0.5723013282
            quad%directions(12)%xmu(1) = 0.2094768248
            quad%directions(13)%xmu(1) = 0.8982071434
            quad%directions(14)%xmu(1) = 0.6575332647
            quad%directions(15)%xmu(1) = 0.2406738786
            quad%directions(16)%xmu(1) = 0.9583213889
            quad%directions(17)%xmu(1) = 0.7015399466
            quad%directions(18)%xmu(1) = 0.2567814421
            quad%directions(1)%xmu(2) = 0.0494735304
            quad%directions(2)%xmu(2) = 0.1351641990
            quad%directions(3)%xmu(2) = 0.1846377296
            quad%directions(4)%xmu(2) = 0.1105893489
            quad%directions(5)%xmu(2) = 0.3021357201
            quad%directions(6)%xmu(2) = 0.4127250690
            quad%directions(7)%xmu(2) = 0.1651683039
            quad%directions(8)%xmu(2) = 0.4512481980
            quad%directions(9)%xmu(2) = 0.6164165018
            quad%directions(10)%xmu(2) = 0.2094768248
            quad%directions(11)%xmu(2) = 0.5723013282
            quad%directions(12)%xmu(2) = 0.7817781530
            quad%directions(13)%xmu(2) = 0.2406738786
            quad%directions(14)%xmu(2) = 0.6575332647
            quad%directions(15)%xmu(2) = 0.8982071434
            quad%directions(16)%xmu(2) = 0.2567814421
            quad%directions(17)%xmu(2) = 0.7015399466
            quad%directions(18)%xmu(2) = 0.9583213889
            quad%directions(1)%xmu(3) = 0.9815606342
            quad%directions(2)%xmu(3) = 0.9815606342
            quad%directions(3)%xmu(3) = 0.9815606342
            quad%directions(4)%xmu(3) = 0.9041172564
            quad%directions(5)%xmu(3) = 0.9041172564
            quad%directions(6)%xmu(3) = 0.9041172564
            quad%directions(7)%xmu(3) = 0.7699026742
            quad%directions(8)%xmu(3) = 0.7699026742
            quad%directions(9)%xmu(3) = 0.7699026742
            quad%directions(10)%xmu(3) = 0.5873179543
            quad%directions(11)%xmu(3) = 0.5873179543
            quad%directions(12)%xmu(3) = 0.5873179543
            quad%directions(13)%xmu(3) = 0.3678314990
            quad%directions(14)%xmu(3) = 0.3678314990
            quad%directions(15)%xmu(3) = 0.3678314990
            quad%directions(16)%xmu(3) = 0.1252334085
            quad%directions(17)%xmu(3) = 0.1252334085
            quad%directions(18)%xmu(3) = 0.1252334085
            quad%directions(1)%wmu = 0.0157251121
            quad%directions(2)%wmu = 0.0157251121
            quad%directions(3)%wmu = 0.0157251121
            quad%directions(4)%wmu = 0.0356464420
            quad%directions(5)%wmu = 0.0356464420
            quad%directions(6)%wmu = 0.0356464420
            quad%directions(7)%wmu = 0.0533594429
            quad%directions(8)%wmu = 0.0533594429
            quad%directions(9)%wmu = 0.0533594429
            quad%directions(10)%wmu = 0.0677224755
            quad%directions(11)%wmu = 0.0677224755
            quad%directions(12)%wmu = 0.0677224755
            quad%directions(13)%wmu = 0.0778308455
            quad%directions(14)%wmu = 0.0778308455
            quad%directions(15)%wmu = 0.0778308455
            quad%directions(16)%wmu = 0.0830490153
            quad%directions(17)%wmu = 0.0830490153
            quad%directions(18)%wmu = 0.0830490153
        else if(ns%state%sn == 16) then
            quad%directions(1)%xmu(1) = 0.1402615760
            quad%directions(2)%xmu(1) = 0.1026786000
            quad%directions(3)%xmu(1) = 0.0375829761
            quad%directions(4)%xmu(1) = 0.3171092351
            quad%directions(5)%xmu(1) = 0.2321400716
            quad%directions(6)%xmu(1) = 0.0849691636
            quad%directions(7)%xmu(1) = 0.4836218252
            quad%directions(8)%xmu(1) = 0.3540357478
            quad%directions(9)%xmu(1) = 0.1295860776
            quad%directions(10)%xmu(1) = 0.6329314604
            quad%directions(11)%xmu(1) = 0.4633379867
            quad%directions(12)%xmu(1) = 0.1695934738
            quad%directions(13)%xmu(1) = 0.7594836797
            quad%directions(14)%xmu(1) = 0.5559806412
            quad%directions(15)%xmu(1) = 0.2035030387
            quad%directions(16)%xmu(1) = 0.8586535310
            quad%directions(17)%xmu(1) = 0.6285780108
            quad%directions(18)%xmu(1) = 0.2300755202
            quad%directions(19)%xmu(1) = 0.9268356232
            quad%directions(20)%xmu(1) = 0.6784907664
            quad%directions(21)%xmu(1) = 0.2483448567
            quad%directions(22)%xmu(1) = 0.9615560538
            quad%directions(23)%xmu(1) = 0.7039078857
            quad%directions(24)%xmu(1) = 0.2576481679
            quad%directions(1)%xmu(2) = 0.0375829761
            quad%directions(2)%xmu(2) = 0.1026786000
            quad%directions(3)%xmu(2) = 0.1402615760
            quad%directions(4)%xmu(2) = 0.0849691636
            quad%directions(5)%xmu(2) = 0.2321400716
            quad%directions(6)%xmu(2) = 0.3171092351
            quad%directions(7)%xmu(2) = 0.1295860776
            quad%directions(8)%xmu(2) = 0.3540357478
            quad%directions(9)%xmu(2) = 0.4836218252
            quad%directions(10)%xmu(2) = 0.1695934738   
            quad%directions(11)%xmu(2) = 0.4633379867   
            quad%directions(12)%xmu(2) = 0.6329314604   
            quad%directions(13)%xmu(2) = 0.2035030387   
            quad%directions(14)%xmu(2) = 0.5559806412   
            quad%directions(15)%xmu(2) = 0.7594836797   
            quad%directions(16)%xmu(2) = 0.2300755202   
            quad%directions(17)%xmu(2) = 0.6285780108   
            quad%directions(18)%xmu(2) = 0.8586535310   
            quad%directions(19)%xmu(2) = 0.2483448567   
            quad%directions(20)%xmu(2) = 0.6784907664   
            quad%directions(21)%xmu(2) = 0.9268356232   
            quad%directions(22)%xmu(2) = 0.2576481679   
            quad%directions(23)%xmu(2) = 0.7039078857   
            quad%directions(24)%xmu(2) = 0.9615560538
            quad%directions(1)%xmu(3) = 0.9894009350
            quad%directions(2)%xmu(3) = 0.9894009350
            quad%directions(3)%xmu(3) = 0.9894009350
            quad%directions(4)%xmu(3) = 0.9445750231
            quad%directions(5)%xmu(3) = 0.9445750231
            quad%directions(6)%xmu(3) = 0.9445750231
            quad%directions(7)%xmu(3) = 0.8656312024
            quad%directions(8)%xmu(3) = 0.8656312024
            quad%directions(9)%xmu(3) = 0.8656312024
            quad%directions(10)%xmu(3) = 0.7554044084
            quad%directions(11)%xmu(3) = 0.7554044084
            quad%directions(12)%xmu(3) = 0.7554044084
            quad%directions(13)%xmu(3) = 0.6178762444
            quad%directions(14)%xmu(3) = 0.6178762444
            quad%directions(15)%xmu(3) = 0.6178762444
            quad%directions(16)%xmu(3) = 0.4580167777
            quad%directions(17)%xmu(3) = 0.4580167777
            quad%directions(18)%xmu(3) = 0.4580167777
            quad%directions(19)%xmu(3) = 0.2816035508
            quad%directions(20)%xmu(3) = 0.2816035508
            quad%directions(21)%xmu(3) = 0.2816035508
            quad%directions(22)%xmu(3) = 0.0950125098
            quad%directions(23)%xmu(3) = 0.0950125098
            quad%directions(24)%xmu(3) = 0.0950125098
            quad%directions(1)%wmu = 0.0090508198
            quad%directions(2)%wmu = 0.0090508198
            quad%directions(3)%wmu = 0.0090508198
            quad%directions(4)%wmu = 0.0207511747
            quad%directions(5)%wmu = 0.0207511747
            quad%directions(6)%wmu = 0.0207511747
            quad%directions(7)%wmu = 0.0317195039
            quad%directions(8)%wmu = 0.0317195039
            quad%directions(9)%wmu = 0.0317195039
            quad%directions(10)%wmu = 0.0415429904
            quad%directions(11)%wmu = 0.0415429904
            quad%directions(12)%wmu = 0.0415429904
            quad%directions(13)%wmu = 0.0498653296
            quad%directions(14)%wmu = 0.0498653296
            quad%directions(15)%wmu = 0.0498653296
            quad%directions(16)%wmu = 0.0563855065
            quad%directions(17)%wmu = 0.0563855065
            quad%directions(18)%wmu = 0.0563855065
            quad%directions(19)%wmu = 0.0608678051
            quad%directions(20)%wmu = 0.0608678051
            quad%directions(21)%wmu = 0.0608678051
            quad%directions(22)%wmu = 0.0631502035
            quad%directions(23)%wmu = 0.0631502035
            quad%directions(24)%wmu = 0.0631502035
        end if           
    
    end subroutine Guass_chebyshev_symmetry_DNTR
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Guass_chebyshev_symmetry (quad, per_octant)
        
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        
        real(8), allocatable  :: x_1d(:)
        real(8), allocatable  :: w_1d(:)
        real(8)  :: theta
        
        integer  :: per_level
        integer  :: is, is_1d
        integer  :: i, j
        integer  :: i_allocate
        
        allocate(x_1d(ns%state%sn), stat=i_allocate)
        allocate(w_1d(ns%state%sn), stat=i_allocate)
        
        call Gauss_Ledrengre (-1.0, 1.0, x_1d, w_1d)
        
        ! assign to 3D
        is = 0
        do i = 1, ns%state%sn/2
            per_level = 3
            
            do j = 1, per_level
                is = is + 1
!                is_1d = i + ns%state%sn/2
                is_1d = ns%state%sn + 1 - i
                
                quad%directions(is)%wmu = w_1d(is_1d) / per_level
                quad%directions(is)%xmu(3) = x_1d(is_1d)
                
                theta = (PI/2.0) * ((2*j-1)/6.0)
                quad%directions(is)%xmu(1) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * COS(theta)
                quad%directions(is)%xmu(2) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * SIN(theta)
            end do
        end do
        
        if (allocated(x_1d))        deallocate(x_1d)
        if (allocated(w_1d))        deallocate(w_1d)
    
    end subroutine Guass_chebyshev_symmetry
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Level_symmetry_DNTR (quad, per_octant)
        
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        
        if(ns%state%sn == 2) then
            quad%directions(1)%xmu(1) = 0.5773503
            quad%directions(1)%xmu(2) = 0.5773503
            quad%directions(1)%xmu(3) = 0.5773503
            quad%directions(1)%wmu = 1.0
        else if(ns%state%sn == 4) then
            quad%directions(1)%xmu(1) = 0.3500212
            quad%directions(2)%xmu(1) = 0.3500212
            quad%directions(3)%xmu(1) = 0.8688903
            quad%directions(1)%xmu(2) = 0.3500212
            quad%directions(2)%xmu(2) = 0.8688903
            quad%directions(3)%xmu(2) = 0.3500212
            quad%directions(1)%xmu(3) = 0.8688903
            quad%directions(2)%xmu(3) = 0.3500212
            quad%directions(3)%xmu(3) = 0.3500213
            quad%directions(1)%wmu = 0.3333333
            quad%directions(2)%wmu = 0.3333333
            quad%directions(3)%wmu = 0.3333333
        else if(ns%state%sn == 6) then
            quad%directions(1)%xmu(1) = 0.266636
            quad%directions(2)%xmu(1) = 0.266636
            quad%directions(3)%xmu(1) = 0.266636
            quad%directions(4)%xmu(1) = 0.681508
            quad%directions(5)%xmu(1) = 0.681508
            quad%directions(6)%xmu(1) = 0.926181
            quad%directions(1)%xmu(2) = 0.266636
            quad%directions(2)%xmu(2) = 0.681508
            quad%directions(3)%xmu(2) = 0.926181
            quad%directions(4)%xmu(2) = 0.266636
            quad%directions(5)%xmu(2) = 0.681508
            quad%directions(6)%xmu(2) = 0.266636
            quad%directions(1)%xmu(3) = 0.926181
            quad%directions(2)%xmu(3) = 0.681508
            quad%directions(3)%xmu(3) = 0.266636
            quad%directions(4)%xmu(3) = 0.681508
            quad%directions(5)%xmu(3) = 0.266636
            quad%directions(6)%xmu(3) = 0.266636
            quad%directions(1)%wmu = 0.176126
            quad%directions(2)%wmu = 0.157207
            quad%directions(3)%wmu = 0.176126
            quad%directions(4)%wmu = 0.157207
            quad%directions(5)%wmu = 0.157207
            quad%directions(6)%wmu = 0.176126
        else if(ns%state%sn == 8) then
            quad%directions(1)%xmu(1) = 0.218218
            quad%directions(2)%xmu(1) = 0.218218
            quad%directions(3)%xmu(1) = 0.218218
            quad%directions(4)%xmu(1) = 0.218218
            quad%directions(5)%xmu(1) = 0.577350
            quad%directions(6)%xmu(1) = 0.577350
            quad%directions(7)%xmu(1) = 0.577350
            quad%directions(8)%xmu(1) = 0.786796
            quad%directions(9)%xmu(1) = 0.786796
            quad%directions(10)%xmu(1) = 0.951190
            quad%directions(1)%xmu(2) = 0.218218
            quad%directions(2)%xmu(2) = 0.577350
            quad%directions(3)%xmu(2) = 0.786796
            quad%directions(4)%xmu(2) = 0.951190
            quad%directions(5)%xmu(2) = 0.218218
            quad%directions(6)%xmu(2) = 0.577350
            quad%directions(7)%xmu(2) = 0.786796
            quad%directions(8)%xmu(2) = 0.218218
            quad%directions(9)%xmu(2) = 0.577350
            quad%directions(10)%xmu(2) = 0.218218
            quad%directions(1)%xmu(3) = 0.951190
            quad%directions(2)%xmu(3) = 0.786796
            quad%directions(3)%xmu(3) = 0.577350
            quad%directions(4)%xmu(3) = 0.218218
            quad%directions(5)%xmu(3) = 0.786796
            quad%directions(6)%xmu(3) = 0.577350
            quad%directions(7)%xmu(3) = 0.218218
            quad%directions(8)%xmu(3) = 0.577350
            quad%directions(9)%xmu(3) = 0.218218
            quad%directions(10)%xmu(3) = 0.218218
            quad%directions(1)%wmu = 0.1209877
            quad%directions(2)%wmu = 0.0907407
            quad%directions(3)%wmu = 0.0907407
            quad%directions(4)%wmu = 0.1209877
            quad%directions(5)%wmu = 0.0907407
            quad%directions(6)%wmu = 0.0925927
            quad%directions(7)%wmu = 0.0907407
            quad%directions(8)%wmu = 0.0907407
            quad%directions(9)%wmu = 0.0907407
            quad%directions(10)%wmu = 0.1209877
        else if(ns%state%sn == 10) then
            quad%directions(1)%xmu(1) = 0.96225  
            quad%directions(2)%xmu(1) = 0.83887  
            quad%directions(3)%xmu(1) = 0.69389  
            quad%directions(4)%xmu(1) = 0.50918  
            quad%directions(5)%xmu(1) = 0.19245  
            quad%directions(6)%xmu(1) = 0.83887  
            quad%directions(7)%xmu(1) = 0.69389  
            quad%directions(8)%xmu(1) = 0.50918  
            quad%directions(9)%xmu(1) = 0.19245  
            quad%directions(10)%xmu(1) = 0.69389 
            quad%directions(11)%xmu(1) = 0.50918 
            quad%directions(12)%xmu(1) = 0.19245 
            quad%directions(13)%xmu(1) = 0.50918 
            quad%directions(14)%xmu(1) = 0.19245 
            quad%directions(15)%xmu(1) = 0.19245 
            quad%directions(1)%xmu(2) = 0.19245 
            quad%directions(2)%xmu(2) = 0.19245 
            quad%directions(3)%xmu(2) = 0.19245 
            quad%directions(4)%xmu(2) = 0.19245 
            quad%directions(5)%xmu(2) = 0.19245 
            quad%directions(6)%xmu(2) = 0.50918 
            quad%directions(7)%xmu(2) = 0.50918 
            quad%directions(8)%xmu(2) = 0.50918 
            quad%directions(9)%xmu(2) = 0.50918 
            quad%directions(10)%xmu(2) = 0.69389
            quad%directions(11)%xmu(2) = 0.69389
            quad%directions(12)%xmu(2) = 0.69389
            quad%directions(13)%xmu(2) = 0.83887
            quad%directions(14)%xmu(2) = 0.83887
            quad%directions(15)%xmu(2) = 0.96225
            quad%directions(1)%xmu(3) = 0.19245 
            quad%directions(2)%xmu(3) = 0.50918 
            quad%directions(3)%xmu(3) = 0.69389 
            quad%directions(4)%xmu(3) = 0.83887 
            quad%directions(5)%xmu(3) = 0.96225 
            quad%directions(6)%xmu(3) = 0.19245 
            quad%directions(7)%xmu(3) = 0.50918 
            quad%directions(8)%xmu(3) = 0.69389 
            quad%directions(9)%xmu(3) = 0.83887 
            quad%directions(10)%xmu(3) = 0.19245
            quad%directions(11)%xmu(3) = 0.50918
            quad%directions(12)%xmu(3) = 0.69389
            quad%directions(13)%xmu(3) = 0.19245
            quad%directions(14)%xmu(3) = 0.50918
            quad%directions(15)%xmu(3) = 0.19245
            quad%directions(1)%wmu = 0.092004 
            quad%directions(2)%wmu = 0.069872 
            quad%directions(3)%wmu = 0.049384 
            quad%directions(4)%wmu = 0.069872 
            quad%directions(5)%wmu = 0.092004 
            quad%directions(6)%wmu = 0.069872 
            quad%directions(7)%wmu = 0.052204 
            quad%directions(8)%wmu = 0.052204 
            quad%directions(9)%wmu = 0.069872 
            quad%directions(10)%wmu = 0.049384
            quad%directions(11)%wmu = 0.052204
            quad%directions(12)%wmu = 0.049384
            quad%directions(13)%wmu = 0.069872
            quad%directions(14)%wmu = 0.069872
            quad%directions(15)%wmu = 0.092004
        else if(ns%state%sn == 12) then
            quad%directions(1)%xmu(1) = 0.96922  
            quad%directions(2)%xmu(1) = 0.87039 
            quad%directions(3)%xmu(1) = 0.75879 
            quad%directions(4)%xmu(1) = 0.62765 
            quad%directions(5)%xmu(1) = 0.46057 
            quad%directions(6)%xmu(1) = 0.17408 
            quad%directions(7)%xmu(1) = 0.87039 
            quad%directions(8)%xmu(1) = 0.75879 
            quad%directions(9)%xmu(1) = 0.62765 
            quad%directions(10)%xmu(1) = 0.46057
            quad%directions(11)%xmu(1) = 0.17408
            quad%directions(12)%xmu(1) = 0.75879
            quad%directions(13)%xmu(1) = 0.62765
            quad%directions(14)%xmu(1) = 0.46057
            quad%directions(15)%xmu(1) = 0.17408
            quad%directions(16)%xmu(1) = 0.62765
            quad%directions(17)%xmu(1) = 0.46057
            quad%directions(18)%xmu(1) = 0.17408
            quad%directions(19)%xmu(1) = 0.46057
            quad%directions(20)%xmu(1) = 0.17408
            quad%directions(21)%xmu(1) = 0.17408
            quad%directions(1)%xmu(2) = 0.17408 
            quad%directions(2)%xmu(2) = 0.17408 
            quad%directions(3)%xmu(2) = 0.17408 
            quad%directions(4)%xmu(2) = 0.17408 
            quad%directions(5)%xmu(2) = 0.17408 
            quad%directions(6)%xmu(2) = 0.17408 
            quad%directions(7)%xmu(2) = 0.46057 
            quad%directions(8)%xmu(2) = 0.46057 
            quad%directions(9)%xmu(2) = 0.46057 
            quad%directions(10)%xmu(2) = 0.46057
            quad%directions(11)%xmu(2) = 0.46057
            quad%directions(12)%xmu(2) = 0.62765
            quad%directions(13)%xmu(2) = 0.62765
            quad%directions(14)%xmu(2) = 0.62765
            quad%directions(15)%xmu(2) = 0.62765
            quad%directions(16)%xmu(2) = 0.75879
            quad%directions(17)%xmu(2) = 0.75879
            quad%directions(18)%xmu(2) = 0.75879
            quad%directions(19)%xmu(2) = 0.87039
            quad%directions(20)%xmu(2) = 0.87039
            quad%directions(21)%xmu(2) = 0.96922
            quad%directions(1)%xmu(3) = 0.17408 
            quad%directions(2)%xmu(3) = 0.46057 
            quad%directions(3)%xmu(3) = 0.62765 
            quad%directions(4)%xmu(3) = 0.75879 
            quad%directions(5)%xmu(3) = 0.87039 
            quad%directions(6)%xmu(3) = 0.96922 
            quad%directions(7)%xmu(3) = 0.17408 
            quad%directions(8)%xmu(3) = 0.46057 
            quad%directions(9)%xmu(3) = 0.62765 
            quad%directions(10)%xmu(3) = 0.75879
            quad%directions(11)%xmu(3) = 0.87039
            quad%directions(12)%xmu(3) = 0.17408
            quad%directions(13)%xmu(3) = 0.46057
            quad%directions(14)%xmu(3) = 0.62765
            quad%directions(15)%xmu(3) = 0.75879
            quad%directions(16)%xmu(3) = 0.17408
            quad%directions(17)%xmu(3) = 0.46057
            quad%directions(18)%xmu(3) = 0.62765
            quad%directions(19)%xmu(3) = 0.17408
            quad%directions(20)%xmu(3) = 0.46057
            quad%directions(21)%xmu(3) = 0.17408
            quad%directions(1)%wmu = 0.0745720 
            quad%directions(2)%wmu = 0.0550400  
            quad%directions(3)%wmu = 0.0392124 
            quad%directions(4)%wmu = 0.0392124  
            quad%directions(5)%wmu = 0.0550400
            quad%directions(6)%wmu = 0.0745720
            quad%directions(7)%wmu = 0.0550400
            quad%directions(8)%wmu = 0.0395648  
            quad%directions(9)%wmu = 0.0306880   
            quad%directions(10)%wmu = 0.0395648 
            quad%directions(11)%wmu = 0.0550400
            quad%directions(12)%wmu = 0.0392124 
            quad%directions(13)%wmu = 0.0306880
            quad%directions(14)%wmu = 0.0306880
            quad%directions(15)%wmu = 0.0392124 
            quad%directions(16)%wmu = 0.0392124 
            quad%directions(17)%wmu = 0.0395648 
            quad%directions(18)%wmu = 0.0392124 
            quad%directions(19)%wmu = 0.0550400
            quad%directions(20)%wmu = 0.0550400
            quad%directions(21)%wmu = 0.0745720
        else if(ns%state%sn == 14) then
            quad%directions(1)%xmu(1) = 0.97402 
            quad%directions(2)%xmu(1) = 0.89156 
            quad%directions(3)%xmu(1) = 0.80064 
            quad%directions(4)%xmu(1) = 0.69798 
            quad%directions(5)%xmu(1) = 0.57735 
            quad%directions(6)%xmu(1) = 0.42366 
            quad%directions(7)%xmu(1) = 0.16013 
            quad%directions(8)%xmu(1) = 0.89156 
            quad%directions(9)%xmu(1) = 0.80064 
            quad%directions(10)%xmu(1) = 0.69798
            quad%directions(11)%xmu(1) = 0.57735
            quad%directions(12)%xmu(1) = 0.42366
            quad%directions(13)%xmu(1) = 0.16013
            quad%directions(14)%xmu(1) = 0.80064
            quad%directions(15)%xmu(1) = 0.69798
            quad%directions(16)%xmu(1) = 0.57735
            quad%directions(17)%xmu(1) = 0.42366
            quad%directions(18)%xmu(1) = 0.16013
            quad%directions(19)%xmu(1) = 0.69798
            quad%directions(20)%xmu(1) = 0.57735
            quad%directions(21)%xmu(1) = 0.42366
            quad%directions(22)%xmu(1) = 0.16013
            quad%directions(23)%xmu(1) = 0.57735
            quad%directions(24)%xmu(1) = 0.42366
            quad%directions(25)%xmu(1) = 0.16013
            quad%directions(26)%xmu(1) = 0.42366
            quad%directions(27)%xmu(1) = 0.16013
            quad%directions(28)%xmu(1) = 0.16013
            quad%directions(1)%xmu(2) = 0.16013 
            quad%directions(2)%xmu(2) = 0.16013 
            quad%directions(3)%xmu(2) = 0.16013 
            quad%directions(4)%xmu(2) = 0.16013 
            quad%directions(5)%xmu(2) = 0.16013 
            quad%directions(6)%xmu(2) = 0.16013 
            quad%directions(7)%xmu(2) = 0.16013 
            quad%directions(8)%xmu(2) = 0.42366 
            quad%directions(9)%xmu(2) = 0.42366 
            quad%directions(10)%xmu(2) = 0.42366
            quad%directions(11)%xmu(2) = 0.42366
            quad%directions(12)%xmu(2) = 0.42366
            quad%directions(13)%xmu(2) = 0.42366
            quad%directions(14)%xmu(2) = 0.57735
            quad%directions(15)%xmu(2) = 0.57735
            quad%directions(16)%xmu(2) = 0.57735
            quad%directions(17)%xmu(2) = 0.57735
            quad%directions(18)%xmu(2) = 0.57735
            quad%directions(19)%xmu(2) = 0.69798
            quad%directions(20)%xmu(2) = 0.69798
            quad%directions(21)%xmu(2) = 0.69798
            quad%directions(22)%xmu(2) = 0.69798
            quad%directions(23)%xmu(2) = 0.80064
            quad%directions(24)%xmu(2) = 0.80064
            quad%directions(25)%xmu(2) = 0.80064
            quad%directions(26)%xmu(2) = 0.89156
            quad%directions(27)%xmu(2) = 0.89156
            quad%directions(28)%xmu(2) = 0.97402
            quad%directions(1)%xmu(3) = 0.16013 
            quad%directions(2)%xmu(3) = 0.42366 
            quad%directions(3)%xmu(3) = 0.57735 
            quad%directions(4)%xmu(3) = 0.69798 
            quad%directions(5)%xmu(3) = 0.80064 
            quad%directions(6)%xmu(3) = 0.89156 
            quad%directions(7)%xmu(3) = 0.97402 
            quad%directions(8)%xmu(3) = 0.16013 
            quad%directions(9)%xmu(3) = 0.42366 
            quad%directions(10)%xmu(3) = 0.57735
            quad%directions(11)%xmu(3) = 0.69798
            quad%directions(12)%xmu(3) = 0.80064
            quad%directions(13)%xmu(3) = 0.89156
            quad%directions(14)%xmu(3) = 0.16013
            quad%directions(15)%xmu(3) = 0.42366
            quad%directions(16)%xmu(3) = 0.57735
            quad%directions(17)%xmu(3) = 0.69798
            quad%directions(18)%xmu(3) = 0.80064
            quad%directions(19)%xmu(3) = 0.16013
            quad%directions(20)%xmu(3) = 0.42366
            quad%directions(21)%xmu(3) = 0.57735
            quad%directions(22)%xmu(3) = 0.69798
            quad%directions(23)%xmu(3) = 0.16013
            quad%directions(24)%xmu(3) = 0.42366
            quad%directions(25)%xmu(3) = 0.57735
            quad%directions(26)%xmu(3) = 0.16013
            quad%directions(27)%xmu(3) = 0.42366
            quad%directions(28)%xmu(3) = 0.16013
            quad%directions(1)%wmu = 0.064028  
            quad%directions(2)%wmu = 0.041780   
            quad%directions(3)%wmu = 0.039020   
            quad%directions(4)%wmu = 0.0220072 
            quad%directions(5)%wmu = 0.039020   
            quad%directions(6)%wmu = 0.041780   
            quad%directions(7)%wmu = 0.064028  
            quad%directions(8)%wmu = 0.041780   
            quad%directions(9)%wmu = 0.0340836 
            quad%directions(10)%wmu = 0.022066 
            quad%directions(11)%wmu = 0.022066 
            quad%directions(12)%wmu = 0.0340836
            quad%directions(13)%wmu = 0.041780  
            quad%directions(14)%wmu = 0.039020  
            quad%directions(15)%wmu = 0.022066 
            quad%directions(16)%wmu = 0.0224612
            quad%directions(17)%wmu = 0.022066 
            quad%directions(18)%wmu = 0.039020  
            quad%directions(19)%wmu = 0.0220072
            quad%directions(20)%wmu = 0.022066 
            quad%directions(21)%wmu = 0.022066 
            quad%directions(22)%wmu = 0.0220072
            quad%directions(23)%wmu = 0.039020  
            quad%directions(24)%wmu = 0.0340836
            quad%directions(25)%wmu = 0.039020  
            quad%directions(26)%wmu = 0.041780  
            quad%directions(27)%wmu = 0.041780  
            quad%directions(28)%wmu = 0.064028
        else if(ns%state%sn == 16) then
            quad%directions(1)%xmu(1) = 0.97753      
            quad%directions(2)%xmu(1) = 0.90676      
            quad%directions(3)%xmu(1) = 0.82999      
            quad%directions(4)%xmu(1) = 0.74536      
            quad%directions(5)%xmu(1) = 0.64979      
            quad%directions(6)%xmu(1) = 0.53748      
            quad%directions(7)%xmu(1) = 0.39441      
            quad%directions(8)%xmu(1) = 0.14907      
            quad%directions(9)%xmu(1) = 0.90676      
            quad%directions(10)%xmu(1) = 0.82999     
            quad%directions(11)%xmu(1) = 0.74536     
            quad%directions(12)%xmu(1) = 0.64979     
            quad%directions(13)%xmu(1) = 0.53748     
            quad%directions(14)%xmu(1) = 0.39441     
            quad%directions(15)%xmu(1) = 0.14907     
            quad%directions(16)%xmu(1) = 0.82999     
            quad%directions(17)%xmu(1) = 0.74536     
            quad%directions(18)%xmu(1) = 0.64979     
            quad%directions(19)%xmu(1) = 0.53748     
            quad%directions(20)%xmu(1) = 0.39441     
            quad%directions(21)%xmu(1) = 0.14907     
            quad%directions(22)%xmu(1) = 0.74536     
            quad%directions(23)%xmu(1) = 0.64979     
            quad%directions(24)%xmu(1) = 0.53748     
            quad%directions(25)%xmu(1) = 0.39441     
            quad%directions(26)%xmu(1) = 0.14907     
            quad%directions(27)%xmu(1) = 0.64979     
            quad%directions(28)%xmu(1) = 0.53748     
            quad%directions(29)%xmu(1) = 0.39441     
            quad%directions(30)%xmu(1) = 0.14907     
            quad%directions(31)%xmu(1) = 0.53748     
            quad%directions(32)%xmu(1) = 0.39441     
            quad%directions(33)%xmu(1) = 0.14907     
            quad%directions(34)%xmu(1) = 0.39441     
            quad%directions(35)%xmu(1) = 0.14907     
            quad%directions(36)%xmu(1) = 0.14907     
            quad%directions(1)%xmu(2) = 0.14907      
            quad%directions(2)%xmu(2) = 0.14907      
            quad%directions(3)%xmu(2) = 0.14907      
            quad%directions(4)%xmu(2) = 0.14907      
            quad%directions(5)%xmu(2) = 0.14907      
            quad%directions(6)%xmu(2) = 0.14907      
            quad%directions(7)%xmu(2) = 0.14907      
            quad%directions(8)%xmu(2) = 0.14907      
            quad%directions(9)%xmu(2) = 0.39441      
            quad%directions(10)%xmu(2) = 0.39441     
            quad%directions(11)%xmu(2) = 0.39441     
            quad%directions(12)%xmu(2) = 0.39441     
            quad%directions(13)%xmu(2) = 0.39441     
            quad%directions(14)%xmu(2) = 0.39441     
            quad%directions(15)%xmu(2) = 0.39441     
            quad%directions(16)%xmu(2) = 0.53748     
            quad%directions(17)%xmu(2) = 0.53748     
            quad%directions(18)%xmu(2) = 0.53748     
            quad%directions(19)%xmu(2) = 0.53748     
            quad%directions(20)%xmu(2) = 0.53748     
            quad%directions(21)%xmu(2) = 0.53748     
            quad%directions(22)%xmu(2) = 0.64979     
            quad%directions(23)%xmu(2) = 0.64979     
            quad%directions(24)%xmu(2) = 0.64979     
            quad%directions(25)%xmu(2) = 0.64979     
            quad%directions(26)%xmu(2) = 0.64979     
            quad%directions(27)%xmu(2) = 0.74536     
            quad%directions(28)%xmu(2) = 0.74536     
            quad%directions(29)%xmu(2) = 0.74536     
            quad%directions(30)%xmu(2) = 0.74536     
            quad%directions(31)%xmu(2) = 0.82999     
            quad%directions(32)%xmu(2) = 0.82999     
            quad%directions(33)%xmu(2) = 0.82999     
            quad%directions(34)%xmu(2) = 0.90676     
            quad%directions(35)%xmu(2) = 0.90676     
            quad%directions(36)%xmu(2) = 0.97753     
            quad%directions(1)%xmu(3) = 0.14907      
            quad%directions(2)%xmu(3) = 0.39441      
            quad%directions(3)%xmu(3) = 0.53748      
            quad%directions(4)%xmu(3) = 0.64979      
            quad%directions(5)%xmu(3) = 0.74536      
            quad%directions(6)%xmu(3) = 0.82999      
            quad%directions(7)%xmu(3) = 0.90676      
            quad%directions(8)%xmu(3) = 0.97753      
            quad%directions(9)%xmu(3) = 0.14907      
            quad%directions(10)%xmu(3) = 0.39441     
            quad%directions(11)%xmu(3) = 0.53748     
            quad%directions(12)%xmu(3) = 0.64979     
            quad%directions(13)%xmu(3) = 0.74536     
            quad%directions(14)%xmu(3) = 0.82999     
            quad%directions(15)%xmu(3) = 0.90676     
            quad%directions(16)%xmu(3) = 0.14907     
            quad%directions(17)%xmu(3) = 0.39441     
            quad%directions(18)%xmu(3) = 0.53748     
            quad%directions(19)%xmu(3) = 0.64979     
            quad%directions(20)%xmu(3) = 0.74536     
            quad%directions(21)%xmu(3) = 0.82999     
            quad%directions(22)%xmu(3) = 0.14907     
            quad%directions(23)%xmu(3) = 0.39441     
            quad%directions(24)%xmu(3) = 0.53748     
            quad%directions(25)%xmu(3) = 0.64979     
            quad%directions(26)%xmu(3) = 0.74536     
            quad%directions(27)%xmu(3) = 0.14907     
            quad%directions(28)%xmu(3) = 0.39441     
            quad%directions(29)%xmu(3) = 0.53748     
            quad%directions(30)%xmu(3) = 0.64979     
            quad%directions(31)%xmu(3) = 0.14907     
            quad%directions(32)%xmu(3) = 0.39441     
            quad%directions(33)%xmu(3) = 0.53748     
            quad%directions(34)%xmu(3) = 0.14907     
            quad%directions(35)%xmu(3) = 0.39441     
            quad%directions(36)%xmu(3) = 0.14907     
            quad%directions(1)%wmu = 0.054344       
            quad%directions(2)%wmu = 0.0390724      
            quad%directions(3)%wmu = 0.0258952      
            quad%directions(4)%wmu = 0.0258536      
            quad%directions(5)%wmu = 0.0258536      
            quad%directions(6)%wmu = 0.0258952      
            quad%directions(7)%wmu = 0.0390724      
            quad%directions(8)%wmu = 0.054344       
            quad%directions(9)%wmu = 0.0390724      
            quad%directions(10)%wmu = 0.020156      
            quad%directions(11)%wmu = 0.0284496     
            quad%directions(12)%wmu = 0.0057524     
            quad%directions(13)%wmu = 0.0284496     
            quad%directions(14)%wmu = 0.020156      
            quad%directions(15)%wmu = 0.0390724     
            quad%directions(16)%wmu = 0.0258952     
            quad%directions(17)%wmu = 0.0284496     
            quad%directions(18)%wmu = 0.0145368     
            quad%directions(19)%wmu = 0.0145368     
            quad%directions(20)%wmu = 0.0284496     
            quad%directions(21)%wmu = 0.0258952     
            quad%directions(22)%wmu = 0.0258536     
            quad%directions(23)%wmu = 0.0057524     
            quad%directions(24)%wmu = 0.0145368     
            quad%directions(25)%wmu = 0.0057524     
            quad%directions(26)%wmu = 0.0258536     
            quad%directions(27)%wmu = 0.0258536     
            quad%directions(28)%wmu = 0.0284496     
            quad%directions(29)%wmu = 0.0284496     
            quad%directions(30)%wmu = 0.0258536     
            quad%directions(31)%wmu = 0.0258952     
            quad%directions(32)%wmu = 0.020156      
            quad%directions(33)%wmu = 0.0258952     
            quad%directions(34)%wmu = 0.0390724     
            quad%directions(35)%wmu = 0.0390724     
            quad%directions(36)%wmu = 0.054344      
        end if

    end subroutine Level_symmetry_DNTR
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Level_symmetry_Hebert (quad, per_octant)
    
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        
        real(8)  :: level(10, 10)  = 0.0
        real(8)  :: weight(10, 12) = 0.0
        
        integer  :: weight_index, x_index, y_index, z_index
        integer  :: per_level
        integer  :: row                                                         ! which row to be used
        integer  :: i, j, is
        real(8)  :: add

        ! error, excess limitation
        if (ns%state%sn > 20)  then
            
        end if
        
        level(1:10, 1)   = [ 0.5773503, 0.3500212, 0.2666354, 0.2182178, 0.1893227, 0.1672308, 0.1519951, 0.1389747, 0.1300795, 0.1206339 ]
        weight(1,  1:1)  = [ 1.0000000 ]
        weight(2,  1:1)  = [ 0.3333333 ]
        weight(3,  1:2)  = [ 0.1761262, 0.1572071 ]
        weight(4,  1:3)  = [ 0.1209876, 0.0907408, 0.0925925 ]
        weight(5,  1:4)  = [ 0.0893043, 0.0725281, 0.0450455, 0.0539274 ]
        weight(6,  1:5)  = [ 0.0707734, 0.0558760, 0.0373436, 0.0502654, 0.0258553 ]
        weight(7,  1:7)  = [ 0.0580031, 0.0488943, 0.0228095, 0.0393955, 0.0380920, 0.0258382, 0.0082759 ]
        weight(8,  1:8)  = [ 0.0489967, 0.0413235, 0.0203158, 0.0265468, 0.0378883, 0.0135404, 0.0326129, 0.0103825 ]
        weight(9,  1:10) = [ 0.0426910, 0.0370806, 0.0139198, 0.0297556, 0.0100159, 0.0306093, 0.0160431, 0.0197011, 0.0011939, 0.0158226 ]
        weight(10, 1:12) = [ 0.0370368, 0.0332667, 0.0112149, 0.0244882, 0.0136043, 0.0317999, 0.0068205, 0.0307489, 0.0000038, 0.0056654, 0.0045702, 0.0281998 ]
        
        row = ns%state%sn / 2
        
        do i = 2, ns%state%sn/2
            add = (i-1.0)*2.0*(1.0-3*level(row, 1)**2) / (ns%state%sn-2.0)
            level(row, i) = SQRT(level(row, 1)**2 + add)
        end do
        
        ! assign to 3D
        is = 0
        z_index = 0
        do i = 1, ns%state%sn / 2
            per_level = ns%state%sn/2 + 1 - i
            z_index = z_index + 1
            
            x_index = 0
            y_index = per_level + 1
            do j = 1, per_level
                is = is + 1
                x_index = x_index + 1
                y_index = y_index - 1
                
                call Get_weight_index (x_index, y_index, z_index, weight_index)
                quad%directions(is)%wmu = weight(row, weight_index)
                quad%directions(is)%xmu(1) = level(row, x_index)
                quad%directions(is)%xmu(2) = level(row, y_index)
                quad%directions(is)%xmu(3) = level(row, z_index)
            end do
        end do
        
    end subroutine Level_symmetry_Hebert
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Guass_Equal_PENTRAN (quad, per_octant)
        
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        
        real(8), allocatable  :: x_1d(:)
        real(8), allocatable  :: w_1d(:)
        real(8)  :: theta
        
        integer  :: per_level
        integer  :: is, is_1d
        integer  :: i, j
        integer  :: i_allocate
        
        allocate(x_1d(ns%state%sn), stat=i_allocate)
        allocate(w_1d(ns%state%sn), stat=i_allocate)
        
        call Gauss_Ledrengre (-1.0, 1.0, x_1d, w_1d)
        
        ! assign to 3D
        is = 0
        do i = 1, ns%state%sn / 2
            per_level = ns%state%sn/2 + 1 - i
            
            theta = (PI/2) / (2*per_level)
            do j = 1, per_level
                is = is + 1
                is_1d = i + ns%state%sn/2
                
                quad%directions(is)%wmu = w_1d(is_1d) / per_level
                quad%directions(is)%xmu(3) = x_1d(is_1d)
                quad%directions(is)%xmu(1) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * COS(theta)
                quad%directions(is)%xmu(2) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * SIN(theta)
                
                theta = theta + (PI/2)/(per_level)
            end do
        end do
        
        if (allocated(x_1d))        deallocate(x_1d)
        if (allocated(w_1d))        deallocate(w_1d)
    
    end subroutine Guass_Equal_PENTRAN
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Guass_ChebyShev_PENTRAN (quad, per_octant)
        
        type(QuadratureSet), intent(in out)  :: quad
        integer, intent(in)  :: per_octant
        
        real(8), allocatable  :: x_1d(:)
        real(8), allocatable  :: w_1d(:)
        real(8)  :: theta
        
        integer  :: per_level
        integer  :: is, is_1d
        integer  :: i, j
        integer  :: i_allocate
        
        allocate(x_1d(ns%state%sn), stat=i_allocate)
        allocate(w_1d(ns%state%sn), stat=i_allocate)
        
        call Gauss_Ledrengre (-1.0, 1.0, x_1d, w_1d)
        
        ! assign to 3D
        is = 0
        do i = 1, ns%state%sn / 2
            per_level = ns%state%sn/2 + 1 - i
            
            do j = 1, per_level
                is = is + 1
                is_1d = i + ns%state%sn/2
                
                quad%directions(is)%wmu = w_1d(is_1d) / per_level
                quad%directions(is)%xmu(3) = x_1d(is_1d)
                
                theta = (PI/2) * (2*j-1) / (2*per_level)
                quad%directions(is)%xmu(1) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * COS(theta)
                quad%directions(is)%xmu(2) = SQRT(1.0 - quad%directions(is)%xmu(3)**2) * SIN(theta)
            end do
        end do
        
        if (allocated(x_1d))        deallocate(x_1d)
        if (allocated(w_1d))        deallocate(w_1d)
    
    end subroutine Guass_ChebyShev_PENTRAN
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! given the lower and upper limits of integration x1 and x2, this routine returns arrays x and w
    ! of length n containing the abscissas and weights of the gauss-legendre n-point quadrature
    ! formula. the parameter eps is the relative precision. note that internal computations are
    ! done in double precision.
    !===============================================================================================
    subroutine Gauss_Ledrengre (x1, x2, x, w)
        
        real(4), intent(in) :: x1, x2
        real(8), dimension(:), intent(in out) :: x, w
        
        integer, parameter :: maxit = 10
        real(8), parameter :: eps = 3.0D-14
        integer :: its, j, m, n
        real(8) :: xl, xm
        real(8), dimension((SIZE(x)+1)/2) :: p1, p2, p3, pp, z, z1, arth
        logical, dimension((SIZE(x)+1)/2) :: unfinished
        
        n = SIZE(x)
        
        ! the roots are symmetric in the interval, so we only have to find half of them.
        m = (n + 1) / 2
        xm = 0.5 * (x2 + x1)    
        xl = 0.5 * (x2 - x1)
        
        ! newton¡¯s method carried out simultaneously on the roots.
        ! initial approximations to the roots.
        do j = 1, m
            arth(j) = j
        end do
        
        z = COS(PI * (arth-0.25) / (n+0.5)) 
        unfinished = .TRUE.
        do its = 1, maxit          
            where(unfinished)       
                p1 = 1.0
                p2 = 0.0
            end where
            
            ! loop up the recurrence relation to get the legendre polynomials evaluated at z.
            do j = 1, n                
                where (unfinished)
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3) / j
                end where
            end do
            
            ! p1 now contains the desired legendre polynomials. we next compute pp, the derivatives,
            ! by a standard relation involving also p2, the polynomials of one lower order.
            where (unfinished)
                pp = n * (z*p1-p2) / (z*z-1.0)
                z1 = z
                ! newton¡¯s method.
                z = z1 - p1 / pp                                                  
                unfinished = (ABS(z-z1) > eps)
            end where
            if (.NOT. ANY(unfinished)) exit
        end do
        
        if (its == maxit+1)  then
        end if
        
        !scale the root to the desired interval, and put in its symmetric counterpart.
        x(1:m) = xm - xl * z 
        x(n:n-m+1:-1) = xm + xl * z 
        
        w(1:m) = 2.0 * xl / ((1.0-z**2)*pp**2) 
        w(n:n-m+1:-1) = w(1:m)
        
    end subroutine Gauss_Ledrengre
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_weight_index(x, y, z, weight)
        
        integer, intent(in)  :: x
        integer, intent(in)  :: y
        integer, intent(in)  :: z
        integer, intent(in out)  :: weight
        
        integer, save  :: count = 0
        integer, save  :: symmetry(20, 3) = 0
        integer  :: current(3), tmp(3)
        integer  :: i, j 
        logical  :: is_equal
        
        current(1) = x
        current(2) = y
        current(3) = z
        
        if (count == 0)  then
            count = count + 1
            weight = count
            symmetry(weight, :) = current(:)

        else
            do i = 1, count
                tmp = symmetry(i, :)
                is_equal = Is_index_equal(tmp, current)
                
                if (is_equal)  then
                    weight = i
                    exit
                end if
            end do
            
            if (.NOT. is_equal)  then
                count = count + 1
                weight = count
                symmetry(weight, :) = current(:)
            end if
        end if
    
    end subroutine Get_weight_index
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Is_index_equal(left, right)  result(is_true)
        
        integer, intent(in)  :: left(3)
        integer, intent(in)  :: right(3)
        logical  :: is_true
        
        integer  :: tmp_left(3), tmp_right(3)
        integer  :: tmp
        integer  :: i, j
        
        is_true = .TRUE.
        tmp_left = left
        tmp_right = right
        
        do i = 1, SIZE(tmp_left)
            do j = i, SIZE(tmp_left)
                if (tmp_left(j) < tmp_left(i))  then
                    tmp = tmp_left(i)
                    tmp_left(i) = tmp_left(j)
                    tmp_left(j) = tmp
                end if
                if (tmp_right(j) < tmp_right(i))  then
                    tmp = tmp_right(i)
                    tmp_right(i) = tmp_right(j)
                    tmp_right(j) = tmp
                end if
            end do
        end do
        
        do i = 1, SIZE(tmp_left)
            if (tmp_left(i) /= tmp_right(i))  then
                is_true = .FALSE.
                exit
            end if
        end do
    
    end function Is_index_equal
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_QuadratureSet_symmetry (this, ns)
        
        class(QuadratureSet), intent(in out)  :: this
        type(SteadyState), intent(in)         :: ns 
        
        integer  :: i, k, is, iss 
        
        integer  :: SN90                                                        ! number of direction per SN90
        integer  :: mapping(12)                                                 ! symmetry mapping for 60 degree
        integer  :: ik                                                          ! scope according to the value of mapping
        integer  :: jk                                                          ! scope according to the position of mapping
        integer  :: k1
        
        integer  :: octant_z, octant_x, octant_y
        integer  :: is_y, is_x, is_z
        
        ! symmetry with X, Y, Z 
        SN90 = ns%deduce%direction / 8
        do is = 1, ns%deduce%direction
            k = ((is-1) - MOD(is-1, SN90))/SN90 + 1
            
            ! sweep for octant
            !  lower:    ^(y)       upper:     ^(y)                             K=8: octant_z=4, octant_x=6, octant_y=7  
            !            |                     |                                K=7: octant_z=3, octant_x=5, octant_y=8  
            !        <6> | <5>             <2> | <1>                            K=6: octant_z=2, octant_x=8, octant_y=5  
            !       ----------->(x)       ----------->(x)                       K=5: octant_z=1, octant_x=7, octant_y=6  
            !        <8> | <7>             <4> | <3>                            K=4: octant_z=8, octant_x=2, octant_y=3  
            !            |                     |                                K=3: octant_z=7, octant_x=1, octant_y=4  
            !                                                                   K=2: octant_z=6, octant_x=4, octant_y=1  
            !                                                                   K=1: octant_z=5, octant_x=3, octant_y=2  
            
            ! is_y  is direction symmetry with y coordinate
            ! is_x  is direction symmetry with x coordinate
            ! is_z  is direction symmetry with z coordinate
            octant_z = k + (-1)**(k/5)*4                                        ! octant symmetry with z
            octant_x = k + (-1)**((k+(-k/5)*4)/3)*2                             !                      x
            octant_y = k + (-1)**(1-k+k/2*2)                                    !                      y

            iss = is - (k-1)*SN90
            this%is_y(is) = (octant_y-1)*SN90 + iss     
            this%is_x(is) = (octant_x-1)*SN90 + iss     
            this%is_z(is) = (octant_z-1)*SN90 + iss     
        end do 
        
        ! symmetry with 60's
        !  lower:                             upper:                   
        !             \  8  /                             \  2  /      
        !           9  \   /  7                         3  \   /   1   
        !               \ /                                 \ /        
        !       -------------------                 -------------------
        !               / \                                 / \        
        !          10  /   \  12                        4  /   \   6   
        !             /  11 \                             /  5  \      
        !  
        !  (1 3 5 4 2 6) --> (5 1 3 6 2 4) 
        if ((ns%flag%is_60degree) .AND. (ns%flag%is_bevel_edge)) then 
            mapping = 0
            if (ns%flag%n_theta == 60 .OR. ns%flag%n_theta == 240)  then            ! symmetry for 60 degree
                mapping(1)  = 6
                mapping(2)  = 3
                mapping(7)  = 5
                mapping(8)  = 4
                mapping(9)  = 10
                mapping(12) = 11
            else if (ns%flag%n_theta == 120 .OR. ns%flag%n_theta == 300)  then      ! symmetry for 120 degree
                mapping(1)  = 11
                mapping(2)  = 10
                mapping(3)  = 4
                mapping(6)  = 5
                mapping(7)  = 12
                mapping(8)  = 9
            else if (ns%flag%n_theta == 180)  then                                  ! symmetry for 180 degree
                mapping(1)  = 7
                mapping(2)  = 8
                mapping(3)  = 9
                mapping(4)  = 10
                mapping(5)  = 11
                mapping(6)  = 12
            end if 
            
            ! get the 60 degree mapping table by (value, position)
            do k = 1, 12
                k1 = mapping(k)
                if (k1 /= 0)  mapping(k1) = k
            end do
            
            do k = 1, 12
                k1 = mapping(k)
                
                ! ik: scope according to the value of mapping
                if (1<=k1 .AND. k1<=3)  then
                    ik = 1
                else if (4<=k1 .AND. k1<=6)  then
                    ik = 2
                else if (7<=k1 .AND. k1<=9)  then
                    ik = 3
                else if (10<=k1 .AND. k1<=12)  then
                    ik = 4
                end if
                
                ! jk: scope according to the position of mapping
                if (k<=3)  then
                    jk = 1
                else if (k<=6)  then
                    jk = 2
                else if (k<=9)  then
                    jk = 3
                else
                    jk = 4
                end if
                
                do i = 1, ns%state%sn/2
                    is = (jk-1)*SN90 - 3*(jk-1) + 3*(i-1) + k                   
                    this%is_symmetry(is) = (ik-1)*SN90 - 3*(ik-1) + 3*(i-1) + k1     
                end do
            end do
            
            do i = (ns%deduce%direction/2 + 1), ns%deduce%direction
                this%is_symmetry(i) = this%is_symmetry(i-ns%deduce%direction/2) + ns%deduce%direction/2
            end do
        end if 
                
    end subroutine Set_QuadratureSet_symmetry
    
    !$
    !===============================================================================================
    ! set SN sweep order 
    !===============================================================================================
    subroutine Set_QuadratureSet_order (this, ns)
        
        class(QuadratureSet), intent(in out)  :: this
        type(SteadyState), intent(in)         :: ns 
        
        integer, allocatable  :: tmp_(:)
        integer  :: i, j, k 
        integer  :: iline, ioctant, isec
        integer  :: iU, iB 
        integer  :: SN90, SN60
        integer  :: ibeg, iend, idx 
        integer  :: i_allocate
        
        integer  :: fullOrder(8)
        integer  :: degreeU(6), degreeB(6)
        integer  :: secline(6,2)
        integer  :: lineidx(12)
        
        fullOrder = [1, 7, 4, 6, 3, 8, 2, 5]
        degreeU = [1, 3, 5, 4, 2 ,6]
        degreeB = [5, 1, 3, 6, 4, 2]
        
        SN90 = ns%deduce%direction / 8
        do i = 1, SIZE(this%is_order)
            this%is_order(i) = i 
        end do 
        
        ! degree-SYM
        if ((ns%flag%is_60degree) .and. (ns%flag%is_bevel_edge)) then
            secline = RESHAPE([1, 2, 3, 6, 5, 4, 10, 11, 12, 9, 8, 7], [6,2], order=[2,1])
            lineidx = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
            SN60 = ns%state%sn/2 * 2 
            allocate(tmp_(SN60), stat=i_allocate)
            
            do i = 1, SIZE(degreeU)
                ! 
                isec = degreeU(i)
                do k = 1, ns%state%sn/2
                    iline = secline(isec, 1); ioctant = lineidx(iline);
                    tmp_(k) = iline + (k-1)*3 + (ioctant-1)*SN90 - (ioctant-1)*3
                    iline = secline(isec, 2); ioctant = lineidx(iline);
                    tmp_(k+ns%state%sn/2) = iline + (k-1)*3 + (ioctant-1)*SN90 - (ioctant-1)*3
                end do 
                ibeg = (i-1)*SN60*2 + 1
                iend = (i-1)*SN60*2 + SN60
                this%is_order(ibeg: iend) = tmp_
                
                ! 
                isec = degreeB(i)
                do k = 1, ns%state%sn/2
                    iline = secline(isec, 1); ioctant = lineidx(iline);
                    tmp_(k) = iline + (k-1)*3 + (ioctant-1)*SN90 - (ioctant-1)*3
                    iline = secline(isec, 2); ioctant = lineidx(iline);
                    tmp_(k+ns%state%sn/2) = iline + (k-1)*3 + (ioctant-1)*SN90 - (ioctant-1)*3
                end do 
                ibeg = (i-1)*SN60*2 + SN60 + 1
                iend = (i-1)*SN60*2 + SN60 + SN60
                this%is_order(ibeg: iend) = tmp_ + SN90*4
            end do 
            if (allocated(tmp_))    deallocate(tmp_)
            
        ! full-SYM
        else
            allocate(tmp_(SN90), stat=i_allocate)
            do i = 1, SN90
                tmp_(i) = i 
            end do 
            
            do idx = 1, 8
                ibeg = fullOrder(idx); this%is_order((idx-1)*SN90+1 : idx*SN90) = tmp_ + (ibeg-1)*SN90;
            end do 
            
            if (allocated(tmp_))     deallocate(tmp_)
        end if 
        
    end subroutine Set_QuadratureSet_order
    
end module quadrature_header
