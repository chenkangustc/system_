!$
!===================================================================================================
!
!   class for mesh sweeping
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          Sweeping                               
!
!===================================================================================================
module sweeping_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use geometry_header,        only : Meshing, Geometry, Boundary
    use quadrature_header,      only : QuadratureSet
    
    implicit none 
    private
    public  :: Sweeping
    
    ! --------------------------------------------------------------------------
    ! type for sweeping order
    type  Sweeping
        integer, public, allocatable   :: order(:, :)                           ! sweeping order per directions
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_Sweeping
        procedure, public  :: set => Set_Sweeping
        procedure, public  :: get_order => Get_sweeping_order
        procedure, public  :: clean =>  Free_Sweeping
    end type Sweeping
    
    ! --------------------------------------------------------------------------
    ! private the real funciton name 
    private  :: Allocate_Sweeping, Set_Sweeping, Get_sweeping_order, Free_Sweeping
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_Sweeping (this)
    
        class(Sweeping), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%order(ns%state%nodal, ns%deduce%direction), stat=i_allocate)
        
        this%order = INT_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + INT_BYTE * SIZE(this%order)
    
    end subroutine Allocate_Sweeping
    
    !$
    !===============================================================================================
    ! finalizer for class of Sweeping
    !===============================================================================================
    subroutine Free_Sweeping (this)
        
        class(Sweeping), intent(in out)  :: this
        
        if (allocated(this%order))         deallocate(this%order)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_Sweeping
    
    !$
    !===============================================================================================
    ! set nodal sweep order for all directions
    !===============================================================================================
    subroutine Set_Sweeping (this, mesh, geom, bound, quad)
        
        class(Sweeping), intent(in out) :: this
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        type(Boundary), intent(in)      :: bound
        type(QuadratureSet), intent(in) :: quad
        
        ! local variables
        integer, parameter  :: RIGHT     =  1                                   !  1, known
        integer, parameter  :: WRONG     = -1                                   ! -1, unknown
        integer             :: i, j, is
        integer             :: i1, i2, i3
        integer             :: j_start, j_end, j_step, isp
        integer             :: is_given(3, ns%state%nodal)                      ! is boundary condition known ? (=1, YES; =-1, NO)
        
        ! ----------------------------------------------------------------------
        do is = 1, ns%deduce%direction
            ! edge j of nodal i is on boundary and is incoming, is_given(j,i)=RIGHT
            do i = 1, ns%state%nodal
                do j = 1, 3
                    is_given(j,i) = bound%nodal(j,i)
!                    if (ABS(bound%nodal(j,i)-bound%INNER) >= EPS_ZERO)  then
                    if (mesh%nearby_nodal(j,i) == i)  then
                        is_given(j,i) = RIGHT
                    end if
                    ! exclude outgoing
                    if (is_given(j,i)==RIGHT .and. quad%directions(is)%projection(j,i)>REAL_ZERO) then
                        is_given(j,i) = WRONG
                    end if
                end do
            end do
            
!            ! negative sweepping in pluse Y hemisphere, positive sweepping in minuns Y hemisphere 
!            if ((ns%deduce%direction/4<is .and. is<=ns%deduce%direction/2) .or. (3*ns%deduce%direction/4<is .and. is<=ns%deduce%direction)) then  
                j_start = ns%state%nodal
                j_end   =  1
                j_step  = -1
!            else
!                j_start = 1
!                j_end   = ns%state%nodal
!                j_step  = 1
!            end if
            
            ! sweep index
            isp = 1   
            do while (isp <= ns%state%nodal)
                nodal_sweep: do j = j_start, j_end, j_step
                    do i1 = 1, 3
                        i2 = i1 + (-2)**(i1/3)
                        i3 = i2 + (-2)**(i2/3)
                        ! one incoming surface: i1--incoming, i2,i3--outgoing 
                        if (quad%directions(is)%projection(i1,j)<REAL_ZERO .and. is_given(i1,j)==RIGHT .and.         &
                            &  quad%directions(is)%projection(i2,j)>REAL_ZERO .and. is_given(i2,j)/=RIGHT .and.      &
                            &  quad%directions(is)%projection(i3,j)>REAL_ZERO .and. is_given(i3,j)/=RIGHT)  then
                            this%order(isp,is) = j
                            
                            ! i2 and i3 become known
                            is_given(i2,j) = RIGHT
                            is_given(i3,j) = RIGHT
                            is_given(mesh%localID(i2,j), mesh%nearby_nodal(i2,j)) = RIGHT
                            is_given(mesh%localID(i3,j), mesh%nearby_nodal(i3,j)) = RIGHT
                            isp = isp + 1
                            exit nodal_sweep
                        end if
                        
                        ! two incoming surface, i1 i2--incoming, i3--outgoing
                        if (quad%directions(is)%projection(i1,j)<REAL_ZERO .and. is_given(i1,j)==RIGHT .and.         &
                            &  quad%directions(is)%projection(i2,j)<REAL_ZERO .and. is_given(i2,j)==RIGHT .and.      &
                            &  quad%directions(is)%projection(i3,j)>REAL_ZERO .and. is_given(i3,j)/=RIGHT)  then
                            this%order(isp,is) = j
                            
                            ! i3 become known
                            is_given(i3,j) = RIGHT
                            is_given(mesh%localID(i3,j),mesh%nearby_nodal(i3,j)) = RIGHT
                            isp = isp + 1
                            exit nodal_sweep
                        end if
                    end do
                    
                    ! sweep to last nodal, still failed
                    if (j == j_end)  then
                        call this_error%set (INFO_LIST_GEOMETRY, 'bad mesh, can not get sweep order')  
                        call this_error%print (FILES%MAIN)
                        call this_error%print (OUTPUT_UNIT)
                    end if
                    
                end do nodal_sweep
            end do
            
        end do
        
    end subroutine Set_Sweeping

    !$
    !===============================================================================================
    ! get sweeping order for specific direction
    !===============================================================================================
    subroutine Get_sweeping_order(this, is, a_order)
    
        class(Sweeping), intent(in out)  :: this
        integer, intent(in)  :: is                                              ! direction
        integer, intent(in out)  :: a_order(:)                                  ! order for this direction
        
        a_order = this%order(:, is)
    
    end subroutine Get_sweeping_order
        
end module sweeping_header
