!$
!===================================================================================================
!
!   module for coefficient parameter used by point kinetics method
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          PKLocalState
!                               DensityVector
!                               RichardsonMatrix
!
!===================================================================================================
module coefficient_pk_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    implicit none 
    private
    public  :: PKLocalState, DensityVector, RichardsonMatrix
    public  :: pk_state, last_vector, current_vector, richard
    
    ! --------------------------------------------------------------------------
    ! state parameter for all the program
    type  PKLocalState
        integer                      :: n_Richard    = 12                       ! size of one row of the Richardson matrix
        real(KREAL)              :: tolerance    = 1.0E-12                  ! 12 digit percision
        character(len=MAX_WORD_LEN)  :: method       = 'VODE'                   ! 'VODE' / 'RBDF'
    end type PKLocalState

    ! type for storing density information
    type  DensityVector
        real(KREAL), allocatable   :: vector(:)
    contains
        procedure, public  :: alloc => Alloc_DensityVector
        procedure, public  :: clean => Free_DensityVector
    end type DensityVector
    
    ! type for richardson matrix
    type  RichardsonMatrix
        type(DensityVector), allocatable  :: matrix(:, :)
    contains
        procedure, public  :: alloc => Alloc_RichardsonMatrix
        procedure, public  :: clean => Free_RichardsonMatrix
    end type RichardsonMatrix
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Alloc_DensityVector, Free_DensityVector
    private  :: Alloc_RichardsonMatrix, Free_RichardsonMatrix
    
    ! --------------------------------------------------------------------------
    ! global parameter for pk solving
    type(PKLocalState)      :: pk_state
    type(DensityVector)     :: last_vector
    type(DensityVector)     :: current_vector
    type(RichardsonMatrix)  :: richard
        
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Alloc_DensityVector (this)
        
        class(DensityVector), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%vector(nt%state%dg+1), stat=i_allocate)
        
        this%vector = REAL_ZERO
    
    end subroutine Alloc_DensityVector
    
    !$
    !===============================================================================================
    ! finalizer for class of DensityVector
    !===============================================================================================
    subroutine Free_DensityVector (this)
    
        class(DensityVector), intent(in out)  :: this
        
        if (allocated(this%vector))             deallocate(this%vector)
    
    end subroutine Free_DensityVector

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Alloc_RichardsonMatrix(this, pk_state)
        
        class(RichardsonMatrix), intent(in out)  :: this
        type(PKLocalState), intent(in)  :: pk_state
        
        integer  :: i_allocate
        integer  :: i,j
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%matrix(pk_state%n_Richard, pk_state%n_Richard), stat=i_allocate)
        
        do i = 1, SIZE(this%matrix, dim=1)
            do j = 1, SIZE(this%matrix, dim=2)
                call this%matrix(i,j)%alloc ()
            end do
        end do
    
    end subroutine Alloc_RichardsonMatrix
    
    !$
    !===============================================================================================
    ! finalizer for class of RichardsonMatrix 
    !===============================================================================================
    subroutine Free_RichardsonMatrix(this)
    
        class(RichardsonMatrix), intent(in out)  :: this
        
        integer  :: i,j
        
        if (allocated(this%matrix))  then
            do i = 1, SIZE(this%matrix, dim=1)
                do j = 1, SIZE(this%matrix, dim=2)
                    call this%matrix(i,j)%clean ()
                end do
            end do
            
            ! free itself
            deallocate(this%matrix)
        end if
        
    end subroutine Free_RichardsonMatrix
    
end module coefficient_pk_header
