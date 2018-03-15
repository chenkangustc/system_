!$
!===================================================================================================
!
!   class for perturbation worth calculation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          CRWorth
!
!===================================================================================================
module worth_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    private
    public  :: CRWorth
    
    ! type for CR-worth calculation
    type  CRWorth
        integer, public  :: pos_type                        = -1                ! input type 
        integer, public  :: bank                            = 1                 ! which bank 
        real(KREAL), public, allocatable  :: pos_input(:)
        real(KREAL), public  :: pos_beg                     = 0.0
        real(KREAL), public  :: pos_end                     = 0.0
        real(KREAL), public  :: pos_step                    = 1.0
        integer  ::  BY_INP  = 0
        integer  ::  BY_STEP = 1 
    contains
        procedure, public  :: alloc => Allocate_CRWorth
        procedure, public  :: clean => Free_CRWorth
        procedure, public  :: getNstep => Get_CRWorth_Nstep 
        procedure, public  :: getistep => Get_CRWorth_istep 
    end type CRWorth
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_CRWorth(this, n)
        
        class(CRWorth), intent(in out)  :: this
        integer, intent(in)   :: n
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%pos_input(n), stat=i_allocate)
        this%pos_input = 0.0 
    
    end subroutine Allocate_CRWorth
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_CRWorth(this)
        
        class(CRWorth), intent(in out)  :: this
        
        if (allocated(this%pos_input))      deallocate(this%pos_input)
    
    end subroutine Free_CRWorth
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_CRWorth_Nstep(this)  result(Nstep)
        
        class(CRWorth), intent(in)  :: this
        integer  :: Nstep
        
        if (this%pos_type == this%BY_INP)  then
            Nstep = SIZE(this%pos_input)
        else if (this%pos_type == this%BY_STEP)  then
            Nstep = CEILING((this%pos_end-this%pos_beg)/this%pos_step)
        end if 
    
    end function Get_CRWorth_Nstep
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    function Get_CRWorth_istep(this, ipert)  result(istep)
        
        class(CRWorth), intent(in)  :: this
        integer, intent(in)  :: ipert 
        integer  :: istep
        
        if (this%pos_type == this%BY_INP)  then
            istep = this%pos_input(ipert)
        else if (this%pos_type == this%BY_STEP)  then
            istep = MIN(this%pos_beg + (ipert-1)*this%pos_step, this%pos_end)
        end if 
    
    end function Get_CRWorth_istep
    
end module worth_header 
