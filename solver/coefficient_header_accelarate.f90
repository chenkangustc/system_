!$
!===================================================================================================
!
!   class for accelaration parameter used in out-in iteration
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          LWExtrapolation
!                               SORIterationInner
!                               FSPScalingFactor
!
!===================================================================================================
module coefficient_header_accelarate
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use iteration_header,           only : IterationFlux
        
    implicit none 
    private
    public  :: LWExtrapolation, SORIterationInner, FSPScalingFactor
    
    ! --------------------------------------------------------------------------
    ! type for LW extrapolation for inner iteration
    type  LWExtrapolation
        real(KREAL)  :: rho               = 0.0                                 ! extrapolation factor
        integer      :: inner_limit       = 6                                   ! from which inner index to LW
        real(KREAL)  :: epsilon_old
        real(KREAL)  :: sigma_old
        real(KREAL)  :: delta_old
        real(KREAL)  :: epsilon_new
        real(KREAL)  :: sigma_new
        real(KREAL)  :: delta_new
    contains
        procedure, public  :: set_LW => Set_LWExtrapolation
        procedure, public  :: LW => Perform_LWExtrapolation
    end type LWExtrapolation
    
    ! type for flux distribution during iteration
    type  SORIterationInner
        real(KREAL), public               :: scaling_factor = 1.2               ! 0.4 or 0.7 is OK!
        real(KREAL), public, allocatable  :: nodal(:, :, :)
    contains
        procedure, public  :: alloc => Alloc_SORIterationInner
        procedure, public  :: clean => Free_SORIterationInner
    end type SORIterationInner
    
    ! type for scaling factor method for FSP iteration
    type  FSPScalingFactor
        integer, public             :: lower_limit = 1
        real(KREAL), public         :: source  = 0.0
        real(KREAL), public         :: q_old   = 0.0 
        real(KREAL), public         :: q_new   = 0.0
        real(KREAL), public         :: factor  = 1.0
    end type FSPScalingFactor
    
contains  
    !$
    !===============================================================================================
    ! generate LW extrapolation parameter
    !===============================================================================================
    subroutine Set_LWExtrapolation (this, flux, count_in, ig)
    
        class(LWExtrapolation), intent(in out)  :: this
        type(IterationFlux), intent(in out)     :: flux
        integer, intent(in)                     :: count_in
        integer, intent(in)                     :: ig
        
        ! local variables
        real(KREAL)  :: error_1, error_2, delta_max
        integer          :: il
        
        this%rho  = 0.0

        error_2   = 0.0
        do il = 1, ns%deduce%nodal_total
            error_1 = flux%info%moment(0,il,ig) - flux%info%old(il,ig)
            if (ABS(flux%info%moment(0,il,ig)) > EPS_ZERO)  then
                error_1 = error_1 / flux%info%moment(0,il,ig)
            end if
            error_2 = MAX(error_2, ABS(error_1))
        end do   
        this%epsilon_new = error_2
        
        if (count_in == 1) then
            this%epsilon_old = 1.0
            this%sigma_old   = 1.0
        end if
        this%sigma_new = this%epsilon_new / this%epsilon_old
        this%delta_new = ABS(1.0 - this%sigma_new/this%sigma_old)
        
        if (count_in >= this%inner_limit) then
            delta_max = MAX(this%delta_new, this%delta_old)
            if (this%sigma_new < 1.0 .AND. this%sigma_old < 1.0 .AND. delta_max < 0.5) then
                if(ABS(1.0-this%sigma_new) >= EPS_ZERO)  this%rho = MIN(50.0, this%sigma_new/(1.0-this%sigma_new))
            else
                this%rho = 0.0
            end if
        end if
        
        ! transfer new value and save
        this%epsilon_old = this%epsilon_new
        this%sigma_old   = this%sigma_new
        this%delta_old   = this%delta_new
    
    end subroutine Set_LWExtrapolation

    !$
    !===============================================================================================
    ! perform LW extrapolation
    !===============================================================================================
    subroutine Perform_LWExtrapolation (this, flux, count_in, ig)
    
        class(LWExtrapolation), intent(in)  :: this                             ! extrapolation factor
        type(IterationFlux), intent(in out) :: flux
        integer, intent(in)                 :: count_in                         ! count for inner iteration
        integer, intent(in)                 :: ig
        
        integer  :: il
        
        if (count_in >= this%inner_limit+1) then
            do il = 1, ns%deduce%nodal_total
                flux%info%moment(0,il,ig) = flux%info%moment(0,il,ig) + (flux%info%moment(0,il,ig)-flux%info%old(il,ig)) * this%rho
            end do
        end if
    
    end subroutine Perform_LWExtrapolation
     
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_SORIterationInner (this)
        
        class(SORIterationInner), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allcoated status first
        call this%clean ()
        
        allocate(this%nodal(ns%state%nodal, ns%state%layer, ns%deduce%direction), stat=i_allocate)
        
        this%nodal    = REAL_ZERO
    
    end subroutine Alloc_SORIterationInner
    
    !$
    !===============================================================================================
    ! finalizer for class of SORIterationInner
    !===============================================================================================
    subroutine Free_SORIterationInner (this)
        
        class(SORIterationInner), intent(in out)  :: this
        
        if (allocated(this%nodal))          deallocate(this%nodal)
    
    end subroutine Free_SORIterationInner
        
end module coefficient_header_accelarate
