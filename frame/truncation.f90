!$
!===================================================================================================
!
!   module of numerical truncation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    truncate_negative
!                               truncate_zero
!
!   Public type lists:          No
!
!===================================================================================================
module truncation
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: truncate_negative, truncate_zero
    
    interface  truncate_negative
        module procedure  truncate_negative_value
        module procedure  truncate_negative_value_1D
    end interface truncate_negative
    
    interface truncate_zero
        module procedure  truncat_zero_value
        module procedure  truncat_zero_value_1D
    end interface truncate_zero
    
contains
    !$
    !===============================================================================================
    ! value >= EPS_ZERO
    !===============================================================================================
    pure subroutine truncate_negative_value (value)
        
        real(KREAL), intent(in out)  :: value
        
        if (value < EPS_ZERO)  then
            value = EPS_ZERO
        end if
    
    end subroutine truncate_negative_value
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    pure subroutine truncate_negative_value_1D (value)
        
        real(KREAL), intent(in out)  :: value(:)
        
        where (value < EPS_ZERO)  
            value = EPS_ZERO
        end where
    
    end subroutine truncate_negative_value_1D
    
    !$
    !===============================================================================================
    ! abs(value) >= EPS_ZERO
    !===============================================================================================
    subroutine truncat_zero_value (value)
    
        real(KREAL), intent(in out)  :: value
        
        if (REAL_ZERO <= value .and. value < EPS_ZERO)  then
            value = EPS_ZERO
        else if (-EPS_ZERO < value .and. value <= REAL_ZERO) then
            value = -EPS_ZERO
        end if
    
    end subroutine truncat_zero_value
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine truncat_zero_value_1D (value)
        
        real(KREAL), intent(in out)  :: value(:)

        where (REAL_ZERO <= value .and. value < EPS_ZERO)
            value = EPS_ZERO
        else where (-EPS_ZERO < value .and. value <= REAL_ZERO)
            value = -EPS_ZERO
        end where
        
    end subroutine truncat_zero_value_1D
    
end module truncation
