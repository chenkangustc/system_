!$
!===================================================================================================
!
!   class for adjoint iteration parameter
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          AdjointIteration
!
!===================================================================================================
module adjoint_header
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use geometry_header,        only : Meshing
    use material_header,        only : CrossSection, ExternalSource
    
    implicit none 
    private
    public  :: AdjointIteration
    
    ! --------------------------------------------------------------------------
    ! type for fission source term when adjoint calculation
    type  AdjointIteration
        logical  :: is_eigenvalue  = .TRUE.                                     ! is adjoint problem eigenvalue type ?
        real(KREAL), allocatable  :: q_moments(:, :)
        real(KREAL), allocatable  :: q_old(:)
        real(KREAL), allocatable  :: sigma_f(:)
        real(KREAL), allocatable  :: source(:, :, :)                            ! save for adjoint source term
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_AdjointIteration
        procedure, public  :: clean =>  Free_AdjointIteration

        procedure, public  :: homo => Set_adjoint_homogeneous
        procedure, public  :: one => Set_adjoint_one
        procedure, public  :: neutron => Set_adjoint_fission_neutron
        procedure, public  :: energy => Set_adjoint_fission_energy
        procedure, public  :: actual => Set_adjoint_actual_source
    end type AdjointIteration
    
    ! private the real function name
    private  :: Allocate_AdjointIteration, Free_AdjointIteration
    private  :: Set_adjoint_one, Set_adjoint_fission_neutron, Set_adjoint_fission_energy, Set_adjoint_actual_source
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_AdjointIteration (this)
        
        class(AdjointIteration), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%q_moments(0:8, ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%q_old(ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%sigma_f(ns%deduce%nodal_total), stat=i_allocate)
        allocate(this%source(ns%state%zone, ns%state%layer, ns%state%ng), stat=i_allocate)
        
        this%q_moments      = REAL_ZERO
        this%q_old          = REAL_ZERO
        this%sigma_f        = REAL_ZERO
        this%source         = REAL_ZERO
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%q_moments)
        this%memory = this%memory + REAL_BYTE * SIZE(this%q_old)
        this%memory = this%memory + REAL_BYTE * SIZE(this%sigma_f)
        this%memory = this%memory + REAL_BYTE * SIZE(this%source)
    
    end subroutine Allocate_AdjointIteration
    
    !$
    !===============================================================================================
    ! finalizer for class of AdjointIteration
    !===============================================================================================
    subroutine  Free_AdjointIteration (this)
        
        class(AdjointIteration)  :: this
        
        if (allocated(this%q_moments))          deallocate(this%q_moments)
        if (allocated(this%q_old))              deallocate(this%q_old)
        if (allocated(this%sigma_f))            deallocate(this%sigma_f)
        if (allocated(this%source))             deallocate(this%source)
        
        this%memory = REAL_ZERO
    
    end subroutine Free_AdjointIteration
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_adjoint_homogeneous (this)
    
        class(AdjointIteration), intent(in out)  :: this
        
        this%source = REAL_ZERO

        this%is_eigenvalue = .TRUE.        
    
    end subroutine Set_adjoint_homogeneous
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_adjoint_one (this)
        
        class(AdjointIteration), intent(in out)  :: this
        
        this%source = REAL_ONE
        
        this%is_eigenvalue = .FALSE.
        
    end subroutine Set_adjoint_one
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_adjoint_fission_neutron (this, xsec)
    
        class(AdjointIteration), intent(in out)  :: this
        type(CrossSection), intent(in)  :: xsec
        
        integer  :: iz, ia
        
        do iz = 1, ns%state%zone
            do ia = 1, ns%state%layer
                this%source(iz, ia, :) = xsec%matrixs(iz, ia)%sigma_f_nu(:)
            end do
        end do
        
        this%is_eigenvalue = .FALSE.
        
    end subroutine Set_adjoint_fission_neutron
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_adjoint_fission_energy (this, xsec)
    
        class(AdjointIteration), intent(in out)  :: this
        type(CrossSection), intent(in)  :: xsec
        
        integer  :: iz, ia
        
        do iz = 1, ns%state%zone
            do ia = 1, ns%state%layer
                this%source(iz, ia, :) = xsec%matrixs(iz, ia)%sigma_f_kappa(:)
            end do
        end do
        
        this%is_eigenvalue = .FALSE.
        
    end subroutine Set_adjoint_fission_energy
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_adjoint_actual_source (this, Q_ext)
    
        class(AdjointIteration), intent(in out)  :: this
        type(ExternalSource), intent(in)         :: Q_ext
        
        integer  :: iz, ia
        
        do iz = 1, ns%state%zone
            do ia = 1, ns%state%layer
                this%source(iz, ia, :) = Q_ext%matrixs(iz, ia)%intensity(:)
            end do
        end do
        
        this%is_eigenvalue = .FALSE.
        
    end subroutine Set_adjoint_actual_source
    
end module adjoint_header
