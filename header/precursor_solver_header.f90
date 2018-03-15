!$
!===================================================================================================
!
!   class for delayed neutron precursor solver
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    None
!
!   Public type lists:          PrecursorSolver      
!
!===================================================================================================
module precursor_solver_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use timestep_header,            only : TimeStepInfo
    use geometry_header,            only : Meshing, Geometry
    use material_header,            only : CrossSection, KineticsParameter
    use contain_header,             only : GroupsFlux, DistributionParameter, TimeListParameter
    
    implicit none
    private
    public  :: PrecursorSolver
    
    ! --------------------------------------------------------------------------
    ! type for precursor solver
    type  PrecursorSolver
        real(KREAL), public, allocatable      :: Lprecursor(:, :, :)
        real(KREAL), public, allocatable      :: Rprecursor(:, :, :)
        real(KREAL), public, allocatable      :: Lflux(:, :, :)
        real(KREAL), public, allocatable      :: Rflux(:, :, :)
        real(KREAL)      :: left
        real(KREAL)      :: right
        real(KREAL)      :: pace
        integer          :: ng
        integer          :: dg
        integer          :: nr
        integer          :: na
        integer          :: method = 1                                          ! 0/1
    contains    
        procedure, public  :: alloc => Allocate_PrecursorSolver
        procedure, public  :: clean => Free_PrecursorSolver
        procedure, public  :: init => Initial_Precursor
        procedure, public  :: advance => Advance_Precursor
    end type PrecursorSolver

    ! private the real function name
    private  :: Allocate_PrecursorSolver, Free_PrecursorSolver
    private  :: Initial_Precursor, Advance_Precursor
    
    integer, parameter  :: FISSION_RATE_CONSTANT = 0
    integer, parameter  :: FISSION_RATE_LINEAR = 1

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_PrecursorSolver (this)
        
        class(PrecursorSolver), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        this%ng = ns%state%ng
        this%nr = ns%state%nodal
        this%na = ns%state%layer
        this%dg = nt%state%dg
        
        allocate(this%Lprecursor(this%dg, this%nr, this%na), stat=i_allocate)
        allocate(this%Rprecursor(this%dg, this%nr, this%na), stat=i_allocate)
        allocate(this%Lflux(this%nr, this%na, this%ng), stat=i_allocate)
        allocate(this%Rflux(this%nr, this%na, this%ng), stat=i_allocate)

        this%Lprecursor = REAL_ZERO
        this%Rprecursor = REAL_ZERO
        this%Lflux = REAL_ZERO
        this%Rflux = REAL_ZERO
        
    end subroutine Allocate_PrecursorSolver

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_PrecursorSolver (this)
        
        class(PrecursorSolver), intent(in out)  :: this
        
        if (allocated(this%Lprecursor))         deallocate(this%Lprecursor)
        if (allocated(this%Rprecursor))         deallocate(this%Rprecursor)
        if (allocated(this%Lflux))              deallocate(this%Lflux)
        if (allocated(this%Rflux))              deallocate(this%Rflux)
    
    end subroutine Free_PrecursorSolver
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Initial_Precursor (this, mesh, geom, xsec, param, flux, dist_dnps, timelist)
        
        class(PrecursorSolver), intent(in out)      :: this
        type(Meshing), intent(in)                   :: mesh
        type(Geometry), intent(in)                  :: geom 
        type(CrossSection), intent(in)              :: xsec
        type(KineticsParameter), intent(in)         :: param
        type(GroupsFlux), intent(in)                :: flux
        type(DistributionParameter), intent(in out) :: dist_dnps(:)
        type(TimeListParameter), intent(in out)     :: timelist
        
        real(KREAL)      :: Rrate
        integer          :: ig, ia, ir, ip, iz
        
        do ig = 1, SIZE(flux%ngs)
            this%Rflux(:, :, ig) = flux%ngs(ig)%scalar(:, :)
        end do
        this%Rprecursor = 0.0D0
        do ia = 1, this%na
            do ir = 1, this%nr
                iz = mesh%zone(ir)
                Rrate = 0.0D0
                do ig = 1, this%ng
                    Rrate = Rrate + this%Rflux(ir, ia, ig)*xsec%matrixs(iz, ia)%sigma_f_nu(ig)
                end do
                
                do ip = 1, this%dg
                    associate(beta => param%matrixs(iz,ia)%beta(ip), lambda => param%matrixs(iz,ia)%lambda(ip))
                    if (lambda > EPS_ZERO)  then
                        this%Rprecursor(ip, ir, ia) = Rrate * (beta/lambda)
                    else 
                        this%Rprecursor(ip, ir, ia) = 0.0D0
                    end if
                    end associate
                end do
            end do
        end do
        
        ! update for the previous time point
        this%Lflux = this%Rflux
        this%Lprecursor = this%Rprecursor
        
        ! return precursor
        do ip = 1, this%dg
            dist_dnps(ip)%matrix(:, :) = this%Rprecursor(ip, :, :)
        end do
        
        timelist%precursor = 0.0D0
        do ia = 1, this%na
            do ir = 1, this%nr
                do ip = 1, this%dg
                    timelist%precursor(ip) = timelist%precursor(ip) + dist_dnps(ip)%matrix(ir,ia)*geom%area(ir)*geom%height(ia)
                end do
            end do
        end do

    end subroutine Initial_Precursor
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Advance_Precursor (this, mesh, xsec, param, flux, dist_dnps, a_step)
        
        class(PrecursorSolver), intent(in out)      :: this
        type(Meshing), intent(in)                   :: mesh
        type(CrossSection), intent(in)              :: xsec
        type(KineticsParameter), intent(in)         :: param
        type(GroupsFlux), intent(in)                :: flux
        type(DistributionParameter), intent(in out) :: dist_dnps(:)
        type(TimeStepInfo), intent(in)              :: a_step
        
        real(KREAL)  :: Rrate, Lrate
        real(KREAL)  :: Rfactor, Lfactor
        integer      :: ig, ia, ir, ip, iz
        
        do ig = 1, SIZE(flux%ngs)
            this%Rflux(:, :, ig) = flux%ngs(ig)%scalar(:, :)
        end do
        this%Rprecursor = 0.0D0
        this%left = a_step%left 
        this%right = a_step%right
        this%pace = this%right - this%left
        
        select case(this%method)
        case (FISSION_RATE_CONSTANT)                                            ! constant approximation
            do ia = 1, this%na
                do ir = 1, this%nr
                    iz = mesh%zone(ir)
                    Rrate = 0.0D0
                    do ig = 1, this%ng
                        Rrate = Rrate + this%Rflux(ir, ia, ig)*xsec%matrixs(iz,ia)%sigma_f_nu(ig)
                    end do
                    
                    do ip = 1, this%dg
                        associate(beta => param%matrixs(iz,ia)%beta(ip), lambda => param%matrixs(iz,ia)%lambda(ip))
                        if (lambda > EPS_ZERO)  then
                            Rfactor = 1.0D0 / (1.0D0 + this%pace*lambda)
                            this%Rprecursor(ip, ir, ia) = this%Lprecursor(ip, ir, ia)*Rfactor + Rrate*this%pace*beta*Rfactor
                        else
                            this%Rprecursor(ip, ir, ia) = 0.0D0
                        end if
                        end associate
                    end do
                end do
            end do

        case (FISSION_RATE_LINEAR)                                              ! linear approximation
            do ia = 1, this%na
                do ir = 1, this%nr
                    iz = mesh%zone(ir)
                    Rrate = 0.0D0; Lrate = 0.0D0
                    do ig = 1, this%ng
                        Rrate = Rrate + this%Rflux(ir, ia, ig)*xsec%matrixs(iz,ia)%sigma_f_nu(ig)
                        Lrate = Lrate + this%Lflux(ir, ia, ig)*xsec%matrixs(iz,ia)%sigma_f_nu(ig)
                    end do
                    
                    do ip = 1, this%dg
                        associate(beta => param%matrixs(iz,ia)%beta(ip), lambda => param%matrixs(iz,ia)%lambda(ip))
                        if (lambda > EPS_ZERO)  then
                            Rfactor = (beta/lambda/this%pace) * (this%pace - (1.0D0/lambda)*(1.0D0-EXP(-lambda*this%pace)))
                            Lfactor = (1.0D0-EXP(-lambda*this%pace)) * (beta/lambda) - Rfactor
                            this%Rprecursor(ip, ir, ia) = this%Lprecursor(ip, ir, ia)*EXP(-lambda*this%pace) + Lfactor*Lrate + Rfactor*Rrate
                        else 
                            this%Rprecursor(ip, ir, ia) = 0.0D0
                        end if
                        end associate
                    end do
                end do
            end do
            
        end select 
        
        ! update for the previous time point
        this%Lflux = this%Rflux
        this%Lprecursor = this%Rprecursor
        
        ! return precursor
        do ip = 1, this%dg
            dist_dnps(ip)%matrix(:, :) = this%Rprecursor(ip, :, :)
        end do        
    
    end subroutine Advance_Precursor

end module precursor_solver_header
