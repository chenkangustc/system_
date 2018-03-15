!$
!===================================================================================================
!
!   class for transient feedback model used in the LRA benchmark
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          LRAmodel
!
!===================================================================================================
module LRAmodel_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use exception_header,           only : ErrorCollector
    use state_header,               only : SteadyState
    use geometry_header,            only : Meshing, Geometry
    use material_header,            only : CrossSection, CrossSectionInput, Material, cross_section_info_tp
    use contain_header,             only : GroupsFlux
    
    implicit none 
    private
    public  :: LRAmodel
    
    ! --------------------------------------------------------------------------
    ! type for feedback model for LRA benchmark
    type  LRAmodel
        real(KREAL), public               :: tf0             = 300.0D0          ! initial fuel temperature
        real(KREAL), public               :: tf_avg          = 300.0D0          ! fuel average temperature 
        real(KREAL), public, allocatable  :: tf(:, :)                           ! real time temperature per zone per layer
        real(KREAL), public, allocatable  :: sigma_a0(:)                        ! keep initial value to perform doppler effect
                                                                                
        real(KREAL), public               :: nu              = 2.43D0           ! neutron number generation per fission
        real(KREAL), public               :: kappa           = 3.204D-11        ! power transfer coefficient
        real(KREAL), public               :: adiabatic       = 1.1954D0         ! adiabatic feedback coefficient (1.1954 = 3.830E-11 / 3.204E-11)
        real(KREAL), public               :: doppler         = 3.034D-03        ! doppler feedback coefficient
    contains
        procedure, public  :: alloc => Allocate_LRAmodel
        procedure, public  :: clean => Free_LRAmodel
        procedure, public  :: update_tf => Update_LRAmodel_tf
        procedure, public  :: update_xsec => Update_LRAmodel_xsec
    end type LRAmodel
    
    ! private the real function name
    private  :: Allocate_LRAmodel, Free_LRAmodel, Update_LRAmodel_tf, Update_LRAmodel_xsec
    
    type(ErrorCollector)  :: a_error
    integer, parameter    :: i_GROUP = 1
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_LRAmodel (this, ns, mat_info, xsec_inp)
        
        class(LRAmodel), intent(in out)     :: this
        type(SteadyState), intent(in)       :: ns
        type(Material), intent(in)      :: mat_info
        type(CrossSectionInput), intent(in) :: xsec_inp
        
        integer  :: i_allocate
        integer  :: ia, iz, im
        
        ! check allocated status first
        call this%clean ()
        
        if (ns%feedback%is_model .and. TRIM(ADJUSTL(ns%feedback%model_name)) == 'LRA')  then
            if (ns%state%ng /= 2)  then
                call a_error%set (INFO_LIST_FEEDBACK, 'the LRA model suit for 2G problem only')
                call a_error%print (OUTPUT_UNIT)
            end if
            
            allocate(this%tf(ns%state%zone, ns%state%layer), stat=i_allocate)
            allocate(this%sigma_a0(ns%state%mat), stat=i_allocate)
            
            ! get the absorption cross section (first group)
            do im = 1, ns%state%mat
                this%sigma_a0(im) = xsec_inp%mats(im)%sigma_t(i_GROUP) - SUM(xsec_inp%mats(im)%sigma_s(i_GROUP, :, 1))
            end do
                        
            this%tf = this%tf0
            this%tf_avg = this%tf0
        end if
        
    end subroutine Allocate_LRAmodel
    
    !$
    !===============================================================================================
    ! finalizer for class of LRAmodel
    !===============================================================================================
    subroutine Free_LRAmodel (this)
        
        class(LRAmodel), intent(in out)  :: this
        
        if (allocated(this%tf))                 deallocate(this%tf)
        if (allocated(this%sigma_a0))           deallocate(this%sigma_a0)
    
    end subroutine Free_LRAmodel 
    
    !$
    !===============================================================================================
    ! update nodal temperature by adiabatic model
    !===============================================================================================
    subroutine Update_LRAmodel_tf (this, ns, mesh, geom, mat_info, xsec, flux, t_step, Tfuel)
        
        class(LRAmodel), intent(in out)    :: this
        type(SteadyState), intent(in)      :: ns
        type(Meshing), intent(in)          :: mesh
        type(Geometry), intent(in)         :: geom
        type(Material), intent(in)         :: mat_info
        type(CrossSection), intent(in)     :: xsec
        type(GroupsFlux), intent(in)       :: flux
        real(KREAL), intent(in)            :: t_step                            ! step length
        real(KREAL), intent(in out), optional  :: Tfuel(:, :)
        
        real(KREAL)  :: alpha
        real(KREAL)  :: delta
        real(KREAL)  :: total, divide
        integer      :: ir, ia, ig, iz
        
        alpha = this%adiabatic * this%kappa
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                delta = 0.0D0
                do ig = 1, ns%state%ng
                    delta = delta + t_step * alpha * xsec%matrixs(iz,ia)%sigma_f_nu(ig) * flux%ngs(ig)%scalar(ir,ia) / this%nu
                end do
                this%tf(iz, ia) = this%tf(iz, ia) + delta
            end do
        end do
        
        ! volume average
        total = 0.0D0; divide = 0.0D0
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                if (mat_info%mask_core(iz, ia))  then
                    total = total + this%tf(iz, ia) * geom%height(ia) * geom%zone_area(iz)
                    divide = divide + geom%height(ia) * geom%zone_area(iz)
                end if
            end do
        end do
        this%tf_avg = total / divide
        
        if (PRESENT(Tfuel))  then
            Tfuel = this%tf
        end if
    
    end subroutine Update_LRAmodel_tf
    
    !$
    !===============================================================================================
    ! update zone xsec by doppler feedback
    !===============================================================================================
    subroutine Update_LRAmodel_xsec (this, ns, iz, ia, im, a_xsec)
    
        class(LRAmodel), intent(in out)  :: this
        type(SteadyState), intent(in)    :: ns
        integer, intent(in)              :: iz
        integer, intent(in)              :: ia
        integer, intent(in)              :: im
        type(cross_section_info_tp), intent(in out)  :: a_xsec
        
        real(KREAL)  :: sigma_a
        real(KREAL)  :: factor
        
        ! perform doppler effect
        factor = 1.0D0 + this%doppler * (SQRT(this%tf(iz,ia)) - SQRT(this%tf0))
        sigma_a = this%sigma_a0(im) * factor
        
        ! reconstruct total cross section by sigma_a
        a_xsec%sigma_t(i_GROUP) = sigma_a + SUM(a_xsec%sigma_s(i_GROUP, :, 1))
        
    end subroutine Update_LRAmodel_xsec
    
end module LRAmodel_header
