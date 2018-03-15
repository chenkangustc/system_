!$
!===================================================================================================
!
!    this module is for Xenon & Samarium parameter definition;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          XenonSamariumInput
!                               XenonSamarium
!
!===================================================================================================
module xesm_header

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use material_header,        only : CrossSection, Material
    
    implicit none
    private
    public  :: XenonSamariumInput, XenonSamarium
    
    ! --------------------------------------------------------------------------
    ! type for Xe & Sm input parameter
    type XenonSamariumInput
        real(KREAL), public  :: lambda_I                                    ! decay constants of I  (Iodine)
        real(KREAL), public  :: lambda_Xe                                   ! decay constants of Xe (Xenon)
        real(KREAL), public  :: lambda_Pm                                   ! decay constants of Pm (Promethium)
        real(KREAL), public  :: lambda_Sm                                   ! decay constants of Sm (Samarium)
        
        real(KREAL), public, allocatable  :: gamma_I(:, :)                  ! fission product yield of I
        real(KREAL), public, allocatable  :: gamma_Xe(:, :)                 ! fission product yield of Xe
        real(KREAL), public, allocatable  :: gamma_Pm(:, :)                 ! fission product yield of Pm
        real(KREAL), public, allocatable  :: gamma_Sm(:, :)                 ! fission product yield of Sm
    contains
        procedure, public  :: alloc => Alloc_XenonSamariumInput
        procedure, public  :: clean => Free_XenonSamariumInput
    end type XenonSamariumInput
    
    ! type for Xe & Sm
    type, private  :: nuclide_info_tp
        real(KREAL), public  :: lambda                                      ! decay constants
        real(KREAL), public  :: density                                     ! nuclide density for current time point
        real(KREAL), public  :: last_density                                ! nuclide density for last time point
        
        real(KREAL), public, allocatable  :: gamma(:)                       ! fission product yield
        real(KREAL), public, allocatable  :: sigma_a(:)                     ! micro absorption xsec for current time point 
        real(KREAL), public, allocatable  :: last_sigma_a(:)                ! micro absorption xsec for last time point
    contains
        procedure, public  :: alloc => Alloc_nuclide_info_tp
        procedure, public  :: clean => Free_nuclide_info_tp
        procedure, public  :: transit => Transit_nuclide_info_tp
    end type  nuclide_info_tp
    
    type  XenonSamarium
        type(nuclide_info_tp), public, allocatable  :: Io(:, :)                 ! Iodine Number Density
        type(nuclide_info_tp), public, allocatable  :: Xe(:, :)                 ! Xenon Number Density
        type(nuclide_info_tp), public, allocatable  :: Pm(:, :)                 ! Promethium Number Density
        type(nuclide_info_tp), public, allocatable  :: Sm(:, :)                 ! Samarium Number Density
    contains
        procedure, public  :: alloc => Alloc_XenonSamarium
        procedure, public  :: clean => Free_XenonSamarium
        procedure, public  :: set => Set_XenonSamarium
        procedure, public  :: steady => Steady_XenonSamarium
        procedure, public  :: transient => Transient_XenonSamarium
        procedure, public  :: update_xsec => Update_XenonSamarium_xsec
        procedure, public  :: transit => Transit_XenonSamarium
    end type XenonSamarium
    
    ! --------------------------------------------------------------------------
    private  :: Alloc_XenonSamariumInput, Free_XenonSamariumInput
    
    private  :: Alloc_nuclide_info_tp, Free_nuclide_info_tp, Transit_nuclide_info_tp
    private  :: Alloc_XenonSamarium, Free_XenonSamarium, Set_XenonSamarium
    private  :: Steady_XenonSamarium, Transient_XenonSamarium, Update_XenonSamarium_xsec, Transit_XenonSamarium
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_XenonSamariumInput (this)
        
        class(XenonSamariumInput), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%gamma_I(ns%state%ng, ns%state%mat), stat=i_allocate)
        allocate(this%gamma_Xe(ns%state%ng, ns%state%mat), stat=i_allocate)
        allocate(this%gamma_Pm(ns%state%ng, ns%state%mat), stat=i_allocate)
        allocate(this%gamma_Sm(ns%state%ng, ns%state%mat), stat=i_allocate)
        
        ! default value, from PARCS
        this%lambda_I  = 0.287500E-4
        this%lambda_Xe = 0.209167E-4
        this%lambda_Pm = 0.355568E-8
        this%lambda_Sm = 0.0
        
        this%gamma_I  = 0.06386
        this%gamma_Xe = 0.00228
        this%gamma_Pm = 0.01130
        this%gamma_Sm = 0.0
    
    end subroutine Alloc_XenonSamariumInput
    
    !$
    !===============================================================================================
    ! finalizer for class of XenonSamariumInput
    !===============================================================================================
    subroutine Free_XenonSamariumInput (this)
    
        class(XenonSamariumInput), intent(in out)  :: this
        
        if (allocated(this%gamma_I))            deallocate(this%gamma_I)
        if (allocated(this%gamma_Xe))           deallocate(this%gamma_Xe)
        if (allocated(this%gamma_Pm))           deallocate(this%gamma_Pm)
        if (allocated(this%gamma_Sm))           deallocate(this%gamma_Sm)
    
    end subroutine Free_XenonSamariumInput
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_nuclide_info_tp (this)
        
        class(nuclide_info_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%gamma(ns%state%ng), stat=i_allocate)
        allocate(this%sigma_a(ns%state%ng), stat=i_allocate)
        allocate(this%last_sigma_a(ns%state%ng), stat=i_allocate)
        
        this%lambda = 0.0
        this%gamma  = 0.0

        this%density  = 1.0
        this%sigma_a  = 0.0

        this%last_density = 1.0
        this%last_sigma_a = 0.0
    
    end subroutine Alloc_nuclide_info_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of nuclide_info_tp
    !===============================================================================================
    subroutine Free_nuclide_info_tp (this)
            
        class(nuclide_info_tp), intent(in out)  :: this
        
        if (allocated(this%gamma))              deallocate(this%gamma)
        if (allocated(this%sigma_a))            deallocate(this%gamma)
        if (allocated(this%last_sigma_a))       deallocate(this%last_sigma_a)
        
    end subroutine Free_nuclide_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_nuclide_info_tp (this)
        
        class(nuclide_info_tp), intent(in out)  :: this
    
        this%last_density = this%density
        this%last_sigma_a = this%sigma_a
        
        this%density = 0.0
        this%sigma_a = 0.0
        
    end subroutine Transit_nuclide_info_tp

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_XenonSamarium (this)
        
        class(XenonSamarium), intent(in out)  :: this

        integer  :: i, j
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%Io(ns%state%zone, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%Io, dim=1)
            do j = 1, SIZE(this%Io, dim=2)
                call this%Io(i, j)%alloc ()
            end do
        end do
    
        allocate(this%Xe(ns%state%zone, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%Xe, dim=1)
            do j = 1, SIZE(this%Xe, dim=2)
                call this%Xe(i, j)%alloc ()
            end do
        end do
    
        allocate(this%Pm(ns%state%zone, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%Pm, dim=1)
            do j = 1, SIZE(this%Pm, dim=2)
                call this%Pm(i, j)%alloc ()
            end do
        end do
    
        allocate(this%Sm(ns%state%zone, ns%state%layer), stat=i_allocate)
        do i = 1, SIZE(this%Sm, dim=1)
            do j = 1, SIZE(this%Sm, dim=2)
                call this%Sm(i, j)%alloc ()
            end do
        end do
    
    end subroutine Alloc_XenonSamarium

    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_XenonSamarium (this)
        
        class(XenonSamarium), intent(in out)  :: this
        integer  :: i, j
        
        if (allocated(this%Io))  then
            do i = 1, SIZE(this%Io, dim=1)
                do j = 1, SIZE(this%Io, dim=2)
                    call this%Io(i, j)%clean ()
                end do
            end do
            ! free itself
            deallocate(this%Io)
        end if
    
        if (allocated(this%Xe))  then
            do i = 1, SIZE(this%Xe, dim=1)
                do j = 1, SIZE(this%Xe, dim=2)
                    call this%Xe(i, j)%clean ()
                end do
            end do
            ! free itself
            deallocate(this%Xe)
        end if
    
        if (allocated(this%Pm))  then
            do i = 1, SIZE(this%Pm, dim=1)
                do j = 1, SIZE(this%Pm, dim=2)
                    call this%Pm(i, j)%clean ()
                end do
            end do
            ! free itself
            deallocate(this%Pm)
        end if
    
        if (allocated(this%Sm))  then
            do i = 1, SIZE(this%Sm, dim=1)
                do j = 1, SIZE(this%Sm, dim=2)
                    call this%Sm(i, j)%clean ()
                end do
            end do
            ! free itself
            deallocate(this%Sm)
        end if
    
    end subroutine Free_XenonSamarium

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_XenonSamarium (this, xesm_inp, mat_info)
        
        class(XenonSamarium), intent(in out)  :: this
        type(XenonSamariumInput), intent(in)  :: xesm_inp
        type(Material), intent(in)            :: mat_info
        
        integer  :: iz, ia, im
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                im = mat_info%loading(iz, ia)
                
                this%Io(iz, ia)%lambda = xesm_inp%lambda_I
                this%Io(iz, ia)%gamma(:) = xesm_inp%gamma_I(:, im)
                
                this%Xe(iz, ia)%lambda = xesm_inp%lambda_Xe
                this%Xe(iz, ia)%gamma(:) = xesm_inp%gamma_Xe(:, im)
                
                this%Pm(iz, ia)%lambda = xesm_inp%lambda_Pm
                this%Pm(iz, ia)%gamma(:) = xesm_inp%gamma_Pm(:, im)
                
                this%Sm(iz, ia)%lambda = xesm_inp%lambda_Sm
                this%Sm(iz, ia)%gamma(:) = xesm_inp%gamma_Sm(:, im)
            end do
        end do
    
    end subroutine Set_XenonSamarium
    
    !$
    !===============================================================================================
    ! get steady-state number densities for Iodine/Xenon and Promethium/Samarium
    !===============================================================================================
    subroutine Steady_XenonSamarium (this, xsec, flux)
        
        class(XenonSamarium), intent(in out)  :: this
        type(CrossSection), intent(in)        :: xsec
        real(KREAL), intent(in)           :: flux(:, :, :)
        
        real(KREAL)  :: rate, micro
        integer          :: iz, ia, ig
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Io(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Io(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                this%Io(iz, ia)%density = rate / this%Io(iz, ia)%lambda
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Xe(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Xe(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                this%Xe(iz, ia)%density = (rate + this%Io(iz, ia)%lambda*this%Io(iz, ia)%density) / (this%Xe(iz, ia)%lambda + micro)
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Pm(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Pm(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                this%Pm(iz, ia)%density = rate / this%Pm(iz, ia)%lambda
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Sm(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Sm(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                this%Sm(iz, ia)%density = (rate + this%Pm(iz, ia)%lambda*this%Pm(iz, ia)%density) / (this%Sm(iz, ia)%lambda + micro)
            end do
        end do
        
    end subroutine Steady_XenonSamarium
    
    !$
    !===============================================================================================
    ! get transient number densities for Iodine/Xenon and Promethium/Samarium
    !===============================================================================================
    subroutine Transient_XenonSamarium (this, xsec, flux, step_length)
        
        class(XenonSamarium), intent(in out)  :: this
        type(CrossSection), intent(in)        :: xsec
        real(KREAL), intent(in)           :: flux(:, :, :)
        real(KREAL), intent(in)           :: step_length
        
        real(KREAL)  :: rate, micro, delta
        integer          :: iz, ia, ig
      
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Io(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Io(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                delta = rate - this%Io(iz, ia)%lambda*this%Io(iz, ia)%last_density
                this%Io(iz, ia)%density = this%Io(iz, ia)%last_density +  delta*step_length
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Xe(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Xe(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                delta = rate + this%Io(iz, ia)%lambda*this%Io(iz, ia)%last_density - this%Xe(iz, ia)%lambda*this%Xe(iz, ia)%last_density - micro*this%Xe(iz, ia)%last_density
                this%Xe(iz, ia)%density = this%Xe(iz, ia)%last_density +  delta*step_length
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Pm(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Pm(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                delta = rate - this%Pm(iz, ia)%lambda*this%Pm(iz, ia)%last_density
                this%Pm(iz, ia)%density = this%Pm(iz, ia)%last_density +  delta*step_length
                
                rate = 0.0
                micro = 0.0
                do ig = 1, ns%state%ng
                    rate = rate + this%Sm(iz, ia)%gamma(ig) * xsec%matrixs(iz, ia)%sigma_f(ig) * flux(iz, ia, ig)
                    micro = micro + this%Sm(iz, ia)%sigma_a(ig) * flux(iz, ia, ig)
                end do
                delta = rate + this%Pm(iz, ia)%lambda*this%Pm(iz, ia)%last_density - this%Sm(iz, ia)%lambda*this%Sm(iz, ia)%last_density - micro*this%Sm(iz, ia)%last_density
                this%Sm(iz, ia)%density = this%Sm(iz, ia)%last_density +  delta*step_length
            end do
        end do
    
    end subroutine Transient_XenonSamarium
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Update_XenonSamarium_xsec (this, xsec)
        
        class(XenonSamarium), intent(in out)  :: this
        type(CrossSection), intent(in out)    :: xsec
    
        real(KREAL)  :: sub, add
        integer  :: iz, ia, ig, ic
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                do ig = 1, ns%state%ng
                    sub = 0.0
                    sub = sub + this%Io(iz,ia)%last_sigma_a(ig) * this%Io(iz,ia)%last_density
                    sub = sub + this%Xe(iz,ia)%last_sigma_a(ig) * this%Xe(iz,ia)%last_density
                    sub = sub + this%Pm(iz,ia)%last_sigma_a(ig) * this%Pm(iz,ia)%last_density
                    sub = sub + this%Sm(iz,ia)%last_sigma_a(ig) * this%Sm(iz,ia)%last_density
                    
                    add = 0.0
                    add = add + this%Io(iz,ia)%sigma_a(ig) * this%Io(iz,ia)%density
                    add = add + this%Xe(iz,ia)%sigma_a(ig) * this%Xe(iz,ia)%density
                    add = add + this%Pm(iz,ia)%sigma_a(ig) * this%Pm(iz,ia)%density
                    add = add + this%Sm(iz,ia)%sigma_a(ig) * this%Sm(iz,ia)%density
                    xsec%matrixs(iz, ia)%sigma_t(ig) = xsec%matrixs(iz, ia)%sigma_t(ig) - sub + add
                end do
            end do
        end do
    
    end subroutine Update_XenonSamarium_xsec
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Transit_XenonSamarium (this)
        
        class(XenonSamarium), intent(in out)  :: this
        
        integer  :: iz, ia
        
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                call this%Io(iz, ia)%transit ()
                call this%Xe(iz, ia)%transit ()
                call this%Pm(iz, ia)%transit ()
                call this%Sm(iz, ia)%transit ()
            end do
        end do
            
    end subroutine Transit_XenonSamarium
    
end module xesm_header
