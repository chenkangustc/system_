!$
!===================================================================================================
!
!   class for transverse coefficient during iteration
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          SurfaceCoefficient
!                               NodalCoefficient
!
!===================================================================================================
module coefficient_header_transverse

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV

    use geometry_header,        only : Meshing, Geometry
    use material_header,        only : CrossSection
    use quadrature_header,      only : QuadratureSet
    
    implicit none 
    private
    public  :: SurfaceCoefficient, NodalCoefficient
    
    ! --------------------------------------------------------------------------
    ! type for order coefficient for radial--integrated surface flux
    type, private  :: rad_surf_tp
        real(KREAL), public, allocatable  :: first(:, :, :)
        real(KREAL), public, allocatable  :: second(:, :, :, :)
        real(KREAL), public, allocatable  :: third(:, :, :, :)
    contains
        procedure, private  :: alloc => Allocate_rad_surf_tp
        procedure, private  :: clean => Free_rad_surf_tp
    end type rad_surf_tp
    
    ! type for order coefficient for radial--integrated nodal flux
    type, private  :: rad_nodal_tp
        real(KREAL), public, allocatable  :: first(:, :, :, :)
        real(KREAL), public, allocatable  :: second(:, :, :, :, :)
        real(KREAL), public, allocatable  :: third(:, :, :, :, :)
        real(KREAL), public, allocatable  :: fourth(:, :, :, :)
    contains
        procedure, private  :: alloc => Allocate_rad_nodal_tp
        procedure, private  :: clean => Free_rad_nodal_tp
    end type rad_nodal_tp
    
    ! type for order coefficient for axial--integrated surface flux 
    type, private  :: axi_surf_tp
        real(KREAL), public, allocatable  :: first(:, :)
        real(KREAL), public, allocatable  :: second(:, :, :)
        real(KREAL), public, allocatable  :: third(:, :)
    contains
        procedure, private  :: alloc => Allocate_axi_surf_tp
        procedure, private  :: clean => Free_axi_surf_tp
    end type axi_surf_tp
    
    ! type for order coefficient for axial--integrated nodal flux
    type, private  :: axi_nodal_tp 
        real(KREAL), public, allocatable  :: first(:, :, :)
        real(KREAL), public, allocatable  :: second(:, :, :, :)
        real(KREAL), public, allocatable  :: third(:, :, :)
    contains
        procedure, private  :: alloc => Allocate_axi_nodal_tp
        procedure, private  :: clean => Free_axi_nodal_tp
    end type axi_nodal_tp
    
    ! --------------------------------------------------------------------------
    ! type for coefficient of surface
    type SurfaceCoefficient
        type(rad_surf_tp), public  :: rad
        type(axi_surf_tp), public  :: axi
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_SurfaceCoefficient
        procedure, public  :: set => Calcu_SurfaceCoefficient2
        procedure, public  :: clean =>  Free_SurfaceCoefficient
    end type SurfaceCoefficient
    
    ! type for coefficient of nodal average
    type NodalCoefficient
        type(rad_nodal_tp), public  :: rad
        type(axi_nodal_tp), public  :: axi
        real(KREAL), public  :: memory = REAL_ZERO
    contains
        procedure, public  :: alloc => Allocate_NodalCoefficient
        procedure, public  :: set => Calcu_NodalCoefficient2
        procedure, public  :: clean =>  Free_NodalCoefficient
    end type NodalCoefficient
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_rad_surf_tp (this)
        
        class(rad_surf_tp), intent(in out)  :: this
        integer  :: i_allocate
        integer, parameter  :: N_SURFACE = 3
        
        ! check allocated status first
        call this%clean ()
        
        ! 3:(surface per radial), 0-2(expansion order)
        allocate (this%first(N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate (this%second(2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate (this%third(0:2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
    
        this%first  = REAL_ZERO
        this%second = REAL_ZERO
        this%third  = REAL_ZERO
        
    end subroutine Allocate_rad_surf_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of rad_surf_tp
    !===============================================================================================
    subroutine Free_rad_surf_tp (this)
    
        class(rad_surf_tp), intent(in out)  :: this 
    
        if (allocated(this%first))      deallocate (this%first)
        if (allocated(this%second))     deallocate (this%second)
        if (allocated(this%third))      deallocate (this%third)
        
    end subroutine Free_rad_surf_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_rad_nodal_tp (this)
        
        class(rad_nodal_tp), intent(in out)  :: this
        integer  :: i_allocate
        integer, parameter  :: N_SURFACE = 3
        
        ! check allocated status first
        call this%clean ()
    
        allocate(this%first(2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%second(2, 2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%third(0:2, 2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%fourth(2, N_SURFACE, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        
        this%first  = REAL_ZERO
        this%second = REAL_ZERO
        this%third  = REAL_ZERO
        this%fourth = REAL_ZERO
        
    end subroutine Allocate_rad_nodal_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of rad_nodal_tp
    !===============================================================================================
    subroutine Free_rad_nodal_tp (this)
    
        class(rad_nodal_tp), intent(in out)  :: this 
    
        if (allocated(this%first))      deallocate (this%first)
        if (allocated(this%second))     deallocate (this%second)
        if (allocated(this%third))      deallocate (this%third)
        if (allocated(this%fourth))     deallocate (this%fourth)
        
    end subroutine Free_rad_nodal_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_axi_surf_tp (this)
        
        class(axi_surf_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%first(ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%second(2, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%third(ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
    
        this%first  = REAL_ZERO
        this%second = REAL_ZERO
        this%third  = REAL_ZERO
    
    end subroutine Allocate_axi_surf_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of axi_surf_tp
    !===============================================================================================
    subroutine Free_axi_surf_tp (this)
    
        class(axi_surf_tp), intent(in out)  :: this 
    
        if (allocated(this%first))      deallocate (this%first)
        if (allocated(this%second))     deallocate (this%second)
        if (allocated(this%third))      deallocate (this%third)
        
    end subroutine Free_axi_surf_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_axi_nodal_tp (this)
        
        class(axi_nodal_tp), intent(in out)  :: this
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%first(2, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%second(2, 2, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        allocate(this%third(2, ns%deduce%nodal_total, ns%deduce%direction), stat=i_allocate)
        
        this%first  = REAL_ZERO
        this%second = REAL_ZERO
        this%third  = REAL_ZERO
    
    end subroutine Allocate_axi_nodal_tp
    
    !$
    !===============================================================================================
    ! finalizer for class of axi_surf_tp
    !===============================================================================================
    subroutine Free_axi_nodal_tp (this)
    
        class(axi_nodal_tp), intent(in out)  :: this 
    
        if (allocated(this%first))      deallocate (this%first)
        if (allocated(this%second))     deallocate (this%second)
        if (allocated(this%third))      deallocate (this%third)
        
    end subroutine Free_axi_nodal_tp

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_SurfaceCoefficient (this)
        
        class(SurfaceCoefficient), intent(in out)  :: this
        
        ! check allocated status
        call this%clean ()
        
        call this%rad%alloc ()
        call this%axi%alloc ()

        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%first)
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%second)
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%third)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%first)
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%second)
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%third)
        
    end subroutine Allocate_SurfaceCoefficient
    
    !$
    !===============================================================================================
    ! finalizer for class of SurfaceCoefficient
    !===============================================================================================
    subroutine Free_SurfaceCoefficient (this)
        
        class(SurfaceCoefficient), intent(in out)  :: this
        
        call this%rad%clean ()
        call this%axi%clean ()
        
        this%memory = REAL_ZERO
    
    end subroutine Free_SurfaceCoefficient
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Allocate_NodalCoefficient (this)
        
        class(NodalCoefficient), intent(in out)  :: this
        
        ! check allocated status
        call this%clean ()
        
        call this%rad%alloc ()
        call this%axi%alloc ()
        
        this%memory = REAL_ZERO
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%first)
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%second)
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%third)
        this%memory = this%memory + REAL_BYTE * SIZE(this%rad%fourth)
        
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%first)
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%second)
        this%memory = this%memory + REAL_BYTE * SIZE(this%axi%third)
    
    end subroutine Allocate_NodalCoefficient
    
    !$
    !===============================================================================================
    ! finalizer for class of NodalCoefficient
    !===============================================================================================
    subroutine Free_NodalCoefficient (this)
        
        class(NodalCoefficient), intent(in out)  :: this
        
        call this%rad%clean ()
        call this%axi%clean ()
        
        this%memory = REAL_ZERO
    
    end subroutine Free_NodalCoefficient
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! calculate coefficitent on surface (radial and axial)
    !===============================================================================================
    subroutine Calcu_SurfaceCoefficient (this, ig, geom, mesh, xsec_iter, quad)
        
        class(SurfaceCoefficient), intent(in out)  :: this
        integer, intent(in)             :: ig                                   ! which energy group to be calculated ?
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        type(CrossSection), intent(in)  :: xsec_iter
        type(QuadratureSet), intent(in) :: quad
        
        real(KDOUBLE)  :: rad_f(0:2, 0:2), rad_g(0:2, 0:2)                      ! refer p51 of the dissertation
        real(KDOUBLE)  :: rad_p(0:2), rad_r(0:2)
        real(KDOUBLE)  :: axi_f(0:2, 0:2), axi_h(0:2), axi_p(0:2), axi_t
        
        real(KREAL)        :: dz                                                ! layer height
        real(KREAL)        :: ttc                                               ! macroscopic total cross section
        real(KREAL)        :: ux, uz                                            ! projection of specific direction
        integer            :: il, ill                                           ! layer index
        integer            :: k, j, j1, j2                                      ! surface index
        integer            :: ml, ia, ir, iz, iia, k1
        integer            :: is
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(is, ia, ml, dz, ir, il, iz, ttc, k, j, ux, j1, j2, k1, uz,    & 
        !$omp &  rad_f, rad_g, rad_p, rad_r, axi_f, axi_h, axi_p, axi_t)
        do is = 1, ns%deduce%direction
            do ia = 1, ns%state%layer
                ml = (ia-1) * ns%state%nodal
                dz = geom%height(ia)
                do ir = 1, ns%state%nodal
                    il = ir + ml
                    iz = mesh%zone(ir)
                    ttc = xsec_iter%matrixs(iz, ia)%sigma_t(ig)
                    
                    k = 0                                                       ! count for surface which ux>0
                    do j = 1, 3                                                 ! sweep for surface
                        ux = quad%directions(is)%projection(j, ir)
                        if (ux > 0.0D0) then                         
                            k = k + 1
                            if (k == 1)  j1 = j
                            if (k == 2)  j2 = j
                            call Coef_rad_plus (ux, ttc, rad_f, rad_g, rad_p, rad_r)
                            this%rad%first(j, il,is) = rad_p(0) / rad_f(0,0)
                            do k1 = 1, 2
                                this%rad%second(k1, j, il,is) = rad_p(k1) - rad_f(0,k1) * this%rad%first(j,il,is)
                            end do
                            do k1 = 0, 2
                                this%rad%third(k1, j, il,is) = rad_r(k1) - rad_g(k1,0) * this%rad%first(j,il,is)
                            end do
                        end if
                    end do
                    
                    uz = quad%directions(is)%xmu(3)
                    call Coef_axial (dz, ttc, uz, axi_f, axi_h, axi_p, axi_t)
                    this%axi%first(il,is) = axi_f(0, 0) / axi_p(0)
                    this%axi%third(il,is) = (axi_f(0, 0) + axi_h(0)*axi_p(0) - axi_f(0, 0)*axi_t) / axi_p(0)
                    this%axi%second(1, il,is) = axi_f(0,1) - axi_p(1)*this%axi%first(il,is)
                    this%axi%second(2, il,is) = axi_f(0,2) - axi_p(2)*this%axi%first(il,is)
                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel 
    
    end subroutine Calcu_SurfaceCoefficient
    
    !$
    !===============================================================================================
    ! calculate coeffcient in nodal (radiao and axial)
    !===============================================================================================
    subroutine Calcu_NodalCoefficient(this, ig, geom, mesh, xsec_iter, quad)
    
        class(NodalCoefficient), intent(in out)  :: this
        integer, intent(in)  :: ig                                              ! which energy group to be calculated ?
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        type(CrossSection), intent(in)  :: xsec_iter
        type(QuadratureSet), intent(in) :: quad
        
        real(KDOUBLE)  :: rad_f(0:2, 0:2), rad_g(0:2, 0:2), rad_h(0:2)          ! refer p51 of the dissertation
        real(KDOUBLE)  :: rad_p(0:2), rad_r(0:2)
        real(KDOUBLE)  :: axi_f(0:2, 0:2), axi_h(0:2), axi_p(0:2), axi_t
        
        real(KREAL)        :: dz                                                ! layer height
        real(KREAL)        :: ttc                                               ! macroscopic total cross section
        real(KREAL)        :: ux, uz                                            ! projection of specific direction
        integer            :: il, ill                                           ! layer index
        integer            :: k, j, j1, j2                                      ! surface index
        integer            :: ml, ia, ir, iia, iz, k1, k2
        integer            :: is
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(is, ia, ml, dz, ir, il, iz, ttc, k, j, ux, j1, j2, k1, k2, uz,   & 
        !$omp &  rad_f, rad_g, rad_p, rad_r, rad_h, axi_f, axi_h, axi_p, axi_t) 
        do is = 1, ns%deduce%direction
            do ia = 1, ns%state%layer
                ml = (ia-1) * ns%state%nodal
                dz = geom%height(ia)
                do ir = 1, ns%state%nodal
                    il = ir + ml
                    iz = mesh%zone(ir)
                    ttc = xsec_iter%matrixs(iz, ia)%sigma_t(ig)
                    
                    k = 0                                                       ! count for surface which ux>0
                    surface: do j = 1, 3                                        ! sweep for surface
                        ux = quad%directions(is)%projection(j, ir)
                        if (ux > 0.0D0)  then                         
                            k = k + 1
                            if (k == 1)  j1 = j
                            if (k == 2)  j2 = j
                            call Coef_rad_plus (ux, ttc, rad_f, rad_g, rad_p, rad_r)
                            do k1 = 1, 2
                                this%rad%first(k1, j, il,is) = rad_f(k1, 0) / rad_f(0, 0)
                                do k2 = 1, 2
                                    this%rad%second(k2, k1, j, il,is) = rad_f(k1, k2) - rad_f(0, k2)*this%rad%first(k1, j, il,is)
                                end do
                                do k2 = 0, 2
                                    this%rad%third(k2, k1, j, il,is) = rad_g(k2, k1) - rad_g(k2, 0)*this%rad%first(k1, j, il,is)
                                end do
                            end do
                        else
                            call Coef_rad_minus (ux, ttc, rad_f, rad_g, rad_h)
                            do k1 = 1, 2
                                this%rad%first(k1, j, il,is) = rad_f(k1, 0) / rad_f(0, 0)
                                this%rad%fourth(k1, j, il,is) = rad_h(k1) - rad_h(0)*this%rad%first(k1, j, il,is)
                                do k2 = 1, 2
                                    this%rad%second(k2, k1, j, il,is) = rad_f(k1, k2) - rad_f(0, k2)*this%rad%first(k1, j, il,is)
                                end do
                                do k2 = 0, 2
                                    this%rad%third(k2, k1, j, il,is) = rad_g(k2, k1) - rad_g(k2, 0)*this%rad%first(k1, j, il,is)
                                end do
                            end do
                        end if
                    end do surface
                    
                    uz = quad%directions(is)%xmu(3)
                    call Coef_axial (dz, ttc, uz, axi_f, axi_h, axi_p, axi_t)
                    do k1 = 1, 2
                        this%axi%first(k1, il,is) = axi_f(k1, 0) / axi_f(0, 0)
                        this%axi%third(k1, il,is) = axi_h(k1) - axi_h(0)*this%axi%first(k1, il,is)
                        this%axi%second(1, k1, il,is) = axi_f(k1, 1) - axi_f(0, 1)*this%axi%first(k1, il,is)
                        this%axi%second(2, k1, il,is) = axi_f(k1, 2) - axi_f(0, 2)*this%axi%first(k1, il,is)
                    end do
                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel 
    
    end subroutine Calcu_NodalCoefficient
    
    !$
    !===============================================================================================
    ! calculate coefficitent on surface (radial and axial)
    !   is_radial =.TRUE./.FALSE., perform radial calculate ?
    !   is_axial  =.TRUE./.FALSE., perform axial calculate ?   
    !   if sigma_t is the same, do not calculte radial coefficient
    !   id sigma_t and height is the same, do not calculate axial coefficient
    !===============================================================================================
    subroutine Calcu_SurfaceCoefficient2 (this, ig, geom, mesh, xsec_iter, quad)
        
        class(SurfaceCoefficient), intent(in out)  :: this
        integer, intent(in)             :: ig                                   ! which energy group to be calculated ?
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        type(CrossSection), intent(in)  :: xsec_iter
        type(QuadratureSet), intent(in) :: quad
        
        ! local variables
        real(KDOUBLE)  :: rad_f(0:2, 0:2), rad_g(0:2, 0:2)                      ! refer p51 of the dissertation
        real(KDOUBLE)  :: rad_p(0:2), rad_r(0:2)
        real(KDOUBLE)  :: axi_f(0:2, 0:2), axi_h(0:2), axi_p(0:2), axi_t
        
        logical            :: is_radial                                         ! whether perform radial coefficient calculation ?
        logical            :: is_axial                                          ! whether perform axial coefficient calculation ?
        integer            :: equal_radial                                      ! equal layer index for radial coefficient
        integer            :: equal_axial                                       ! equal layer index for axial coefficient
        
        real(KREAL)        :: dz, dza                                           ! layer height
        real(KREAL)        :: ttc, ttca                                         ! macroscopic total cross section
        real(KREAL)        :: ux, uz                                            ! projection of specific direction
        integer            :: il, ill                                           ! layer index
        integer            :: k, j, j1, j2                                      ! surface index
        integer            :: ml, ia, ir, iz, iia, k1
        integer            :: is
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(is, is_radial, is_axial, equal_radial, equal_axial,  &
        !$omp &  ia, ml, dz, ir, il, iz, ttc, iia, ttca, dza, k, j, ux, j1, j2, k1, ill, uz,              & 
        !$omp &  rad_f, rad_g, rad_p, rad_r, axi_f, axi_h, axi_p, axi_t) 
        do is = 1, ns%deduce%direction
            equal_radial = 1
            do ia = 1 , ns%state%layer
                ml = (ia-1) * ns%state%nodal
                dz = geom%height(ia)
                do ir = 1, ns%state%nodal
                    il = ir + ml
                    iz = mesh%zone(ir)
                    ttc = xsec_iter%matrixs(iz, ia)%sigma_t(ig)
                    ! whether perform coefficient calculationg(radial and axial) for a nodal ?
                    ! is NOT, get the layer index of the equal nodal
                    
                    ! the first layer should be calculated
                    if (ia == 1) then      
                        is_radial = .TRUE.
                        is_axial  = .TRUE.
                    else
                        do iia = 1, ia-1
                            ttca = xsec_iter%matrixs(iz, iia)%sigma_t(ig)
                            if (ABS(ttca - ttc) < EPS_EQUAL)  then
                                is_radial = .FALSE.
                                equal_radial = iia
                                exit
                            else
                                is_radial = .TRUE.
                            end if
                        end do
                        do iia = equal_radial, ia-1
                            ttca = xsec_iter%matrixs(iz, iia)%sigma_t(ig)
                            dza = geom%height(iia)
                            if (ABS(ttca- ttc)<EPS_EQUAL .and. ABS(dza - dz)<EPS_EQUAL)  then
                                is_axial = .FALSE.
                                equal_axial = iia
                                exit
                            else
                                is_axial = .TRUE.
                            end if
                        end do
                    end if
            
                    ! --------------------------------------------------------------
                    ! if calculate coefficient for radial, invoke subroutine to get reslut
                    if (is_radial)  then
                        k = 0                                                   ! count for surface which ux>0
                        do j = 1, 3                                             ! sweep for surface
                            ux = quad%directions(is)%projection(j, ir)
                            if (ux > 0.0D0) then                         
                                k = k + 1
                                if (k == 1)  j1 = j
                                if (k == 2)  j2 = j
                                call Coef_rad_plus (ux, ttc, rad_f, rad_g, rad_p, rad_r)
                                
                                this%rad%first(j, il,is) = rad_p(0) / rad_f(0,0)
                                do k1 = 1, 2
                                    this%rad%second(k1, j, il,is) = rad_p(k1) - rad_f(0,k1) * this%rad%first(j,il,is)
                                end do
                                do k1 = 0, 2
                                    this%rad%third(k1, j, il,is) = rad_r(k1) - rad_g(k1,0) * this%rad%first(j,il,is)
                                end do
                            end if
                        end do
                    ! if do NOT calculate, just equal ill layer to il
                    else
                        ill = (equal_radial-1) * ns%state%nodal + ir
                        this%rad%first(:, il,is) = this%rad%first(:, ill,is)
                        this%rad%second(:, :, il,is) = this%rad%second(:, :, ill,is)
                        this%rad%third(:, :, il,is) = this%rad%third(:, :, ill,is)
                    end if
                    
                    if (is_axial)  then
                        uz = quad%directions(is)%xmu(3)
                        call Coef_axial (dz, ttc, uz, axi_f, axi_h, axi_p, axi_t)
                        
                        this%axi%first(il,is) = axi_f(0, 0) / axi_p(0)
                        this%axi%third(il,is) = (axi_f(0, 0) + axi_h(0)*axi_p(0) - axi_f(0, 0)*axi_t) / axi_p(0)
                        this%axi%second(1, il,is) = axi_f(0,1) - axi_p(1)*this%axi%first(il,is)
                        this%axi%second(2, il,is) = axi_f(0,2) - axi_p(2)*this%axi%first(il,is)
                    else
                        ill = (equal_axial-1) * ns%state%nodal + ir
                        this%axi%first(il,is) = this%axi%first(ill,is)
                        this%axi%third(il,is) = this%axi%third(ill,is)
                        this%axi%second(:, il,is)  = this%axi%second(:, ill,is)
                    end if
                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel 
    
    end subroutine Calcu_SurfaceCoefficient2
    
    !$
    !===============================================================================================
    ! calculate coeffcient in nodal (radiao and axial)
    !   is_radial =.TRUE./.FALSE., perform radial calculate ?
    !   is_axial  =.TRUE./.FALSE., perform axial calculate ?   
    !   if sigma_t is the same, do not calculte radial coefficient
    !   id sigma_t and height is the same, do not calculate axial coefficient
    !===============================================================================================
    subroutine Calcu_NodalCoefficient2(this, ig, geom, mesh, xsec_iter, quad)
    
        class(NodalCoefficient), intent(in out)  :: this
        integer, intent(in)  :: ig                                              ! which energy group to be calculated ?
        type(Meshing), intent(in)       :: mesh
        type(Geometry), intent(in)      :: geom
        type(CrossSection), intent(in)  :: xsec_iter
        type(QuadratureSet), intent(in) :: quad
        
        ! local variables
        real(KDOUBLE)  :: rad_f(0:2, 0:2), rad_g(0:2, 0:2), rad_h(0:2)          ! refer p51 of the dissertation
        real(KDOUBLE)  :: rad_p(0:2), rad_r(0:2)
        real(KDOUBLE)  :: axi_f(0:2, 0:2), axi_h(0:2), axi_p(0:2), axi_t
        
        logical            :: is_radial                                         ! whether perform radial coefficient calculation ?
        logical            :: is_axial                                          ! whether perform axial coefficient calculation ?
        integer            :: equal_radial                                      ! equal layer index for radial coefficient
        integer            :: equal_axial                                       ! equal layer index for axial coefficient
        
        real(KREAL)        :: dz, dza                                           ! layer height
        real(KREAL)        :: ttc, ttca                                         ! macroscopic total cross section
        real(KREAL)        :: ux, uz                                            ! projection of specific direction
        integer            :: il, ill                                           ! layer index
        integer            :: k, j, j1, j2                                      ! surface index
        integer            :: ml, ia, ir, iia, iz, k1, k2
        integer            :: is
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(is, is_radial, is_axial, equal_radial, equal_axial,      &
        !$omp &  ia, ml, dz, ir, il, iz, ttc, iia, ttca, dza, k, j, ux, j1, j2, k1, k2, ill, uz,              & 
        !$omp &  rad_f, rad_g, rad_p, rad_r, rad_h, axi_f, axi_h, axi_p, axi_t)                                
        do is = 1, ns%deduce%direction
            equal_radial = 1
            do ia = 1 , ns%state%layer
                ml = (ia-1) * ns%state%nodal
                dz = geom%height(ia)
                do ir = 1, ns%state%nodal
                    il = ir + ml
                    iz = mesh%zone(ir)
                    ttc = xsec_iter%matrixs(iz, ia)%sigma_t(ig)
                    ! whether perform coefficient calculationg(radial and axial) for a nodal ?
                    ! is NOT, get the layer index of the equal nodal
                    
                    ! the first layer should be calculated
                    if (ia == 1)  then      
                        is_radial = .TRUE.
                        is_axial  = .TRUE.
                    else
                        do iia = 1, ia-1
                            ttca = xsec_iter%matrixs(iz, iia)%sigma_t(ig)
                            if (ABS(ttca - ttc) < EPS_EQUAL)  then
                                is_radial = .FALSE.
                                equal_radial = iia
                                exit
                            else
                                is_radial = .TRUE.
                            end if
                        end do
                        do iia = equal_radial, ia-1
                            ttca = xsec_iter%matrixs(iz, iia)%sigma_t(ig)
                            dza = geom%height(iia)
                            if (ABS(ttca- ttc)<EPS_EQUAL .and. ABS(dza - dz)<EPS_EQUAL)  then
                                is_axial = .FALSE.
                                equal_axial = iia
                                exit
                            else
                                is_axial = .TRUE.
                            end if
                        end do
                    end if
            
                    ! --------------------------------------------------------------
                    ! if calculate coefficient for radial, invoke subroutine to get reslut
                    if (is_radial)  then
                        k = 0                                                   ! count for surface which ux>0
                        surface: do j = 1, 3                                    ! sweep for surface
                            ux = quad%directions(is)%projection(j, ir)
                            if (ux > 0.0D0)  then                         
                                k = k + 1
                                if (k == 1)  j1 = j
                                if (k == 2)  j2 = j
                                call Coef_rad_plus (ux, ttc, rad_f, rad_g, rad_p, rad_r)
                                
                                do k1 = 1, 2
                                    this%rad%first(k1, j, il,is) = rad_f(k1, 0) / rad_f(0, 0)
                                    do k2 = 1, 2
                                        this%rad%second(k2, k1, j, il,is) = rad_f(k1, k2) - rad_f(0, k2)*this%rad%first(k1, j, il,is)
                                    end do
                                    do k2 = 0, 2
                                        this%rad%third(k2, k1, j, il,is) = rad_g(k2, k1) - rad_g(k2, 0)*this%rad%first(k1, j, il,is)
                                    end do
                                end do
                            else
                                call Coef_rad_minus (ux, ttc, rad_f, rad_g, rad_h)
                                
                                do k1 = 1, 2
                                    this%rad%first(k1, j, il,is) = rad_f(k1, 0) / rad_f(0, 0)
                                    this%rad%fourth(k1, j, il,is) = rad_h(k1) - rad_h(0)*this%rad%first(k1, j, il,is)
                                    do k2 = 1, 2
                                        this%rad%second(k2, k1, j, il,is) = rad_f(k1, k2) - rad_f(0, k2)*this%rad%first(k1, j, il,is)
                                    end do
                                    do k2 = 0, 2
                                        this%rad%third(k2, k1, j, il,is) = rad_g(k2, k1) - rad_g(k2, 0)*this%rad%first(k1, j, il,is)
                                    end do
                                end do
                            end if
                        end do surface
                    ! if do NOT calculate, just equal ill layer to il
                    else
                        ill = (equal_radial-1) * ns%state%nodal + ir
                        this%rad%first(:, :, il,is) = this%rad%first(:, :, ill,is)
                        this%rad%second(:, :, :,il,is) = this%rad%second(:, :, :, ill,is)
                        this%rad%third(:, :,:, il,is) = this%rad%third(:, :, :, ill,is)
                        this%rad%fourth(:, :, il,is) = this%rad%fourth(:, :, ill,is)
                    end if
                    
                    ! --------------------------------------------------------------
                    ! if calculate coefficient for axial, ivoke subroutine to get result
                    if (is_axial)  then
                        uz = quad%directions(is)%xmu(3)
                        call Coef_axial (dz, ttc, uz, axi_f, axi_h, axi_p, axi_t)
                        
                        do k1 = 1, 2
                            this%axi%first(k1, il,is) = axi_f(k1, 0) / axi_f(0, 0)
                            this%axi%third(k1, il,is) = axi_h(k1) - axi_h(0)*this%axi%first(k1, il,is)
                            this%axi%second(1, k1, il,is) = axi_f(k1, 1) - axi_f(0, 1)*this%axi%first(k1, il,is)
                            this%axi%second(2, k1, il,is) = axi_f(k1, 2) - axi_f(0, 2)*this%axi%first(k1, il,is)
                        end do
                    ! if do NOT calculate, just equal ill layer to il
                    else
                        ill = (equal_axial-1) * ns%state%nodal + ir
                        this%axi%first(:, il,is) = this%axi%first(:, ill,is)
                        this%axi%third(:, il,is) = this%axi%third(:, ill,is)
                        this%axi%second(:, :, il,is) = this%axi%second(:,:,ill,is)
                    end if
                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel 
    
    end subroutine Calcu_NodalCoefficient2
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    pure subroutine Coef_rad_plus (ux, tc, F, G, P, R)
    
        real(KDOUBLE), intent(in out)  :: F(0:2,0:2), G(0:2,0:2), P(0:2), R(0:2)  ! F(j,i)
        real(KREAL), intent(in)      :: ux
        real(KREAL), intent(in)      :: tc                                  ! macro scopic total cross section
        ! page 113 of the dissertation
        ! a=>F  b=>G 
        ! c=>P  d=>R
        ! Fi,j => F(i,j); Gi,j => G(j,i)
        
        real(KDOUBLE) :: aa, ea, a1, a2, a3, a4, a5, a6
       
        if ( tc < EPS_ZERO ) then
            aa = 1.0D0 / ux
                   
            F = aa
            G = aa
            P = aa
            R = aa
        else
            aa = tc / ux
            ea = EXP(-aa)
            a1 = 1 / (tc*aa)
            a2 = a1 / aa
            a3 = a2 / aa
            a4 = a3 / aa
            a5 = a4 / aa
            a6 = a5 / aa
        
            F(0,0) = (aa*(aa-2)+2*(1-ea))*a2
            F(0,1) = ((2*(2+ea)-aa)*aa+6*(ea-1))*2.0/3.0*a3
            F(0,2) = ((((3-ea)*3-aa)*aa-12*(3+2*ea))*aa+60*(1-ea))*0.2*a4
            F(1,0) = -F(0,1)*9
            F(1,1) = (((4*(1+2*ea)+aa*(aa-4))*aa+24*(1+2*ea))*aa+72*(ea-1))*a4
            F(1,2) = ((((3*(2-ea)-aa)*aa-3*(3+11*ea))*aa-12.0*(4+11*ea))*aa+180*(1-ea))*1.2*a5
            F(2,0) = F(0,2)*100
            F(2,1) = -100/9.0*F(1,2)
            F(2,2) = (((((aa*(aa-6)+6*(1-3*ea))*aa+96*(1-3*ea))*aa-144*(2+13*ea))*aa-1440*(1+4*ea))*aa+7200*(1-ea))*a6
        
            G(0,0) = -2*(-1+aa+ea)*a1
            G(0,1) = -(6*(ea-1)+aa*(-aa+2*(ea+2)))*6*a2
            G(0,2) = -(((aa+3*(ea-3))*aa+12*(2*ea+3))*aa+60*(ea-1))*20*a3
            G(1,0) = -(((1+2*ea)*2+aa)*(-aa)+6*(1-ea))/3.0*a2
            G(1,1) = -(((2*aa-5-4*ea)*aa-18*ea)*aa+18*(1-ea))*2*a3
            G(1,2) = -(((((5-2*ea)*3-2*aa)*aa-3*(15+19*ea))*aa+12*(1-16*ea))*aa+180*(1-ea))*20/3.0*a4
            G(2,0) = -(((aa+4*ea-1)*aa+6*(1+2*ea))*aa+18*(ea-1))*2/9.0*a3
            G(2,1) = -((((-5*aa+8*(1+2*ea))*aa+12*(1+8*ea))*aa+216*ea)*aa+216*(ea-1))/3.0*a4
            G(2,2) = -(((((2*aa+6*(ea-2))*aa+3*(9+22*ea))*aa+3*(5+97*ea))*aa+36*(16*ea-1))*aa+540*(ea-1))*40/9.0*a5
        
            P(0) = (aa-1+ea)*a1
            P(1) = ((aa-2*(2+ea))*aa+6*(1-ea))/3.0*a2
            P(2) = (((aa+3*(ea-3))*aa+12*(3+2*ea))*aa+60*(ea-1))*0.1*a3
               
            R(0) = -(1-ea)/tc
            R(1) = -(aa*(2*ea+1)+3*(ea-1))/3.0*a1
            R(2) = ((aa*(4*ea-1)+6*(2*ea+1))*aa+18*(ea-1))/9.0*a2
        end if 

    end subroutine Coef_rad_plus
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    pure subroutine Coef_rad_minus(ux, tc, F, G, H)    

        real(KDOUBLE), intent(in out)  :: F(0:2,0:2), G(0:2,0:2), H(0:2)
        real(KREAL), intent(in)      :: ux, tc
        ! a=>F  b=>G  c=>H
        ! Fi,j => F(i,j); Gi,j => G(j,i)

        ! use precision number double than usual
        real(KDOUBLE) :: aa, ea, a1, a2, a3, a4, a5, a6

        if ( tc < EPS_ZERO ) then
            aa = 1.0D0 / abs(ux)
            
            F = aa
            G = aa
            H = 1.0D0
        else
            aa = tc / ABS(ux) 
            ea = EXP(-aa)
            a1 = 1 / (aa*tc)
            a2 = a1 / aa
            a3 = a2 / aa
            a4 = a3 / aa
            a5 = a4 / aa
            a6 = a5 / aa
        
            F(0,0) = ((aa+2*ea)*aa+2*(ea-1))*a2
            F(0,1) = (((2*ea+1)*2+aa*ea)*2*aa+12*(ea-1))/3.0*a3
            F(0,2) = ((((3*ea-1)*3+aa*ea)*aa+12*(2+3*ea))*aa+60*(ea-1))*0.2*a4
            F(1,0) = -(((5*ea-2)*2+aa*(4*ea+3))*6*aa+36*(ea-1))*a3
            F(1,1) = (((aa*(aa-8*ea)-4*(11*ea+4))*aa+24*(1-4*ea))*aa+72*(1-ea))*a4
            F(1,2) = -((((2*aa*ea+3*(7*ea-2))*aa+3*(33*ea+13))*aa+12*(19*ea-4))*aa+180*(ea-1))*1.2*a5
            F(2,0) = -(((aa*(2-9*ea)-3*(15*ea+7))*aa+12*(3-8*ea))*aa                    &
                &  +60*(1-ea))*20*a4
            F(2,1) = -((((aa*(5-9*ea)-18*(1+4*ea))*aa-6*(3+43*ea))*aa                   &
                &  +24*(4-19*ea))*aa+360*(1-ea))*20/3.0*a5
            F(2,2) = ((((((aa+18*ea)*aa+18*(13*ea-3))*aa+24*(ea*59+9))*aa               &
                &  +144*(33*ea+2))*aa+1440*(6*ea-1))*aa+7200*(ea-1))*a6
        
            G(0,0) = 2*(1-aa-ea)*a1
            G(0,1) = ((aa+2*(1+2*ea))*aa+6*(ea-1))*6*a2
            G(0,2) = (((-aa+3*(1-3*ea))*aa-12*(2+3*ea))*aa+60*(1-ea))*20*a3
            G(1,0) = ((aa-2*(2+ea))*aa+6*(1-ea))/3.0*a2
            G(1,1) = (((-2*aa+5+4*ea)*aa+18*ea)*aa+18*(ea-1))*2*a3
            G(1,2) = ((((2*aa-9*ea)*aa-3*(5+21*ea))*aa-12*(1+14*ea))*aa+180*(1-ea))*20/3.0*a4
            G(2,0) = (((-aa+4-ea)*aa-6*(2+ea))*aa+18*(1-ea))*2/9.0*a3
            G(2,1) = ((((5*aa+4*(2*ea-5))*aa+12*(4+5*ea))*aa+216*ea)*aa+216*(ea-1))/3.0*a4
            G(2,2) = (((((-4*aa+3*(5-3*ea))*aa-6*(4+15*ea))*aa-6*(5+73*ea))*aa-72*(1+14*ea))*aa+1080*(1-ea))*20/9.0*a5
        
            H(0) = 2*(1-ea)/aa
            H(1) = 12*(aa*(2*ea+1)+3*(ea-1))/(aa*aa)
            H(2) = -60*(((3*ea-1)*aa+4*(3*ea+2))*aa+20*(ea-1))/aa**3
        end if 
        
    end subroutine Coef_rad_minus

    !$
    !===============================================================================================
    !
    !===============================================================================================
    pure subroutine Coef_axial(dr, tc, uz, F, H, P, T)
        real(KDOUBLE), intent(in out)  :: F(0:2,0:2), H(0:2), P(0:2), T      
        real(KREAL), intent(in)      :: dr, tc, uz
        ! a=>F  b=>H  c=>P
        ! Fi,j => F(i,j)
        ! Hi/(dr/ABS(uz)) => H(i)
        
        ! use precision number double than usual
        real(KDOUBLE) :: bb,eb,b1,b2,b3,b4,b5
        real(KDOUBLE) :: k
        
        if ( tc < EPS_ZERO ) then
            k = uz / ABS(uz)
            bb = dr / ABS(uz)
                   
            F = bb
            P = bb
            H = 1.0D0
            T = 1.0D0
        else
            k = uz / ABS(uz)
            bb = dr*tc / ABS(uz)
            eb = EXP(-bb)
            b1 = 1.0D0 / (bb*tc)
            b2 = b1 / bb
            b3 = b2 / bb
            b4 = b3 / bb
            b5 = b4 / bb
        
            F(0,0) = (-1.0D0+eb+bb)*b1
            F(0,1) = k*(2.0D0*(1.0D0-eb)-bb*(1.0D0+eb))*0.5D0*b2
            F(0,2) = (bb*(bb*(eb-1.0D0)+6.0D0*(1.0D0+eb))+12.0D0*(eb-1.0D0))/6.0D0*b3
            F(1,0) = -F(0,1)*12.0D0
            F(1,1) = ((bb*(bb-3.0D0*(1.0D0+eb))-12.0D0*eb)*bb+12.0D0*(1.0D0-eb))*b3
            F(1,2) = k*((((eb-1.0D0)*bb+4.0D0*(1.0D0+2.0D0*eb))*bb+24.0D0*eb)*bb+24.0D0*(eb-1.0D0))*b4
            F(2,0) = F(0,2)*180.0D0
            F(2,1) = -F(1,2)*15.0D0
            F(2,2) = (((((bb+5.0D0*(eb-1.0D0))*bb+60.0D0*eb)*bb+60.0D0*(1.0D0+5.0D0*eb))*bb+720.0D0*eb)*bb+720.0D0*(eb-1.0D0))*b5
        
            H(0) = (1.0D0-eb)/bb
            H(1) = k*(2.0D0*(1.0D0-eb)-bb*(1.0D0+eb))*6.0D0/bb/bb
            H(2) = ((bb*(1.0D0-eb)-6.0D0*(1.0D0+eb))*bb+12.0D0*(1.0D0-eb))*30.0D0/bb/bb/bb
               
            P(0) = (1.0D0-eb)/tc
            P(1) = k*(2.0D0*(eb-1.0D0)+bb*(1.0D0+eb))*0.5D0*b1
            P(2) = ((bb*(1.0D0-eb)-6.0D0*(1.0D0+eb))*bb+12.0D0*(1.0D0-eb))/6.0D0*b2
        
            T = eb
        end if 
       
    end subroutine Coef_axial 

end module coefficient_header_transverse
