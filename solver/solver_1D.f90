!$
!===================================================================================================
!
!   class for one dimensional solver
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    onedsolver_meshrefining
!                               onedsolver_spacialsweep
!                               onedsolver_diffusionfdm
!                               onedsolver_snfdm
!                               onedsolver_fluxcondensation
!                               onedsolver_multiplication
!                               
!   Public type lists:          No
!
!===================================================================================================
module solver_1D
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: onedsolver_meshrefining, onedsolver_spacialsweep
    public  :: onedsolver_diffusionfdm, onedsolver_snfdm
    public  :: onedsolver_fluxcondensation, onedsolver_multiplication
    
    public  :: onedsolver_type
    
    ! --------------------------------------------------------------------------
    type  onedsolver_type
        logical                       :: defined = .FALSE.                      ! index of the object array allocation state
        integer                       :: number_matregion                       ! number of total material regions
        integer,allocatable           :: index_matregion(:)                     ! corresponding material index 
        real(KREAL), allocatable  :: size_matregion(:)                      ! size of each mesh
        real(KREAL), allocatable  :: diffcoef_matregion(:, :)               ! diffusion coefficient of each material region
        real(KREAL), allocatable  :: remvxs_matregion(:, :)                 ! removal xs of each material region
        real(KREAL), allocatable  :: toxs_matregion(:, :)                   ! total xs of each material region     
        real(KREAL), allocatable  :: scatxs_matregion(:, :, :)              ! scattering xs of each material region
        real(KREAL), allocatable  :: nuproxs_matregion(:, :)                ! neutrion production xs of each material region
    end type onedsolver_type

    integer                         :: number_group                             ! number of energy groups
    integer                         :: order_sn                                 ! order of sn
    integer                         :: polar_angle                              ! number of polar angles for sn calculation
    integer                         :: MODe                                     ! calculation MODe   
    integer                         :: number_outputregion                      ! number of flux condensation regions
    real(KREAL), allocatable    :: boundary_condition(:)                    ! boundary condtions : 0/1 vaccum/reflection
    real(KREAL), allocatable    :: kappa(:)                                 ! fission spectrum
    real(KREAL), allocatable    :: flux_condensed(:, :)                     ! condensed flux values of each condensation region
    real(KREAL), allocatable    :: flux_finemesh(:, :)                      ! fine mesh flux to be solved
    real(KREAL), allocatable    :: leakage_finemesh(:, :)                   ! transverse leakage item of the fine mesh
    
contains 
    !$
    !===============================================================================================
    ! refine meshes from the provided nodal informations and correspond sizes and material xss
    !===============================================================================================
    subroutine onedsolver_meshrefining (coarseonedsolver_type,fineonedsolver_type,number_group)    

        type(onedsolver_type), intent(in) :: coarseonedsolver_type                  ! data structure for the coarse nodal mesh
        type(onedsolver_type), intent(in out) :: fineonedsolver_type                ! data structure for the fine mesh 
        integer,intent(in) :: number_group

        integer  :: i,j,k,jj,ig1,ig2,nn
        integer  :: totalfinemesh
    
        totalfinemesh = 0
        nn = 1
        
        ! calculate the number of fine meshes
        do i = 1, coarseonedsolver_type%number_matregion
            totalfinemesh = totalfinemesh + CEILING(coarseonedsolver_type%size_matregion(i))*nn
        end do
    
        call onedsolver_define (fineonedsolver_type,totalfinemesh,number_group)
        
        k = 1
        j = 1
        jj = 1
        do i = 1, coarseonedsolver_type%number_matregion
            do j = 1, CEILING(coarseonedsolver_type%size_matregion(i))*nn
            
                ! material index of each fine mesh
                fineonedsolver_type%index_matregion(jj)=coarseonedsolver_type%index_matregion(i)
                ! sizes of fine meshes
                fineonedsolver_type%size_matregion(jj)=coarseonedsolver_type%size_matregion(i)/CEILING(coarseonedsolver_type%size_matregion(i))/nn
                ! xs informations
                do ig1 = 1, number_group
                    fineonedsolver_type%diffcoef_matregion(jj,ig1)=coarseonedsolver_type%diffcoef_matregion(i,ig1)
                    fineonedsolver_type%remvxs_matregion(jj,ig1)=coarseonedsolver_type%remvxs_matregion(i,ig1)  
                    fineonedsolver_type%toxs_matregion(jj,ig1)=coarseonedsolver_type%toxs_matregion(i,ig1)  
                    fineonedsolver_type%nuproxs_matregion(jj,ig1)=coarseonedsolver_type%nuproxs_matregion(i,ig1) 
                    do ig2 = 1, number_group      
                        fineonedsolver_type%scatxs_matregion(jj,ig1,ig2)=coarseonedsolver_type%scatxs_matregion(i,ig1,ig2)
                    end do
                end do
                jj=jj+1
            end do
            k=k+CEILING(coarseonedsolver_type%size_matregion(i))
        end do      
    
    end subroutine onedsolver_meshrefining
    
    !$
    !===============================================================================================
    ! spacial sweep scheme for the fdm
    !===============================================================================================
    subroutine onedsolver_spacialsweep (fineonedsolver_type,flux1d,source,number_meshes,ig, boundary_condition)   

        type(onedsolver_type),intent(in) :: fineonedsolver_type
        real(KREAL), intent(in)      :: source,boundary_condition           ! values of sources
        integer,intent(in)               :: number_meshes,ig                    ! number of meshes,bc,group index
        real(KREAL), intent(out)     :: flux1d                              ! 1-dimensional flux distribution
        dimension                        :: flux1d(number_meshes), source(number_meshes),boundary_condition(2)
        
        ! local 
        real(KREAL), allocatable :: coefficient_matrix(:), output_matrix(:), a(:,:),b(:,:),source2(:),coefficient_matrix2(:),cc(:),dd(:)
        real(KREAL) :: d1,d2,d3,h1,h2,h3,sigr,beta
        integer         :: i,j,k,l
    
        ! initializations
        allocate(coefficient_matrix(3*number_meshes-2),coefficient_matrix2(3*number_meshes-2))
        
        flux1d             = 0.0 
        coefficient_matrix = 0.0 
        d1                 = 0.0
        d2                 = 0.0
        d3                 = 0.0
        h1                 = 0.0
        h2                 = 0.0
        h3                 = 0.0
        sigr               = 0.0
        
        if (number_meshes <= 2) then
            write(*,10)
            stop
            
        else
            ! left boundary. the first element of the matrix
            d1 = fineonedsolver_type%diffcoef_matregion(1,ig)
            d2 = fineonedsolver_type%diffcoef_matregion(2,ig)
            h1 = fineonedsolver_type%size_matregion(1)
            h2 = fineonedsolver_type%size_matregion(2)
            beta = boundary_condition(1)
            
            coefficient_matrix(1) = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d1*(1.0-beta) /                            &
                &   (4.0*d1*(1.0+beta)+h1*(1.0-beta))+fineonedsolver_type%remvxs_matregion(1,ig)*h1
            coefficient_matrix(2) = -2.0*(d2*d1/(d2*h1+d1*h2))
            
            flux1d(1) = source(1)*h1
            do i = 3, 3*(number_meshes-1)-1
                k = i/3 + 1
                d1 = fineonedsolver_type%diffcoef_matregion(k-1,ig) 
                d2 = fineonedsolver_type%diffcoef_matregion(k,ig)
                d3 = fineonedsolver_type%diffcoef_matregion(k+1,ig)
                h1 = fineonedsolver_type%size_matregion(k-1)  
                h2 = fineonedsolver_type%size_matregion(k)      
                h3 = fineonedsolver_type%size_matregion(k+1)              
                j = MOD(i,3)
                select case(j)
                case(0)
                    coefficient_matrix(i) = -2.0*(d1*d2/(d1*h2+d2*h1))
                case(1)
                    coefficient_matrix(i) = 2.0*(d1*d2/(d1*h2+d2*h1)+d2*d3/(d2*h3+d3*h2))+                   &
                        &   fineonedsolver_type%remvxs_matregion(k,ig)*h2
                case(2)
                    coefficient_matrix(i) = -2.0*(d2*d3/(d2*h3+d3*h2))
                end select
                flux1d(k) = source(k)*h2
            end do
            
            ! the second last element of the matrix
            d1 = fineonedsolver_type%diffcoef_matregion(number_meshes-1,ig)  
            d2 = fineonedsolver_type%diffcoef_matregion(number_meshes,ig)
            h1 = fineonedsolver_type%size_matregion(number_meshes-1)    
            h2 = fineonedsolver_type%size_matregion(number_meshes)       
            coefficient_matrix(3*(number_meshes-1)) = -2.0*(d1*d2/(d1*h2+d2*h1))
            
            ! right boundary. the last element of the matrix
            beta = boundary_condition(2)
            coefficient_matrix(3*number_meshes-2) = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d2*(1.0-beta)/               &
                &   (4.0*d2*(1.0+beta)+h2*(1.0-beta))+fineonedsolver_type%remvxs_matregion(number_meshes,ig)*h2
            flux1d(number_meshes) = source(number_meshes)*h2
            
            open(2,file = "coeff.out")
            write(2,*)
            write(2,*) coefficient_matrix
            write(2,*)
            write(2,*) flux1d
            coefficient_matrix2 = coefficient_matrix
            call onedsolver_chasingmethod(coefficient_matrix,number_meshes,3*number_meshes-2,flux1d,l)     
        end if   
    
        allocate(a(number_meshes,number_meshes),source2(number_meshes))
        
        a = 0.0
        source2 = 0.0
        d1 = fineonedsolver_type%diffcoef_matregion(1,ig)
        d2 = fineonedsolver_type%diffcoef_matregion(2,ig)
        h1 = fineonedsolver_type%size_matregion(1)
        h2 = fineonedsolver_type%size_matregion(2)
        beta = boundary_condition(1)
        
        a(1,1) = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d1*(1.0-beta)/                                                   &
            &   (4.0*d1*(1.0+beta)+h1*(1.0-beta))+fineonedsolver_type%remvxs_matregion(1,ig)*h1
        a(1,2) = -2.0*(d2*d1/(d2*h1+d1*h2))
        
        source2(1) = source(1)*h2
        do i = 2,number_meshes-1
            d1 = fineonedsolver_type%diffcoef_matregion(i-1,ig) 
            d2 = fineonedsolver_type%diffcoef_matregion(i,ig)
            d3 = fineonedsolver_type%diffcoef_matregion(i+1,ig)
            h1 = fineonedsolver_type%size_matregion(i-1)  
            h2 = fineonedsolver_type%size_matregion(i)      
            h3 = fineonedsolver_type%size_matregion(i+1)  
            a(i,i-1) = -2.0*(d1*d2/(d1*h2+d2*h1))
            a(i,i) = 2.0*(d1*d2/(d1*h2+d2*h1)+d2*d3/(d2*h3+d3*h2)) + fineonedsolver_type%remvxs_matregion(i,ig)*h2
            a(i,i+1) = -2.0*(d2*d3/(d2*h3+d3*h2))
            source2(i) = source(i)*h2
        end do
        
        d1 = fineonedsolver_type%diffcoef_matregion(number_meshes-1,ig)  
        d2 = fineonedsolver_type%diffcoef_matregion(number_meshes,ig)
        h1 = fineonedsolver_type%size_matregion(number_meshes-1)    
        h2 = fineonedsolver_type%size_matregion(number_meshes)       
        a(number_meshes,number_meshes-1) = -2.0*(d1*d2/(d1*h2+d2*h1))
        
        ! right boundary. the last element of the matrix
        beta = boundary_condition(2)
        a(number_meshes,number_meshes) = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d2*(1.0-beta)/                           &
            &   (4.0*d2*(1.0+beta)+h2*(1.0-beta)) + fineonedsolver_type%remvxs_matregion(number_meshes,ig)*h2
            
        open(2,file = "a.out")
        write(2,*)
        do i = 1,number_meshes
            write(2,"(<number_meshes>(tr2,es12.5))") (a(i,j),j = 1,number_meshes)
        end do
        write(2,*)
10      format("the mesh number of the problem can not be less than 5!")       
    end subroutine onedsolver_spacialsweep

    !$
    !===============================================================================================
    ! one-dimensional diffusion equation solver with the finite-difference-method
    !===============================================================================================
    subroutine onedsolver_diffusionfdm (fineonedsolver_type,flux_finemesh,leakage_finemesh,         &
        &   boundary_condition,kappa,err,mMAX,keffective,number_meshes,number_group,MODe)
 
        type(onedsolver_type),intent(in) :: fineonedsolver_type                     ! fine mesh relevant informations
        real(KREAL), intent(in)      :: leakage_finemesh,kappa,err,boundary_condition
        integer,intent(in)               :: number_group,number_meshes,mMAX,MODe

        real(KREAL), intent(out)    :: flux_finemesh
        real(KREAL), intent(in out) :: keffective
        
        dimension :: boundary_condition(2),kappa(number_group),                 &
            &   flux_finemesh(number_meshes,number_group) ,                     &
            &   leakage_finemesh(number_meshes,number_group),err(2)             ! flux, keff error criteria
        
        ! local 
        real(KREAL), allocatable :: element_matrix(:),q(:),fission_rate(:),scat_source(:), q_fission(:,:),q_fissionold(:,:),flux_old(:,:),km(:,:)
        integer :: m,l,nout,ig,ig1,im,i
        real(KREAL) :: error_MAX,error_MIN,vac,fqold,fqnew,eigv,error1,error2,fMAX,yita,kmean,kerr
        
        !initializations
        allocate(q(number_meshes),q_fission(number_meshes,number_group),                           &
            &   fission_rate(number_meshes),scat_source(number_meshes),                            &
            &   q_fissionold(number_meshes,number_group),flux_old(number_meshes,number_group),     &
            &   km(number_meshes,number_group))
        
        flux_finemesh    = 1.0
        q                = 0.0 
        q_fission        = 1.0
        q_fissionold     = 1.0
        flux_old         = 1.0
        fission_rate     = 0.0
        scat_source      = 0.0
        vac              = 0.0
        fqold            = 0.0
        fqnew            = 0.0
        eigv             = 1.0
        error1           = 0.0
        error2           = 0.0
        nout             = 0
        fMAX             = 0.0
        km               = 0.0
        kmean            = 0.0
        
        ! outer iteration
        outer : do
            nout = nout + 1
            
            ! group sweepings  
            itgroup: do ig = 1, number_group
                q = 0.0
                scat_source = 0.0 
                do im = 1, number_meshes  
                    ! add fission source to total source
                    q(im)=q(im)+q_fission(im,ig)
                    
                    ! add scattering source to total source
                    do ig1 = 1, number_group
                        if (ig1 /= ig) then
                            scat_source(im) = scat_source(im) + flux_finemesh(im,ig1) * fineonedsolver_type%scatxs_matregion(im,ig1,ig)
                        end if
                    end do
                    q(im) = q(im) + scat_source(im)
                    
                    ! add transverse leakage to total source
                    q(im) = q(im) + leakage_finemesh(im,ig)
                end do 
                
                ! spacial sweeping
                call onedsolver_spacialsweep(fineonedsolver_type,flux_finemesh(:,ig),q, number_meshes,ig,boundary_condition)
            end do itgroup
            
            !judge the criteria of keffective
            fqold     = 0.0
            fqnew     = 0.0
            error_MAX = 0.0
            error_MIN = 2.0
            fission_rate= 0.0
            
            do im = 1, number_meshes
                do ig = 1, number_group
                    fission_rate(im)=fission_rate(im)+fineonedsolver_type%nuproxs_matregion(im,ig)*flux_finemesh(im,ig) 
                    q_fission(im,ig)=fission_rate(im)
                end do
            end do
            do im = 1, number_meshes
                do ig = 1, number_group
                    fqold = fqold+q_fissionold(im,ig)*fineonedsolver_type%size_matregion(im)
                    fqnew = fqnew+q_fission(im,ig)*fineonedsolver_type%size_matregion(im)
                    if (q_fissionold(im,ig) /=0 ) then 
                        vac = ABS(q_fission(im,ig)/q_fissionold(im,ig))
                        error_MAX = MAX(error_MAX,vac)
                        error_MIN = MIN(error_MIN,vac)
                    end if
                end do
            end do
            
            ! ------------------------------------------------------------------
            if (MODe == 1) then
                eigv = keffective*fqnew/fqold
                q_fissionold = q_fission
                error1 = ABS(eigv-keffective)/eigv
                keffective = eigv 
        
                error_MAX = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        if (flux_old(im,ig) /= 0) then 
                            vac=ABS((flux_finemesh(im,ig)-flux_old(im,ig))/flux_finemesh(im,ig))
                            error_MAX = MAX(error_MAX,vac)
                        end if
                    end do
                end do
                flux_old = flux_finemesh
                error2 = error_MAX
                
            ! do not change keff, use km to judge convergence
            else if (MODe == 2) then
                km = flux_finemesh/flux_old
                kmean = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        kmean = kmean + km(im,ig)
                    end do
                end do
                kmean = kmean/number_meshes/number_group
                kerr = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        kerr = kerr+(km(im,ig)-kmean)*(km(im,ig)-kmean)
                    end do
                end do
                kerr = SQRT(kerr)/number_meshes/number_group
                error2 = kerr
                flux_old = flux_finemesh
                
            !convergence criteria for flux(fixed source problem)
            else if (MODe == 0) then
                error_MAX = 0.0
                error_MIN = 2.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        if (flux_old(im,ig) /= 0) then 
                            vac = ABS((flux_finemesh(im,ig)-flux_old(im,ig))/flux_finemesh(im,ig))
                            error_MAX = MAX(error_MAX,vac)
                            error_MIN = MIN(error_MIN,vac)
                        end if
                    end do
                end do
                flux_old = flux_finemesh
                error1 = 0.0
                error2 = error_MAX
            end if
            
            fission_rate = 0.0
            do im = 1, number_meshes
                do ig1 = 1, number_group
                    fission_rate(im)=fission_rate(im)+fineonedsolver_type%nuproxs_matregion(im,ig1)*flux_finemesh(im,ig1) 
                end do
                do ig = 1, number_group
                    q_fission(im,ig)=fission_rate(im)*kappa(ig)/keffective
                end do
            end do

            if ((error1 <= err(1)).and.(error2 <= err(2)).and.(nout <= mMAX))then
                exit outer
            else
                cycle outer
            end if
        end do outer
        
        fMAX = 0.0
        do ig = 1, number_group
            do im = 1, number_meshes
                fMAX=fMAX+flux_finemesh(im,ig)
            end do
        end do
        fMAX = fMAX/number_meshes/number_group
        do ig = 1, number_group
            do im = 1, number_meshes
                flux_finemesh(im,ig)=flux_finemesh(im,ig)/fMAX
            end do
        end do
    end subroutine onedsolver_diffusionfdm
    
    !$
    !===============================================================================================
    ! one-dimensional sn equation solver with the finite-difference-method
    !===============================================================================================
    subroutine onedsolver_snfdm(fineonedsolver_type,flux_finemesh,leakage_finemesh,                             &
        &   boundary_condition,kappa,err,mMAX,keffective,number_meshes, number_group,polar_angle,MODe)

        type(onedsolver_type),intent(in)    :: fineonedsolver_type                                      ! fine mesh relevant informations
        real(KREAL), intent(in)         :: leakage_finemesh,kappa,err,boundary_condition
        integer,intent(in)                  :: number_group,number_meshes,mMAX,polar_angle,MODe
        real(KREAL), intent(out)        :: flux_finemesh
        real(KREAL), intent(in out)     :: keffective
        dimension :: boundary_condition(2),kappa(number_group), flux_finemesh(number_meshes,number_group),      &
            &   leakage_finemesh(number_meshes,number_group),err(2)                                 !flux, keff error criteria
        
        ! local 
        real(KREAL), allocatable :: element_matrix(:),q(:),fission_rate(:),angular_flux(:,:,:),             &
            &   scat_source(:),q_fission(:,:),q_fissionold(:,:),flux_old(:,:),                                  &
            &   boundary_flux(:,:),flux2(:),km(:,:),angular_flux_old(:,:,:)
        integer :: n,m,l,nout,ig,ig1,im,is,mMAX2,nin 
        real(KREAL) :: error_MAX,error_MIN,vac,fqold,fqnew,eigv,error1,error2,fMAX,surfaceflux_in,          &
            &   surfaceflux_out,h,toxs,err2,error_MAXin,kmean,kerr
        real(KREAL), allocatable :: xmu(:)                                  ! mu
        real(KREAL), allocatable :: wmu(:)                                  ! weighting factors
    
        ! initializations
        allocate(q(number_meshes),q_fission(number_meshes,number_group),fission_rate(number_meshes),            &
            &   scat_source(number_meshes),q_fissionold(number_meshes,number_group),                            &
            &   flux_old(number_meshes,number_group),xmu(polar_angle),wmu(polar_angle),                         &
            &   angular_flux(polar_angle,number_meshes,number_group),flux2(number_meshes),                      &
            &   boundary_flux(polar_angle,number_group),km(number_meshes,number_group),                         &
            &   angular_flux_old(polar_angle,number_meshes,number_group))
    
        flux_finemesh    = 1.0
        flux_old         = 1.0
        q                = 0.0 
        q_fission        = 1.0
        q_fissionold     = 1.0
        fission_rate     = 0.0
        scat_source      = 0.0
        vac              = 0.0
        fqold            = 0.0
        fqnew            = 0.0
        eigv             = 1.0
        error1           = 0.0
        error2           = 0.0
        nout             = 0
        nin              = 0
        fMAX             = 0.0 
        surfaceflux_in   = 0.0
        surfaceflux_out  = 0.0
        xmu              = 0.0 
        wmu              = 0.0
        angular_flux     = 0.1 
        flux2            = 0.0
        boundary_flux    = 1.0
        h                = 0.0
        toxs             = 0.0
        error_MAXin      = 0.0
        angular_flux_old = 0.0
        
        mMAX2            = 5                                                    ! MAXimum number of inner iterations
        err2             = 1E-6                                                 !inner iteration convergence criteria
        
        open(1,file="1dsolver.out")
        
        ! select gauss quadrature sets according to polar angle number 
        call onedsolver_gaussquadraturesets (polar_angle,xmu,wmu)
    
        ! begin outer iteration
        outer: do
            nout = nout + 1
            ! update fission rate
            fission_rate = 0.0 
            do im = 1, number_meshes
                do ig1 = 1, number_group
                    fission_rate(im)=fission_rate(im)+fineonedsolver_type%nuproxs_matregion(im,ig1)* flux_finemesh(im,ig1) 
                end do
                do ig1 = 1, number_group
                    q_fission(im,ig1)=fission_rate(im)*kappa(ig1)/keffective/2.0
                end do
            end do 
            
            ! group sweepings
            itgroup: do ig = 1, number_group
                nin = 0
                inner: do 
                    nin = nin + 1
                    scat_source = 0.0 
                    q           = 0.0
                    do im = 1, number_meshes  
                        ! add fission source and leakage source to total source
                        q(im) = q(im)+q_fission(im,ig)+leakage_finemesh(im,ig)/2.0
                    end do    
                    do im = 1, number_meshes
                        ! add scattering source to total source
                        do ig1 = 1, number_group
                            scat_source(im)=scat_source(im)+flux_finemesh(im,ig1)*fineonedsolver_type%scatxs_matregion(im,ig1,ig)
                        end do
                        q(im) = q(im)+scat_source(im)/2.0
                        ! add transverse leakage to total source
                        q(im)=q(im)
                    end do 
                    
                    ! angular sweeping
                    ! from left to right
                    do is = polar_angle/2+1, polar_angle
                        surfaceflux_in = boundary_flux(is,ig)
                        ! spacial sweeping
                        do im = 1, number_meshes
                            h = fineonedsolver_type%size_matregion(im)
                            toxs = fineonedsolver_type%toxs_matregion(im,ig)
                            ! angular flux of mesh im, group ig
                            angular_flux(is,im,ig) = (h*q(im)+2.0*xmu(is)*surfaceflux_in)/(h*toxs+2.0*xmu(is))
                            surfaceflux_out = 2*angular_flux(is,im,ig)-surfaceflux_in
                            surfaceflux_in = surfaceflux_out  
                        end do
                        ! boundary conditions
                        boundary_flux(polar_angle+1-is,ig) = boundary_condition(2)*surfaceflux_in
                    end do
    
                    !from right to left
                    do is = 1, polar_angle/2
                        surfaceflux_in = boundary_flux(is,ig)
                        !spacial sweeping
                        do im = number_meshes, 1, -1
                            h = fineonedsolver_type%size_matregion(im)
                            toxs = fineonedsolver_type%toxs_matregion(im,ig)
                            !angular flux of mesh im, group ig
                            angular_flux(is,im,ig) = (h*q(im)-2.0*xmu(is)*surfaceflux_in)/(h*toxs-2.0*xmu(is))
                            surfaceflux_out = 2*angular_flux(is,im,ig)-surfaceflux_in
                            surfaceflux_in = surfaceflux_out   
                        end do
                        ! boundary conditions
                        boundary_flux(polar_angle+1-is,ig) = boundary_condition(1)*surfaceflux_in
                    end do 
    
                    flux2 = 0.0
                    !calculate scalar flux 
                    do im = 1, number_meshes
                        do is = 1, polar_angle
                            flux2(im)=flux2(im)+wmu(is)*angular_flux(is,im,ig)
                        end do
                    end do 
    
                    error_MAXin = 0.0
                    do im=1, number_meshes
                        if (flux2(im) /= 0) then 
                            vac = ABS((flux2(im)-flux_finemesh(im,ig))/flux2(im))
                            error_MAXin = MAX(error_MAXin,vac)
                        end if
                    end do
    
                    if(nin==mMAX2 .or. error_MAXin < err2)then
                        exit inner
                    end if
                    flux_finemesh(:,ig) = flux2
                    flux2 = 0.0
                end do inner
            end do itgroup
    
            ! judge the criteria of keffective
            fqold     = 0.0
            fqnew     = 0.0
            error_MAX = 0.0
            error_MIN = 2.0
            fission_rate= 0.0
            
            do im = 1, number_meshes
                do ig = 1, number_group
                    fission_rate(im) = fission_rate(im)+fineonedsolver_type%nuproxs_matregion(im,ig)*flux_finemesh(im,ig) 
                    q_fission(im,ig) = fission_rate(im)
                end do
            end do
            do im = 1, number_meshes
                do ig = 1, number_group
                    fqold = fqold+q_fissionold(im,ig)*fineonedsolver_type%size_matregion(im)
                    fqnew = fqnew+q_fission(im,ig)*fineonedsolver_type%size_matregion(im)
                    if (q_fissionold(im,ig) /= 0) then 
                        vac = ABS(q_fission(im,ig)/q_fissionold(im,ig))
                        error_MAX = MAX(error_MAX,vac)
                        error_MIN = MIN(error_MIN,vac)
                    end if
                end do
            end do
            
            ! ------------------------------------------------------------------
            if (MODe == 1) then
                eigv = keffective*fqnew/fqold
                q_fissionold = q_fission
                error1 = ABS(eigv-keffective)/eigv
                keffective = eigv 
    
                error_MAX = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        if (flux_old(im,ig) /= 0) then 
                            vac = ABS((flux_finemesh(im,ig)-flux_old(im,ig))/flux_finemesh(im,ig))
                            error_MAX = MAX(error_MAX,vac)
                        end if
                    end do
                end do
                flux_old = flux_finemesh
                error2 = error_MAX
                
            ! do not change keff, use km to judge convergence
            else if (MODe == 2) then
                km = flux_finemesh/flux_old
                kmean = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        kmean = kmean+km(im,ig)
                    end do
                end do
                kmean = kmean/number_meshes/number_group
                kerr = 0.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        kerr=kerr+(km(im,ig)-kmean)*(km(im,ig)-kmean)
                    end do
                end do
                kerr = SQRT(kerr)/number_meshes/number_group
                error2 = kerr
                flux_old = flux_finemesh
            
            else if(MODe==0)then
                ! convergence criteria for flux(fixed source problem)
                error_MAX = 0.0
                error_MIN = 2.0
                do im = 1, number_meshes
                    do ig = 1, number_group
                        if (flux_old(im,ig)/=0) then 
                            vac = ABS((flux_finemesh(im,ig)-flux_old(im,ig))/flux_finemesh(im,ig))
                            error_MAX = MAX(error_MAX,vac)
                            error_MIN = MIN(error_MIN,vac)
                        end if
                    end do
                end do
                flux_old = flux_finemesh
                error1 = 0.0
                error2 = error_MAX
            end if
            
            if ((error1 <= err(1)).and.(error2 <= err(2)).and.(nout <= mMAX)) then
                exit outer
            else
                cycle outer
            end if
        end do outer  

        close(1)            
    end subroutine onedsolver_snfdm          

    !$
    !===============================================================================================
    ! condense volume-integrated fine mesh flux distribution to the specified node
    !===============================================================================================
    subroutine onedsolver_fluxcondensation(fineonedsolver_type,flux_finemesh,number_outputregion,  &
            &   index_outputregion,number_group,number_finemesh,flux_condensed)

        type(onedsolver_type),intent(in)    :: fineonedsolver_type
        real(KREAL), intent(in)         :: flux_finemesh
        integer,intent(in)                  :: number_outputregion,number_group,index_outputregion,number_finemesh
        real(KREAL), intent(in out)     :: flux_condensed  
        dimension   :: flux_finemesh(number_finemesh,number_group),             &
            &   flux_condensed(number_outputregion,number_group),               &
            &   index_outputregion(number_outputregion)
    
        !local 
        integer          :: i,j,im,ig
        real(KREAL)  :: volume 
    
        flux_condensed = 0.0
        do i = 1, number_outputregion
            do im = 1, number_finemesh
                j = fineonedsolver_type%index_matregion(im)
                if (j == index_outputregion(i)) then
                    volume = fineonedsolver_type%size_matregion(im)
                    do ig = 1, number_group
                        flux_condensed(i,ig) = flux_condensed(i,ig)+flux_finemesh(im,ig)*volume
                    end do
                end if
            end do
        end do
    
    end subroutine onedsolver_fluxcondensation   

    !$
    !===============================================================================================
    ! spacial sweep scheme for the fdm
    !===============================================================================================
    subroutine onedsolver_multiplication (fineonedsolver_type,kappa,keffective,xin,xout,number_meshes,number_group,boundary_condition)   

        type(onedsolver_type),intent(in) :: fineonedsolver_type
        real(KREAL), intent(in)      :: boundary_condition,kappa,keffective ! values of sources
        integer,intent(in)               :: number_meshes,number_group          ! number of meshes,bc,number of groups
        real(KREAL)                  :: xin, xout                           !output vector
        dimension :: xin(number_meshes*number_group), xout(number_meshes*number_group), boundary_condition(2),kappa(number_group)
        
        !local 
        real(KREAL) :: bprime(number_meshes*number_group)
        real(KREAL) :: d1,d2,d3,h1,h2,h3,sigr,beta,a1,a2,a3,temp,co_loss
        integer i,j,k,l,im,im2,ig1,ig
        real(KREAL), allocatable :: aa(:,:),bb(:,:)
        ! initializations
        bprime  = 0.0
        d1      = 0.0
        d2      = 0.0
        d3      = 0.0
        a1      = 0.0
        a2      = 0.0
        a3      = 0.0
        h1      = 0.0
        h2      = 0.0
        h3      = 0.0
        sigr    = 0.0
        temp    = 0.0
        co_loss = 0.0

        allocate(aa(number_meshes*number_group,number_meshes*number_group),bb(number_meshes*number_group,number_meshes*number_group))
        
        aa = 0.0
        bb = 0.0
        
        if (number_meshes <= 2) then
            write(*,10)
            stop
        else
            do i = 1, number_group*number_meshes
                ig = (i-1) / number_meshes+1                                    ! group index
                im = MOD(i,number_meshes)                                       ! mesh index
                if (im /= 0) then
                    im = im
                else if (im == 0) then
                    im = number_meshes
                end if
                
                if (i == number_group*number_meshes) then
                    write(*,*)
                end if
                
                if (im == 1) then
                    ! left boundary. ax vector for group ig, mesh im
                    d1 = fineonedsolver_type%diffcoef_matregion(1,ig)
                    d2 = fineonedsolver_type%diffcoef_matregion(2,ig)
                    h1 = fineonedsolver_type%size_matregion(1)
                    h2 = fineonedsolver_type%size_matregion(2)
                    beta = boundary_condition(1)
                    a2 = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d1*(1.0-beta)/(4.0*d1*(1.0+beta)+h1*(1.0-beta))+ fineonedsolver_type%remvxs_matregion(1,ig)*h1
                    a3 = -2.0*(d2*d1/(d2*h1+d1*h2))
                    bprime(i) = a2*xin(i)+a3*xin(i+1)
                    aa(i,i) = a2
                    aa(i,i+1) = a3
                    temp = 0.0
                    
                    do ig1 = 1, number_group
                        if (ig1 == ig) then
                            ! mesh im, self contribution of group ig
                            im2 = (ig1-1)*number_meshes+im    !index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1))
                            bb(i,im2) = co_loss
                            temp = temp+xin(im2)*co_loss*h2
                            ! mesh im , out-of-group contributions 
                        else if (ig1 /= ig) then
                            im2 = (ig1-1)*number_meshes+im                      ! index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1)+ fineonedsolver_type%scatxs_matregion(im,ig1,ig))
                            bb(i,im2) = co_loss
                            temp = temp+xin(im2)*co_loss*h2
                        end if
                    end do
                    bprime(i) = bprime(i)+temp
                    
                else if (im/=1 .and. im/=number_meshes) then
                    d1 = fineonedsolver_type%diffcoef_matregion(im-1,ig) 
                    d2 = fineonedsolver_type%diffcoef_matregion(im,ig)
                    d3 = fineonedsolver_type%diffcoef_matregion(im+1,ig)
                    h1 = fineonedsolver_type%size_matregion(im-1)  
                    h2 = fineonedsolver_type%size_matregion(im)      
                    h3 = fineonedsolver_type%size_matregion(im+1)              
                    a1 = -2.0*(d1*d2/(d1*h2+d2*h1))
                    a2 = 2.0*(d1*d2/(d1*h2+d2*h1)+d2*d3/(d2*h3+d3*h2))+fineonedsolver_type%remvxs_matregion(im,ig)*h2
                    a3 = -2.0*(d2*d3/(d2*h3+d3*h2))
                    aa(i,i-1) = a1
                    aa(i,i) = a2
                    aa(i,i+1) = a3
                    bprime(i) = a1*xin(i-1)+a2*xin(i)+a3*xin(i+1)
                    temp = 0.0
                    
                    do ig1 = 1, number_group
                        if (ig1 == ig) then
                            ! mesh im, self contribution of group ig
                            im2 = (ig1-1)*number_meshes+im    !index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1))
                            bb(i,im2) = co_loss
                            temp = temp+xin(im2)*co_loss*h2
                            ! mesh im, out-of-group contributions 
                        else if (ig1 /= ig) then
                            im2 = (ig1-1)*number_meshes+im                      ! index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1)+ fineonedsolver_type%scatxs_matregion(im,ig1,ig))
                            bb(i,im2) = co_loss
                            temp = temp+xin(im2)*co_loss*h2
                        end if
                    end do
                    bprime(i) = bprime(i)+temp
                    
                else
                    ! right boundary
                    d1 = fineonedsolver_type%diffcoef_matregion(number_meshes-1,ig)  
                    d2 = fineonedsolver_type%diffcoef_matregion(number_meshes,ig)
                    h1 = fineonedsolver_type%size_matregion(number_meshes-1)    
                    h2 = fineonedsolver_type%size_matregion(number_meshes)   
                    beta = boundary_condition(2)    
                    a1 = -2.0*(d1*d2/(d1*h2+d2*h1))
                    a2 = 2.0*d2*d1/(d2*h1+d1*h2)+2.0*d2*(1.0-beta)/(4.0*d2*(1.0+beta)+h2*(1.0-beta))+fineonedsolver_type%remvxs_matregion(number_meshes,ig)*h2
                    aa(i,i-1) = a1
                    aa(i,i) = a2
                    bprime(i) = a1*xin(i-1)+a2*xin(i)
                    temp = 0.0
                    
                    do ig1 = 1, number_group
                        if (ig1 == ig) then
                            ! mesh im, self contribution of group ig
                            im2 = (ig1-1)*number_meshes+im                      ! index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1))
                            temp = temp+xin(im2)*co_loss*h2
                            bb(i,im2) = co_loss
                            ! mesh im , out-of-group contributions 
                        else if (ig1 /= ig) then
                            im2 = (ig1-1)*number_meshes+im                      ! index in xin
                            co_loss = -(kappa(ig)/keffective*fineonedsolver_type%nuproxs_matregion(im,ig1)+ fineonedsolver_type%scatxs_matregion(im,ig1,ig))
                            temp = temp+xin(im2)*co_loss*h2
                            bb(i,im2) = co_loss
                        end if
                    end do
                    bprime(i) = bprime(i)+temp
                end if
            end do
            xout = bprime
        end if   
    
10      format("the mesh number of the problem can not be less than 3!")       
    end subroutine onedsolver_multiplication
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine onedsolver_define (useronedsolver_type, number_matregion, number_group)

        type(onedsolver_type), intent(in out) :: useronedsolver_type                ! the user provided data set_exponent
        integer, intent(in)  :: number_matregion                                ! number of material regions
        integer, intent(in)  :: number_group                                    ! number of material groups

        integer  :: ios, i, j, k  
    
        if (.NOT. useronedsolver_type%defined) then 
            if(number_matregion <= 0)then
            
            end if
            if(number_group <= 0)then

            end if
    
            ! allocate the object
            allocate(useronedsolver_type%index_matregion(number_matregion), stat = ios)                        
            allocate(useronedsolver_type%size_matregion(number_matregion), stat = ios)                         
            allocate(useronedsolver_type%diffcoef_matregion(number_matregion,number_group), stat = ios)        
            allocate(useronedsolver_type%remvxs_matregion(number_matregion,number_group), stat = ios)          
            allocate(useronedsolver_type%toxs_matregion(number_matregion,number_group), stat = ios)            
            allocate(useronedsolver_type%nuproxs_matregion(number_matregion,number_group), stat = ios)             
            allocate(useronedsolver_type%scatxs_matregion(number_matregion,number_group,number_group), stat = ios)
    
            ! initialize the object
            useronedsolver_type%defined = .true.
            useronedsolver_type%number_matregion = number_matregion
            
            do i = 1, number_matregion
                useronedsolver_type%index_matregion(i) = 0 
                useronedsolver_type%size_matregion(i)  = 0.0
                do j = 1, number_group
                    useronedsolver_type%diffcoef_matregion(i,j) = 0.0
                    useronedsolver_type%remvxs_matregion(i,j)   = 0.0
                    useronedsolver_type%toxs_matregion(i,j)     = 0.0
                    useronedsolver_type%nuproxs_matregion(i,j)  = 0.0
                    do k = 1, number_group
                        useronedsolver_type%scatxs_matregion(i,j,k) = 0.0
                    end do
                end do
            end do
        else 

        end if
   
    end subroutine onedsolver_define

    !$
    !===============================================================================================
    ! void object onedsolver_type
    !===============================================================================================
    subroutine onedsolver_void (useronedsolver_type)   

        type(onedsolver_type), intent(in out) :: useronedsolver_type                ! the user provided data set
        integer ios
    
        if(useronedsolver_type%defined) then
            deallocate(useronedsolver_type%index_matregion,            &
                &   useronedsolver_type%size_matregion,                &
                &   useronedsolver_type%diffcoef_matregion,            &
                &   useronedsolver_type%remvxs_matregion,              &
                &   useronedsolver_type%toxs_matregion,                &
                &   useronedsolver_type%nuproxs_matregion,             &
                &   useronedsolver_type%scatxs_matregion, stat =ios)
                
            if (ios /= 0) then 

            end if
        end if
    
        ! redefine values
        useronedsolver_type%defined  = .FALSE.
        useronedsolver_type%number_matregion = 0 

    end subroutine onedsolver_void
    
    !$
    !===============================================================================================
    ! the chasing method solver for the tridiagonal matrix
    !===============================================================================================
    subroutine onedsolver_chasingmethod (b,n,m,d,l)
 
        real(KREAL), intent(in out) :: b
        integer,intent(in)              :: n,m
        integer,intent(out)             :: l
        real(KREAL), intent(in out) :: d
        dimension b(m),d(n)
        
        integer k,j
    
        l = 1
        if (m /= (3*n-2)) then
            l = -1
            write(*,10)
            return
        end if
        do k = 1, n-1
            j = 3*k-2
            if (ABS(b(j))+1.0 == 1.0) then
                l = 0
                write(*,10)
                return
            end if
            b(j+1) = b(j+1)/b(j)
            d(k) = d(k)/b(j)
            b(j+3) = b(j+3)-b(j+2)*b(j+1)
            d(k+1) = d(k+1)-b(j+2)*d(k)
        end do
        if (ABS(b(3*n-2))+1.0 == 1.0) then
            l = 0
            write(*,10)
            return
        end if
        
        d(n) = d(n) / b(3*n-2)
        do k = n-1, 1, -1
            d(k) = d(k)-b(3*k-1)*d(k+1)
        end do
        
        return
10      format(1x,'  err  ')
    end subroutine onedsolver_chasingmethod 
    
    !$
    !===============================================================================================
    ! select angular cosine and corresponding weighting factors according to the polar angle number
    !===============================================================================================
    subroutine onedsolver_gaussquadraturesets(polar_angle,xmu,wmu)

        integer,intent(in)           :: polar_angle
        real(KREAL), intent(out) :: xmu(polar_angle)
        real(KREAL), intent(out) :: wmu(polar_angle)
    
        select case (polar_angle)
        case (2)   
            xmu(1) = -0.577350269189626
            wmu(1) = 1.0
            xmu(2) = -xmu(1) 
            wmu(2) = wmu(1)
        case (4)   
            xmu(1) = -0.861136311594053  
            wmu(1) = 0.347854845137454
            xmu(2) = -0.339981043584856
            wmu(2) = 0.652145154862546
            xmu(3) = -xmu(2)
            wmu(3) = wmu(2)
            xmu(4) = -xmu(1)
            wmu(4) = wmu(1)
        case (6)   
            xmu(1) = -0.932469514203152
            wmu(1) = 0.17132449237917
            xmu(2) = -0.661209386466265
            wmu(2) = 0.360761573048139
            xmu(3) = -0.238619186083197
            wmu(3) = 0.467913934572691
            xmu(4) = -xmu(3)
            wmu(4) = wmu(3)
            xmu(5) = -xmu(2)
            wmu(5) = wmu(2)
            xmu(6) = -xmu(1)
            wmu(6) = wmu(1)
        case (8)   
            xmu(1) = -0.960289856497536
            wmu(1) = 0.101228536290376
            xmu(2) = -0.796666477413627
            wmu(2) = 0.222381034453374
            xmu(3) = -0.525532409916329
            wmu(3) = 0.313706645877887
            xmu(4) = -0.18343464249565
            wmu(4) = 0.362683783378362
            xmu(5) = -xmu(4)
            wmu(5) = wmu(4)
            xmu(6) = -xmu(3)
            wmu(6) = wmu(3)
            xmu(7) = -xmu(2)
            wmu(7) = wmu(2)
            xmu(8) = -xmu(1)
            wmu(8) = wmu(1)
        case (10)   
            xmu(1) = -0.973906528517172
            wmu(1) = 0.066671344308688
            xmu(2) = -0.865063366688985
            wmu(2) = 0.149451349150581
            xmu(3) = -0.679409568299024
            wmu(3) = 0.219086362515982
            xmu(4) = -0.433395394129247
            wmu(4) = 0.269266719309996
            xmu(5) = -0.148874338981631
            wmu(5) = 0.295524224714753
            xmu(6) = -xmu(5)
            wmu(6) = wmu(5)
            xmu(7) = -xmu(4)
            wmu(7) = wmu(4)
            xmu(8) = -xmu(3)
            wmu(8) = wmu(3)
            xmu(9) = -xmu(2)
            wmu(9) = wmu(2)
            xmu(10) = -xmu(1)
            wmu(10) = wmu(1)
        case (12)   
            xmu(1) = -0.981560634246719
            wmu(1) = 0.047175336386512
            xmu(2) = -0.904117256370475
            wmu(2) = 0.106939325995318
            xmu(3) = -0.769902674194305
            wmu(3) = 0.160078328543346
            xmu(4) = -0.587317954286617
            wmu(4) = 0.203167426723066
            xmu(5) = -0.36783149899818
            wmu(5) = 0.233492536538355
            xmu(6) = -0.125233408511469
            wmu(6) = 0.249147045813403
            xmu(7) = -xmu(6)
            wmu(7) = wmu(6)
            xmu(8) = -xmu(5)
            wmu(8) = wmu(5)
            xmu(9) = -xmu(4)
            wmu(9) = wmu(4)
            xmu(10) = -xmu(3)
            wmu(10) = wmu(3)
            xmu(11) = -xmu(2)
            wmu(11) = wmu(2)
            xmu(12) = -xmu(1)
            wmu(12) = wmu(1)
        case (14)   
            xmu(1) = -0.986283808696812
            wmu(1) = 0.035119460331752
            xmu(2) = -0.928434883663573
            wmu(2) = 0.08015808715976
            xmu(3) = -0.827201315069765
            wmu(3) = 0.121518570687903
            xmu(4) = -0.687292904811685
            wmu(4) = 0.157203167158193
            xmu(5) = -0.515248636358154
            wmu(5) = 0.185538397477938
            xmu(6) = -0.31911236892789
            wmu(6) = 0.205198463721296
            xmu(7) = -0.108054948707344
            wmu(7) = 0.215263853463158
            xmu(8) = -xmu(7)
            wmu(8) = wmu(7)
            xmu(9) = -xmu(6)
            wmu(9) = wmu(6)
            xmu(10) = -xmu(5)
            wmu(10) = wmu(5)
            xmu(11) = -xmu(4)
            wmu(11) = wmu(4)
            xmu(12) = -xmu(3)
            wmu(12) = wmu(3)
            xmu(13) = -xmu(2)
            wmu(13) = wmu(2)
            xmu(14) = -xmu(1)
            wmu(14) = wmu(1)
        case (16)   
            xmu(1) = -0.98940093499165
            wmu(1) = 0.027152459411754
            xmu(2) = -0.944575023073233
            wmu(2) = 0.062253523938648
            xmu(3) = -0.865631202387832
            wmu(3) = 0.095158511682493
            xmu(4) = -0.755404408355003
            wmu(4) = 0.124628971255534
            xmu(5) = -0.617876244402644
            wmu(5) = 0.149595988816577
            xmu(6) = -0.458016777657227
            wmu(6) = 0.169156519395003
            xmu(7) = -0.281603550779259
            wmu(7) = 0.182603415044924
            xmu(8) = -0.095012509837637
            wmu(8) = 0.189450610455068
            xmu(9) = -xmu(8)
            wmu(9) = wmu(8)
            xmu(10) = -xmu(7)
            wmu(10) = wmu(7)
            xmu(11) = -xmu(6)
            wmu(11) = wmu(6)
            xmu(12) = -xmu(5)
            wmu(12) = wmu(5)
            xmu(13) = -xmu(4)
            wmu(13) = wmu(4)
            xmu(14) = -xmu(3)
            wmu(14) = wmu(3)
            xmu(15) = -xmu(2)
            wmu(15) = wmu(2)
            xmu(16) = -xmu(1)
            wmu(16) = wmu(1)
        case (18)   
            xmu(1) = -0.99156516842093
            wmu(1) = 0.021616013526483
            xmu(2) = -0.955823949571397
            wmu(2) = 0.04971454889497
            xmu(3) = -0.892602466497556
            wmu(3) = 0.076425730254889
            xmu(4) = -0.803704958972523
            wmu(4) = 0.100942044106287
            xmu(5) = -0.691687043060353
            wmu(5) = 0.122555206711478
            xmu(6) = -0.559770831073947
            wmu(6) = 0.140642914670651
            xmu(7) = -0.411751161462843
            wmu(7) = 0.154684675126265
            xmu(8) = -0.251886225691505
            wmu(8) = 0.164276483745833
            xmu(9) = -0.084775013041735
            wmu(9) = 0.169142382963143
            xmu(10) = -xmu(9)
            wmu(10) = wmu(9)
            xmu(11) = -xmu(8)
            wmu(11) = wmu(8)
            xmu(12) = -xmu(7)
            wmu(12) = wmu(7)
            xmu(13) = -xmu(6)
            wmu(13) = wmu(6)
            xmu(14) = -xmu(5)
            wmu(14) = wmu(5)
            xmu(15) = -xmu(4)
            wmu(15) = wmu(4)
            xmu(16) = -xmu(3)
            wmu(16) = wmu(3)
            xmu(17) = -xmu(2)
            wmu(17) = wmu(2)
            xmu(18) = -xmu(1)
            wmu(18) = wmu(1)
        case (20)   
            xmu(1) = -0.993128599185095
            wmu(1) = 0.017614007139152
            xmu(2) = -0.963971927277914
            wmu(2) = 0.040601429800387
            xmu(3) = -0.912234428251326
            wmu(3) = 0.062672048334109
            xmu(4) = -0.839116971822219
            wmu(4) = 0.083276741576705
            xmu(5) = -0.746331906460151
            wmu(5) = 0.10193011981724
            xmu(6) = -0.636053680726515
            wmu(6) = 0.118194531961518
            xmu(7) = -0.510867001950827
            wmu(7) = 0.131688638449177
            xmu(8) = -0.37370608871542
            wmu(8) = 0.142096109318382
            xmu(9) = -0.227785851141645
            wmu(9) = 0.149172986472604
            xmu(10) = -0.076526521133497
            wmu(10) = 0.152753387130726
            xmu(11) = -xmu(10)
            wmu(11) = wmu(10)
            xmu(12) = -xmu(9)
            wmu(12) = wmu(9)
            xmu(13) = -xmu(8)
            wmu(13) = wmu(8)
            xmu(14) = -xmu(7)
            wmu(14) = wmu(7)
            xmu(15) = -xmu(6)
            wmu(15) = wmu(6)
            xmu(16) = -xmu(5)
            wmu(16) = wmu(5)
            xmu(17) = -xmu(4)
            wmu(17) = wmu(4)
            xmu(18) = -xmu(3)
            wmu(18) = wmu(3)
            xmu(19) = -xmu(2)
            wmu(19) = wmu(2)
            xmu(20) = -xmu(1)
            wmu(20) = wmu(1)
        case (24)   
            xmu(1) = -0.995187219997021
            wmu(1) = 0.012341229799987
            xmu(2) = -0.974728555971309
            wmu(2) = 0.028531388628934
            xmu(3) = -0.938274552002733
            wmu(3) = 0.04427743881742
            xmu(4) = -0.886415527004401
            wmu(4) = 0.059298584915437
            xmu(5) = -0.820001985973903
            wmu(5) = 0.07334648141108
            xmu(6) = -0.740124191578554
            wmu(6) = 0.086190161531953
            xmu(7) = -0.648093651936975
            wmu(7) = 0.097618652104114
            xmu(8) = -0.545421471388839
            wmu(8) = 0.107444270115966
            xmu(9) = -0.433793507626045
            wmu(9) = 0.115505668053726
            xmu(10) = -0.315042679696163
            wmu(10) = 0.121670472927803
            xmu(11) = -0.191118867473616
            wmu(11) = 0.125837456346828
            xmu(12) = -0.064056892862605
            wmu(12) = 0.127938195346752
            xmu(13) = -xmu(12)
            wmu(13) = wmu(12)
            xmu(14) = -xmu(11)
            wmu(14) = wmu(11)
            xmu(15) = -xmu(10)
            wmu(15) = wmu(10)
            xmu(16) = -xmu(9)
            wmu(16) = wmu(9)
            xmu(17) = -xmu(8)
            wmu(17) = wmu(8)
            xmu(18) = -xmu(7)
            wmu(18) = wmu(7)
            xmu(19) = -xmu(6)
            wmu(19) = wmu(6)
            xmu(20) = -xmu(5)
            wmu(20) = wmu(5)
            xmu(21) = -xmu(4)
            wmu(21) = wmu(4)
            xmu(22) = -xmu(3)
            wmu(22) = wmu(3)
            xmu(23) = -xmu(2)
            wmu(23) = wmu(2)
            xmu(24) = -xmu(1)
            wmu(24) = wmu(1)
        case (28)   
            xmu(1) = -0.996442497573954
            wmu(1) = 0.009124282593094
            xmu(2) = -0.981303165370873
            wmu(2) = 0.021132112592771
            xmu(3) = -0.954259280628938
            wmu(3) = 0.032901427782304
            xmu(4) = -0.915633026392132
            wmu(4) = 0.044272934759004
            xmu(5) = -0.865892522574395
            wmu(5) = 0.055107345675717
            xmu(6) = -0.805641370917179
            wmu(6) = 0.065272923966999
            xmu(7) = -0.735610878013632
            wmu(7) = 0.074646214234569
            xmu(8) = -0.656651094038865
            wmu(8) = 0.083113417228901
            xmu(9) = -0.569720471811402
            wmu(9) = 0.090571744393033
            xmu(10) = -0.475874224955118
            wmu(10) = 0.09693065799793
            xmu(11) = -0.376251516089079
            wmu(11) = 0.102112967578061
            xmu(12) = -0.272061627635178
            wmu(12) = 0.106055765922846
            xmu(13) = -0.16456928213338
            wmu(13) = 0.108711192258294
            xmu(14) = -0.055079289884034
            wmu(14) = 0.110047013016475
            xmu(15) = -xmu(14)
            wmu(15) = wmu(14)
            xmu(16) = -xmu(13)
            wmu(16) = wmu(13)
            xmu(17) = -xmu(12)
            wmu(17) = wmu(12)
            xmu(18) = -xmu(11)
            wmu(18) = wmu(11)
            xmu(19) = -xmu(10)
            wmu(19) = wmu(10)
            xmu(20) = -xmu(9)
            wmu(20) = wmu(9)
            xmu(21) = -xmu(8)
            wmu(21) = wmu(8)
            xmu(22) = -xmu(7)
            wmu(22) = wmu(7)
            xmu(23) = -xmu(6)
            wmu(23) = wmu(6)
            xmu(24) = -xmu(5)
            wmu(24) = wmu(5)
            xmu(25) = -xmu(4)
            wmu(25) = wmu(4)
            xmu(26) = -xmu(3)
            wmu(26) = wmu(3)
            xmu(27) = -xmu(2)
            wmu(27) = wmu(2)
            xmu(28) = -xmu(1)
            wmu(28) = wmu(1)
        case (32)
            xmu(1) = -0.997263861849482
            wmu(1) = 0.00701861000947
            xmu(2) = -0.985611511545268
            wmu(2) = 0.016274394730906
            xmu(3) = -0.964762255587506
            wmu(3) = 0.025392065309262
            xmu(4) = -0.93490607593774
            wmu(4) = 0.034273862913021
            xmu(5) = -0.896321155766052
            wmu(5) = 0.042835898022227
            xmu(6) = -0.84936761373257
            wmu(6) = 0.050998059262376
            xmu(7) = -0.794483795967942
            wmu(7) = 0.058684093478536
            xmu(8) = -0.73218211874029
            wmu(8) = 0.065822222776362
            xmu(9) = -0.663044266930215
            wmu(9) = 0.072345794108849
            xmu(10) = -0.587715757240762
            wmu(10) = 0.07819389578707
            xmu(11) = -0.506899908932229
            wmu(11) = 0.083311924226947
            xmu(12) = -0.421351276130635
            wmu(12) = 0.087652093004404
            xmu(13) = -0.331868602282128
            wmu(13) = 0.091173878695764
            xmu(14) = -0.239287362252137
            wmu(14) = 0.093844399080805
            xmu(15) = -0.144471961582796
            wmu(15) = 0.095638720079275
            xmu(16) = -0.048307665687738
            wmu(16) = 0.096540088514728
            xmu(17) = -xmu(16)
            wmu(17) = wmu(16)
            xmu(18) = -xmu(15)
            wmu(18) = wmu(15)
            xmu(19) = -xmu(14)
            wmu(19) = wmu(14)
            xmu(20) = -xmu(13)
            wmu(20) = wmu(13)
            xmu(21) = -xmu(12)
            wmu(21) = wmu(12)
            xmu(22) = -xmu(11)
            wmu(22) = wmu(11)
            xmu(23) = -xmu(10)
            wmu(23) = wmu(10)
            xmu(24) = -xmu(9)
            wmu(24) = wmu(9)
            xmu(25) = -xmu(8)
            wmu(25) = wmu(8)
            xmu(26) = -xmu(7)
            wmu(26) = wmu(7)
            xmu(27) = -xmu(6)
            wmu(27) = wmu(6)
            xmu(28) = -xmu(5)
            wmu(28) = wmu(5)
            xmu(29) = -xmu(4)
            wmu(29) = wmu(4)
            xmu(30) = -xmu(3)
            wmu(30) = wmu(3)
            xmu(31) = -xmu(2)
            wmu(31) = wmu(2)
            xmu(32) = -xmu(1)
            wmu(32) = wmu(1)
        case (64)
            xmu(1) = -0.999305041735772
            wmu(1) = 0.001783280721696
            xmu(2) = -0.996340116771955
            wmu(2) = 0.004147033260562
            xmu(3) = -0.991013371476744
            wmu(3) = 0.006504457968978
            xmu(4) = -0.983336253884626
            wmu(4) = 0.008846759826364
            xmu(5) = -0.973326827789911
            wmu(5) = 0.011168139460131
            xmu(6) = -0.961008799652054
            wmu(6) = 0.013463047896719
            xmu(7) = -0.946411374858403
            wmu(7) = 0.015726030476025
            xmu(8) = -0.92956917213194
            wmu(8) = 0.017951715775697
            xmu(9) = -0.910522137078503
            wmu(9) = 0.02013482315353
            xmu(10) = -0.889315445995114
            wmu(10) = 0.022270173808383
            xmu(11) = -0.865999398154093
            wmu(11) = 0.024352702568711
            xmu(12) = -0.84062929625258
            wmu(12) = 0.026377469715055
            xmu(13) = -0.813265315122798
            wmu(13) = 0.028339672614259
            xmu(14) = -0.783972358943341
            wmu(14) = 0.030234657072402
            xmu(15) = -0.752819907260532
            wmu(15) = 0.032057928354852
            xmu(16) = -0.719881850171611
            wmu(16) = 0.033805161837142
            xmu(17) = -0.685236313054233
            wmu(17) = 0.035472213256882
            xmu(18) = -0.648965471254657
            wmu(18) = 0.03705512854024
            xmu(19) = -0.611155355172393
            wmu(19) = 0.038550153178616
            xmu(20) = -0.571895646202634
            wmu(20) = 0.03995374113272
            xmu(21) = -0.531279464019894
            wmu(21) = 0.041262563242624
            xmu(22) = -0.489403145707053
            wmu(22) = 0.042473515123654
            xmu(23) = -0.446366017253464
            wmu(23) = 0.043583724529323
            xmu(24) = -0.402270157963992
            wmu(24) = 0.044590558163757
            xmu(25) = -0.357220158337668
            wmu(25) = 0.045491627927418
            xmu(26) = -0.311322871990211
            wmu(26) = 0.046284796581314
            xmu(27) = -0.264687162208767
            wmu(27) = 0.04696818281621
            xmu(28) = -0.217423643740007
            wmu(28) = 0.04754016571483
            xmu(29) = -0.169644420423993
            wmu(29) = 0.047999388596458
            xmu(30) = -0.121462819296121
            wmu(30) = 0.048344762234803
            xmu(31) = -0.072993121787799
            wmu(31) = 0.048575467441503
            xmu(32) = -0.024350292663424
            wmu(32) = 0.04869095700914
            xmu(33) = -xmu(32)
            wmu(33) = wmu(32)
            xmu(34) = -xmu(31)
            wmu(34) = wmu(31)
            xmu(35) = -xmu(30)
            wmu(35) = wmu(30)
            xmu(36) = -xmu(29)
            wmu(36) = wmu(29)
            xmu(37) = -xmu(28)
            wmu(37) = wmu(28)
            xmu(38) = -xmu(27)
            wmu(38) = wmu(27)
            xmu(39) = -xmu(26)
            wmu(39) = wmu(26)
            xmu(40) = -xmu(25)
            wmu(40) = wmu(25)
            xmu(41) = -xmu(24)
            wmu(41) = wmu(24)
            xmu(42) = -xmu(23)
            wmu(42) = wmu(23)
            xmu(43) = -xmu(22)
            wmu(43) = wmu(22)
            xmu(44) = -xmu(21)
            wmu(44) = wmu(21)
            xmu(45) = -xmu(20)
            wmu(45) = wmu(20)
            xmu(46) = -xmu(19)
            wmu(46) = wmu(19)
            xmu(47) = -xmu(18)
            wmu(47) = wmu(18)
            xmu(48) = -xmu(17)
            wmu(48) = wmu(17)
            xmu(49) = -xmu(16)
            wmu(49) = wmu(16)
            xmu(50) = -xmu(15)
            wmu(50) = wmu(15)
            xmu(51) = -xmu(14)
            wmu(51) = wmu(14)
            xmu(52) = -xmu(13)
            wmu(52) = wmu(13)
            xmu(53) = -xmu(12)
            wmu(53) = wmu(12)
            xmu(54) = -xmu(11)
            wmu(54) = wmu(11)
            xmu(55) = -xmu(10)
            wmu(55) = wmu(10)
            xmu(56) = -xmu(9)
            wmu(56) = wmu(9)
            xmu(57) = -xmu(8)
            wmu(57) = wmu(8)
            xmu(58) = -xmu(7)
            wmu(58) = wmu(7)
            xmu(59) = -xmu(6)
            wmu(59) = wmu(6)
            xmu(60) = -xmu(5)
            wmu(60) = wmu(5)
            xmu(61) = -xmu(4)
            wmu(61) = wmu(4)
            xmu(62) = -xmu(3)
            wmu(62) = wmu(3)
            xmu(63) = -xmu(2)
            wmu(63) = wmu(2)
            xmu(64) = -xmu(1)
            wmu(64) = wmu(1)
        end select
        
    end subroutine onedsolver_gaussquadraturesets

end module solver_1D
