!$
!===================================================================================================
!
!   module for iteraion calculation kernel
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Generate_source_moments
!                               Add_in_group_source
!                               Get_average_flux
!                               Add_in_group_source_moments
!                               Re_generate_flux_moments
!                               Update_source_moment
!                               Get_forward_error
!
!   Public type lists:          No
!
!===================================================================================================
module iteration_process
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use omp_lib
    
    use truncation
    
    use global
    use coefficient_iteration
    
    use iteration_header,               only : IterationCriterion
    use vector_operation
    
    implicit none
    private
    public  :: Generate_source_moments, Add_in_group_source, Get_average_flux   
    public  :: Re_generate_flux_moments, Add_in_group_source_moments
    public  :: Update_source_moment, Get_forward_error

contains
    !$
    !===============================================================================================
    ! generate in-group source moments, out-group source and out-group source moments
    !===============================================================================================
    subroutine Generate_source_moments (ig, is_adjoint)

        ! intent parameters
        integer, intent(in)  :: ig
        logical, intent(in)  :: is_adjoint
        
        integer  :: ia, ir, il, id, ml, iig, is, ic1, ic2, iz, ie
        real(KREAL)  :: q_out_group                                         ! hold for out-group moment
        real(KREAL)  :: summary                                             ! hold for anosotropic scatter source 
        
        ! out-of-group source and source moments, due to fission, up-scatter source and down-scatter source
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = ml + ir
                do id = 0, 8
                    ! out-of-group source due to fission and external source
                    if (.NOT. is_adjoint)  then
                        if (id == 0) then
                            iter_q%info%out_group_moment(id,il,1) = (iter_q%fission%moment(id,il)*xsec_iter%matrixs(iz,ia)%chi_steady(ig)/iter_count%eigenvalue)    &
                            &   + Q_ext%iter_scalar(ir,ia)%intensity(ig)*iter_count%coeff_FSP
                        else
                            iter_q%info%out_group_moment(id,il,1) = iter_q%fission%moment(id,il)*xsec_iter%matrixs(iz,ia)%chi_steady(ig)/iter_count%eigenvalue
                        end if
                    ! @adjoint
                    else
                        if (id == 0) then
                            iter_q%info%out_group_moment(id,il,1) = (iter_adjoint%q_moments(id,il)*xsec_iter%matrixs(iz,ia)%sigma_f_nu(ig)/iter_count%eigenvalue)   &
                            &   + Q_ext%iter_scalar(ir,ia)%intensity(ig)*iter_count%coeff_FSP
                        else
                            iter_q%info%out_group_moment(id,il,1) = iter_adjoint%q_moments(id,il)*xsec_iter%matrixs(iz,ia)%sigma_f_nu(ig)/iter_count%eigenvalue
                        end if
                    end if 
                    
                    ! out-of-group source due to scatter
                    if (ns%state%ng /= 1)  then
                        if (ig == ns%state%ng)  then
                            do iig = 1, ns%state%ng-1
                                iter_q%info%out_group_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,iig)*xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,1)
                            end do
                        else if (ig == 1)  then
                            do iig = 2, ns%state%ng
                                iter_q%info%out_group_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,iig)*xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,1)
                            end do
                        else
                            do iig = ig+1, ns%state%ng
                                iter_q%info%out_group_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,iig)*xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,1)
                            end do
                            do iig = 1, ig-1
                                iter_q%info%out_group_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,iig)*xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,1)
                            end do
                        end if
                    end if
                end do
            end do
        end do
        
        ! out-of-group source and source moments due to anisotropic scatter, refer p59 of the dissertation
        if (ns%state%scat_order /= 0)  then
            !!!   !$omp parallel do default(shared) private(is, ia, ir, id, ic1, ic2,     &
            !!!   !$omp &  ml, iz, il, q_out_group, summary)                              &
            !!!   !$omp &  schedule(dynamic) ordered
            do is = ns%deduce%direction, 1, -1
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do id = 0, 8
                            q_out_group = iter_q%info%out_group_moment(id, il, 1)
                            
                            if (ns%state%ng /= 1)  then
                                if (ig == ns%state%ng)  then
                                    do iig = 1, ns%state%ng-1
                                        do ic1 = 1, ns%state%scat_order
                                            summary = 0.0D0
                                            do ic2 = 1, ic1
                                                summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_cos(ic2,ic1,id,il)               &
                                                    &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_sin(ic2,ic1,id,il)
                                            end do
                                            q_out_group = q_out_group + (flux_scat%ngs(iig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is) + 2*summary)   &
                                                &  * xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,ic1+1)
                                        end do
                                    end do
                                else if (ig == 1)  then
                                    do iig = 2, ns%state%ng
                                        do ic1 = 1, ns%state%scat_order
                                            summary = 0.0D0
                                            do ic2 = 1, ic1
                                                summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_cos(ic2,ic1,id,il)               &
                                                    &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_sin(ic2,ic1,id,il)
                                            end do
                                            q_out_group = q_out_group + (flux_scat%ngs(iig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is)+2*summary)     &
                                                &  * xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,ic1+1)
                                        end do
                                    end do
                                else
                                    do iig = ig + 1, ns%state%ng
                                        do ic1 = 1, ns%state%scat_order
                                            summary = 0.0D0
                                            do ic2 = 1, ic1
                                                summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_cos(ic2,ic1,id,il)               &
                                                    &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_sin(ic2,ic1,id,il)
                                            end do
                                            q_out_group = q_out_group + (flux_scat%ngs(iig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is)+2*summary)     &
                                                &  * xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,ic1+1)
                                        end do
                                    end do
                                    do iig = 1, ig-1
                                        do ic1 = 1, ns%state%scat_order
                                            summary = 0.0D0
                                            do ic2 = 1, ic1
                                                summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_cos(ic2,ic1,id,il)               &
                                                    &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(iig)%aniso_sin(ic2,ic1,id,il)
                                            end do
                                            q_out_group = q_out_group + (flux_scat%ngs(iig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is)+2*summary)     &
                                                &  * xsec_iter%matrixs(iz,ia)%sigma_s(iig,ig,ic1+1)
                                        end do
                                    end do
                                end if
                            end if
                            
                            iter_q%info%out_group_moment(id,il,is) = q_out_group
                        end do
                    end do
                end do
            end do
            !!!   !$omp end parallel do
            
        end if
        
        ! high order moment for in-group source, add to out-of-group source
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = ml + ir
                do id = 1, 8
                    iter_q%info%total_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
                end do
            end do
        end do
        
        ! high order moment for in-group-source when anisotropic scatter, add to out-of-group source      
        if (ns%state%scat_order /= 0)  then
            !!!   !$omp parallel do default(shared) private(is, ia, ir, id,   &
            !!!   !$omp &  ml, iz, il)                                        
            do is = 2, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do id = 1, 8
                            iter_q%info%total_moment(id,il,is) = iter_q%info%out_group_moment(id,il,is) + iter_flux%info%moment(id,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
                        end do
                    end do
                end do
            end do
            !!!   !$omp end parallel do
            
            !!!   !$omp parallel do default(shared) private(is, ia, ir, id, ic1, ic2,     &
            !!!   !$omp &  ml, iz, il, summary)                                           
            do is = 1, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do id = 1, 8
                            do ic1 = 1, ns%state%scat_order
                                summary = 0.0
                                do ic2 = 1, ic1
                                    summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_cos(ic2,ic1,id,il)                                                            &
                                        &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_sin(ic2,ic1,id,il)
                                end do
                                iter_q%info%total_moment(id,il,is) = iter_q%info%total_moment(id,il,is) + (flux_scat%ngs(ig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is)+2*summary)    &
                                    &  * xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,ic1+1)
                            end do
                        end do
                    end do
                end do
            end do
            !!!   !$omp end parallel do
            
        end if
    
    end subroutine Generate_source_moments

    !$
    !===============================================================================================
    ! add in-group source to generate total source
    !===============================================================================================
    subroutine Add_in_group_source (ig)
        
        integer, intent(in)  :: ig
        
        integer  :: ia, ir, il, ml, is, ic1, ic2, iz
        real(KREAL)  :: summary                                             ! hold for anosotropic scatter source
        
        ! in-group-source when isotropic, so only the first direction
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = ml + ir
                iter_q%info%total_moment(0,il,1) = iter_q%info%out_group_moment(0,il,1) + iter_flux%info%moment(0,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
            end do
        end do
        
        ! in-group-source when anisotropic scatters
        if (ns%state%scat_order /= 0)  then
            do is = 2, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        iter_q%info%total_moment(0,il,is) = iter_q%info%out_group_moment(0,il,is) + iter_flux%info%moment(0,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
                    end do
                end do
            end do
            
            do is = 1, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do ic1 = 1, ns%state%scat_order
                            summary = 0.0
                            do ic2 = 1, ic1
                                summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_cos(ic2,ic1,0,il)                                                         &
                                    &  +coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_sin(ic2,ic1,0,il)
                            end do
                            iter_q%info%total_moment(0,il,is) = iter_q%info%total_moment(0,il,is) + (flux_scat%ngs(ig)%aniso_zero(ic1,0,il)*coeff_source%pl_zero(ic1,is)+2*summary)   &
                                &  * xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,ic1+1)
                        end do
                    end do
                end do
            end do
            
            ! re-set zero order to zero
            call flux_scat%value_zero (ig)
            
        end if
        
    end subroutine Add_in_group_source
    
    !$
    !===============================================================================================
    ! generate average flux per nodal
    !===============================================================================================
    subroutine Get_average_flux (ig, is_adjoint)
        
        ! intent parameters
        integer, intent(in)  :: ig
        logical, intent(in)  :: is_adjoint
        
        ! loop index
        integer  :: k, ia, ir, ic1, ic2, ml, il, iy, is, is0, iz
        
        ! ----------------------------------------------------------------------
        ! local variables
        real(KREAL)  :: sp11, sp12, sp13, sp14, sp15                            ! store temporary coefficient for expansion term
        real(KREAL)  :: sp21, sp22, sp23, sp24, sp25                            
        real(KREAL)  :: sp31, sp32, sp33, sp34, sp35                            
        real(KREAL)  :: sp41, sp45                                              
        real(KREAL)  :: ps01, ps02                                              
        real(KREAL)  :: ps11, ps12                                              
        real(KREAL)  :: ps21, ps22                                              
        real(KREAL)  :: psi1, psi2, psi3                                        
        real(KREAL)  :: ra1dm, ra2dm                                            
        real(KREAL)  :: deltaz                                                  
                                                                                
        real(KREAL)  :: rr1, rr2, ss1, ss2                                      
                                                                                    
        real(KREAL)  :: u1, u2, u3, uz                                          ! projection
        real(KREAL)  :: hlayer                                                  ! height of specific layer
        real(KREAL)  :: sigt                                                    ! total xsec_iter of specific nodal
                                                                                    
        real(KREAL)  :: surf_1, surf_2, surf_3                                  ! three surface flux for nodal
        real(KREAL)  :: surf_z1, surf_z2                                        ! two axial surface for nodal, up and low
        real(KREAL)  :: denominator                                             ! temporary variable to get flux
        real(KREAL)  :: avgflux                                                 ! nodal average flux
        real(KREAL)  :: weighted_flux                                           ! use for anisotropic flux
        
        integer  :: SN90
        integer  :: is_y, is_x, is_z                                            ! symmetry direction for specific direction
        integer  :: ia_start, ia_end, ia_step                                   ! axial sweep
        integer  :: iss, isa, ks                                                ! 
        integer  :: j1, j2, j3, impinge                                         ! for nodal eage
        integer  :: i_count
        integer  :: ix0, iy0 
        
        real(KREAL)  :: SOR_factor

        ! ----------------------------------------------------------------------    
        SN90 = ns%deduce%direction / 8
        
        if (iter_count%in == 1)  then
            SOR_factor = 1.0D0 
        else 
            SOR_factor = accele_SOR%scaling_factor                              ! set to 1 to quit SOR 
        end if
        accele_SOR%nodal(:, :, :) = iter_flux%dist%nodal(:, :, :, ig)
            
        ! for the implement of openMP, collpase the direction sweep into one loop;
        !   sweep for octant + sweep for direction => sweep for direction
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(k, ia_start, ia_end, ia_step,                        &
        !$omp &  is0, is, iss, is_y, is_x, is_z, uz, ia, ir, ml, j2, j3, iy, i_count,                       &
        !$omp &  hlayer, impinge, iz, sigt, il, isa, rr1, ss1, rr2, ss2, deltaz, ks, u1, u2, u3, surf_1, surf_2, surf_3, psi1, psi2, psi3,        &
        !$omp &  sp12, sp14, sp21, sp22, sp23, sp24, sp32, sp34, denominator, avgflux, surf_z1, surf_z2, ix0, iy0,          &
        !$omp &  j1, ps01, ps02, ps11, ps12, ps21, ps22, ra1dm, ra2dm,                                              &
        !$omp &  sp11, sp13, sp15, sp25, sp31, sp33, sp35, sp41, sp45, weighted_flux ,ic1, ic2)                      
            
        ! sweep for direction per octant
        direction: do is0 = ns%deduce%direction, 1, -1
            is = quad%is_order(is0)
            k = ((is-1) - MOD(is-1, SN90))/SN90 + 1
            
            ia_start = (ns%state%layer)**(k/5)             
            ia_end = ns%state%layer / ia_start               
            ia_step  = 1 - 2*(k/5) 
            
            ! is_y  is direction symmetry with y coordinate
            ! is_x  is direction symmetry with x coordinate
            ! is_z  is direction symmetry with z coordinate
            is_x = quad%is_x(is)
            is_y = quad%is_y(is)
            is_z = quad%is_z(is)
            
            uz = quad%directions(is)%xmu(3)
            
            ! boundary treatment
            ! ------------------------------------------------------------------
            ! important NOTE: 
            !   when sweep for boundary treatment, "ml+ir"  means the nodal order
            !   but when sweep for calculation, "ml+impinge" is the true nodal order
            ! ------------------------------------------------------------------
            ! up to down in axial(octant 5-8) or. down to up(octant 1-4)
            boundary: do ia = ia_start, ia_end, ia_step        
                ml = (ia-1) * ns%state%nodal
                i_count = 0
                do ir = 1, ns%state%nodal
                    do iy = 1, 3
                        if (ABS(bound%nodal(iy,ir)-bound%INNER) < EPS_ZERO)  cycle
                        
                        ! boundary condition is known
                        i_count = i_count + 1
                        if (quad%directions(is)%projection(iy,ir) < 0.0)  then
                            j2 = iy + (-2)**(iy/3)
                            j3 = j2 + (-2)**(j2/3)
                            
                            ! parallel with coordinate, simply income equals outcome or zero
                            if (.NOT. ns%flag%is_bevel_edge)  then
                                select case(bound%segmentID(iy,ir))
                                case(1)
                                    iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_y,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                case(2)
                                    iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_x,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                case(3)
                                    iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_y,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                case default
                                    iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_x,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                    iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                end select
                            ! not parallel with coordinate
                            else
                                ! 'vaccum' boundary
                                iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_y,ig)*bound%nodal(iy,ir)
                                iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                
                                ! not symmetry
                                if (.NOT. ns%flag%is_60degree)  then
                                    if (bound%segmentID(iy,ir) == (ns%state%segment-1))  then
                                        iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_y,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_y,ig)*bound%nodal(iy,ir)
                                    else if (bound%segmentID(iy,ir) == ns%state%segment)  then
                                        iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_x,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                    end if
                                ! symmetry
                                else
                                    if (bound%segmentID(iy,ir) == (ns%state%segment-1))  then
                                        iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,quad%is_symmetry(is),ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,quad%is_symmetry(is),ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,quad%is_symmetry(is),ig)*bound%nodal(iy,ir)
                                    else if (bound%segmentID(iy,ir) == ns%state%segment)  then
                                        iter_flux%dist%surface(iy,ir,ia,is) = iter_flux%dist%rad_surf(i_count,ia,is_x,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j2,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j2,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                        iter_flux%dist%point(mesh%point(j3,ir),ia,is,ig) = iter_flux%dist%point(mesh%point(j3,ir),ia,is_x,ig)*bound%nodal(iy,ir)
                                    end if
                                end if
                            end if
                        end if
                    end do
                    
                    ! starting layer: the lower surface when K>4 , the upper surface when K<=4
                    if (ia == ia_start)  then
                        iter_flux%dist%surface(k/5+4,ir,ia,is) = iter_flux%dist%axi_surf(ir,1+k/5,is_z,ig) * bound%axial(1+k/5)
                    end if 
                end do
            end do boundary
            
            ! -----------------------------------------------------------------------------------------------------
            ! sweep for nodal
            flux: do ia = ia_start, ia_end, ia_step
                hlayer = geom%height(ia)
                ml = (ia-1) * ns%state%nodal
                do ir = 1, ns%state%nodal
                    impinge = sweep%order(ir,is)
                    iz = mesh%zone(impinge)
                    sigt = xsec_iter%matrixs(iz,ia)%sigma_t(ig)
                    il = ml + impinge
                    if (ns%state%scat_order == 0) then
                        isa = 1
                    else
                        isa = is
                    end if
                    
                    ! equation 55 at page 51, second term of the RHS
                    rr1 = coeff_surface%axi%second(1,il,is); rr2 = coeff_surface%axi%second(2,il,is);
                    ss1 = iter_q%info%total_moment(7,il,isa); ss2 = iter_q%info%total_moment(8,il,isa); 
!                    deltaz = (coeff_surface%axi%second(1,il,is)*iter_q%info%total_moment(7,il,isa) + coeff_surface%axi%second(2,il,is)*iter_q%info%total_moment(8,il,isa))
                    deltaz = rr1*ss1 + rr2*ss2 
                    
                    ! solve the right higher order moment in radial plant
                    ! sweep for surface, we known the information of income surface
                    ks = 0
                    do iy = 1, 3
                        ! count for unknown surfaces
                        if (quad%directions(is)%projection(iy,impinge) > 0.0)  then
                            ks = ks + 1
                            if (ks == 1)  j1 = iy
                            if (ks == 2)  j2 = iy
                        end if
                    end do
                    
                    ! one surface unknown, refer p51 of the dissertation, (J2,J3 income surfaces)
                    if (ks == 1)  then
                        j2 = j1 + (-2)**(j1/3)
                        j3 = j2 + (-2)**(j2/3)      
                        ! solve average angular flux, surface angular flux and transfer value
                        u1 = quad%directions(is)%projection(j1, impinge)
                        u2 = quad%directions(is)%projection(j2, impinge)
                        u3 = quad%directions(is)%projection(j3, impinge)
                        ! average angular flux of : two radial and one axial income surfaces
                        surf_2 = iter_flux%dist%surface(j2, impinge, ia, is)
                        surf_3 = iter_flux%dist%surface(j3, impinge, ia, is)
                        surf_z2 = iter_flux%dist%surface(4+k/5, impinge, ia, is)
                        
                        if (sigt < EPS_ZERO) then            ! vacuum region
                            !!! flat approximation
                            avgflux = (-2*u2*surf_2 - 2*u3*surf_3 + ABS(uz)/hlayer*surf_z2)/(2*u1 + ABS(uz)/hlayer)
                            surf_z1 = avgflux
                            surf_1 = avgflux
                        else
                            ! angular flux of three points of the triangle
                            psi1 = iter_flux%dist%point(mesh%point(j1,impinge), ia, is, ig)
                            psi2 = iter_flux%dist%point(mesh%point(j2,impinge), ia, is, ig)
                            psi3 = iter_flux%dist%point(mesh%point(j3,impinge), ia, is, ig)
                            
                            sp12 = coeff_surface%rad%first(j1,il,is)
                            sp14 = - coeff_surface%rad%second(1,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)          &
                                &  - coeff_surface%rad%second(2,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)            &
                                &  - ((4*surf_2-psi1)*u2+(4*surf_3-psi1)*u3)/3.0*coeff_surface%rad%third(0,j1,il,is)        &
                                &  - ((psi3-surf_2)*u2+(psi2-surf_3)*u3)*2*coeff_surface%rad%third(1,j1,il,is)              &
                                &  - ((psi1+psi3-2*surf_2)*u2+(psi1+psi2-2*surf_3)*u3)*3*coeff_surface%rad%third(2,j1,il,is)      
                            sp21 = 2 * u1
                            sp22 = sigt
                            sp23 = ABS(uz) / hlayer
                            sp24 = -2*(u2*surf_2+u3*surf_3) + sp23*surf_z2 + iter_q%info%total_moment(0,il,isa)
                            sp32 = 1 / coeff_surface%axi%first(il,is)
                            sp34 = ((coeff_surface%axi%third(il,is)-coeff_surface%axi%first(il,is))*surf_z2+deltaz) * sp32
                            
                            ! denonimator of average nodal flux
                            denominator = sp22 + sp32*sp23 + sp21*sp12
                            if (ABS(denominator) >= EPS_ZERO)  then
                                avgflux = ((sp24+sp23*sp34)+sp14*sp21) / denominator                        ! update average nodal flux
                                surf_1 = (sp12*(sp24+sp23*sp34)-sp14*(sp22+sp32*sp23)) / denominator        ! update angular flux for outgoing surface in radial
                                surf_z1 = (sp32*(sp24+sp14*sp21)-(sp22+sp12*sp21)*sp34) / denominator       ! update angular flux for outgoing surface in axial
                            else
                                avgflux = 0.0D0
                                surf_1 = 0.0D0
                                surf_z1 = 0.0D0
                            end if
                        end if 
                        
                        ! negative flux fix-up
                        if (.NOT. is_adjoint)  then
                            call Truncate_negative (avgflux)
                            call Truncate_negative (surf_1)
                            call Truncate_negative (surf_z1)
                        else 
                            call Truncate_negative (avgflux)
                            call Truncate_negative (surf_1)
                            call Truncate_negative (surf_z1)
                        end if
                        
                        iter_flux%dist%nodal(impinge,ia,is,ig) = avgflux * SOR_factor + (1.0-SOR_factor) * accele_SOR%nodal(impinge,ia,is)
                        iter_flux%dist%surface(j1,impinge,ia,is) = surf_1 
                        iter_flux%dist%surface(5-k/5,impinge,ia,is) = surf_z1 
                        
                        if (mesh%nearby_nodal(j1,impinge) /= impinge)  then
                            ix0 = mesh%localID(j1,impinge)
                            iy0 = mesh%nearby_nodal(j1,impinge)
                        
                            iter_flux%dist%surface(ix0,iy0,ia,is) = surf_1 
                        end if 
                        
                        if (ia/=ia_end)  then
                            iter_flux%dist%surface(4+k/5,impinge,ia+ia_step,is) = surf_z1 
                        end if 
                        
                    ! two surface unknown, refer p52 of the dissertation, (J3 income surfaces£¬J1,J2 outcome)
                    else
                        j3 = j2 + (-2)**(j2/3)
                        if (ABS(j1-j2) == 2)  then
                            j1 = 3
                            j2 = 1
                            j3 = 2
                        end if
                        ! solve average angular flux,surface angular flux and transfer value
                        ! j3 in radial, zi in axial is known
                        ! sp2x, sp3x: old value for unkown surface in radial
                        surf_3 = iter_flux%dist%surface(j3,impinge,ia,is)
                        surf_z2 = iter_flux%dist%surface(4+k/5,impinge,ia,is)
                        
                        u1 = quad%directions(is)%projection(j1,impinge)
                        u2 = quad%directions(is)%projection(j2,impinge)
                        u3 = quad%directions(is)%projection(j3,impinge)
                        
                        if (sigt < EPS_ZERO)  then
                            avgflux = (-2*u3*surf_3 + ABS(uz)/hlayer*surf_z2)/(2*u1 + 2*u2 + ABS(uz)/hlayer)
                            surf_z1 = avgflux
                            surf_1 = avgflux
                            surf_2 = avgflux
                            psi3 = (5*(surf_1+surf_2)-surf_3)/3.0 - 2*avgflux
                        else 
                        
                            ps01 = coeff_surface%rad%third(0,j1,il,is)
                            ps11 = coeff_surface%rad%third(1,j1,il,is)
                            ps21 = coeff_surface%rad%third(2,j1,il,is)
                            ps02 = coeff_surface%rad%third(0,j2,il,is)
                            ps12 = coeff_surface%rad%third(1,j2,il,is)
                            ps22 = coeff_surface%rad%third(2,j2,il,is)
                            
                            psi1 = iter_flux%dist%point(mesh%point(j1,impinge),ia,is,ig)
                            psi2 = iter_flux%dist%point(mesh%point(j2,impinge),ia,is,ig)
                            ra1dm = u2*(2/3.0*ps11 + ps21)
                            ra2dm = u1*(2/3.0*ps12 + ps22)
                            
                            sp11 = sigt
                            sp12 = 2 * u1
                            sp13 = 2 * u2
                            sp14 = ABS(uz) / hlayer
                            sp15 = iter_q%info%total_moment(0,il,isa) - 2*u3*surf_3 + sp14*surf_z2
                            
                            sp21 = 6*ra1dm - coeff_surface%rad%first(j1,il,is)
                            sp22 = 1 - ra1dm*5
                            sp23 = u2*(ps21-4/3.0*(ps01+ps11))
                            sp25 = coeff_surface%rad%second(1,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)     &
                                &  + coeff_surface%rad%second(2,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)     &
                                &  + surf_3*(u3*(4/3.0*ps01-2*ps11-6*ps21)-ra1dm)                                 &
                                &  + psi1*(-ps01/3.0+ps21*3)*(u2+u3) + psi2*u3*(2*ps11+3*ps21)
                            
                            sp31 = 6*ra2dm - coeff_surface%rad%first(j2,il,is)
                            sp32 = u1*(ps22-4/3.0*(ps02+ps12))
                            sp33 = 1 - ra2dm*5
                            sp34 = 0.0
                            sp35 = coeff_surface%rad%second(1,j2,il,is)*iter_q%info%total_moment(2*j2-1,il,isa)     &
                                &  + coeff_surface%rad%second(2,j2,il,is)*iter_q%info%total_moment(2*j2,il,isa)     &
                                &  + surf_3*(u3*(4/3.0*ps02-2*ps12-6*ps22)-ra2dm)                                 &
                                &  + psi2*(-ps02/3.0+ps22*3)*(u1+u3) + psi1*u3*(2*ps12+3*ps22)
                            
                            sp41 = 1 / coeff_surface%axi%first(il,is)
                            sp45 = ((coeff_surface%axi%third(il,is)-coeff_surface%axi%first(il,is))*surf_z2+deltaz)*sp41
                            
                            denominator = (sp41*sp34+sp31)*(sp23*sp12-sp13*sp22) + sp21*(sp13*sp32-sp33*sp12)   &
                                & + (sp41*sp14+sp11)*(sp22*sp33-sp32*sp23)
                                
                            if (ABS(denominator) >= EPS_ZERO)  then
                                avgflux = ((sp22*sp33-sp32*sp23)*(sp45*sp14+sp15)+(sp35+sp34*sp45)         &
                                    &  *(sp23*sp12-sp13*sp22)+sp25*(sp13*sp32-sp33*sp12)) / denominator
                                surf_1 = (sp21*(sp13*sp35-sp15*sp33)+(sp25*sp33-sp23*sp35)             &
                                    &  *(sp14*sp41+sp11)+(sp41*sp34+sp31)*(sp23*sp15-sp13*sp25)       &
                                    &  +sp45*(sp34*(sp21*sp13-sp11*sp23)+sp14*(sp23*sp31              &
                                    &  -sp33*sp21))) / denominator  
                                surf_2 = (sp41*(sp22*(sp35*sp14-sp15*sp34)+sp12*(sp25*sp34)            &
                                    &  +sp32*(-sp14*sp25))+(sp35+sp45*sp34)*(sp11*sp22-sp12*sp21)     &
                                    &  +sp25*(sp12*sp31-sp11*sp32)+(sp32*sp21-sp22*sp31)              &
                                    &  *(sp15+sp14*sp45)) / denominator
                                surf_z1 = (sp41*(sp23*(sp35*sp12-sp32*sp15)+sp25*(sp32*sp13-sp33        &
                                    &  *sp12)+sp22*(sp15*sp33-sp13*sp35))+sp45*(sp12*(sp33*sp21       &
                                    &  -sp23*sp31)+sp13*(sp22*sp31-sp32*sp21)+sp11*(sp32*sp23         &
                                    &  -sp33*sp22)))/ denominator
                            else
                                avgflux = 0.0
                                surf_1 = 0.0
                                surf_2 = 0.0
                                surf_z1 = 0.0
                            end if
                            psi3 = (5*(surf_1+surf_2)-surf_3)/3.0 - 2*avgflux
                        end if 
                        
                        ! negative flux fix-up
                        if (.NOT. is_adjoint)  then
                            call Truncate_negative (avgflux)
                            call Truncate_negative (surf_1)
                            call Truncate_negative (surf_2)
                            call Truncate_negative (surf_z1)
                        else 
                            call Truncate_negative (avgflux)
                            call Truncate_negative (surf_1)
                            call Truncate_negative (surf_2)
                            call Truncate_negative (surf_z1)
                        end if
                        
                        
                        if (.NOT. is_adjoint)  then
                            call Truncate_negative (psi3)
                        else 
                            call Truncate_negative (psi3)
                        end if
                        
                        iter_flux%dist%nodal(impinge,ia,is,ig) = avgflux * SOR_factor + (1.0-SOR_factor) * accele_SOR%nodal(impinge,ia,is)
                                
                        iter_flux%dist%surface(j1,impinge,ia,is) = surf_1 
                        iter_flux%dist%surface(j2,impinge,ia,is) = surf_2 
                        iter_flux%dist%surface(5-k/5,impinge,ia,is) = surf_z1 
                        iter_flux%dist%point(mesh%point(j3,impinge),ia,is,ig) = psi3 
                        
                        if (mesh%nearby_nodal(j1,impinge) /= impinge)  then
                            ix0 = mesh%localID(j1,impinge)
                            iy0 = mesh%nearby_nodal(j1,impinge)
                            iter_flux%dist%surface(ix0,iy0,ia,is) = surf_1
                        end if
                        
                        if (mesh%nearby_nodal(j2,impinge) /= impinge)  then 
                            ix0 = mesh%localID(j2,impinge)
                            iy0 = mesh%nearby_nodal(j2,impinge)
                            iter_flux%dist%surface(ix0,iy0,ia,is) = surf_2
                        end if
                        
                        if (ia/=ia_end)  then
                            iter_flux%dist%surface(4+k/5,impinge, ia+ia_step,is) = surf_z1
                        end if
                    end if
                end do
                
            end do flux
            
            ! ------------------------------------------------------------------
            ! get iteration save for boundary surface
            do ia = ia_start, ia_end, ia_step
                i_count = 0
                
                ! for radial surface
                do ir = 1, ns%state%nodal
                    do iy = 1, 3
                        if (ABS(bound%nodal(iy,ir)-bound%INNER) >= EPS_ZERO)  then
                            i_count = i_count + 1
                            iter_flux%dist%rad_surf(i_count,ia,is,ig) = iter_flux%dist%surface(iy,ir,ia,is)
                        end if
                    end do
                end do
                
                ! for axial surface
                do ir = 1, ns%state%nodal
                    if (ia == ia_start)  then
                       iter_flux%dist%axi_surf(ir,k/5+1,is,ig) = iter_flux%dist%surface(k/5+4,ir,ia,is)
                    end if
                    if (ia == ia_end)  then
                       iter_flux%dist%axi_surf(ir,2-k/5,is,ig) = iter_flux%dist%surface(5-k/5,ir,ia,is)
                    end if
                end do
                
!!!                ! for axial surface
!!!                do ir = 1, ns%state%nodal
!!!                   if (ia == 1)  then
!!!                       iter_flux%dist%axi_surf(ir,1,is,ig) = iter_flux%dist%surface(4,ir,ia,is)
!!!                   end if
!!!                   if (ia == ns%state%layer)  then
!!!                       iter_flux%dist%axi_surf(ir,2,is,ig) = iter_flux%dist%surface(5,ir,ia,is)
!!!                   end if
!!!                end do
            end do
            
        end do direction
        !$omp end do
        !$omp end parallel 
        
        do is = 1, ns%deduce%direction
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    il = (ia-1)*ns%state%nodal + ir
                    weighted_flux = iter_flux%dist%nodal(ir,ia,is,ig) * quad%directions(is)%wmu
                    iter_flux%info%moment(0,il,ig) = iter_flux%info%moment(0,il,ig) + weighted_flux
                    
                    if (ns%state%scat_order /= 0) then
                        do ic1 = 1, ns%state%scat_order
                            flux_scat%ngs(ig)%aniso_zero(ic1,0,il) = flux_scat%ngs(ig)%aniso_zero(ic1,0,il) + weighted_flux*coeff_source%pl_zero(ic1,is)
                            do ic2 = 1, ic1
                                flux_scat%ngs(ig)%aniso_cos(ic2,ic1,0,il) = flux_scat%ngs(ig)%aniso_cos(ic2,ic1,0,il) + weighted_flux*coeff_source%pl_cos(ic2,ic1,is)
                                flux_scat%ngs(ig)%aniso_sin(ic2,ic1,0,il) = flux_scat%ngs(ig)%aniso_sin(ic2,ic1,0,il) + weighted_flux*coeff_source%pl_sin(ic2,ic1,is)
                            end do
                        end do
                    end if
                end do
            end do
        end do
        
    end subroutine Get_average_flux
    
    !$
    !===============================================================================================
    ! add in-group source moments to get total source moments
    !===============================================================================================
    subroutine Add_in_group_source_moments (ig)
        
        ! intent parameters
        integer, intent(in)  :: ig
        
        integer  :: ia, ir, il, id, ml, is, ic1, ic2, iz
        real(KREAL)  :: summary                                             ! hold for anosotropic scatter source 
        
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                il = ml + ir
                do id = 1, 8
                    iter_q%info%total_moment(id,il,1) = iter_q%info%out_group_moment(id,il,1) + iter_flux%info%moment(id,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
                end do
            end do
        end do
        
        ! add in-group-source to total source when anisotropic scatter
        if (ns%state%scat_order /= 0)  then
            do is = 2, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do id = 1, 8
                            iter_q%info%total_moment(id,il,is) = iter_q%info%out_group_moment(id,il,is) + iter_flux%info%moment(id,il,ig)*xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,1)
                        end do
                    end do
                end do
            end do
            
            do is = 1, ns%deduce%direction
                do ia = 1, ns%state%layer
                    ml = (ia-1) * ns%state%nodal
                    do ir = 1, ns%state%nodal
                        iz = mesh%zone(ir)
                        il = ml + ir
                        do id = 1, 8
                            do ic1 = 1, ns%state%scat_order
                                summary = 0.0
                                do ic2 = 1, ic1
                                    summary = summary + coeff_source%pl_cos(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_cos(ic2,ic1,id,il)                                                            &
                                        &  + coeff_source%pl_sin(ic2,ic1,is)*flux_scat%ngs(ig)%aniso_sin(ic2,ic1,id,il)
                                end do
                                iter_q%info%total_moment(id,il,is) = iter_q%info%total_moment(id,il,is) + (flux_scat%ngs(ig)%aniso_zero(ic1,id,il)*coeff_source%pl_zero(ic1,is)+2*summary)    &
                                        &   * xsec_iter%matrixs(iz,ia)%sigma_s(ig,ig,ic1+1)
                            end do
                        end do
                    end do
                end do
            end do
        end if
    
    end subroutine Add_in_group_source_moments
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Re_generate_flux_moments(ig)
        
        integer, intent(in)  :: ig
        
        integer  :: k, ia, ir, id, ic1, ic2, il, is
        
        ! ----------------------------------------------------------------------
        ! local variables
        real(KREAL)  :: a01, a02, a03                                     
        real(KREAL)  :: a11, a12, a13
        real(KREAL)  :: a21, a22, a23
        real(KREAL)  :: psi1, psi2, psi3
       
        real(KREAL)  :: u1, u2, u3, uz
        real(KREAL)  :: surf_1,surf_2,surf_3
        real(KREAL)  :: avgflux
        real(KREAL)  :: SNwt(8)

        real(KREAL)  :: hlayer
        real(KREAL)  :: sigt
        
        integer  :: ml, j1, j2, j3
        integer  :: SN90                                                        ! number of directions per octant
        integer  :: ia_start, ia_end, ia_step                                   ! axial layer sweepping
        integer  :: isa, ks, impinge
        integer  :: iy
    
        ! re-set high-order flux moments to zero
        iter_flux%info%moment(1:8, :, ig) = 0.0D0
        if (ns%state%scat_order /= 0)  then
            call flux_scat%moment_zero (ig)
        end if
        
        ! sweep space-angle mesh for high-order moments 
        SN90 = ns%deduce%direction / 8
        
        !$omp parallel default(shared)
        !$omp do schedule(static, 1) ordered private(is, k, ia_start, ia_end, ia_step, ia, hlayer, ml, ir, impinge, sigt, il,       &
        !$omp &  isa, ks, iy, id, j1, j2, j3, surf_1, surf_2, surf_3, psi1, psi2, psi3, u1, u2, u3,     &
        !$omp &  a01, a11, a21, a02, a12, a22, a03, a13, a23, avgflux, ic1, ic2)                         
        do is = ns%deduce%direction, 1, -1
            k = ((is-1) - MOD(is-1, SN90))/SN90 + 1
            
            ia_start = (ns%state%layer)**(k/5)
            ia_end = ns%state%layer / ia_start
            ia_step = 1 - 2*(k/5)
            
            do ia = ia_start, ia_end, ia_step
                hlayer = geom%height(ia)
                ml = (ia-1) * ns%state%nodal
                do ir = 1, ns%state%nodal
                    impinge = sweep%order(ir,is)
                    sigt = xsec_iter%matrixs(mesh%zone(impinge),ia)%sigma_t(ig)
                    il = ml + impinge
                    if (ns%state%scat_order == 0)  then
                        isa = 1
                    else
                        isa = is
                    end if
                    
                    ! solve the right higher moment in radial plant
                    ! count for unknown surface
                    ks = 0
                    do iy = 1, 3
                        if (quad%directions(is)%projection(iy,impinge) > 0.0D0)  then
                            ks = ks + 1
                            if (ks == 1) j1 = iy
                            if (ks == 2) j2 = iy
                        end if
                    end do
                    
                    ! one surface unknown, (J2,J3 income surfaces)
                    if (ks == 1)  then
                        j2 = j1 + (-2)**(j1/3)
                        j3 = j2 + (-2)**(j2/3)
                        
                        avgflux = iter_flux%dist%nodal(impinge,ia,is,ig)
                        surf_1 = iter_flux%dist%surface(j1,impinge,ia,is)
                        surf_2 = iter_flux%dist%surface(j2,impinge,ia,is)
                        surf_3 = iter_flux%dist%surface(j3,impinge,ia,is)
                        psi1 = iter_flux%dist%point(mesh%point(j1,impinge),ia,is,ig)
                        psi2 = iter_flux%dist%point(mesh%point(j2,impinge),ia,is,ig)
                        psi3 = iter_flux%dist%point(mesh%point(j3,impinge),ia,is,ig)
                        u1 = quad%directions(is)%projection(j1,impinge)
                        u2 = quad%directions(is)%projection(j2,impinge)
                        u3 = quad%directions(is)%projection(j3,impinge)
                        
                        a01 = ((4*surf_2-psi1)*u2+(4*surf_3-psi1)*u3) / 3.0D0
                        a11 = ((psi3-surf_2)*u2+(psi2-surf_3)*u3) * 2.0D0
                        a21 = ((psi1+psi3-2*surf_2)*u2+(psi1+psi2-2*surf_3)*u3) * 3.0D0
                        a02 = ((4*surf_3-psi2)*u3+(4*surf_1-psi2)*u1) / 3.0D0
                        a12 = ((psi1-surf_3)*u3+(psi3-surf_1)*u1) * 2.0D0
                        a22 = ((psi2+psi1-2*surf_3)*u3+(psi2+psi3-2*surf_1)*u1) * 3.0D0
                        a03 = ((4*surf_1-psi3)*u1+(4*surf_2-psi3)*u2) / 3.0D0
                        a13 = ((psi2-surf_1)*u1+(psi1-surf_2)*u2) * 2.0D0
                        a23 = ((psi2+psi3-2*surf_1)*u1+(psi1+psi3-2*surf_2)*u2) * 3.0D0 
                        
                        iter_flux%info%moment_omp(j1*2-1,il,is) = coeff_nodal%rad%first(1,j1,il,is)*avgflux             &
                            &  + coeff_nodal%rad%second(1,1,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,1,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)             &
                            &  + a01*coeff_nodal%rad%third(0,1,j1,il,is) + a11*coeff_nodal%rad%third(1,1,j1,il,is)      &
                            &  + a21*coeff_nodal%rad%third(2,1,j1,il,is)           
                        iter_flux%info%moment_omp(j1*2,il,is) = coeff_nodal%rad%first(2,j1,il,is)*avgflux               &
                            &  + coeff_nodal%rad%second(1,2,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,2,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)             &
                            &  + a01*coeff_nodal%rad%third(0,2,j1,il,is) + a11*coeff_nodal%rad%third(1,2,j1,il,is)      &
                            &  + a21*coeff_nodal%rad%third(2,2,j1,il,is)
                        
                        iter_flux%info%moment_omp(j2*2-1,il,is) = coeff_nodal%rad%first(1,j2,il,is)*avgflux             &
                            &  + coeff_nodal%rad%second(1,1,j2,il,is)*iter_q%info%total_moment(2*j2-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,1,j2,il,is)*iter_q%info%total_moment(2*j2,il,isa)             &
                            &  + a02*coeff_nodal%rad%third(0,1,j2,il,is) + a12*coeff_nodal%rad%third(1,1,j2,il,is)      &
                            &  + a22*coeff_nodal%rad%third(2,1,j2,il,is) + coeff_nodal%rad%fourth(1,j2,il,is)*surf_2
                        iter_flux%info%moment_omp(j2*2,il,is) = coeff_nodal%rad%first(2,j2,il,is)*avgflux               &
                            &  + coeff_nodal%rad%second(1,2,j2,il,is)*iter_q%info%total_moment(2*j2-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,2,j2,il,is)*iter_q%info%total_moment(2*j2,il,isa)             &
                            &  + a02*coeff_nodal%rad%third(0,2,j2,il,is) + a12*coeff_nodal%rad%third(1,2,j2,il,is)      &
                            &  + a22*coeff_nodal%rad%third(2,2,j2,il,is) + coeff_nodal%rad%fourth(2,j2,il,is)*surf_2  
                        
                        iter_flux%info%moment_omp(j3*2-1,il,is) = coeff_nodal%rad%first(1,j3,il,is)*avgflux             &
                            &  + coeff_nodal%rad%second(1,1,j3,il,is)*iter_q%info%total_moment(2*j3-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,1,j3,il,is)*iter_q%info%total_moment(2*j3,il,isa)             &
                            &  + a03*coeff_nodal%rad%third(0,1,j3,il,is) + a13*coeff_nodal%rad%third(1,1,j3,il,is)      &
                            &  + a23*coeff_nodal%rad%third(2,1,j3,il,is) + coeff_nodal%rad%fourth(1,j3,il,is)*surf_3     
                        iter_flux%info%moment_omp(j3*2,il,is) = coeff_nodal%rad%first(2,j3,il,is)*avgflux               &
                            &  + coeff_nodal%rad%second(1,2,j3,il,is)*iter_q%info%total_moment(2*j3-1,il,isa)           &
                            &  + coeff_nodal%rad%second(2,2,j3,il,is)*iter_q%info%total_moment(2*j3,il,isa)             &
                            &  + a03*coeff_nodal%rad%third(0,2,j3,il,is) + a13*coeff_nodal%rad%third(1,2,j3,il,is)      &
                            &  + a23*coeff_nodal%rad%third(2,2,j3,il,is) + coeff_nodal%rad%fourth(2,j3,il,is)*surf_3 
                        
                        iter_flux%info%moment_omp(7,il,is) = coeff_nodal%axi%first(1,il,is)*avgflux                     &
                            &  + coeff_nodal%axi%third(1,il,is)*iter_flux%dist%surface(4+k/5,impinge,ia,is)             &
                            &  + coeff_nodal%axi%second(1,1,il,is)*iter_q%info%total_moment(7,il,isa)                   &
                            &  + coeff_nodal%axi%second(2,1,il,is)*iter_q%info%total_moment(8,il,isa)   
                        iter_flux%info%moment_omp(8,il,is) = coeff_nodal%axi%first(2,il,is)*avgflux                     &
                            &  + coeff_nodal%axi%third(2,il,is)*iter_flux%dist%surface(4+k/5,impinge,ia,is)             &
                            &  + coeff_nodal%axi%second(1,2,il,is)*iter_q%info%total_moment(7,il,isa)                   &
                            &  + coeff_nodal%axi%second(2,2,il,is)*iter_q%info%total_moment(8,il,isa)
                    end if
                    
                    ! two surface unknown, (J3 income surfaces£¬J1,J2 outcome)
                    if (ks == 2)  then      
                        j3 = j2 + (-2)**(j2/3)
                        if (ABS(j1-j2) == 2)  then
                            j1 = 3
                            j2 = 1
                            j3 = 2
                        end if
                        
                        avgflux = iter_flux%dist%nodal(impinge,ia,is,ig)
                        surf_1 = iter_flux%dist%surface(j1,impinge,ia,is)
                        surf_2 = iter_flux%dist%surface(j2,impinge,ia,is)
                        surf_3 = iter_flux%dist%surface(j3,impinge,ia,is)
                        psi1 = iter_flux%dist%point(mesh%point(j1,impinge),ia,is,ig)
                        psi2 = iter_flux%dist%point(mesh%point(j2,impinge),ia,is,ig)
                        psi3 = iter_flux%dist%point(mesh%point(j3,impinge),ia,is,ig)
                        u1 = quad%directions(is)%projection(j1,impinge)
                        u2 = quad%directions(is)%projection(j2,impinge)
                        u3 = quad%directions(is)%projection(j3,impinge)
                        
                        a01 = ((4*surf_2-psi1)*u2+(4*surf_3-psi1)*u3) / 3.0D0
                        a11 = ((psi3-surf_2)*u2+(psi2-surf_3)*u3) * 2.0D0
                        a21 = ((psi1+psi3-2*surf_2)*u2+(psi1+psi2-2*surf_3)*u3) * 3.0D0
                        a02 = ((4*surf_3-psi2)*u3+(4*surf_1-psi2)*u1) / 3.0D0
                        a12 = ((psi1-surf_3)*u3+(psi3-surf_1)*u1) * 2.0D0
                        a22 = ((psi2+psi1-2*surf_3)*u3+(psi2+psi3-2*surf_1)*u1) * 3.0D0 
                        a03 = ((4*surf_1-psi3)*u1+(4*surf_2-psi3)*u2) / 3.0D0
                        a13 = ((psi2-surf_1)*u1+(psi1-surf_2)*u2) * 2.0D0
                        a23 = ((psi2+psi3-2*surf_1)*u1+(psi1+psi3-2*surf_2)*u2) * 3.0D0 
                        
                        iter_flux%info%moment_omp(j1*2-1,il,is) = coeff_nodal%rad%first(1,j1,il,is)*avgflux                     &
                            &  + coeff_nodal%rad%second(1,1,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,1,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)                     &
                            &  + a01*coeff_nodal%rad%third(0,1,j1,il,is) + a11*coeff_nodal%rad%third(1,1,j1,il,is)              &
                            &  + a21*coeff_nodal%rad%third(2,1,j1,il,is)                                    
                        iter_flux%info%moment_omp(j1*2,il,is) = coeff_nodal%rad%first(2,j1,il,is)*avgflux                       &
                            &  + coeff_nodal%rad%second(1,2,j1,il,is)*iter_q%info%total_moment(2*j1-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,2,j1,il,is)*iter_q%info%total_moment(2*j1,il,isa)                     &
                            &  + a01*coeff_nodal%rad%third(0,2,j1,il,is) + a11*coeff_nodal%rad%third(1,2,j1,il,is)              &
                            &  + a21*coeff_nodal%rad%third(2,2,j1,il,is)        
                            
                        iter_flux%info%moment_omp(j2*2-1,il,is) = coeff_nodal%rad%first(1,j2,il,is)*avgflux                     &
                            &  + coeff_nodal%rad%second(1,1,j2,il,is)*iter_q%info%total_moment(2*j2-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,1,j2,il,is)*iter_q%info%total_moment(2*j2,il,isa)                     &
                            &  + a02*coeff_nodal%rad%third(0,1,j2,il,is) + a12*coeff_nodal%rad%third(1,1,j2,il,is)              &
                            &  + a22*coeff_nodal%rad%third(2,1,j2,il,is)                                    
                        iter_flux%info%moment_omp(j2*2,il,is) = coeff_nodal%rad%first(2,j2,il,is)*avgflux                       &
                            &  + coeff_nodal%rad%second(1,2,j2,il,is)*iter_q%info%total_moment(2*j2-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,2,j2,il,is)*iter_q%info%total_moment(2*j2,il,isa)                     &
                            &  + a02*coeff_nodal%rad%third(0,2,j2,il,is) + a12*coeff_nodal%rad%third(1,2,j2,il,is)              &
                            &  + a22*coeff_nodal%rad%third(2,2,j2,il,is)    
                            
                        iter_flux%info%moment_omp(j3*2-1,il,is) = coeff_nodal%rad%first(1,j3,il,is)*avgflux                     &
                            &  + coeff_nodal%rad%second(1,1,j3,il,is)*iter_q%info%total_moment(2*j3-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,1,j3,il,is)*iter_q%info%total_moment(2*j3,il,isa)                     &
                            &  + coeff_nodal%rad%fourth(1,j3,il,is)*surf_3 + a03*coeff_nodal%rad%third(0,1,j3,il,is)            &
                            &  + a13*coeff_nodal%rad%third(1,1,j3,il,is) + a23*coeff_nodal%rad%third(2,1,j3,il,is)             
                        iter_flux%info%moment_omp(j3*2,il,is) = coeff_nodal%rad%first(2,j3,il,is)*avgflux                       &
                            &  + coeff_nodal%rad%second(1,2,j3,il,is)*iter_q%info%total_moment(2*j3-1,il,isa)                   &
                            &  + coeff_nodal%rad%second(2,2,j3,il,is)*iter_q%info%total_moment(2*j3,il,isa)                     &
                            &  + coeff_nodal%rad%fourth(2,j3,il,is)*surf_3 + a03*coeff_nodal%rad%third(0,2,j3,il,is)            &
                            &  + a13*coeff_nodal%rad%third(1,2,j3,il,is) + a23*coeff_nodal%rad%third(2,2,j3,il,is)  
                        
                        iter_flux%info%moment_omp(7,il,is) = coeff_nodal%axi%first(1,il,is)*avgflux                             &
                            &  + coeff_nodal%axi%third(1,il,is)*iter_flux%dist%surface(4+k/5,impinge,ia,is)                     &
                            &  + coeff_nodal%axi%second(1,1,il,is)*iter_q%info%total_moment(7,il,isa)                           &
                            &  + coeff_nodal%axi%second(2,1,il,is)*iter_q%info%total_moment(8,il,isa)                             
                        iter_flux%info%moment_omp(8,il,is) = coeff_nodal%axi%first(2,il,is)*avgflux                             &
                            &  + coeff_nodal%axi%third(2,il,is)*iter_flux%dist%surface(4+k/5,impinge,ia,is)                     &
                            &  + coeff_nodal%axi%second(1,2,il,is)*iter_q%info%total_moment(7,il,isa)                           &
                            &  + coeff_nodal%axi%second(2,2,il,is)*iter_q%info%total_moment(8,il,isa)
                    end if
                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel 
        
        do is = 1, ns%deduce%direction
            do il = 1, ns%deduce%nodal_total
                do id = 1, 8
                    SNwt(id) = iter_flux%info%moment_omp(id,il,is)*quad%directions(is)%wmu
                    iter_flux%info%moment(id,il,ig) = iter_flux%info%moment(id,il,ig) + SNwt(id)
                end do
                if (ns%state%scat_order /= 0)  then
                    do id = 1, 8
                        do ic1 = 1, ns%state%scat_order
                            flux_scat%ngs(ig)%aniso_zero(ic1,id,il) = flux_scat%ngs(ig)%aniso_zero(ic1,id,il) + SNwt(id)*coeff_source%pl_zero(ic1,is)
                            do ic2 = 1, ic1
                                flux_scat%ngs(ig)%aniso_cos(ic2,ic1,id,il) = flux_scat%ngs(ig)%aniso_cos(ic2,ic1,id,il) + SNwt(id)*coeff_source%pl_cos(ic2,ic1,is)
                                flux_scat%ngs(ig)%aniso_sin(ic2,ic1,id,il) = flux_scat%ngs(ig)%aniso_sin(ic2,ic1,id,il) + SNwt(id)*coeff_source%pl_sin(ic2,ic1,is)
                            end do
                        end do
                    end do
                end if
            end do
        end do
        
    end subroutine Re_generate_flux_moments
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_source_moment ()
        
        integer  :: ia, ml, ir, il, iz, id, iig
        
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                il = ml + ir
                iz = mesh%zone(ir)
                do id = 0, 8
                    iter_q%fission%moment(id,il) = 0.0D0
                    do iig = 1, ns%state%ng
                        iter_q%fission%moment(id,il) = iter_q%fission%moment(id,il) + xsec_iter%matrixs(iz,ia)%sigma_f_nu(iig)*iter_flux%info%moment(id,il,iig)
                    end do 
                end do
            end do
        end do 
    
    end subroutine Update_source_moment
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_forward_error (criteria, error_fission_rate, error_flux, fq_new, fq_old)
        
        type(IterationCriterion), intent(in)  :: criteria
        real(KREAL), intent(out)  :: error_fission_rate
        real(KREAL), intent(out)  :: error_flux
        real(KREAL), intent(out)  :: fq_new
        real(KREAL), intent(out)  :: fq_old
        
        real(KREAL)  :: old_rate(ns%state%layer * ns%state%nodal)
        real(KREAL)  :: new_rate(ns%state%layer * ns%state%nodal)
        real(KREAL)  :: division
        integer  :: ia, ml, ir, il
        
        fq_new = 0.0D0
        fq_old = 0.0D0
        new_rate = 0.0D0
        old_rate = 0.0D0
        
        do ia = 1, ns%state%layer
            ml = (ia-1) * ns%state%nodal
            do ir = 1, ns%state%nodal
                il = ml + ir
                new_rate(il) = iter_q%fission%moment(0,il)*geom%area(ir)*geom%height(ia)
                old_rate(il) = iter_q%fission%old(il)*geom%area(ir)*geom%height(ia)
                fq_new = fq_new + new_rate(il)
                fq_old = fq_old + old_rate(il)
            end do
        end do
        
        error_flux = get_vector_error (iter_q%fission%old, iter_q%fission%moment(0,:), criteria%error_type)
        error_fission_rate = get_vector_error (old_rate, new_rate, criteria%error_type)
  
        do il = 1, ns%deduce%nodal_total
            iter_q%fission%old(il) = iter_q%fission%moment(0,il)
        end do 
    
    end subroutine Get_forward_error
    
end module iteration_process
