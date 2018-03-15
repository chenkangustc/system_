!$
!===================================================================================================
!
!   module for xsec transit before iteration
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Transit_xsec_steady
!                               Transit_xsec_adjoint
!                               Transit_xsec_gpt
! 
!                               Transit_xsec_theta
!                               Transit_xsec_pcqs
!
!   Public type lists:          No
!
!===================================================================================================
module transit_to_solver
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    
    use process_theta,                  only : Generate_TFSP
    
    implicit none 
    private
    public  :: Transit_xsec_steady, Transit_xsec_adjoint, Transit_xsec_gpt
    public  :: Transit_xsec_theta, Transit_xsec_pcqs
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_xsec_steady (is_eigen)
        
        logical, intent(in)  :: is_eigen
        
        call Xsec_to_iteration ()
        
        if (is_eigen )  then
            call Zero_source_to_iteration ()
            
        else
            call Steady_source_to_iteration ()
            
            if (ns%method%is_Ks  .and. xsec%is_fission () )  then
                call Q_ext%normal (geom, is_normal=.TRUE.)
            else
                call Q_ext%normal (geom, is_normal=.FALSE.)
            end if
        end if
        
    end subroutine Transit_xsec_steady
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_xsec_adjoint (is_eigen)
    
        logical, intent(in)  :: is_eigen
        
        call Xsec_to_iteration ()
        
        ! transpose scatter xsec
        call xsec_iter%transpose_scat ()
        
        if (is_eigen)  then
            call Zero_source_to_iteration ()
            
        else
            call Adjoint_source_to_iteration ()
            
            if (ns%method%is_Ks  .and. xsec%is_fission () )  then
                call Q_ext%normal (geom, is_normal=.TRUE.)
            else
                call Q_ext%normal (geom, is_normal=.FALSE.)
            end if
        end if
                
    end subroutine Transit_xsec_adjoint
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_xsec_gpt (eigenvalue)
    
        real(KREAL), intent(in)  :: eigenvalue
        
        call Xsec_to_iteration ()
        
        ! transpose scatter xsec
        call xsec_iter%transpose_scat ()
        call xsec_iter%critical (eigenvalue)
        
        call Adjoint_source_to_iteration ()
        
        if (ns%method%is_Ks  .and. xsec%is_fission () )  then
            call Q_ext%normal (geom, is_normal=.TRUE.)
        else
            call Q_ext%normal (geom, is_normal=.FALSE.)
        end if
                
    end subroutine Transit_xsec_gpt
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_xsec_theta (time_pace)
        
        real(KREAL), intent(in)  :: time_pace
        
        call Xsec_to_iteration ()
        
        call Steady_source_to_iteration ()
        
        ! re-arrange cross section to generate the TFSP(time dependent fixed source problem)
        call Generate_TFSP (time_pace)
               
        if (ns%method%is_Ks  .and. xsec%is_fission () )  then
            call Q_ext%normal (geom, is_normal=.TRUE.)
        else
            call Q_ext%normal (geom, is_normal=.FALSE.)
        end if
            
    end subroutine Transit_xsec_theta
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Transit_xsec_pcqs (time_pace)
    
        real(KREAL)  :: time_pace
    
        call Xsec_to_iteration ()
        
        call Steady_source_to_iteration ()
        
        ! re-arrange cross section to generate the TFSP(time dependent fixed source problem)
        call Generate_TFSP (time_pace)
        
        if (ns%method%is_Ks  .and. xsec%is_fission () )  then
            call Q_ext%normal (geom, is_normal=.TRUE.)
        else
            call Q_ext%normal (geom, is_normal=.FALSE.)
        end if
                        
    end subroutine Transit_xsec_pcqs 
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! transit status xsec to iteration xsec
    !===============================================================================================
    subroutine Xsec_to_iteration ()
    
        integer  :: iz, ia
        
        do iz = 1, ns%state%zone
            do ia = 1, ns%state%layer
                xsec_iter%matrixs(iz, ia) = xsec%matrixs(iz, ia)
            end do
        end do
    
    end subroutine Xsec_to_iteration
    
    !$
    !===============================================================================================
    ! transit status source to iteration source
    !===============================================================================================
    subroutine Steady_source_to_iteration ()
        
        integer  :: ir, ia, iz, is, ig

        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                Q_ext%iter_scalar(ir,ia)%intensity = 0.0
                Q_ext%iter_scalar(ir,ia)%intensity = Q_ext%matrixs(iz,ia)%intensity
            end do
        end do
    
    end subroutine Steady_source_to_iteration 
    
    !$
    !===============================================================================================
    ! steady calculation without source
    !===============================================================================================
    subroutine Zero_source_to_iteration ()
        
        integer  :: ir, ia, iz, is, ig

        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                Q_ext%iter_scalar(ir,ia)%intensity = 0.0
            end do
        end do
    
    end subroutine Zero_source_to_iteration 
    
    !$
    !===============================================================================================
    ! transit status source to iteration source
    !===============================================================================================
    subroutine Adjoint_source_to_iteration ()
        
        integer  :: ir, ia, iz, is, ig

        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                iz = mesh%zone(ir)
                Q_ext%iter_scalar(ir,ia)%intensity = 0.0
                do ig = 1, ns%state%ng
                    Q_ext%iter_scalar(ir,ia)%intensity(ig) = iter_adjoint%source(iz,ia,ig)
                end do
            end do
        end do
    
    end subroutine Adjoint_source_to_iteration 
    
end module transit_to_solver
