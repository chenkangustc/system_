!$
!===================================================================================================
!
!   control subroutine for steady calculation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module driver_steady 

    use constants
    use global_state    
    use, intrinsic  :: ISO_FORTRAN_ENV

    use global
    
    use stastics,                   only : stastics_average_value
    use transit_to_solver,          only : Transit_xsec_steady
    use iteration_initialize,       only : Init_iteration_variable
    use iteration_control,          only : Driving_iteration
    use feedback,                   only : Check_feedback_steady, Check_feedback_steady2
    use time_advancing,             only : Get_physics_kinetics_parameter
    
    use output_timelist,            only : Print_timelist
    use output_hdf5,                only : Print_binary_hdf5
    use output_visit,               only : Print_vtk_files
    
    use hdf5_interpolation,         only : Check_interpolation
    
    implicit none
    private
    public  :: Run_steady, Run_search, Run_initial
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_steady ()
        
        logical  :: is_pass
        integer  :: i, j
        integer  :: idx  
        
        idx = 0 
        nested: do 
            idx = idx + 1
            
!!!            call Check_interpolation ()
!!!            stop ('@#')
            
!            if (idx <= 3)  then
!                call xsec%print (900+idx)
!                call param%print (950+idx)
!            end if 
            
            call Transit_xsec_steady (is_eigen=ns%flag%is_eigen)
            call Init_iteration_variable (is_adjoint=.FALSE.)
            call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.FALSE.)
            
            if (.NOT. ns%feedback%is_feedback)  then
                exit nested
            end if
            
            call Check_feedback_steady (is_pass)
            
            if (is_pass)  then
                exit nested
            end if
        end do nested
        
    end subroutine Run_steady
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_search ()
        
        logical  :: is_pass = .FALSE.
        integer  :: ia, iz 
        real(KREAL), allocatable  :: volume(:, :)
        
        ! ----------------------------------------------------------------------
        ! for FSP problem, use the real flux as initial condition
        if (.NOT. ns%flag%is_eigen)  then
            
            ! normal external source
            call Q_ext%adjust (timelist%normal_factor)
            
            ! perform another transport
            call Init_iteration_variable (is_adjoint=.FALSE.)
            call Transit_xsec_steady (is_eigen=.FALSE.)
            call Driving_iteration (is_eigen=.FALSE., is_adjoint=.FALSE., is_transient=.TRUE., tidx=0, ctime=time_step%get_start ())
            
        ! ----------------------------------------------------------------------
        ! for eigenvalue problem
        else 
            if (ns%feedback%is_feedback)  then
                ! search critical boron
                if (nt%flag%is_boron_search .AND. self_link%is_active ("CB"))  then
                    ppm: do 
                        call self_fdbk%CB_search (iter_count%eigenvalue, is_pass)
                        if (is_pass)  then
                            exit
                        end if
                        call Check_feedback_steady (is_pass)
                        call Transit_xsec_steady (is_eigen=ns%flag%is_eigen)
                        call Init_iteration_variable (is_adjoint=.FALSE.)
                        call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.FALSE.)
                    end do ppm
                end if
                
                ! search critical CR
                if (cr_bank%is_search .AND. nt%flag%is_CR_rod)  then
                    rod1: do
                        call cr_bank%rod_search (iter_count%eigenvalue, is_pass)
                        call cr_bank%move (mat_info, cr_bank%init_step)
                        call cr_bank%homo (xsec_inp, param_inp, geom)
                        call cr_bank%map (xsec, param)
                        if (is_pass)  then
                            exit
                        end if
                        call Check_feedback_steady (is_pass)
                        call Transit_xsec_steady (is_eigen=ns%flag%is_eigen)
                        call Init_iteration_variable (is_adjoint=.FALSE.)
                        call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.FALSE.)
                    end do rod1
                end if 
                
                ! once-more, for initical condition, no-under-relation
                call Check_feedback_steady2 (is_pass)
                call Transit_xsec_steady (is_eigen=ns%flag%is_eigen)
                call Init_iteration_variable (is_adjoint=.FALSE.)
                call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.FALSE.)
            else 
                ! search critical CR
                if (cr_bank%is_search .AND. nt%flag%is_CR_rod)  then
                    rod2: do
                        call cr_bank%rod_search (iter_count%eigenvalue, is_pass)
                        call cr_bank%move (mat_info, cr_bank%init_step)
                        call cr_bank%homo (xsec_inp, param_inp, geom)
                        call cr_bank%map (xsec, param)
                        if (is_pass)  then
                            exit
                        end if
                        call Transit_xsec_steady (is_eigen=ns%flag%is_eigen)
                        call Init_iteration_variable (is_adjoint=.FALSE.)
                        call Driving_iteration (is_eigen=ns%flag%is_eigen, is_adjoint=.FALSE.)
                    end do rod2
                else 
                    ! force xsec to critical 
                    iter_count%kcritical = iter_count%eigenvalue
                    call xsec_inp%critical (iter_count%kcritical)
                    call xsec%critical (iter_count%kcritical)
                    if (nt%flag%is_CR_rod)  then
                        call cr_bank%critical (iter_count%kcritical)
                    end if
                end if 
!                ! perform another transport
!                call Init_iteration_variable (is_adjoint=.FALSE.)
!                call Transit_xsec_steady (is_eigen=.TRUE.)
!                call Driving_iteration (is_eigen=.TRUE., is_adjoint=.FALSE., is_transient=.TRUE., tidx=0, ctime=time_step%get_start ())
            end if 
        end if
        
    end subroutine Run_search
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Run_initial ()
        
        logical  :: is_pass = .FALSE.
        integer  :: ia, iz 
        real(KREAL), allocatable  :: volume(:, :)
        
        ! store initial xsec
        do ia = 1, ns%state%layer
            do iz = 1, ns%state%zone
                xsec_init%matrixs(iz,ia) = xsec%matrixs(iz,ia)
            end do
        end do
        
        if (nt%perturb%is_CR_move) then
            call pert_cr%init (cr_bank)
        end if 
        
        ! for SRAC format, perform transfer
        if ((nt%flag%is_transient) .and. (nt%flag%kinetics_type == 'SRAC'))  then
            call Get_physics_kinetics_parameter ()
        end if
        
        ! print the information of steady state into output.timelist
        call timelist%set (geom, mesh, dist_power, self_fdbk, 0.0_KREAL, 0.0_KREAL)
        call Print_timelist (timelist, 0, 0.0_KREAL, FILES%TIMELIST)
        call det%get (ns, geom, mesh, flux_forward)
        call det%pvalue (FILES%DET, 0, 0.0_KREAL)
!        call self_fdbk%print (200)
        
        call Print_binary_hdf5 (is_adjoint=.FALSE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        call Print_vtk_files (is_adjoint=.FALSE., is_transient=.FALSE., tidx=0, ctime=0.0_KREAL)
        
!        ! @ feedback-mapping
!        do ia = ns%state%layer_bottom+1, ns%state%layer - ns%state%layer_top
!            do iz = 1, ns%state%zone
!                if (mat_info%mask_core(iz, ia))  then
!                    write(999, fmt="(1x, 2(I5, TR2), *(ES12.5, TR2))")  ia, iz, self_fdbk%Tf%new(iz,ia), self_fdbk%Tm%new(iz,ia), self_fdbk%Rho_m%new(iz,ia)
!                end if  
!            end do 
!        end do 
        
    end subroutine Run_initial
    
end module driver_steady
