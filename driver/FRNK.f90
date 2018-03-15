!$
!===================================================================================================
!
!Method:
!    The code based on --steady code (3D discrete ordinates nodal transport in triangular mesh by hllu@NECP lib.);
!                      --with PCQM (predictor-corrector quasi-static method);
!                      --with point kinetics approximation;
!    reference:
!        [1] Haoliang Lu, A nodal SN transport method for three-dimensional trangular-z geometry[J], 
!            Nuclear Engineering and Design, 2007, 830-839.
!        [2] Alain Hebert, Applied Reactor Physics[M], 2009, 331-336.
!        [3] Sandra Dulla, et al, The quasi-static method revisited[J], Progress in Nuclear Energy, 2008, 908-920.
!        [4] B. Ganapol, P. Picca et al, Benchmarks for the point kinetics equations[C], M&C 2013.
!
!Revised records:
!    Author                Date                 History                                
!    -------------         --------             ---------------------------------------
!    hllu@NECP lib.        2007.08              steady version
!    xxxx@NECP lib.        2016.08              xxxx
!    -------------         --------             ---------------------------------------
!
!===================================================================================================
program DAISY
    
    use global_state
    use driver_pre_process,                 only : Run_pre_process
    use driver_post_process,                only : Run_post_process

    use driver_steady,                      only : Run_steady, Run_search, Run_initial
    use driver_adjoint,                     only : Run_adjoint
    use driver_perturb,                     only : Run_perturb
    
    use driver_transient_theta,             only : Run_transient_theta
    use driver_transient_pcqs,              only : Run_transient_pcqs
    use driver_transient_pk,                only : Run_transient_pk
    
    ! use thermal driver
    use th_pre_process,                     only : Driving_th_pre_process
    use th_post_process,                    only : Driving_th_post_process
    use th_check_model,                     only : Driving_th_check_model
	
	!IMPC
	!use imp_assm_global
	use imp_driving_pre_process
	!use imp_driving_output
    !test
    use testNK2TH
    
    implicit none
	!local
	integer tNK2TH
    tNK2TH=0
    ! --------------------------------------------------------------------------
	!    call Driving_th_check_model ()
    
    call Run_pre_process ()
    
    call Driving_th_pre_process ()
	
	call Sys_pre_process()
	
	! --------------------------------------------------------------------------
	if (tNK2TH==1) then 
        call driving_testNK2TH()   
    endif
    ! --------------------------------------------------------------------------
    ! perturbation calculation
    if (nt%flag%is_perturb)  then
        call Run_steady ()
        call Run_adjoint ()
        call Run_perturb ()
        
    else 
        ! steady calculation
        if (.NOT. nt%flag%is_perturb .and. .NOT. nt%flag%is_transient)  then
            call Run_steady ()
            if (nt%flag%is_boron_search .OR. nt%flag%is_CR_rod)  then
                call Run_search ()
            end if 
            call Run_initial ()
            if (nt%flag%is_adjoint )  then
                call Run_adjoint ()
            end if
        end if
        

        ! transient calculation
        if (nt%flag%is_transient)  then
            call Run_steady ()
            call Run_search ()
            call Run_initial ()
            call Run_adjoint ()
        
            select case(TRIM(nt%method%scheme))
            case ('THETA', 'FIM')
                call Run_transient_theta ()
                
            case ('PCQS', 'PCQM') 
                call Run_transient_pcqs ()
            
            case ('PK', 'EPKM')
                call Run_transient_pk ()
            end select
        end if
    end if 
    ! --------------------------------------------------------------------------
	
    call Driving_th_post_process ()
    
    call Run_post_process ()

end program DAISY
