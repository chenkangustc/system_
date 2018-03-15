
!****************************************************************************
!
!  PROGRAM: system
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program loop
    use Imp_driving_presys
	use Imp_cal_loop
    use Imp_timer_global
    use Imp_inputcard
    use Imp_driving_syspost
    implicit none
    !local
    integer::i,Nt
    call Timer%set()
    Nt=Timer%Nt
	!call test_loop_hydraulic()
    call driving_presys()
    
    do i=1,Nt,1
        call driving_loop_transient(timer%current)
        call Timer%update()
        write(unit=file_o,fmt=*)timer%current,core%Q/pump1%Qe
    end do
    
    call driving_postsys()
    end program loop

