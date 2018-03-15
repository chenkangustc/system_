module Imp_driving_syspost
	use Imp_loop_global!common data
    use Imp_inputcard!IOfile
	use Imp_driving_presys, only:driving_free_loop
    implicit none
	contains
	subroutine driving_postsys()
        implicit none
        !free
        call driving_free_loop()
        !close
        close(FILE_IN)
        close(FILE_O)
    end subroutine driving_postsys
end module Imp_driving_syspost