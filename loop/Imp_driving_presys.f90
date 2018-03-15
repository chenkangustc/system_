module Imp_driving_presys
	use Imp_loop_global
	use Imp_inputcard
	contains
	subroutine driving_presys()
		implicit none
        !scan & read Ny
		call driving_input_read()
		!alloc
		call driving_alloc_loop()
		!call loop%alloc()
		call driving_init_loop()
		!read after alloc
		!call driving_plain_read()
		!init
		!call loop%init()
	end subroutine driving_presys
	
	subroutine driving_plain_scan()
		implicit none
		
	end subroutine driving_plain_scan
	
	subroutine driving_init_loop()
		implicit none
		call core%init()
		call IHX1%init()
		call PipeRI%init()
		call PipeIP%init()
		call PipePR%init()
		call pump1%init()		
    end subroutine driving_init_loop
	
	subroutine driving_alloc_loop()
		implicit none
		call IHX1%free()
		call IHX1%alloc()
	end subroutine driving_alloc_loop
	
	subroutine driving_free_loop()
		implicit none
		call IHX1%free()
	end subroutine driving_free_loop
end module Imp_driving_presys