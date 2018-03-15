module Imp_inputcard
	use Imp_loop_global
    use,intrinsic::ISO_FORTRAN_ENV
    use constants
	implicit none
	integer::file_unit,file_o
	integer,parameter,private::N_keyword=6
    integer,parameter,private::MAX_REAL_PARAMETER=50
	integer,parameter,private::MAX_INT_PARAMETER=50
    character(len=MAX_WORD_LEN),parameter::FILE_IN='input.case'
    character(len=MAX_WORD_LEN),parameter::FILE_OUT='output.txt'
    character(len=MAX_WORD_LEN)::INP_SECTION(N_keyword)
    
	
	contains
	subroutine driving_input_read()
		implicit none
		!local
		integer::io_error
		real::dummy_real(MAX_REAL_PARAMETER)
		integer::dummy_int(MAX_INT_PARAMETER)
		character(len=MAX_WORD_LEN)::aline
		character(len=MAX_WORD_LEN)::section_name,keyword
		! Variables
		call set_section_keyword()
        open(newunit=file_o,file=FILE_OUT,status='replace',action='write',iostat=io_error)
		open(newunit=file_unit,file=FILE_IN,status='old',action='read',iostat=io_error)      
        !read(unit=file_unit,fmt='(A)',iostat=io_error) aline
		do
			read(unit=file_unit,fmt='(A)',iostat=io_error) aline
            if(io_error==IOSTAT_END) exit
			read(unit=aline,fmt=*,iostat=io_error) section_name       
			! if(is_keyword(INP_SECTION,section_name)) then
			!     backspace(file_unit,iostat=io_error)
			!     exit
			! end if        
			if(is_keyword(INP_SECTION,section_name)) then
				select case(trim(adjustl(section_name)))
					case('pump')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5)
					pump1%I=dummy_real(1)
					pump1%He=dummy_real(2)
					pump1%Qe=dummy_real(3)
					pump1%omegae=dummy_real(4)
					pump1%yita=dummy_real(5)
					
					case('pipePR')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5)
					pipePR%length=dummy_real(1)
					pipePR%r=dummy_real(2)
					pipePR%theta=dummy_real(3)
					pipePR%Q=dummy_real(4)
					pipePR%T=dummy_real(5)
					
					case('reactor')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:7)
					core%length=dummy_real(1)
					core%r=dummy_real(2)
					core%theta=dummy_real(3)
					core%Q=dummy_real(4)
					core%T=dummy_real(5)
					
					case('pipeRI')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5)
					pipeRI%length=dummy_real(1)
					pipeRI%r=dummy_real(2)
					pipeRI%theta=dummy_real(3)
					pipeRI%Q=dummy_real(4)
					pipeRI%T=dummy_real(5)
					
					case('IHX')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5),dummy_int(1:2)
					IHX1%Lsingle=dummy_real(1)
					IHX1%AreaOuter=dummy_real(2)
					IHX1%Qouter=dummy_real(3)
					IHX1%Rtube=dummy_real(4)
					IHX1%Ti=dummy_real(5)
					IHX1%Ntube=dummy_int(1)
					IHX1%N=dummy_int(2)
					
					case('pipeIP')
					read(unit=aline,fmt=*,iostat=io_error) keyword,dummy_real(1:5)
					pipeIP%length=dummy_real(1)
					pipeIP%r=dummy_real(2)
					pipeIP%theta=dummy_real(3)
					pipeIP%Q=dummy_real(4)
					pipeIP%T=dummy_real(5)
				end select
			end if
			!print*,trim(aline)
			print*,trim(adjustl(section_name))
		end do
        close(file_unit)
	end subroutine driving_input_read
	
    subroutine Set_section_keyword()
        implicit none
        INP_SECTION(1:N_keyword)=['pump   ',   &
                                & 'pipePR ',   &
                                & 'reactor',   &        
                                & 'pipeRI ',   &
                                & 'IHX    ',   &
                                & 'pipeIP '    ]
    end subroutine Set_section_keyword
	
	function is_keyword(input,key) result(is_true)
        implicit none
        character(len=MAX_WORD_LEN),intent(in)::input(:)
        character(len=MAX_WORD_LEN),intent(in)::key
        !local
        logical::is_true
        integer::i,list
        list=size(input)
        is_true=.FALSE.
        do i=1,list,1
            if(trim(adjustl(input(i)))==trim(adjustl(key))) then
                is_true=.TRUE.
                exit
            endif
        end do
    end function is_keyword
end module Imp_inputcard