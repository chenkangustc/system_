!$
!===================================================================================================
!
!   scan input file and get necessity parameter, get several global parameter
!   ----------------------------------------------------------------------------
!   |ng, scat_order, dg, hg
!   |is_60symmetry, theta, sn
!   |n_parameter
!   |mat, source
!   |point, nodal, segment, zone
!   |layer
!   |layer_bottom, layer_top
!   |rod, bank
!   |n_assm_geom
!   !nf, nc
!   !n_pert
!   ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_plain_scan
!                               Driving_plain_read
!
!   Public type lists:          No
!
!===================================================================================================
module input_plain
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use geomregular_header,         only : GeomRec, GeomHex
    
    use global
    use th_global
	
	use imp_re_input_global
    
    use string 
    use files_dirs
    use input_keyword
    use input_xsec,                 only : xsec_read_known
    use output_main,                only : Print_header_main
    use output_timelist,            only : Print_header_timelist
    use output_visit 
    
    implicit none 
    private
    public  :: Driving_plain_scan, Driving_plain_read
    
    ! --------------------------------------------------------------------------
    type  Configure
        integer   :: n_matRow   = 0
        integer, allocatable  :: matRow_assign(:)
        integer, allocatable  :: matRow_config(:, :)
        
        integer   :: n_matCol   = 0
        integer, allocatable  :: matCol_assign(:)
        integer, allocatable  :: matCol_config(:, :)
        
        integer   :: n_extqRow   = 0
        integer, allocatable  :: extqRow_assign(:)
        integer, allocatable  :: extqRow_config(:, :)
        
        integer   :: n_extqCol  = 0
        integer, allocatable  :: extqCol_assign(:)
        integer, allocatable  :: extqCol_config(:, :)
    end type Configure
    
    type(Configure)  :: config
    
    type(GeomRec)  :: geom_rec 
    type(GeomHex)  :: geom_hex 
    
    ! --------------------------------------------------------------------------
    ! max parameter number for different kind
    integer, parameter, private  :: MAX_INT_PARAMETER     = 3000                 ! for assembly configure, so this more larger than other kind
    integer, parameter, private  :: MAX_REAL_PARAMETER    = 1000
    integer, parameter, private  :: MAX_CHAR_PARAMETER    = 1000
    integer, parameter, private  :: MAX_LOGICAL_PARAMETER = 1000
    
    ! dummy container 
    integer                         :: dummy_int(MAX_INT_PARAMETER)
    real(KREAL)                     :: dummy_real(MAX_REAL_PARAMETER)
    logical                         :: dummy_log(MAX_LOGICAL_PARAMETER)
    character(len=MAX_WORD_LEN)     :: dummy_char(MAX_CHAR_PARAMETER)
    

    character(len=MAX_WORD_LEN)  :: words(MAX_WORDS)                        ! hold for words after a line splited
    character(len=MAX_LINE_LEN)  :: aline                                  ! hold for a line content
    character(len=MAX_WORD_LEN)  :: keyword
    integer  :: n_word                                                      ! how many words of a line

contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Driving_plain_scan ()
    
        character(len=MAX_WORD_LEN)  :: section_name                            ! hold for section key word
        integer  :: io_error
        
        do 
            read(unit=FILES%CASENAME, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
        
            read(unit=aline, fmt=*, iostat=io_error) section_name
            
            if (Is_keyword(INP_SECTION, section_name))  then
                select case (TRIM(ADJUSTL(section_name)))
                case ('CASENAME:')
                    call Scan_casename (FILES%CASENAME, FILES%MAIN)
                    
                case ('CONTROL:')
                    call Scan_control (FILES%CASENAME, FILES%MAIN)
                    
                case ('METHOD:')
                    call Scan_method (FILES%CASENAME, FILES%MAIN)
                    
                case ('LINK:')
                    call Scan_link (FILES%CASENAME, FILES%MAIN)
                    
                case ('MATERIAL:')
                    call Scan_material (FILES%CASENAME, FILES%MAIN)
                    
                case ('GEOMETRY:')
                    call Scan_geometry (FILES%CASENAME, FILES%MAIN)
                    
                case ('TRANSIENT:')
                    call Scan_transient (FILES%CASENAME, FILES%MAIN)
                    
                case ('PKMODEL:')
                    call Scan_pkmodel (FILES%CASENAME, FILES%MAIN)
                    
                case ('FEEDBACK:')
                    call Scan_feedback (FILES%CASENAME, FILES%MAIN)
                    
                case ('PERTURBATION:')
                    call Scan_perturb (FILES%CASENAME, FILES%MAIN)
                    
                case ('END:')
                    exit 
                    
                end select
            
            ! input is out of pre-define range
            else 
                call this_error%set (INFO_LIST_INPUT, 'section key word input error: '//TRIM(ADJUSTL(section_name)))  
                call this_error%print (FILES%MAIN)
            end if 
        end do
    
    end subroutine Driving_plain_scan
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Driving_plain_read ()
    
        character(len=MAX_WORD_LEN)  :: section_name                            ! hold for section key word
        integer  :: io_error
        
        do 
            read(unit=FILES%CASENAME, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
        
            read(unit=aline, fmt=*, iostat=io_error) section_name
            
            if (Is_keyword(INP_SECTION, section_name))  then
                select case (TRIM(ADJUSTL(section_name)))
                case ('CASENAME:')
                    call Read_casename (FILES%CASENAME, FILES%MAIN)
                    
                case ('CONTROL:')
                    call Read_control (FILES%CASENAME, FILES%MAIN)
                    
                case ('METHOD:')
                    call Read_method (FILES%CASENAME, FILES%MAIN)
                    
                case ('LINK:')
                    call Read_link (FILES%CASENAME, FILES%MAIN)
                    
                case ('MATERIAL:')
                    call Read_material (FILES%CASENAME, FILES%MAIN)
                    
                case ('GEOMETRY:')
                    call Read_geometry (FILES%CASENAME, FILES%MAIN)
                    
                case ('TRANSIENT:')
                    call Read_transient (FILES%CASENAME, FILES%MAIN)
                    
                case ('PKMODEL:')
                    call Read_pkmodel (FILES%CASENAME, FILES%MAIN)
                    
                case ('FEEDBACK:')
                    call Read_feedback (FILES%CASENAME, FILES%MAIN)
                    
                case ('PERTURBATION:')
                    call Read_perturb (FILES%CASENAME, FILES%MAIN) 
                    
                case ('END:')
                    call this_information%set (INFO_LIST_INPUT, 'end of input parse succeed')
                    exit 
                    
                end select
            
            ! input is out of pre-define range
            else 
                call this_error%set (INFO_LIST_INPUT, 'section key word input error: '//TRIM(ADJUSTL(section_name)))  
                call this_error%print (FILES%MAIN)
            end if 
        end do

    end subroutine Driving_plain_read
    
    ! --------------------------------------------------------------------------
    ! scanning input file
    ! --------------------------------------------------------------------------

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_casename (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        integer  :: io_error

        ! re-get the case name
        backspace(file_in, iostat=io_error)
        
        read(unit=file_in, fmt="(A)", iostat=io_error) aline
        read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_char(1)
    
    end subroutine Scan_casename
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_control (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(CONTROL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('PK_MODEL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1) 
                    nt%flag%is_pkmodel = dummy_log(1)
            
                case ('STEADY')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1)
                    
                case ('PPM_SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1)
                
                case ('ROD_SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1), dummy_int(1)
                
                case ('ADJOINT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1) 

                case ('TRANSIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_char(1)
                    nt%flag%is_transient = dummy_log(1)
                    nt%flag%kinetics_type = dummy_char(1)
            
                case ('LINK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_char(1)
                    
                case ('FEEDBACK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    
                case ('PERTURBATION')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    
                case ('DECAY_HEAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    
                case ('ENERGY_GROUP') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%ng = dummy_int(1)
                    ns%state%scat_order = dummy_int(2)
                    
                case ('DELAYD_GROUP', 'DELAYED_GROUP') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    nt%state%dg = dummy_int(1)
                    
                case ('OUTPUT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:n_word-1)
                
                case ('THREAD')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                
                case ('INITVAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:2)
                
                case ('MISC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                
                end select
                    
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word)))
                call this_error%print (file_out)
            end if
        end do
        
    end subroutine Scan_control
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_method (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
                        
            read(unit=aline, fmt=*, iostat=io_error) a_word

            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(METHOD_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('QUDARATURE', 'QUADRATURE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1:2)
                    ns%flag%is_60degree = dummy_log(1)
                    ns%flag%n_theta = dummy_int(1)
                    ns%state%sn = dummy_int(2)
                
                case ('STEADY')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:2), dummy_int(1), dummy_log(3)
                
                case ('TRANSIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1), dummy_log(1), dummy_char(2)
                    
                case ('ERROR_TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                
                case ('ERROR_EIGEN')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:4)
                
                case ('ERROR_FSP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:4)

                end select

            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_method
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_link (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
                        
            read(unit=aline, fmt=*, iostat=io_error) a_word

            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(LINK_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('DATA')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    self_link%n_parameter = dummy_int(1)
                    
                case ('PARAMETER')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1:dummy_int(1))
                
                case ('REFERENCE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:dummy_int(1))

                case ('MAX_LIMIT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:dummy_int(1))

                case ('MIN_LIMIT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:dummy_int(1))

                end select

            else
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_link

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_material (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        integer  :: i_allocate
        
        character(len=MAX_WORD_LEN)  :: tmp_file
        real(KREAL), allocatable  :: tmp_xsec(:)                                ! tmp for xsec input
        integer  :: i, im, ic, ig, iig, ie, ip
        integer  :: scat_matrix
        integer  :: im_start, im_end                                            ! use for one group xsec adapt to material range
        logical  :: is_true
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(MATERIAL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('SCALE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%mat = dummy_int(1)
                    ns%state%source = dummy_int(2)
                    
                case ('IMPORT')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_log(1), dummy_int(1), dummy_char(1)
                    ns%misc%is_import = dummy_log(1)
                    
                case ('TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    ns%misc%is_external = dummy_log(1)
                    
                case ('MAT')
                    backspace(file_in, iostat=io_error)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3)
                    
                    if (.NOT. ns%misc%is_external)  then
                        allocate(tmp_xsec(ns%state%ng), stat=i_allocate)
                        read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                        read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                        read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                        read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                        
                        ! scatter matrix (ig-->iig, iig=1, ns%state%ng)
                        scat_matrix = ns%state%scat_order + 1
                        do ic = 1, scat_matrix
                            do ig = 1, ns%state%ng
                                read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            end do
                        end do
                        deallocate(tmp_xsec)
                        
                    ! check existence
                    else
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        tmp_file = TRIM(DIR_XSEC) // TRIM(ADJUSTL(keyword))
                        inquire(file=tmp_file, exist=is_true)
                        if (.NOT. is_true)  then
                            call this_error%set (INFO_LIST_INPUT, 'input file does not exist: '//TRIM(ADJUSTL(tmp_file))) 
                            call this_error%print (file_out)
                            stop
                        end if
                    end if
                
                case ('MAT_D')
                    backspace(file_in, iostat=io_error)

                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    
                    if (.NOT. ns%misc%is_external)  then
                        select case(TRIM(nt%flag%kinetics_type))
                        case ('NORMAL')
                            allocate(tmp_xsec(ns%state%ng), stat=i_allocate)
                            do ip = 1, nt%state%dg
                                read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            end do
                            read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            deallocate(tmp_xsec)
                            
                            allocate(tmp_xsec(nt%state%dg), stat=i_allocate)
                            read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            deallocate(tmp_xsec)
                            
                        case ('SRAC')
                            allocate(tmp_xsec(ns%state%ng), stat=i_allocate)
                            do ip = 1, nt%state%dg
                                read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            end do
                            read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            deallocate(tmp_xsec)
                            
                            allocate(tmp_xsec(nt%state%dg), stat=i_allocate)
                            do ig = 1, ns%state%ng
                                read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            end do
                            do ig = 1, ns%state%ng
                                read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                            end do
                            deallocate(tmp_xsec)
                        
                        case default
                        end select
                    end if
                    
                case ('Q_EXTERNAL')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    
                    if (.NOT. ns%misc%is_external)  then
                        allocate(tmp_xsec(ns%state%ng), stat=i_allocate)
                        read(unit=file_in, fmt=*, iostat=io_error)  tmp_xsec
                        deallocate(tmp_xsec)
                    
                    ! chech input file existence
                    else
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        tmp_file = TRIM(DIR_XSEC) // TRIM(ADJUSTL(keyword))
                        inquire(file=tmp_file, exist=is_true)
                        if (.NOT. is_true)  then
                            call this_error%set (INFO_LIST_INPUT, 'input file does not exist: '//TRIM(ADJUSTL(tmp_file))) 
                            call this_error%print (file_out)
                            stop
                        end if
                    end if
                    
                case ('DETINFO')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    det%nxs = dummy_int(1)
                    det%ndet = dummy_int(2)
                    
                case ('DETXS')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    
                    if (.NOT. ns%misc%is_external)  then
                        read(unit=file_in, fmt=*, iostat=io_error) dummy_real(1: ns%state%ng)
                    else
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        tmp_file = TRIM(DIR_XSEC) // TRIM(ADJUSTL(keyword))
                        inquire(file=tmp_file, exist=is_true)
                        if (.NOT. is_true)  then
                            call this_error%set (INFO_LIST_INPUT, 'input file does not exist: '//TRIM(ADJUSTL(tmp_file))) 
                            call this_error%print (file_out)
                            stop
                        end if
                    end if 
                    
                case ('DETPOINT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2:4)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_material
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_geometry (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_LINE_LEN)  :: tmp_file, out_file
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        integer  :: i, j
        integer  :: ibeg, iend 
        integer  :: ia
        logical  :: is_true
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(GEOMETRY_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    mesh%meshtype = dummy_int(1)
                    
                ! -------------------------------------------------------------- ansys-mesh
                ! --------------------------------------------------------------
                case ('TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                
                case ('ANSYS_FILE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1:2)
                    call Lower_case (dummy_char(1))
                    call Lower_case (dummy_char(2))
                    mesh%xyfile = dummy_char(1)
                    mesh%iifile = dummy_char(2)
                
                case ('XY_EXPAND')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2)
                 
                case ('SCALE') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:4)
                    ns%state%point = dummy_int(1)
                    ns%state%nodal =  dummy_int(2)
                    ns%state%segment = dummy_int(3)
                    ns%state%zone =  dummy_int(4)
                    
                case ('BC_RADIAL')
                    backspace(file_in, iostat=io_error)
                    do i = 1, ns%state%segment
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1)
                    end do
                    
                ! -------------------------------------------------------------- rec-mesh
                ! --------------------------------------------------------------
                case ('REC_XDIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1)
                    geom_rec%nx = dummy_int(1)
                    geom_rec%is_xdefine = .TRUE.
                    geom_rec%xchar = TRIM(ADJUSTL(aline))
                    if (geom_rec%is_xdefine .AND. geom_rec%is_ydefine)  then
                        call geom_rec%alloc ()
                        read(unit=geom_rec%xchar, fmt=*, iostat=io_error) keyword, dummy_int(1), geom_rec%xdim(1: geom_rec%nx)
                        read(unit=geom_rec%ychar, fmt=*, iostat=io_error) keyword, dummy_int(1), geom_rec%ydim(1: geom_rec%ny)
                    end if 
                    
                case ('REC_YDIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1)
                    geom_rec%ny = dummy_int(1)
                    geom_rec%is_ydefine = .TRUE.
                    geom_rec%ychar = TRIM(ADJUSTL(aline))
                    if (geom_rec%is_xdefine .AND. geom_rec%is_ydefine)  then
                        call geom_rec%alloc ()
                        read(unit=geom_rec%xchar, fmt=*, iostat=io_error) keyword, dummy_int(1), geom_rec%xdim(1: geom_rec%nx)
                        read(unit=geom_rec%ychar, fmt=*, iostat=io_error) keyword, dummy_int(1), geom_rec%ydim(1: geom_rec%ny)
                    end if
                    
                case ('REC_CONF')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword
                    
                    call Concatenate (geom_rec%str0(1), words(2:), n_word-1)
                    do i = 2, SIZE(geom_rec%str0) 
                        read(unit=file_in, fmt="(A)", iostat=io_error)  geom_rec%str0(i)
                    end do 
                    
                case ('REC_MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    geom_rec%meshsize = dummy_int(1)
                    geom_rec%is_mesh_define = .TRUE.
                    
                    if (geom_rec%is_mesh_define .AND. geom_rec%is_bc_define)  then
                        call geom_rec%set (ns, geom, mesh, bound, mesh_vtk, model=0)
                    end if
                    
                case ('REC_BC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:4)
                    geom_rec%bc_xmin = dummy_int(1)
                    geom_rec%bc_xmax = dummy_int(2)
                    geom_rec%bc_ymin = dummy_int(3)
                    geom_rec%bc_ymax = dummy_int(4)
                    geom_rec%is_bc_define = .TRUE.
                    
                    if (geom_rec%is_mesh_define .AND. geom_rec%is_bc_define)  then
                        call geom_rec%set (ns, geom, mesh, bound, mesh_vtk, model=0)
                    end if
                    
                ! -------------------------------------------------------------- hex-mesh
                ! --------------------------------------------------------------
                case ('HEX_DIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1)
                    geom_hex%degree = dummy_int(1)
                    geom_hex%ring = dummy_int(2)
                    geom_hex%pitch = dummy_real(1)
                    call geom_hex%alloc ()
                    
                case ('HEX_CONF')
                    if (n_word == 1)  then
                        do i = 1, SIZE(geom_hex%str0) 
                            read(unit=file_in, fmt="(A)", iostat=io_error)  geom_hex%str0(i)
                        end do 
                        
                    else 
                        call Concatenate (geom_hex%str0(1), words(2:), n_word-1)
                        do i = 2, SIZE(geom_hex%str0) 
                            read(unit=file_in, fmt="(A)", iostat=io_error)  geom_hex%str0(i)
                        end do 
                    end if
                
                case ('HEX_MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    geom_hex%meshsize = dummy_int(1)
                    geom_hex%is_mesh_define = .TRUE.
                    
                    if (geom_hex%is_mesh_define .AND. geom_hex%is_bc_define)  then
                        call geom_hex%set (ns, geom, mesh, bound, mesh_vtk, model=0) 
                    end if 
                    
                case ('HEX_BC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    geom_hex%bc_symline = dummy_int(1)
                    geom_hex%bc_outer = dummy_int(2)
                    geom_hex%is_bc_define = .TRUE.
                    
                    if (geom_hex%is_mesh_define .AND. geom_hex%is_bc_define)  then
                        call geom_hex%set (ns, geom, mesh, bound, mesh_vtk, model=0) 
                    end if 

                ! -------------------------------------------------------------- continue 
                ! --------------------------------------------------------------
                case ('BC_AXIAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2)

                case ('LAYER')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:3)
                    ns%state%layer = dummy_int(1)
                    
                    config%n_matRow = dummy_int(2)
                    config%n_matCol = dummy_int(2)
                    config%n_extqRow = dummy_int(3)
                    config%n_extqCol = dummy_int(3)
                    
                case ('REF_AXIAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%layer_bottom  = dummy_int(1)
                    ns%state%layer_top     = dummy_int(2)
                    
                case ('HEIGHT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1: ns%state%layer)
                
                case ('PLANE_MAT')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%matRow_config))   allocate(config%matRow_config(ns%state%zone, config%n_matRow), stat=i_allocate)
                    if (.NOT. allocated(config%matRow_assign))   allocate(config%matRow_assign(ns%state%layer), stat=i_allocate)
                    
                    do i = 1, config%n_matRow
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%zone)
!                        config%matRow_config(:, dummy_int(1)) = dummy_real(1: ns%state%zone)
                    end do
                    
                case ('ASSIGN_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, config%matRow_assign
                
                case ('PLANE_SOURCE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%extqRow_config))   allocate(config%extqRow_config(ns%state%zone, config%n_extqRow), stat=i_allocate)
                    if (.NOT. allocated(config%extqRow_assign))   allocate(config%extqRow_assign(ns%state%layer), stat=i_allocate)
                    
                    do i = 1, config%n_extqRow
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%zone)
!                        config%extqRow_config(:, dummy_int(1)) = dummy_real(1: ns%state%zone)
                    end do
                
                case ('ASSIGN_SOURCE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, config%extqRow_assign
                    
                case ('FA_TYPE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%matCol_config))   allocate(config%matCol_config(config%n_matCol, ns%state%layer), stat=i_allocate)
                    if (.NOT. allocated(config%matCol_assign))   allocate(config%matCol_assign(ns%state%zone), stat=i_allocate)
                    
                    do i = 1, config%n_matCol
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%layer)
!                        config%matCol_config(dummy_int(1), :) = dummy_real(1: ns%state%layer)
                    end do
                    
                case ('ASSIGN_FA')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, config%matCol_assign 
                    
                case ('EXTQ_TYPE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%extqCol_config))   allocate(config%extqCol_config(config%n_extqCol, ns%state%layer), stat=i_allocate)
                    if (.NOT. allocated(config%extqCol_assign))   allocate(config%extqCol_assign(ns%state%zone), stat=i_allocate)
                    
                    do i = 1, config%n_extqCol
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%layer)
!                        config%extqCol_config(dummy_int(1), :) = dummy_real(1: ns%state%layer)
                    end do
                    
                case ('ASSIGN_EXTQ')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, config%extqCol_assign 
                    
                case ('PRINT_MASK')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)
                
                case ('CR_STEP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:4)
                    
                case ('CR_BANK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    cr_bank%state%n_rod = dummy_int(1)
                    cr_bank%state%n_bank = dummy_int(2)
                    
                case ('CR_CONFIGURE')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)
                    nt%flag%is_CR_rod = .TRUE.
                    
                case ('CR_POSITION')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:cr_bank%state%n_bank)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
        
        if (nt%flag%is_pkmodel)  then
            return 
        end if 
        
        if (mesh%meshtype == mesh%MESH_ANSYS)  then 
            tmp_file = TRIM(mesh%xyfile)
            call Is_file_existence (tmp_file, is_input=.TRUE.)
            call Preprocess_ansys (tmp_file, out_file, dummy_int(1))
            if (ns%state%point /= dummy_int(1))  then
                ns%state%point = dummy_int(1)
                call this_warning%set (INFO_LIST_INPUT, 'input point-# is not correct') 
                call this_warning%print (file_out)
            end if 
            
            tmp_file = TRIM(mesh%iifile)
            call Is_file_existence (tmp_file, is_input=.TRUE.)
            call Preprocess_ansys (tmp_file, out_file, dummy_int(2))
            if (ns%state%nodal /= dummy_int(2))  then
                ns%state%nodal = dummy_int(2)
                call this_warning%set (INFO_LIST_INPUT, 'input nodal-# is not correct') 
                call this_warning%print (file_out)
            end if 
        end if 
    
    end subroutine Scan_geometry
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_transient (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        integer  :: it, i
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(TRANSIENT_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('STEP_TR')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2), dummy_int(1) 
                    
                case ('SECTION')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:2)
                    
                case ('INTERVAL')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:3)
                    
                case ('PERTURB_XSEC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_xsec%n_pert = dummy_int(1)
                    
                case ('SIGMA_S', 'SIGMA_T', 'SIGMA_F_NU')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3)
                    
                case ('PERTURB_SOURCE') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_q%n_pert = dummy_int(1)
                    
                case ('INTENSITY')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3)
                    
                case ('PERTURB_CR')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_cr%n_pert = dummy_int(1)
                    
                case ('CR_MOVE')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:4), dummy_char(1)
                
                case ('CR_TRIP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1:2)
                
                case ('PERTURB_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_mat%n_pert = dummy_int(1)
                    
                case ('MAT')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2), dummy_real(2:4), dummy_char(1)
                                    
                case ('PERTURB_FLOW')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_th%n_flow = dummy_int(1)
                    pert_th%var_flow = 4 
                    
                case ('FLOW')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:6) 
                                    
                case ('PERTURB_TM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    pert_th%n_Tm = dummy_int(1)
                    pert_th%var_Tm = 4 
                    
                case ('TM')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:6)
                                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_transient
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_feedback (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(FEEDBACK_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('FDBK')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_log(1:2), dummy_char(1)
                    
                case ('RELAXATION')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:2)
                    
                case ('EFFTF')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1), dummy_real(1)
                    
                case ('DESIGN')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1), dummy_log(1)
                    
                case ('FLOW')
                    backspace(file_in, iostat=io_error)
                
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:ns%state%zone)
                    
                case ('SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:4)
                    
                case ('GAMMA')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1)
                    
                case ('FQ_LATTICE')
                    backspace(file_in, iostat=io_error)
                
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:ns%state%zone)
                    
                case ('CRITERIA')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:2)
                    
                case ('MESH_SIZE')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1:2)
                    nth%nf = dummy_int(1)
                    nth%nc = dummy_int(2)
                    
                case ('COOLANT')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1)
                    
                case ('GEOM_INFO')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1), dummy_log(1)
                    nth%n_assm_geom = dummy_int(1)
                    
                case ('ASSEMBLY')
                    backspace(file_in, iostat=io_error)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:5)
                    
                case ('PROPERTY_TYPE')
                    backspace(file_in, iostat=io_error)

                    read(unit=file_in, fmt=*, iostat=io_error) keyword, dummy_int(1:4)
                    
                case ('TH_CONFIGURE')
                    backspace(file_in, iostat=io_error)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)

                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_feedback
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_perturb (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        
        integer  :: i, j
        integer  :: ia
        integer  :: n_pane
        
        ! use for material configuration
        integer, allocatable  :: plane_mat_assign(:)                            ! planar index per layer
        integer, allocatable  :: plane_mat_configure(:, :)                      ! material ID per zone per planar
        integer, allocatable  :: plane_mat_ID(:)                                ! contains planar index
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(PERTURBATION_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('PERTURB')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    pert_case%n_pert = dummy_int(1)
                    
                case ('PERTURB_INFO')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_log(1), dummy_int(2)
                    n_pane = dummy_int(2)
                    
                case ('PLANE_MAT')
                    backspace(file_in, iostat=io_error)
                    
                    allocate(plane_mat_ID(n_pane), stat=i_allocate)
                    allocate(plane_mat_configure(ns%state%zone, n_pane), stat=i_allocate)
                    allocate(plane_mat_assign(ns%state%layer), stat=i_allocate)
                    
                    do i = 1, n_pane
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, plane_mat_ID(i), (plane_mat_configure(j, i), j=1,ns%state%zone)
                    end do
                    
                case ('ASSIGN_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, plane_mat_assign
                    
                    if (allocated(plane_mat_ID))            deallocate(plane_mat_ID)
                    if (allocated(plane_mat_configure))     deallocate(plane_mat_configure)
                    if (allocated(plane_mat_assign))        deallocate(plane_mat_assign)
                    
                case ('MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_char(1), dummy_int(2:3)
                    
                case ('CR_WORTH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    worth_cr%pos_type = dummy_int(1)
                    
                case ('CR_POS')
                    if (worth_cr%pos_type == 0)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:n_word-2)
                    else if (worth_cr%pos_type == 1)  then 
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:3)
                    end if 
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_perturb
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Scan_pkmodel (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        
        integer  :: i, j
        integer  :: ia
        integer  :: n_pane
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(PKMODEL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('REACTIVITY')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:5), dummy_char(1)
                    
                case ('COEFFICIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:4)
                    
                case ('BLOCK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1:2)
                    
                case ('PKPARAMETER')
                    if (n_word == 1)  then
                        read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1)
                    else
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1)
                    end if 
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1)
                    read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Scan_pkmodel
    
    ! --------------------------------------------------------------------------
    ! read input file
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_casename (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        integer  :: io_error
        
        ! re-get the case name
        backspace(file_in, iostat=io_error)
        
        read(unit=file_in, fmt="(A)", iostat=io_error) aline
        read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_char(1)
        
        ns%flag%case_title = dummy_char(1)
        
        ! print header
        call Print_header_main (file_out, TRIM(ns%flag%case_title))
        call Print_header_main (OUTPUT_UNIT, TRIM(ns%flag%case_title))
        
        call Print_header_timelist (FILES%TIMELIST, TRIM(ns%flag%case_title))
        
    end subroutine Read_casename
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_control (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(CONTROL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('PK_MODEL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1) 
                    nt%flag%is_pkmodel = dummy_log(1)
                
                case ('STEADY')
                    if (n_word == 3)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1)
                        ns%flag%is_eigen = dummy_log(1)
                        ns%flag%rated_power = dummy_real(1)
                        ns%flag%power_frac = 1.0D0
                    else if (n_word == 4)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1:2)
                        ns%flag%is_eigen = dummy_log(1)
                        ns%flag%rated_power = dummy_real(1)
                        ns%flag%power_frac = dummy_real(2)
                    end if 
                    ns%flag%power_level = ns%flag%rated_power * ns%flag%power_frac
                    
                case ('PPM_SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_log(1), dummy_real(1)
                    nt%flag%is_boron_search = dummy_log(1)
                    self_fdbk%is_CB_search = dummy_log(1)
                    self_fdbk%CB_target = dummy_real(1)
                
                case ('ROD_SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1), dummy_int(1)
                    cr_bank%is_search = dummy_log(1)
                    cr_bank%search_value = dummy_real(1)
                    cr_bank%bank_idx = dummy_int(1)
                    
                case ('ADJOINT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1) 
                    nt%flag%is_adjoint = dummy_log(1)

                case ('TRANSIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_char(1)
                    nt%flag%is_transient = dummy_log(1)
                    nt%flag%kinetics_type = dummy_char(1)
            
                case ('LINK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_char(1)
                    ns%flag%is_link = dummy_log(1)
                    ns%flag%link_type = dummy_char(1)
                    
                case ('FEEDBACK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    ns%feedback%is_feedback = dummy_log(1)
                    ns%feedback%is_nested = .TRUE.
                    
                case ('PERTURBATION')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    nt%flag%is_perturb = dummy_log(1)
                    
                case ('DECAY_HEAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    nt%flag%is_decay_heat = dummy_log(1)
                    
                case ('ENERGY_GROUP') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%ng = dummy_int(1)
                    ns%state%scat_order = dummy_int(2)
                    
                case ('DELAYD_GROUP', 'DELAYED_GROUP') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    nt%state%dg = dummy_int(1)
                    
                case ('OUTPUT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:n_word-1)
                    ns%output%is_HDF5 = dummy_log(1)
                    ns%output%is_vtk = dummy_log(2)
                    if (n_word >= 4)  then
                        ns%output%is_log = dummy_log(3)
                    end if 
                
                case ('THREAD')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    ns%misc%nthread = dummy_int(1)
                
                case ('INITVAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:2)
                    iter_init%is_initValin = dummy_log(1)
                    iter_init%is_initValout = dummy_log(2)
                
                case ('MISC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    ns%misc%is_angwise = dummy_log(1)
                
                end select
                    
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word)))
                call this_error%print (file_out)
            end if
        end do
        
    end subroutine Read_control
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_method (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
                        
            read(unit=aline, fmt=*, iostat=io_error) a_word

            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(METHOD_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('QUDARATURE', 'QUADRATURE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1:2)
                    ns%flag%is_60degree = dummy_log(1)
                    ns%flag%n_theta = dummy_int(1)
                    ns%state%sn = dummy_int(2)
                
                case ('STEADY')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:2), dummy_int(1), dummy_log(3)
                    ns%method%is_LW = dummy_log(1)
                    ns%method%is_upscatter_cycle = dummy_log(2) 
                    ns%method%n_upscatter_cycle = dummy_int(1)
                    ns%method%is_Ks = dummy_log(3)
                
                case ('TRANSIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1), dummy_log(1), dummy_char(2)
                    nt%method%scheme = dummy_char(1)
                    nt%method%is_extrapolation = dummy_log(1)
                    nt%method%extrapolation_scheme = dummy_char(2)
                    
                case ('ERROR_TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    criteria_eigen%error_type = dummy_int(1)
                    criteria_fsp%error_type = dummy_int(2)
                
                case ('ERROR_EIGEN')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:4)
                    criteria_eigen%max_inner = dummy_int(1)
                    criteria_eigen%max_outer = dummy_int(2)
                    criteria_eigen%error_inner_flux = dummy_real(1)
                    criteria_eigen%error_outer_flux = dummy_real(2)
                    criteria_eigen%error_fission_rate = dummy_real(3)
                    criteria_eigen%error_eigen = dummy_real(4)
                
                case ('ERROR_FSP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:4)
                    criteria_fsp%max_inner = dummy_int(1)
                    criteria_fsp%max_outer = dummy_int(2)
                    criteria_fsp%error_inner_flux = dummy_real(1)
                    criteria_fsp%error_outer_flux = dummy_real(2)
                    criteria_fsp%error_fission_rate = dummy_real(3)
                    criteria_fsp%error_eigen = dummy_real(4)
                end select

            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_method
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_link (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
                        
            read(unit=aline, fmt=*, iostat=io_error) a_word

            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(LINK_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('DATA')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    self_link%n_parameter = dummy_int(1)
                
                case ('PARAMETER')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1:self_link%n_parameter)
                    call self_link%set_type (dummy_char(1:self_link%n_parameter))
                        
                case ('REFERENCE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:self_link%n_parameter)
                    call self_link%set_value (dummy_real(1:self_link%n_parameter))
                    
                case ('MAX_LIMIT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:self_link%n_parameter)
                    call self_link%set_max (dummy_real(1:self_link%n_parameter))
                    
                case ('MIN_LIMIT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:self_link%n_parameter)
                    call self_link%set_min (dummy_real(1:self_link%n_parameter))
                    
                end select

            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_link

    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_material (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_LINE_LEN)  :: tmp_line
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        integer  :: unit_xsec                                                   ! file unit for external xsec
        
        character(MAX_WORD_LEN)  :: data_in
        character(MAX_WORD_LEN)  :: data_out
        
        character(len=MAX_WORD_LEN)  :: tmp_file
        integer  :: i, im, ic, ig, iig, ie, ip
        integer  :: im_start, im_end                                            ! use for one group xsec adapt to material range
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(MATERIAL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('SCALE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%mat = dummy_int(1)
                    ns%state%source = dummy_int(2)
                    
                case ('IMPORT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1), dummy_char(1)
                    ns%misc%is_import = dummy_log(1)
                    ns%misc%burnstep = dummy_int(1)
                    ns%misc%file = dummy_char(1)
                    
                case ('TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    ns%misc%is_external = dummy_log(1)
                    
                case ('MAT')
                    backspace(file_in, iostat=io_error)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3) 
                    ! hold for material ID
                    im = dummy_int(1)
                    mat_info%libs(im)%ID = dummy_int(1)
                    mat_info%libs(im)%ID_sub = dummy_int(2)
                    mat_info%libs(im)%ID_add = dummy_int(3)
                    
                    if ((mat_info%libs(im)%ID == mat_info%libs(im)%ID_add) .AND. (mat_info%libs(im)%ID_sub /= mat_info%libs(im)%ID_add))  then
                        mat_info%libs(im)%is_CR = .TRUE.
                    else
                        mat_info%libs(im)%is_CR = .FALSE.
                    end if
                    
                    if (.NOT. ns%misc%is_external)  then
                        call xsec_read_known (file_in, xsec_inp%mats(im))
                    
                    ! xsec from external file
                    else
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        mat_info%libs(im)%file = TRIM(keyword) 
                        call xsec_read_known (DIR_XSEC, keyword, a_xsec=xsec_inp%mats(im))
                    end if
                    
                    ! set sigma_f & sigma_f_kappa
                    call xsec_inp%mats(im)%set (is_nu=.FALSE., is_kappa=.FALSE.)
                
                case ('MAT_D')
                    backspace(file_in, iostat=io_error)

                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ! hold for started material & end material ID
                    im_start = dummy_int(1)
                    im_end = dummy_int(2)
                        
                    if (.NOT. ns%misc%is_external)  then
                        call xsec_read_known (file_in, param_inp%mats(im_start))
                        
                    ! xsec from external file
                    else
                        tmp_file = TRIM(DIR_XSEC) // TRIM(ADJUSTL(mat_info%libs(im_start)%file))
                        call xsec_read_known (DIR_XSEC, mat_info%libs(im_start)%file, a_param=param_inp%mats(im_start))
                    end if
                           
                    ! assign to other material
                    if (im_start < im_end)  then
                        do im = im_start+1, im_end
                            param_inp%mats(im) = param_inp%mats(im_start)
                        end do
                    end if
                    
                case ('Q_EXTERNAL')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    
                    if (.NOT. ns%misc%is_external)  then
                        ! hold for external source ID
                        ie = dummy_int(1)
                        read(unit=file_in, fmt=*, iostat=io_error)  (Q_ext%kinds(ie)%intensity(ig), ig=1, ns%state%ng)
                        
                    ! read from external file
                    else 
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        tmp_file = TRIM(DIR_XSEC)//TRIM(keyword)
                        open(newunit=unit_xsec, file=tmp_file, status='old', access='sequential', form='formatted', action='read', iostat=io_error)
                        
                        do i = 1, ns%state%source
                            read(unit=unit_xsec, fmt=*, iostat=io_error)  ie
                            read(unit=unit_xsec, fmt=*, iostat=io_error)  (Q_ext%kinds(ie)%intensity(ig), ig=1, ns%state%ng)
                        end do
                        
                        close(unit=unit_xsec, status='keep', iostat=io_error)
                    end if
                    
                case ('DETINFO')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    det%nxs = dummy_int(1)
                    det%ndet = dummy_int(2)
                    
                case ('DETXS')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    
                    if (.NOT. ns%misc%is_external)  then
                        im = dummy_int(1)
                        read(unit=file_in, fmt=*, iostat=io_error) det%detxs(:, im)
                    else
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword
                        tmp_file = TRIM(DIR_XSEC)//TRIM(keyword)
                        open(newunit=unit_xsec, file=tmp_file, status='old', access='sequential', form='formatted', action='read', iostat=io_error)
                        
                        do i = 1, det%nxs
                            read(unit=unit_xsec, fmt=*, iostat=io_error)  im
                            read(unit=unit_xsec, fmt=*, iostat=io_error)  det%detxs(:, im)
                        end do
                        
                        close(unit=unit_xsec, status='keep', iostat=io_error)
                    end if 
                    
                case ('DETPOINT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2:4)
                    im = dummy_int(1)
                    det%detpoint(im)%const = dummy_real(1)
                    det%detpoint(im)%idx = dummy_int(1)
                    det%detpoint(im)%ixs = dummy_int(2)
                    det%detpoint(im)%iz = dummy_int(3)
                    det%detpoint(im)%ia = dummy_int(4)
                    if (n_word >= 7)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2:5)
                        det%detpoint(im)%ir = dummy_int(5)
                    end if 
                    if (n_word >= 8)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2:6)
                        det%detpoint(im)%is = dummy_int(6)
                    end if 
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_material
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_geometry (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_LINE_LEN)  :: tmp_file, out_file
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        integer  :: i, j
        integer  :: ia
        integer  :: unit1_, unit2_ 
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(GEOMETRY_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    mesh%meshtype = dummy_int(1)
                    
                ! -------------------------------------------------------------- ansys-mesh
                ! --------------------------------------------------------------
                case ('TYPE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1)
                    ns%flag%is_bevel_edge = dummy_log(1)
                
                case ('ANSYS_FILE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_char(1:2)
                    call Lower_case (dummy_char(1))
                    call Lower_case (dummy_char(2))
                    mesh%xyfile = dummy_char(1)
                    mesh%iifile = dummy_char(2)
                
                case ('XY_EXPAND')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2)
                    geom%dx0 = dummy_real(1)
                    geom%dx1 = dummy_real(2)
                 
                case ('SCALE') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:4)
!                    ns%state%point = dummy_int(1)
!                    ns%state%nodal =  dummy_int(2)
                    ns%state%segment = dummy_int(3)
                    ns%state%zone =  dummy_int(4)
                    
                case ('BC_RADIAL')
                    backspace(file_in, iostat=io_error)
                    
                    ! ----------------------------------------------------------
                    ! NOTE:
                    !   set information, because some parameter dependent on this, 
                    !   especially control rod.
                    ! ----------------------------------------------------------
                    
                    ! geometry input
                    tmp_file = TRIM(mesh%xyfile)
                    call Preprocess_ansys (tmp_file, out_file, dummy_int(1))
                    call Open_file (out_file, .TRUE., unit1_)
                    tmp_file = TRIM(mesh%iifile)
                    call Preprocess_ansys (tmp_file, out_file, dummy_int(2))
                    call Open_file (out_file, .TRUE., unit2_)
                    
                    do i = 1, ns%state%point
                        read(unit1_, fmt=*)  dummy_int(1), dummy_real(1:2)
                        geom%coordinate(1, i) = dummy_real(1)  
                        geom%coordinate(2, i) = dummy_real(2)  
                    end do
                    geom%coordinate = geom%coordinate * (geom%dx1 / geom%dx0)
                    
                    do i = 1, ns%state%nodal
                        read(unit2_, fmt=*)  dummy_int(1:9)
                        mesh%zone(i) = dummy_int(3)
                        mesh%point(1, i) = dummy_int(7)
                        mesh%point(2, i) = dummy_int(8)
                        mesh%point(3, i) = dummy_int(9)
                    end do
                    
                    close(unit=unit1_, status='delete', iostat=io_error)
                    close(unit=unit2_, status='delete', iostat=io_error)
                    
                    call Set_section_1 ()
                                    
                    bound%nodal = bound%VACCUM
                    do i = 1, ns%state%segment
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1)
                        bound%radial(i) = dummy_real(1)
                        call bound%set (mesh, geom, dummy_int(1:3), i)
                    end do
                    
                    call Set_section_2 ()
                    
                ! -------------------------------------------------------------- rec-mesh
                ! --------------------------------------------------------------
                case ('REC_XDIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: geom_rec%nx)
                    geom_rec%nx = dummy_int(1)
                    geom_rec%xdim = dummy_real(1: geom_rec%nx)
                    geom_rec%is_mesh_define = .FALSE.
                    geom_rec%is_bc_define = .FALSE.
                    
                case ('REC_YDIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: geom_rec%ny)
                    geom_rec%ny = dummy_int(1)
                    geom_rec%ydim = dummy_real(1: geom_rec%ny)
                    
                case ('REC_CONF')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword
                    
                    call Concatenate (geom_rec%str0(1), words(2:), n_word-1)
                    do i = 2, SIZE(geom_rec%str0) 
                        read(unit=file_in, fmt="(A)", iostat=io_error)  geom_rec%str0(i)
                    end do 
                    
                case ('REC_MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    geom_rec%meshsize = dummy_int(1)
                    geom_rec%is_mesh_define = .TRUE.
                    
                    if (geom_rec%is_mesh_define .AND. geom_rec%is_bc_define)  then
                        call geom_rec%set (ns, geom, mesh, bound, mesh_vtk, model=1)
                        call Set_section_1 ()
                        call Set_section_2 ()
                    end if
                    
                case ('REC_BC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:4)
                    geom_rec%bc_xmin = dummy_int(1)
                    geom_rec%bc_xmax = dummy_int(2)
                    geom_rec%bc_ymin = dummy_int(3)
                    geom_rec%bc_ymax = dummy_int(4)
                    geom_rec%is_bc_define = .TRUE.
                    
                    if (geom_rec%is_mesh_define .AND. geom_rec%is_bc_define)  then
                        call geom_rec%set (ns, geom, mesh, bound, mesh_vtk, model=1)
                        call Set_section_1 ()
                        call Set_section_2 ()
                    end if
                
                ! -------------------------------------------------------------- hex-mesh
                ! --------------------------------------------------------------
                case ('HEX_DIM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1)
                    geom_hex%degree = dummy_int(1)
                    geom_hex%ring = dummy_int(2)
                    geom_hex%pitch = dummy_real(1)
                    call geom_hex%alloc ()
                
                case ('HEX_CONF')
                    if (n_word == 1)  then
                        do i = 1, SIZE(geom_hex%str0) 
                            read(unit=file_in, fmt="(A)", iostat=io_error)  geom_hex%str0(i)
                        end do 
                        
                    else 
                        call Concatenate (geom_hex%str0(1), words(2:), n_word-1)
                        do i = 2, SIZE(geom_hex%str0) 
                            read(unit=file_in, fmt="(A)", iostat=io_error)  geom_hex%str0(i)
                        end do 
                    end if 
                    
                case ('HEX_MESH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    geom_hex%meshsize = dummy_int(1)
                    geom_hex%is_mesh_define = .TRUE.
                    
                    if (geom_hex%is_mesh_define .AND. geom_hex%is_bc_define)  then
                        call geom_hex%set (ns, geom, mesh, bound, mesh_vtk, model=1) 
                        call Set_section_1 ()
                        call Set_section_2 ()
                    end if 
                    
                case ('HEX_BC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    geom_hex%bc_symline = dummy_int(1)
                    geom_hex%bc_outer = dummy_int(2)
                    geom_hex%is_bc_define = .TRUE.
                    
                    if (geom_hex%is_mesh_define .AND. geom_hex%is_bc_define)  then
                        call geom_hex%set (ns, geom, mesh, bound, mesh_vtk, model=1) 
                        call Set_section_1 ()
                        call Set_section_2 ()
                    end if 
                
                ! -------------------------------------------------------------- continue 
                ! --------------------------------------------------------------
                case ('BC_AXIAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2)
                    bound%axial(1) = dummy_real(1)                              ! bottom boudary
                    bound%axial(2) = dummy_real(2)                              ! top boundary

                case ('LAYER')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:3)
                    ns%state%layer = dummy_int(1)
                    config%n_matRow = dummy_int(2)
                    config%n_matCol = dummy_int(2)
                    config%n_extqRow = dummy_int(3)
                    config%n_extqCol = dummy_int(3)
                 
                case ('REF_AXIAL')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    ns%state%layer_bottom  = dummy_int(1)
                    ns%state%layer_top = dummy_int(2)
                    cr_bank%state%n_bottom = dummy_int(1)
                    cr_bank%state%n_top = dummy_int(2)
                    
                case ('HEIGHT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1: ns%state%layer)
                    geom%height = dummy_real(1: ns%state%layer)

                case ('PLANE_MAT')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%matRow_config))   allocate(config%matRow_config(ns%state%zone, config%n_matRow), stat=i_allocate)
                    if (.NOT. allocated(config%matRow_assign))   allocate(config%matRow_assign(ns%state%layer), stat=i_allocate)
                    
                    do i = 1, config%n_matRow
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%zone)
                        config%matRow_config(:, dummy_int(1)) = dummy_real(1: ns%state%zone)
                    end do

                case ('ASSIGN_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, config%matRow_assign
                    
                    call mat_info%set_row (config%matRow_assign, config%matRow_config)

                case ('PLANE_SOURCE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%extqRow_config))   allocate(config%extqRow_config(ns%state%zone, config%n_extqRow), stat=i_allocate)
                    if (.NOT. allocated(config%extqRow_assign))   allocate(config%extqRow_assign(ns%state%layer), stat=i_allocate)
                    
                    do i = 1, config%n_extqRow
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%zone)
                        config%extqRow_config(:, dummy_int(1)) = dummy_real(1: ns%state%zone)
                    end do
                    
                case ('ASSIGN_SOURCE')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, config%extqRow_assign
                    
                    call Q_ext%set_row (config%extqRow_assign, config%extqRow_config)
                    
                case ('FA_TYPE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%matCol_config))   allocate(config%matCol_config(config%n_matCol, ns%state%layer), stat=i_allocate)
                    if (.NOT. allocated(config%matCol_assign))   allocate(config%matCol_assign(ns%state%zone), stat=i_allocate)
                    
                    do i = 1, config%n_matCol
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%layer)
                        config%matCol_config(dummy_int(1), :) = dummy_real(1: ns%state%layer)
                    end do
                    
                case ('ASSIGN_FA')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, config%matCol_assign 
                    
                    call mat_info%set_col (config%matCol_assign, config%matCol_config)
                    
                case ('EXTQ_TYPE')
                    backspace(file_in, iostat=io_error)
                    
                    if (.NOT. allocated(config%extqCol_config))   allocate(config%extqCol_config(config%n_extqCol, ns%state%layer), stat=i_allocate)
                    if (.NOT. allocated(config%extqCol_assign))   allocate(config%extqCol_assign(ns%state%zone), stat=i_allocate)
                    
                    do i = 1, config%n_extqCol
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1: ns%state%layer)
                        config%extqCol_config(dummy_int(1), :) = dummy_real(1: ns%state%layer)
                    end do
                    
                case ('ASSIGN_EXTQ')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, config%extqCol_assign 
                    
                    call Q_ext%set_col (config%extqCol_assign, config%extqCol_config)
                
                case ('PRINT_MASK')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)
                    mesh_vtk%mapping = dummy_int(1:ns%state%zone)

                case ('CR_STEP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:4)
                    cr_bank%state%step_size = dummy_real(1)
                    cr_bank%state%cr_gap = dummy_real(2)
                    cr_bank%state%min_step = dummy_real(3)
                    cr_bank%state%max_step = dummy_real(4)
                    call cr_bank%alloc (ns, geom)
                    
                case ('CR_BANK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:2)
                    cr_bank%state%n_rod = dummy_int(1)
                    cr_bank%state%n_bank = dummy_int(2)
                    
                case ('CR_CONFIGURE')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)
                    nt%flag%is_CR_rod = .TRUE.
                    call cr_bank%alloc (ns, geom) 
                    call cr_bank%set (geom, mat_info, xsec, param, dummy_int(1:ns%state%zone))
                
                case ('CR_POSITION')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:cr_bank%state%n_bank)
                    call cr_bank%init (dummy_real(1:cr_bank%state%n_bank))
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
        
        ! NOTE--import from output of LAVENDER code 
        if (ns%misc%is_import )  then
            ns%state%mat = ns%state%layer * ns%state%zone
            
            call xsec_inp%clean ()
            call xsec_inp%alloc ()
            call param_inp%clean ()
            call param_inp%alloc ()
            
            if (nt%flag%is_transient)  then
                call xsec_read_known (DIR_XSEC, ns%misc%file, ns%misc%burnstep, xsec_inp, param_inp)
            else
                call xsec_read_known (DIR_XSEC, ns%misc%file, ns%misc%burnstep, xsec_inp)
            end if
            call mat_info%fix ()
        end if
        
        call geom_rec%clean ()
        call geom_hex%clean ()

    end subroutine Read_geometry
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_transient (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        integer  :: it, i
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(TRANSIENT_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('STEP_TR')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:2), dummy_int(1)
                    time_step%start_time = dummy_real(1)
                    time_step%end_time = dummy_real(2)
                    time_step%n_section = dummy_int(1)
                    time_step%n_interval = dummy_int(1)
                    
                case ('SECTION')
                    backspace(file_in, iostat=io_error)
                    call time_step%alloc ()
                    do it = 1, time_step%n_section
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:2)
                        time_step%sections(it)%point = dummy_real(1)
                        time_step%sections(it)%length = dummy_real(2)
                    end do
                    
                    call time_step%set_section ()
                    
                case ('INTERVAL')
                    backspace(file_in, iostat=io_error)
                    call time_step%alloc ()
                    do it = 1, time_step%n_interval 
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:3)
                        time_step%intervals(it)%tend = dummy_real(1)
                        time_step%intervals(it)%dmacro = dummy_real(2)
                        time_step%intervals(it)%dmiddle = dummy_real(3)
                    end do
                    
                    call time_step%set_interval ()
                                        
                case ('PERTURB_XSEC')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_xsec = dummy_log(1)
                    pert_xsec%n_pert = dummy_int(1)
                    
                case ('SIGMA_S', 'SIGMA_T', 'SIGMA_F_NU')
                    backspace(file_in, iostat=io_error)
                    
                    do i = 1, pert_xsec%n_pert
                        read(unit=file_in, fmt="(A)", iostat=io_error)  aline
                        read(unit=aline, fmt=*, iostat=io_error)  keyword
                        
                        select case(TRIM(keyword))
                        case ('SIGMA_S')
                            read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3), dummy_int(4:6)
                            pert_xsec%perts(i)%xs_type = keyword
                            
                            pert_xsec%perts(i)%matID = dummy_int(1)
                            pert_xsec%perts(i)%ng_start = dummy_int(2)
                            pert_xsec%perts(i)%ng_end = dummy_int(3)
                            pert_xsec%perts(i)%time_start = dummy_real(1)
                            pert_xsec%perts(i)%time_end = dummy_real(2)
                            pert_xsec%perts(i)%type = dummy_char(1)
                            pert_xsec%perts(i)%percent = dummy_real(3)
                            
                            pert_xsec%perts(i)%to_ng_start = dummy_int(4)
                            pert_xsec%perts(i)%to_ng_end = dummy_int(5)
                            pert_xsec%perts(i)%scat_index = dummy_int(6)
                            
                        case ('SIGMA_T')    
                            read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3)
                            pert_xsec%perts(i)%xs_type = keyword
                            
                            pert_xsec%perts(i)%matID = dummy_int(1)
                            pert_xsec%perts(i)%ng_start = dummy_int(2)
                            pert_xsec%perts(i)%ng_end = dummy_int(3)
                            pert_xsec%perts(i)%time_start = dummy_real(1)
                            pert_xsec%perts(i)%time_end = dummy_real(2)
                            pert_xsec%perts(i)%type = dummy_char(1)
                            pert_xsec%perts(i)%percent = dummy_real(3)
                        
                        case ('SIGMA_F_NU')
                            read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3)
                            pert_xsec%perts(i)%xs_type = keyword
                            
                            pert_xsec%perts(i)%matID = dummy_int(1)
                            pert_xsec%perts(i)%ng_start = dummy_int(2)
                            pert_xsec%perts(i)%ng_end = dummy_int(3)
                            pert_xsec%perts(i)%time_start = dummy_real(1)
                            pert_xsec%perts(i)%time_end = dummy_real(2)
                            pert_xsec%perts(i)%type = dummy_char(1)
                            pert_xsec%perts(i)%percent = dummy_real(3)
                            
                        end select
                    end do
                    
                case ('PERTURB_SOURCE') 
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_source = dummy_log(1)
                    pert_q%n_pert = dummy_int(1)
                    
                case ('INTENSITY')
                    backspace(file_in, iostat=io_error)
                    
                    do i = 1, pert_q%n_pert
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:2), dummy_char(1), dummy_real(3)
                        pert_q%perts(i)%kindID = dummy_int(1)
                        
                        pert_q%perts(i)%ng_start = dummy_int(2)
                        pert_q%perts(i)%ng_end = dummy_int(3)
                        pert_q%perts(i)%time_start = dummy_real(1)
                        pert_q%perts(i)%time_end = dummy_real(2)

                        pert_q%perts(i)%type = dummy_char(1)
                        pert_q%perts(i)%percent = dummy_real(3)
                    end do
                    
                case ('PERTURB_CR')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_CR_move = dummy_log(1)
                    pert_cr%n_pert = dummy_int(1)
                    
                case ('CR_MOVE')
                    backspace(file_in, iostat=io_error)
                    do i = 1, pert_cr%n_pert
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:4), dummy_char(1)
                        pert_cr%perts(i)%bankID = dummy_int(1)
                        pert_cr%perts(i)%step_start = dummy_real(1)             ! useless 
                        pert_cr%perts(i)%step_end = dummy_real(2)
                        pert_cr%perts(i)%time_start = dummy_real(3)
                        pert_cr%perts(i)%time_end = dummy_real(4)
                        pert_cr%perts(i)%type = dummy_char(1)
                    end do
                    
                case ('CR_TRIP')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1:2)
                    pert_cr%is_trip = dummy_log(1)
                    pert_cr%trip%pstart = dummy_real(1)
                    pert_cr%trip%tdelay = dummy_real(2)
                
                case ('PERTURB_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_mat = dummy_log(1)
                    pert_mat%n_pert = dummy_int(1)
                    
                case ('MAT')
                    backspace(file_in, iostat=io_error)
                    do i = 1, pert_mat%n_pert
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1), dummy_int(2), dummy_real(2:4), dummy_char(1)
                        pert_mat%perts(i)%mat_beg = dummy_int(1)
                        pert_mat%perts(i)%per_beg = dummy_real(1)
                        pert_mat%perts(i)%mat_end = dummy_int(2)
                        pert_mat%perts(i)%per_end = dummy_real(2)
                        pert_mat%perts(i)%time_start = dummy_real(3)
                        pert_mat%perts(i)%time_end = dummy_real(4)
                        pert_mat%perts(i)%type = dummy_char(1)
                    end do
                    
                case ('PERTURB_FLOW')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_flow = dummy_log(1)
                    pert_th%n_flow = dummy_int(1)
                    pert_th%var_flow = 4 
                    
                case ('FLOW')
                    backspace(file_in, iostat=io_error)
                    
                    do i = 1, pert_th%n_flow
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:2), dummy_real(1:6) 
                        pert_th%flow_perts(i)%channelID = dummy_int(1)
                        pert_th%flow_perts(i)%type = dummy_int(2)
                        pert_th%flow_perts(i)%natural = dummy_real(1)
                        pert_th%flow_perts(i)%time_start = dummy_real(2)
                        pert_th%flow_perts(i)%variables = dummy_real(3:6)
                    end do 
                                    
                case ('PERTURB_TM')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_int(1)
                    nt%perturb%is_Tm = dummy_log(1)
                    pert_th%n_Tm = dummy_int(1)
                    pert_th%var_Tm = 4 
                    
                case ('TM')
                    backspace(file_in, iostat=io_error)
                    
                    do i = 1, pert_th%n_Tm                    
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:6)
                        pert_th%Tm_perts(i)%type = dummy_int(1)
                        pert_th%Tm_perts(i)%time_start = dummy_real(1)
                        pert_th%Tm_perts(i)%time_end = dummy_real(2)
                        pert_th%Tm_perts(i)%variables = dummy_real(3:6)
                    end do 
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_transient
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_feedback (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error
        
        logical  :: is_hexagonal
        logical  :: is_square
        integer  :: i_type
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(FEEDBACK_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('FDBK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1:2), dummy_char(1)
                    ns%feedback%is_inner = dummy_log(1)
                    ns%feedback%is_model = dummy_log(2)
                    ns%feedback%model_name = dummy_char(1)
                    self_fdbk%is_th_inner = dummy_log(1)
                    self_fdbk%is_model = dummy_log(2)
                    self_fdbk%model_name = dummy_char(1)
                    
                case ('RELAXATION')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:2)
                    self_fdbk%relax_TH = dummy_real(1)
                    self_fdbk%relax_CB = dummy_real(2)
                    
                case ('EFFTF')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1), dummy_real(1)
                    design%tftype = dummy_int(1)
                    design%tfweight = dummy_real(1)
                    
                case ('DESIGN')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1), dummy_log(1)
                    design%tcoolin = dummy_real(1)
                    design%init_tcoolin = dummy_real(1)
                    design%is_search = dummy_log(1)
                    
                case ('FLOW')
                    backspace(file_in, iostat=io_error)
                
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:ns%state%zone)
                    design%assembly_flow = dummy_real(1:ns%state%zone)
                    
                case ('SEARCH')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:4)
                    design%tmincoolout = dummy_real(1)
                    design%tmaxcoolout = dummy_real(2)
                    design%tmaxcladsurf = dummy_real(3)
                    design%max_velocity = dummy_real(4)
                    
                case ('GAMMA')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1)
                    th_power%gamma = dummy_real(1)
                    
                case ('FQ_LATTICE')
                    backspace(file_in, iostat=io_error)
                
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_real(1:ns%state%zone)
                    th_power%fq_lattice = dummy_real(1:ns%state%zone)

                case ('CRITERIA')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_real(1:2)
                    self_fdbk%Tm%limit = dummy_real(1)
                    self_fdbk%Tf%limit = dummy_real(2)
                    
                case ('MESH_SIZE')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1:2)
                    nth%nf = dummy_int(1)
                    nth%nc = dummy_int(2)
                    
                case ('COOLANT')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1)
                    geom_th%coolant_type = dummy_int(1)

                case ('GEOM_INFO')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1), dummy_log(1)
                    nth%n_assm_geom = dummy_int(1)
                    is_hexagonal = dummy_log(1)
                    
                case ('ASSEMBLY')
                    backspace(file_in, iostat=io_error)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:3), dummy_real(1:5)
                    i_type = dummy_int(1)
                    
                    geom_assm(i_type)%n_pin = dummy_int(2)
                    geom_assm(i_type)%n_fuelpin = dummy_int(3)
                    geom_assm(i_type)%pitch = dummy_real(1)
                    geom_assm(i_type)%rod = dummy_real(2)
                    geom_assm(i_type)%cladth = dummy_real(3)
                    geom_assm(i_type)%bond = dummy_real(4)
                    geom_assm(i_type)%hole = dummy_real(5)
					
					reInputdata%npin = dummy_int(2)
                    reInputdata%nFuelPin = dummy_int(3)
                    reInputdata%pd = dummy_real(1)
                    reInputdata%xf = dummy_real(2)*0.01D0
                    reInputdata%xs = dummy_real(3)*0.01D0
                    reInputdata%xg = dummy_real(4)*0.01D0
                    
                    geom_assm(i_type)%is_hexagonal = is_hexagonal
                    call geom_assm(i_type)%set (nth)
                    
                case ('PROPERTY_TYPE')
                    read(unit=aline, fmt=*, iostat=io_error) keyword, dummy_int(1:4)    
                    i_type = dummy_int(1)
                    
                    geom_th%cladding_type(i_type) = dummy_int(2)
                    geom_th%gap_type(i_type) = dummy_int(3)
                    geom_th%fuel_type(i_type) = dummy_int(4)
                    
                case ('TH_CONFIGURE')
                    backspace(file_in, iostat=io_error)
                
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1:ns%state%zone)
                    
                    geom_th%geom_type(1:ns%state%zone) = dummy_int(1:ns%state%zone)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_feedback
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_perturb (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        integer  :: i_pert, i_mat
        integer  :: n_plane
        
        integer  :: i, j
        integer  :: ia
        
        ! use for material configuration
        integer, allocatable  :: plane_mat_assign(:)                            ! planar index per layer
        integer, allocatable  :: plane_mat_configure(:, :)                      ! material ID per zone per planar
        integer, allocatable  :: plane_mat_ID(:)                                ! contains planar index
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(PERTURBATION_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('PERTURB')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    pert_case%n_pert = dummy_int(1)
                    
                case ('PERTURB_INFO')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_log(1), dummy_int(2)
                    i_pert = dummy_int(1)
                    n_plane = dummy_int(2)
                    pert_case%perts(i_pert)%caseID = dummy_int(1)
                    pert_case%perts(i_pert)%is_new_configure = dummy_log(1)
                    
                case ('PLANE_MAT')
                    backspace(file_in, iostat=io_error)
                    allocate(plane_mat_ID(n_plane), stat=i_allocate)
                    allocate(plane_mat_configure(ns%state%zone, n_plane), stat=i_allocate)
                    allocate(plane_mat_assign(ns%state%layer), stat=i_allocate)
                    do i = 1, n_plane
                        read(unit=file_in, fmt=*, iostat=io_error)  keyword, plane_mat_ID(i), (plane_mat_configure(j, i), j=1,ns%state%zone)
                    end do
                    
                case ('ASSIGN_MAT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, plane_mat_assign
                    call pert_case%set_conf (i_pert, plane_mat_assign, plane_mat_ID, plane_mat_configure)
                    if (allocated(plane_mat_ID))            deallocate(plane_mat_ID)
                    if (allocated(plane_mat_configure))     deallocate(plane_mat_configure)
                    if (allocated(plane_mat_assign))        deallocate(plane_mat_assign)
                    
                case ('MAT')
                    backspace(file_in, iostat=io_error)
                    read(unit=file_in, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_char(1), dummy_int(2:3)
                    call pert_case%set_xsec(i_pert, dummy_char(1), dummy_int(1:3))
                    
                case ('CR_WORTH')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1)
                    worth_cr%pos_type = dummy_int(1)
                    
                case ('CR_POS')
                    if (worth_cr%pos_type == 0)  then
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:n_word-2)
                        call worth_cr%alloc(n_word-2)
                        worth_cr%pos_input = dummy_real(1:n_word-2)
                    else if (worth_cr%pos_type == 1)  then 
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_int(1), dummy_real(1:3)
                        worth_cr%pos_beg = dummy_real(1)
                        worth_cr%pos_end = dummy_real(2)
                        worth_cr%pos_step = dummy_real(3)
                    end if 
                    worth_cr%bank = dummy_int(1)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_perturb
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Read_pkmodel (file_in, file_out)
        
        integer, intent(in)            :: file_in
        integer, intent(in), optional  :: file_out
        
        character(MAX_WORD_LEN)  :: a_word
        integer  :: io_error, i_allocate
        
        integer  :: i, j
        integer  :: ia
        integer  :: n_pane
        
        do 
            read(unit=file_in, fmt="(A)", iostat=io_error)  aline
            
            ! end of file, then exit
            if (io_error == IOSTAT_END)  exit
            
            read(unit=aline, fmt=*, iostat=io_error) a_word
            
            ! reaching to next section, backspace and exit
            if (Is_keyword(INP_SECTION, a_word))  then
                backspace(file_in, iostat=io_error)
                exit
            end if
            
            if (Is_keyword(PKMODEL_CARD, a_word))  then
                call Split_string (aline, words, n_word)
                
                select case (TRIM(ADJUSTL(a_word)))
                case ('REACTIVITY')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:5), dummy_char(1)
                    pk_rho%rho_start = dummy_real(1)
                    pk_rho%rho_end = dummy_real(2)
                    pk_rho%first = dummy_real(3)
                    pk_rho%second = dummy_real(4)
                    pk_rho%third = dummy_real(5)
                    pk_rho%type = dummy_char(1)
                    
                case ('COEFFICIENT')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1:4)
                    call pk_rho%set_coeff (dummy_real(1:4))
                    
                case ('BLOCK')
                    read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_log(1), dummy_real(1:2)
                    design%is_blockage = dummy_log(1)
                    design%block_time = dummy_real(1)
                    design%percentage = dummy_real(2)
                    
                case ('PKPARAMETER')
                    if (n_word == 1)  then
                        read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1)
                    else
                        read(unit=aline, fmt=*, iostat=io_error)  keyword, dummy_real(1)
                    end if 
                    pk_param%generation_time = dummy_real(1)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1: nt%state%dg)
                    pk_param%partial_beta = dummy_real(1: nt%state%dg)
                    
                    read(unit=file_in, fmt=*, iostat=io_error)  dummy_real(1: nt%state%dg)
                    pk_param%partial_lambda = dummy_real(1: nt%state%dg)
                    
                end select
                
            else  
                call this_error%set (INFO_LIST_INPUT, 'key word input error: '//TRIM(ADJUSTL(a_word))) 
                call this_error%print (file_out)
            end if
        end do
    
    end subroutine Read_pkmodel
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_section_1 ()
        
        integer  :: i, j 
    
        ! INPORTANT NOTE: re-range mesh point, then set meshing information
        ! get area per nodal or per zone
        call geom%get_area (mesh)
        call geom%get_zone_area (mesh)
        
        ! get encoding information for mesh
        call mesh%set ()
        call mesh_vtk%set (mesh, geom)
        
!        do i = 1, SIZE(geom%area)
!            write(121, fmt="(1x, I3, TR3, F8.4)")  i, geom%area(i)
!        end do 
!        
!        do i = 1, SIZE(geom%zone_area)
!            write(122, fmt="(1x, I3, TR3, F8.4)")  i, geom%zone_area(i)
!        end do 
            
    end subroutine Set_section_1 
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_section_2 ()
        
        integer  :: i, j 
        
!        call Print_vtk_geometry_2D (150, 'test')
        
        call quad%set_weight ()
        call quad%set_project (mesh, geom)
        call quad%set_sym (ns)
        call quad%set_order (ns)
        
        call bound%inbc (mesh)
        call sweep%set (mesh, geom, bound, quad)
        
        if (ns%flag%is_link)  then
            call self_fdbk%set (self_link)
        end if
    
    end subroutine Set_section_2
    
end module input_plain
