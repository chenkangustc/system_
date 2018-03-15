module hdf5_interpolation
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use timer_header,               only : TotalTimeCounter, CPUTimeCounter
    use material_header,            only : Material, CrossSection, CrossSectionInput, cross_section_info_tp,        &
                                        &   KineticsParameter, KineticsParameterInput, kinetics_parameter_info_tp
    
    use global
    
    use stastics
    use input_xsec 
    use interpolation,              only : array_interpolation, two_variables_interpolation
    
    implicit none 
    private
    public  :: Check_interpolation
    
    type(TotalTimeCounter)  :: a_total
    type(CPUTimeCounter)    :: a_cpu
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Check_interpolation ()
        
        integer  :: ia, iz, im 
        integer  :: iia, iiz 
        integer  :: idx
        character(len=MAX_WORD_LEN)  :: file_in
        
        integer  :: index_00
        integer  :: index_01
        integer  :: index_10
        integer  :: index_11
        integer  :: index_xx
        type(cross_section_info_tp)  :: xsec_00
        type(cross_section_info_tp)  :: xsec_01
        type(cross_section_info_tp)  :: xsec_10
        type(cross_section_info_tp)  :: xsec_11
        type(cross_section_info_tp)  :: xsec_a
        type(cross_section_info_tp)  :: xsec_b
        type(cross_section_info_tp)  :: xsec_xx
        type(cross_section_info_tp)  :: xsec_tmp
        type(kinetics_parameter_info_tp)  :: param_00
        type(kinetics_parameter_info_tp)  :: param_01
        type(kinetics_parameter_info_tp)  :: param_10
        type(kinetics_parameter_info_tp)  :: param_11
        type(kinetics_parameter_info_tp)  :: param_a
        type(kinetics_parameter_info_tp)  :: param_b
        type(kinetics_parameter_info_tp)  :: param_xx
        type(kinetics_parameter_info_tp)  :: param_tmp
        
        type(cross_section_info_tp)        :: a_xsec
        type(kinetics_parameter_info_tp)   :: a_param
        
        real(KREAL)  :: x_min, x, x_max
        real(KREAL)  :: y_min, y, y_max
        
        call a_total%start ()
        call a_cpu%start ()
        
        call xsec_00%alloc (); call param_00%alloc ();
        call xsec_01%alloc (); call param_01%alloc ();
        call xsec_10%alloc (); call param_10%alloc ();
        call xsec_11%alloc (); call param_11%alloc ();
        call xsec_a%alloc (); call param_a%alloc ();
        call xsec_b%alloc (); call param_b%alloc ();
        call xsec_xx%alloc (); call param_xx%alloc ();
        
        call a_xsec%alloc (); call a_param%alloc ();
        
        ! read from fort.999
        do ia = ns%state%layer_bottom+1, ns%state%layer - ns%state%layer_top
            do iz = 1, ns%state%zone
                if (mat_info%mask_core(iz, ia))  then
                    read(999, fmt="(1x, 2(I5, TR2), *(ES12.5, TR2))")  iia, iiz, self_fdbk%Tf%new(iz,ia), self_fdbk%Tm%new(iz,ia), self_fdbk%Rho_m%new(iz,ia)
                end if  
            end do 
        end do 
        
        do idx = 1, 1
            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param)
        end do 
        call a_total%update ()
        call a_cpu%update ()
        
        do idx = 1, 1
            call read_xsec_unknown (self_fdbk, self_link, mat_info, xsec, param)
        end do 
        call a_total%update ()
        call a_cpu%update ()
        
        ! test-a, mesh-center
        if (.TRUE.)  then 
            iz = 2
            ia = 9
            im = 5 
            self_fdbk%Tf%new(iz,ia) = (1200.0D0 + 1500.0D0) / 2.0D0
            self_fdbk%Rho_m%new(iz,ia) = (873.0D0 + 802.5D0) / 2.0D0 
            file_in = TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(EXTEND_H5))
            
            call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(841); call a_param%print(841); 
            
            call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(842); call a_param%print(842); 
            
            call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(843); call a_param%print(843); 
            
            index_00 = 9
            index_01 = 12
            index_10 = 8
            index_11 = 11
            xsec_00 = hdf5_interp%mats(im)%xsecs(index_00) .MULTI. 0.25D0
            xsec_01 = hdf5_interp%mats(im)%xsecs(index_01) .MULTI. 0.25D0
            xsec_10 = hdf5_interp%mats(im)%xsecs(index_10) .MULTI. 0.25D0
            xsec_11 = hdf5_interp%mats(im)%xsecs(index_11) .MULTI. 0.25D0
            param_00 = hdf5_interp%mats(im)%params(index_00) .MULTI. 0.25D0
            param_01 = hdf5_interp%mats(im)%params(index_01) .MULTI. 0.25D0
            param_10 = hdf5_interp%mats(im)%params(index_10) .MULTI. 0.25D0
            param_11 = hdf5_interp%mats(im)%params(index_11) .MULTI. 0.25D0
            
            xsec_tmp = xsec_00 .ADD. xsec_01 .ADD. xsec_10 .ADD. xsec_11
            param_tmp = param_00 .ADD. param_01 .ADD. param_10 .ADD. param_11
            call xsec_tmp%print(840); call param_tmp%print(840); 
        end if 
        
        ! test-b, line-center 
        if (.TRUE.)  then
            iz = 2
            ia = 9
            im = 5 
            self_fdbk%Tf%new(iz,ia) = (1200.0D0 + 1200.0D0) / 2.0D0
            self_fdbk%Rho_m%new(iz,ia) = (873.0D0 + 802.5D0) / 2.0D0 
            file_in = TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(EXTEND_H5))
            
            call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(851); call a_param%print(851); 
            
            call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(852); call a_param%print(852); 
            
            call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(853); call a_param%print(853); 
            
            index_00 = 9
            index_01 = 9
            index_10 = 8
            index_11 = 8
            xsec_00 = hdf5_interp%mats(im)%xsecs(index_00) .MULTI. 0.25D0
            xsec_01 = hdf5_interp%mats(im)%xsecs(index_01) .MULTI. 0.25D0
            xsec_10 = hdf5_interp%mats(im)%xsecs(index_10) .MULTI. 0.25D0
            xsec_11 = hdf5_interp%mats(im)%xsecs(index_11) .MULTI. 0.25D0
            param_00 = hdf5_interp%mats(im)%params(index_00) .MULTI. 0.25D0
            param_01 = hdf5_interp%mats(im)%params(index_01) .MULTI. 0.25D0
            param_10 = hdf5_interp%mats(im)%params(index_10) .MULTI. 0.25D0
            param_11 = hdf5_interp%mats(im)%params(index_11) .MULTI. 0.25D0
            
            xsec_tmp = xsec_00 .ADD. xsec_01 .ADD. xsec_10 .ADD. xsec_11
            param_tmp = param_00 .ADD. param_01 .ADD. param_10 .ADD. param_11
            call xsec_tmp%print(850); call param_tmp%print(850); 
        end if 
        
        ! test-c, point-center 
        if (.TRUE.)  then
            iz = 2
            ia = 9
            im = 5 
            self_fdbk%Tf%new(iz,ia) = (1200.0D0 + 1200.0D0) / 2.0D0
            self_fdbk%Rho_m%new(iz,ia) = (873.0D0 + 873.0D0) / 2.0D0 
            file_in = TRIM(ADJUSTL(mat_info%libs(im)%file)) // TRIM(ADJUSTL(EXTEND_H5))
            
            call xsec_interpolation_lin (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(861); call a_param%print(861); 
            
            call xsec_interpolation_pwl (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(862); call a_param%print(862); 
            
            call xsec_interpolation_lag (self_fdbk, self_link, file_in, iz, ia, im, a_xsec, a_param)
            call a_xsec%print(863); call a_param%print(863); 
            
            index_00 = 8
            index_01 = 8
            index_10 = 8
            index_11 = 8
            xsec_00 = hdf5_interp%mats(im)%xsecs(index_00) .MULTI. 0.25D0
            xsec_01 = hdf5_interp%mats(im)%xsecs(index_01) .MULTI. 0.25D0
            xsec_10 = hdf5_interp%mats(im)%xsecs(index_10) .MULTI. 0.25D0
            xsec_11 = hdf5_interp%mats(im)%xsecs(index_11) .MULTI. 0.25D0
            param_00 = hdf5_interp%mats(im)%params(index_00) .MULTI. 0.25D0
            param_01 = hdf5_interp%mats(im)%params(index_01) .MULTI. 0.25D0
            param_10 = hdf5_interp%mats(im)%params(index_10) .MULTI. 0.25D0
            param_11 = hdf5_interp%mats(im)%params(index_11) .MULTI. 0.25D0
            
            xsec_tmp = xsec_00 .ADD. xsec_01 .ADD. xsec_10 .ADD. xsec_11
            param_tmp = param_00 .ADD. param_01 .ADD. param_10 .ADD. param_11
            call xsec_tmp%print(860); call param_tmp%print(860); 
        end if 
        
    end subroutine Check_interpolation
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine test_interpolation ()
        
        real(KREAL)  :: key(10) = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
        real(KREAL)  :: value(10) = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
        integer      :: scheme = 1 
        real(KREAL)  :: point 
        real(KREAL)  :: output 
        
        integer  :: i 
        
        write(*, *) TRIM(CHAR_SUBMARK)
        do i = 1, 20
            point = 1.92 + 0.0049*i
            output = array_interpolation (value, point, scheme)
            write(*, fmt='(1x, *(A, ES12.5))') 'point = ', point, ';  output = ', output
        end do 
        
        write(*, *) TRIM(CHAR_SUBMARK)
        do i = 1, 20
            point = 1.92 + 0.0049*i
            output = array_interpolation (key, point, scheme)
            write(*, fmt='(1x, *(A, ES12.5))') 'point = ', point, ';  output = ', output
        end do 
        
    end subroutine test_interpolation 

end module hdf5_interpolation
