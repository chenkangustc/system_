!$
!===================================================================================================
!
!   module for thermal calculation pre process
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_th_pre_process
!
!   Public type lists:          No
!
!===================================================================================================
module th_pre_process

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global 
    use th_global
    
    use files_dirs
    
    implicit none
    private 
    public  :: Driving_th_pre_process
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_th_pre_process ()
        
        character(len=MAX_WORD_LEN)  :: file_in
        integer  :: file_unit
        
        integer  :: i_kind                                                      
        integer  :: i_type                                                      ! 0 for water, 1 for REFPROP
        integer  :: option 
        
        ! ----------------------------------------------------------------------      
        ! allocate memory   
        call avg_channel%alloc (nth)
        call hot_channel%alloc (nth)
        
        ! pointer allocated
        call Prepare_property_pointer ()
        
        ! coolant selection
        i_type = MOD(geom_th%coolant_type, 10)
        i_kind = (geom_th%coolant_type - i_type) / 10
        option = 0 
        select case(i_kind)
        case(0)
            select case(i_type)
            case(1)                                                             
                a_coolant => a_coolant_parcs
                call a_coolant%set (i_type, option)
            case(2)                                                             ! water @ 15.5-MPa
                option = 1 
                a_coolant => a_coolant_refprop
                call a_coolant%set (i_type, option)
            case(3)                                                             ! water @ 7.0-Mpa
                option = 2 
                a_coolant => a_coolant_refprop
                call a_coolant%set (i_type, option) 
            end select
            
        case(1)                                                                 ! HLM
            a_coolant => a_coolant_HLM
            call a_coolant%set (i_type, option)
        case default
        end select
        
        ! ----------------------------------------------------------------------
        ! set initial value
        call geom_th%set_height (geom%height)
        
        ! thermal output
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(TH_WARNING)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%TH_WARNING = file_unit

        file_in = TRIM(PREFIX_OUTPUT) // TRIM(TH_HOT)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%TH_HOT = file_unit
    
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(TH_AVERAGE)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%TH_AVERAGE = file_unit
    
        file_in = TRIM(PREFIX_OUTPUT) // TRIM(TH_RBFD)
        call Open_file (file_in, is_input=.FALSE., file_unit=file_unit)
        FILES%TH_RBFD = file_unit
    
    end subroutine Driving_th_pre_process
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Prepare_property_pointer ()
        
        integer  :: i_allocate
        
        if (associated(a_coolant_parcs))         deallocate(a_coolant_parcs, stat=i_allocate)
        if (associated(a_coolant_refprop))       deallocate(a_coolant_refprop, stat=i_allocate)
        if (associated(a_coolant_HLM))           deallocate(a_coolant_HLM, stat=i_allocate)
        if (associated(a_clad_Zr))               deallocate(a_clad_Zr, stat=i_allocate)
        if (associated(a_clad_steels))           deallocate(a_clad_steels, stat=i_allocate)
        if (associated(a_gap_gas))               deallocate(a_gap_gas, stat=i_allocate)
        if (associated(a_fuel_metallic))         deallocate(a_fuel_metallic, stat=i_allocate)
        if (associated(a_fuel_ceramic))          deallocate(a_fuel_ceramic, stat=i_allocate)
        
        allocate(a_coolant_parcs, stat=i_allocate)
        allocate(a_coolant_refprop, stat=i_allocate)
        allocate(a_coolant_HLM, stat=i_allocate)
        allocate(a_clad_Zr, stat=i_allocate)
        allocate(a_clad_steels, stat=i_allocate)
        allocate(a_gap_gas, stat=i_allocate)
        allocate(a_fuel_metallic, stat=i_allocate)
        allocate(a_fuel_ceramic, stat=i_allocate)
    
    end subroutine Prepare_property_pointer
    
end module th_pre_process
