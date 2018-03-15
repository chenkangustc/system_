!$
!===================================================================================================
!
!   module for thermal calculation pre process
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_th_check_model 
!
!   Public type lists:          No
!
!===================================================================================================
module th_check_property 

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use th_global
    
    use string
    use files_dirs
    
    implicit none
    private
    public  :: Driving_th_check_property

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Driving_th_check_property ()
        
        integer  :: fuel_metallic_idx(3) = [1, 2, 3]
        integer  :: fuel_ceramic_idx(9) = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        integer  :: clad_Zr_idx(6) = [1, 2, 3, 4, 5, 6]
        integer  :: clad_steels_idx(3) = [1, 2, 3]
        integer  :: coolant_HLM_idx(4) = [1, 2, 3, 4]
        integer  :: gap_idx(5) = [1, 2, 3, 4, 5]
        
        character(len=200)  :: fname
        real(KREAL)  :: tmin, tmax
        real(KREAL)  :: tset
        real(KREAL)  :: hset
        integer      :: nstep
        integer      :: unit_        
        integer      :: idx, i
        
        call Prepare_property_pointer ()
        
        ! check fuel
        tmin = 300.0
        tmax = 1700.0
        do i = 1, SIZE(fuel_metallic_idx)
            idx = fuel_metallic_idx(i)
            fname = "outp/fuel_metallic_" // Int_to_string(idx, 1)
            a_fuel => a_fuel_metallic
            call a_fuel%set (idx)
            call Check_fuelrod (fname, tmin, tmax)
        end do
       
        tmin = 300.0
        tmax = 3000.0
        do i = 1, SIZE(fuel_ceramic_idx)
            idx = fuel_ceramic_idx(i)
            fname = "outp/fuel_ceramic_" // Int_to_string(idx, 1)
            a_fuel => a_fuel_ceramic
            call a_fuel%set (idx)
            call Check_fuelrod (fname, tmin, tmax)
        end do
        
        ! check clad
        tmin = 300.0
        tmax = 2100.0
        do i = 1, SIZE(clad_Zr_idx)
            idx = clad_Zr_idx(i)
            fname = "outp/clad_Zr_" // Int_to_string(idx, 1)
            a_clad => a_clad_Zr
            call a_clad%set (idx)
            call Check_clad (fname, tmin, tmax)
        end do
       
        tmin = 300.0
        tmax = 1600.0
        do i = 1, SIZE(clad_steels_idx)
            idx = clad_steels_idx(i)
            fname = "outp/clad_steels_" // Int_to_string(idx, 1)
            a_clad => a_clad_steels
            call a_clad%set (idx)
            call Check_clad (fname, tmin, tmax)
        end do
        
        ! check coolant
        tmin = 300.0
        tmax = 1300.0
        do i = 1, SIZE(coolant_HLM_idx)
            idx = coolant_HLM_idx(i)
            fname = "outp/coolant_HLM_" // Int_to_string(idx, 1)
            a_coolant => a_coolant_HLM
            call a_coolant%set (idx, option=0)
            call Check_coolant (fname, tmin, tmax)
        end do
        
        tmin = 550.0
        tmax = 615.0
        fname = "outp/coolant_water_parcs"
        a_coolant => a_coolant_parcs
        call a_coolant%set (idx, option=0)
        call Check_coolant (fname, tmin, tmax)
        
!        fname = "outp/coolant_water_refprop"
!        a_coolant => a_coolant_refprop
!        call a_coolant%set (idx, option=0)
!        call Check_coolant (fname, tmin, tmax)
        
        ! check coolant enthalpy
        tmin = 300.0
        tmax = 1300.0
        do i = 1, SIZE(coolant_HLM_idx)
            idx = coolant_HLM_idx(i)
            fname = "outp/coolantH_HLM_" // Int_to_string(idx, 1)
            a_coolant => a_coolant_HLM
            call a_coolant%set (idx, option=0)
            call Check_enthalpy (fname, tmin, tmax)
        end do
        
        tmin = 550.0
        tmax = 615.0
        fname = "outp/coolantH_water_parcs"
        a_coolant => a_coolant_parcs
        call a_coolant%set (idx, option=0)
        call Check_enthalpy (fname, tmin, tmax)
        
!        fname = "outp/coolantH_water_refprop"
!        a_coolant => a_coolant_refprop
!        call a_coolant%set (idx, option=0)
!        call Check_enthalpy (fname, tmin, tmax)
        
        ! check gap
        tmin = 300.0
        tmax = 1700.0
        do i = 1, SIZE(gap_idx)
            idx = gap_idx(i)
            fname = "outp/gap_" // Int_to_string(idx, 1)
            a_gap => a_gap_gas
            call a_gap%set (idx)
            call Check_gap (fname, tmin, tmax)
        end do
        
        ! ----------------------------------------------------------------------
        ! compared with parcs
        open(newunit=unit_, file='test-FUEL')
        tmin = 300.0
        tmax = 3000.0
        nstep = NINT(tmax-tmin)
        a_fuel => a_fuel_ceramic
        call a_fuel%set (4)
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,             &
                &   a_fuel%get_conductivity(tset), a_fuel%get_capacity(tset)
        end do
        close(unit=unit_)

        open(newunit=unit_, file='test-CLAD')
        tmin = 300.0
        tmax = 2100.0
        nstep = NINT(tmax-tmin)
        a_clad => a_clad_Zr
        call a_clad%set (4)
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,             &
                &   a_clad%get_conductivity(tset), a_clad%get_capacity(tset)
        end do
        close(unit=unit_)
        
        open(newunit=unit_, file='test-WATER')
        tmin = 550.0
        tmax = 615.0
        nstep = NINT(tmax-tmin)
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            hset = a_coolant%get_enthalpy (tset)
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,             &
                &   a_coolant%get_conductivity(tset), a_coolant%get_density(tset), a_coolant%get_temperature (hset), a_coolant%get_viscosity(tset)
        end do
        close(unit=unit_)
        
        
        call Clean_property_pointer ()
        
    end subroutine Driving_th_check_property

    !$
    !===============================================================================================
    ! check fuel property
    !===============================================================================================
    subroutine Check_fuelrod (fname, tmin, tmax)
        
        character(len=200), intent(in)  :: fname
        real(KREAL), intent(in)         :: tmin
        real(KREAL), intent(in)         :: tmax
        
        ! local
        real(KREAL)  :: tset
        integer      :: nstep
        integer      :: unit_
        integer      :: i
        
        call Open_file (fname, .FALSE., unit_)
        
        nstep = NINT(tmax-tmin)
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,             &
                &   a_fuel%get_density(tset),                                   &
                &   a_fuel%get_capacity(tset),                                  &
                &   a_fuel%get_density(tset)*a_fuel%get_capacity(tset),         &
                &   a_fuel%get_conductivity(tset)
        end do
        
        call Close_file (unit_, .TRUE.)
    
    end subroutine Check_fuelrod
        
    !$
    !===============================================================================================
    ! check clad property
    !===============================================================================================
    subroutine Check_clad (fname, tmin, tmax)
        
        character(len=200), intent(in)  :: fname
        real(KREAL), intent(in)         :: tmin
        real(KREAL), intent(in)         :: tmax
        
        ! local
        real(KREAL)  :: tset
        integer      :: nstep
        integer      :: unit_
        integer      :: i
        
        call Open_file (fname, .FALSE., unit_)
        
        nstep = NINT(tmax-tmin)
        
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,             &
                &   a_clad%get_density(tset),                                   &
                &   a_clad%get_capacity(tset),                                  &
                &   a_clad%get_density(tset)*a_clad%get_capacity(tset),         &
                &   a_clad%get_conductivity(tset)
        end do
        
        call Close_file (unit_, .TRUE.)
    
    end subroutine Check_clad
            
    !$
    !===============================================================================================
    ! check coolant property
    !===============================================================================================
    subroutine Check_coolant (fname, tmin, tmax)
        
        character(len=200), intent(in)  :: fname
        real(KREAL), intent(in)         :: tmin
        real(KREAL), intent(in)         :: tmax
        
        ! local
        real(KREAL)  :: tset
        integer      :: nstep
        integer      :: unit_
        integer      :: i
        
        call Open_file (fname, .FALSE., unit_)
        
        nstep = NINT(tmax-tmin)
        
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,               &
                &   a_coolant%get_density(tset),                                &
                &   a_coolant%get_capacity(tset),                               &
                &   a_coolant%get_density(tset)*a_coolant%get_capacity(tset),   &
                &   a_coolant%get_conductivity(tset),                           &
                &   a_coolant%get_viscosity(tset)
        end do
        
        call Close_file (unit_, .TRUE.)
    
    end subroutine Check_coolant
        
    !$
    !===============================================================================================
    ! check temperature & enthalpy
    !===============================================================================================
    subroutine Check_enthalpy (fname, tmin, tmax)
        
        character(len=200), intent(in)  :: fname
        real(KREAL), intent(in)         :: tmin
        real(KREAL), intent(in)         :: tmax
        
        ! local
        real(KREAL)  :: tset, tget, terror
        real(KREAL)  :: enthalpy
        integer      :: nstep
        integer      :: unit_
        integer      :: i
        
        call Open_file (fname, .FALSE., unit_)
        
        nstep = NINT(tmax-tmin)
        
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            enthalpy = a_coolant%get_enthalpy (tset)
            tget = a_coolant%get_temperature (enthalpy)
            terror = 100.0 * (tget - tset) / tset
            write(unit=unit_, fmt="(1x, I5, TR3, *(ES11.4, TR3))") i, tset, tget, terror
        end do
        
        call Close_file (unit_, .TRUE.)
    
    end subroutine Check_enthalpy
    
    !$
    !===============================================================================================
    ! check gas property
    !===============================================================================================
    subroutine Check_gap (fname, tmin, tmax)
        
        character(len=200), intent(in)  :: fname
        real(KREAL), intent(in)         :: tmin
        real(KREAL), intent(in)         :: tmax
        
        ! local
        real(KREAL)  :: tset
        integer      :: nstep
        integer      :: unit_
        integer      :: i
        
        call Open_file (fname, .FALSE., unit_)
        
        nstep = NINT(tmax-tmin)
        
        do i = 1, nstep+1
            tset = tmin + (i-1) * (tmax-tmin) / nstep
            write(unit=unit_, fmt="(1x, I5, TR3, F8.2, TR3, *(ES11.4, TR3))") i, tset,               &
                &   a_gap%get_transfer(tset, 0.003951_KREAL, 0.000059_KREAL, .FALSE.)
        end do
        
        call Close_file (unit_, .TRUE.)
    
    end subroutine Check_gap
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Prepare_property_pointer ()
        
        integer  :: i_allocate
        
        if (associated(a_coolant_parcs))         deallocate(a_coolant_parcs, stat=i_allocate)
!        if (associated(a_coolant_refprop))       deallocate(a_coolant_refprop, stat=i_allocate)
        if (associated(a_coolant_HLM))           deallocate(a_coolant_HLM, stat=i_allocate)
        if (associated(a_clad_Zr))               deallocate(a_clad_Zr, stat=i_allocate)
        if (associated(a_clad_steels))           deallocate(a_clad_steels, stat=i_allocate)
        if (associated(a_gap_gas))               deallocate(a_gap_gas, stat=i_allocate)
        if (associated(a_fuel_metallic))         deallocate(a_fuel_metallic, stat=i_allocate)
        if (associated(a_fuel_ceramic))          deallocate(a_fuel_ceramic, stat=i_allocate)
        
        allocate(a_coolant_parcs, stat=i_allocate)
!        allocate(a_coolant_refprop, stat=i_allocate)
        allocate(a_coolant_HLM, stat=i_allocate)
        allocate(a_clad_Zr, stat=i_allocate)
        allocate(a_clad_steels, stat=i_allocate)
        allocate(a_gap_gas, stat=i_allocate)
        allocate(a_fuel_metallic, stat=i_allocate)
        allocate(a_fuel_ceramic, stat=i_allocate)
    
    end subroutine Prepare_property_pointer
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Clean_property_pointer ()
        
        integer  :: i_allocate
        
        if (associated(a_coolant_parcs))         deallocate(a_coolant_parcs, stat=i_allocate)
!        if (associated(a_coolant_refprop))       deallocate(a_coolant_refprop, stat=i_allocate)
        if (associated(a_coolant_HLM))           deallocate(a_coolant_HLM, stat=i_allocate)
        if (associated(a_clad_Zr))               deallocate(a_clad_Zr, stat=i_allocate)
        if (associated(a_clad_steels))           deallocate(a_clad_steels, stat=i_allocate)
        if (associated(a_gap_gas))               deallocate(a_gap_gas, stat=i_allocate)
        if (associated(a_fuel_metallic))         deallocate(a_fuel_metallic, stat=i_allocate)
        if (associated(a_fuel_ceramic))          deallocate(a_fuel_ceramic, stat=i_allocate)
        
    end subroutine Clean_property_pointer
    
end module th_check_property
