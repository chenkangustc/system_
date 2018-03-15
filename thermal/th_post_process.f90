!$
!===================================================================================================
!
!   module for thermal calculation post process
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Driving_th_post_process
!
!   Public type lists:          No
!
!===================================================================================================
module th_post_process

    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV

    use files_dirs
    use th_global
    
    implicit none
    private 
    public  :: Driving_th_post_process

contains
    !$
    !===============================================================================================
    ! free memory in thermal calculation
    !===============================================================================================
    subroutine Driving_th_post_process ()
        
        logical  :: is_keep
        integer  :: i
        
        if (.FALSE.)  then
            do i = 1, 20
                call avg_channel%print_rod (997, nth, ir=i)
            end do 
        end if 
        
        ! destory pointer
        call Clean_property_pointer ()
        
        call design%clean ()
        call th_power%clean ()
        call avg_channel%clean ()
        call hot_channel%clean ()
        
        ! assembly information
        call geom_th%clean ()

        if (allocated(geom_assm))  then
            do i = 1, SIZE(geom_assm)
                call geom_assm(i)%clean ()
            end do
            deallocate(geom_assm)
        end if
        
        ! close thermal output
        is_keep = (ns%feedback%is_feedback) .and. (ns%feedback%is_inner)
        
        call Close_file (FILES%TH_WARNING, is_keep=is_keep)
        call Close_file (FILES%TH_HOT, is_keep=is_keep)
        call Close_file (FILES%TH_AVERAGE, is_keep=is_keep)
        call Close_file (FILES%TH_RBFD, is_keep=is_keep)
    
    end subroutine Driving_th_post_process
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Clean_property_pointer ()
        
        integer  :: i_allocate
        
        if (associated(a_coolant_parcs))         deallocate(a_coolant_parcs, stat=i_allocate)
        if (associated(a_coolant_refprop))       deallocate(a_coolant_refprop, stat=i_allocate)
        if (associated(a_coolant_HLM))           deallocate(a_coolant_HLM, stat=i_allocate)
        if (associated(a_clad_Zr))               deallocate(a_clad_Zr, stat=i_allocate)
        if (associated(a_clad_steels))           deallocate(a_clad_steels, stat=i_allocate)
        if (associated(a_gap_gas))               deallocate(a_gap_gas, stat=i_allocate)
        if (associated(a_fuel_metallic))         deallocate(a_fuel_metallic, stat=i_allocate)
        if (associated(a_fuel_ceramic))          deallocate(a_fuel_ceramic, stat=i_allocate)
        
    end subroutine Clean_property_pointer
    
end module th_post_process
