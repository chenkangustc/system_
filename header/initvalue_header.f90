!$
!===================================================================================================
!
!   class for initial value of iteration
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          InitialValue
!
!===================================================================================================
module initvalue_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use snap_interface,         only : snap_array_read, snap_array_write
    use iteration_header,       only : IterationSource, IterationFlux, IterationCounter, AnisotropicScatterFlux
    
    implicit none
    private
    public  :: InitialValue
    
    ! type for initial value calculation
    type  InitialValue
        logical, public                      :: is_initValin  = .FALSE.         ! import initial value 
        logical, public                      :: is_initValout = .FALSE.         ! export initial value 
        character(len=MAX_WORD_LEN), public      :: fluxfile      = 'initFlux.bin'
    contains
        procedure, public  :: import => Import_InitialValue
        procedure, public  :: export => Export_InitialValue
    end type InitialValue
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Import_InitialValue(this, iter_q, iter_flux, iter_count, flux_scat)
        
        class(InitialValue), intent(in out)  :: this
        type(IterationSource), intent(in out)  :: iter_q
        type(IterationFlux), intent(in out)    :: iter_flux
        type(IterationCounter), intent(in out) :: iter_count
        type(AnisotropicScatterFlux), intent(in out)     :: flux_scat
        
        integer  :: unit_, io_error
        logical  :: is_true
        
        if (.NOT. this%is_initValin)  then
            return
        else
            inquire(file=this%fluxfile, exist=is_true)
            
            if (is_true)  then
                open(newunit=unit_, file=this%fluxfile, access='sequential', form='unformatted', status='old', action='read', iostat=io_error)
                call snap_array_read(unit_, iter_q%fission%moment)
                call snap_array_read(unit_, iter_q%fission%old)
!                call snap_array_read(unit_, iter_q%info%total_moment)
!                call snap_array_read(unit_, iter_q%info%out_group_moment) 
                
                call snap_array_read(unit_, iter_flux%info%moment)
                call snap_array_read(unit_, iter_flux%info%moment_omp)
                call snap_array_read(unit_, iter_flux%info%old)
                
                call snap_array_read(unit_, iter_flux%dist%surface)
                call snap_array_read(unit_, iter_flux%dist%nodal)
                call snap_array_read(unit_, iter_flux%dist%point)
                call snap_array_read(unit_, iter_flux%dist%axi_surf)
                call snap_array_read(unit_, iter_flux%dist%rad_surf)
                
                read(unit_)  iter_count%eigenvalue, iter_count%ks, iter_count%ksub, iter_count%coeff_FSP, iter_count%kcritical
                close(unit=unit_, status='keep', iostat=io_error)
            end if 
            
            this%is_initValin = .FALSE.
        end if 
        
    end subroutine Import_InitialValue
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Export_InitialValue(this, iter_q, iter_flux, iter_count, flux_scat)
        
        class(InitialValue), intent(in out)  :: this
        type(IterationSource), intent(in out)  :: iter_q
        type(IterationFlux), intent(in out)    :: iter_flux
        type(IterationCounter), intent(in out) :: iter_count
        type(AnisotropicScatterFlux), intent(in out)     :: flux_scat
        
        integer  :: unit_, io_error
        
        if (.NOT. this%is_initValout)  then
            return
        else
            open(newunit=unit_, file=this%fluxfile, access='sequential', form='unformatted', status='replace', action='write', iostat=io_error)
            call snap_array_write(unit_, iter_q%fission%moment)
            call snap_array_write(unit_, iter_q%fission%old)
!            call snap_array_write(unit_, iter_q%info%total_moment)
!            call snap_array_write(unit_, iter_q%info%out_group_moment)
            
            call snap_array_write(unit_, iter_flux%info%moment)
            call snap_array_write(unit_, iter_flux%info%moment_omp)
            call snap_array_write(unit_, iter_flux%info%old)
            
            call snap_array_write(unit_, iter_flux%dist%surface)
            call snap_array_write(unit_, iter_flux%dist%nodal)
            call snap_array_write(unit_, iter_flux%dist%point)
            call snap_array_write(unit_, iter_flux%dist%axi_surf)
            call snap_array_write(unit_, iter_flux%dist%rad_surf)
            
            write(unit_)  iter_count%eigenvalue, iter_count%ks, iter_count%ksub, iter_count%coeff_FSP, iter_count%kcritical
            close(unit=unit_, status='keep', iostat=io_error)
            
            this%is_initValout = .FALSE.
        end if 
        
    end subroutine Export_InitialValue
    
end module initvalue_header 
