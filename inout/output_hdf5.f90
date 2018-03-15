!$
!===================================================================================================
!
!   module of subroutines to output distribution information to HDF5 file
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    
!                               
!                               
!
!   Public type lists:          No
!
!===================================================================================================
module output_hdf5
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use global
    use th_global
    use HDF5
    use H5LT
    use hdf5_interface
    
    implicit none 
    private
    public  :: Print_binary_hdf5
    
    character(len=MAX_WORD_LEN), parameter   :: HDF5_FILENAME = 'output.h5'
    logical, save                            :: is_hdf5_exist = .FALSE.         ! NOTE: save
    
    integer  :: hdferr
    integer  :: rank
    integer  :: dim_scale(1:10)
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_binary_hdf5 (is_adjoint, is_transient, tidx, ctime)
        
        logical, intent(in)  :: is_adjoint
        logical, intent(in)  :: is_transient
        integer, intent(in)  :: tidx
        real(KREAL), intent(in) :: ctime
        
        character(len=MAX_WORD_LEN)  :: groupname
        character(len=MAX_WORD_LEN)  :: dataname
        character(len=MAX_WORD_LEN)  :: tmp
        integer(HID_T)  :: file
        integer(HID_T)  :: casegroup, typegroup, nggroup, dggroup
        
        ! container for data
        real(KREAL)  :: hdf5_nodal_layer(ns%state%nodal, ns%state%layer)
        real(KREAL)  :: hdf5_nodal(ns%state%nodal)
        real(KREAL)  :: hdf5_zone_layer(ns%state%zone, ns%state%layer)
        real(KREAL)  :: hdf5_zone(ns%state%zone)
        real(KREAL)  :: hdf5_layer(ns%state%layer)
        
        integer  :: i_allocate
        integer  :: io_error
        integer  :: ig, is, ip
        
        ! do not output
        if (.NOT. ns%output%is_HDF5)  then
            return
        end if
        
        ! initialize fortran interface
        call h5open_f(hdferr)

        ! create a new file using the default properties
        if (.NOT. is_hdf5_exist)  then
            call hdf5_file_create(TRIM(ADJUSTL(HDF5_FILENAME)), file)
            is_hdf5_exist = .TRUE.
        else 
            call hdf5_file_open(TRIM(ADJUSTL(HDF5_FILENAME)), file, 'w')
        end if
        
        ! ----------------------------------------------------------------------
        ! for steady stat forward
        if (.NOT. is_adjoint .and. .NOT. is_transient)  then
            groupname = 'steady'
            call hdf5_open_group(file, groupname, casegroup)

            call hdf5_set_attr(file, groupname, 'nodal', ns%state%nodal)
            call hdf5_set_attr(file, groupname, 'layer', ns%state%layer)
            call hdf5_set_attr(file, groupname, 'zone', ns%state%zone)
            call hdf5_set_attr(file, groupname, 'point', ns%state%point)
            call hdf5_set_attr(file, groupname, 'ng', ns%state%ng)
            
            ! angular flux
            groupname = 'flux_angular'
            call hdf5_open_group(casegroup, groupname, typegroup)
            call hdf5_set_attr(casegroup, groupname, 'direction', ns%deduce%direction)
            do ig = 1, ns%state%ng
                if (ig <= 9)  then
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_00' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else if (ig <= 99)  then
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_0' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                end if
                call hdf5_open_group(typegroup, groupname, nggroup)
                
                dataname = 'all_direction'
                rank = 3
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                dim_scale(3) = ns%deduce%direction
                call hdf5_write_data(nggroup, dataname, flux_forward%ngs(ig)%angular, dim_scale(1:rank))
                
                do is = 1, ns%deduce%direction
                    if (is <= 9)  then
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_00' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    else if (is <= 99)  then
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_0' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    else
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    end if
                    
                    rank = 2
                    dim_scale(1) = ns%state%nodal
                    dim_scale(2) = ns%state%layer
                    call hdf5_write_data(nggroup, dataname, flux_forward%ngs(ig)%angular(:, :, is), dim_scale(1:rank))
                end do
                
                call hdf5_close_group(nggroup)
                
            end do
            call hdf5_close_group(typegroup)
            
            ! scalar flux
            groupname = 'flux_scalar'
            call hdf5_open_group(casegroup, groupname, typegroup)
            do ig = 1, ns%state%ng
                if (ig <= 9)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_00' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else if (ig <= 99)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_0' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                end if
                rank = 2
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                call hdf5_write_data(typegroup, dataname, flux_forward%ngs(ig)%scalar, dim_scale(1:rank))
            end do
            call hdf5_close_group(typegroup)
            
            ! flux mapping
            groupname = 'mapping_flux'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_flux%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_flux%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_flux%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_flux%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_flux%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            ! power mapping
            groupname = 'mapping_power'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_power%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_power%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_power%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_power%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_power%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            ! fission rate mapping
            groupname = 'mapping_fission_rate'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_fission_rate%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_fission_rate%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_fission_rate%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_fission_rate%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_fission_rate%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            if (ns%feedback%is_feedback .and. ns%feedback%is_inner)  then
                call Print_thermal_hdf5 (casegroup)
            end if
            
            ! close steady
            call hdf5_close_group(casegroup)
        end if
        
        ! ----------------------------------------------------------------------
        ! for steady stat adjoint
        if (is_adjoint )  then
            groupname = 'adjoint'
            call hdf5_open_group(file, groupname, casegroup)

            call hdf5_set_attr(file, groupname, 'nodal', ns%state%nodal)
            call hdf5_set_attr(file, groupname, 'layer', ns%state%layer)
            call hdf5_set_attr(file, groupname, 'zone', ns%state%zone)
            call hdf5_set_attr(file, groupname, 'point', ns%state%point)
            call hdf5_set_attr(file, groupname, 'ng', ns%state%ng)
            
            ! angular flux
            groupname = 'flux_angular'
            call hdf5_open_group(casegroup, groupname, typegroup)
            call hdf5_set_attr(casegroup, groupname, 'direction', ns%deduce%direction)
            
            do ig = 1, ns%state%ng
                if (ig <= 9)  then
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_00' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else if (ig <= 99)  then
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_0' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else
                    write(unit=tmp, fmt=*) ig
                    groupname = 'group_' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                end if
                
                call hdf5_open_group(typegroup, groupname, nggroup)
                
                dataname = 'all_direction'
                rank = 3
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                dim_scale(3) = ns%deduce%direction
                call hdf5_write_data(nggroup, dataname, flux_adjoint%ngs(ig)%angular, dim_scale(1:rank))
                
                do is = 1, ns%deduce%direction
                    if (is <= 9)  then
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_00' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    else if (is <= 99)  then
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_0' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    else
                        write(unit=tmp, fmt=*) is
                        dataname = 'direction_' // TRIM(ADJUSTL(tmp))
                        tmp = ' '
                    end if
                    
                    rank = 2
                    dim_scale(1) = ns%state%nodal
                    dim_scale(2) = ns%state%layer
                    call hdf5_write_data(nggroup, dataname, flux_adjoint%ngs(ig)%angular(:, :, is), dim_scale(1:rank))
                end do
                
                call hdf5_close_group(nggroup)
                
            end do
            call hdf5_close_group(typegroup)
            
            ! scalar flux
            groupname = 'flux_scalar'
            call hdf5_open_group(casegroup, groupname, typegroup)
            
            do ig = 1, ns%state%ng
                if (ig <= 9)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_00' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else if (ig <= 99)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_0' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                end if
                
                rank = 2
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                call hdf5_write_data(typegroup, dataname, flux_adjoint%ngs(ig)%scalar, dim_scale(1:rank))
                
            end do
            call hdf5_close_group(typegroup)
            
            ! close adjoint
            call hdf5_close_group(casegroup)
        end if
        
        ! ----------------------------------------------------------------------
        ! for transient stat
        if (is_transient )  then
            if (tidx <= 9)  then
                write(unit=tmp, fmt=*) tidx
                groupname = 'index_000' // TRIM(ADJUSTL(tmp))
                tmp = ' '
            else if (tidx <= 99)  then
                write(unit=tmp, fmt=*) tidx
                groupname = 'index_00' // TRIM(ADJUSTL(tmp))
                tmp = ' '
            else if (tidx <= 999)  then
                write(unit=tmp, fmt=*) tidx
                groupname = 'index_0' // TRIM(ADJUSTL(tmp))
                tmp = ' '
            else 
                write(unit=tmp, fmt=*) tidx
                groupname = 'index_' // TRIM(ADJUSTL(tmp))
                tmp = ' '
            end if
        
            call hdf5_open_group(file, groupname, casegroup)

            call hdf5_set_attr(file, groupname, 'nodal', ns%state%nodal)
            call hdf5_set_attr(file, groupname, 'layer', ns%state%layer)
            call hdf5_set_attr(file, groupname, 'zone', ns%state%zone)
            call hdf5_set_attr(file, groupname, 'point', ns%state%point)
            call hdf5_set_attr(file, groupname, 'ng', ns%state%ng)
            
            call hdf5_set_attr(file, groupname, 'dg', nt%state%dg)
            call hdf5_set_attr(file, groupname, 'index', tidx)
            call hdf5_set_attr(file, groupname, 'time', ctime)
            
            ! scalar flux
            groupname = 'flux_scalar'
            call hdf5_open_group(casegroup, groupname, typegroup)
            do ig = 1, ns%state%ng
                if (ig <= 9)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_00' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else if (ig <= 99)  then
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_0' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                else
                    write(unit=tmp, fmt=*) ig
                    dataname = 'group_' // TRIM(ADJUSTL(tmp))
                    tmp = ' '
                end if
                rank = 2
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                call hdf5_write_data(typegroup, dataname, flux_forward%ngs(ig)%scalar, dim_scale(1:rank))
            end do
            call hdf5_close_group(typegroup)
            
            ! flux mapping
            groupname = 'mapping_flux'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_flux%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_flux%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_flux%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_flux%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_flux%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            ! power mapping
            groupname = 'mapping_power'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_power%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_power%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_power%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_power%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_power%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            ! fission rate mapping
            groupname = 'mapping_fission_rate'
            call hdf5_open_group(casegroup, groupname, typegroup)
                
            dataname = 'nodal_layer' 
            hdf5_nodal_layer = dist_fission_rate%matrix
            rank = 2
            dim_scale(1) = ns%state%nodal
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
            
            dataname = 'nodal_only'
            call dist_fission_rate%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
            rank = 1
            dim_scale(1) = ns%state%nodal
            call hdf5_write_data(typegroup, dataname, hdf5_nodal, dim_scale(1:rank))
            
            dataname = 'zone_layer'
            call dist_fission_rate%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
            rank = 2
            dim_scale(1) = ns%state%zone
            dim_scale(2) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
            
            dataname = 'zone_only'
            call dist_fission_rate%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
            rank = 1
            dim_scale(1) = ns%state%zone
            call hdf5_write_data(typegroup, dataname, hdf5_zone, dim_scale(1:rank))
            
            dataname = 'layer'
            call dist_fission_rate%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
            rank = 1
            dim_scale(1) = ns%state%layer
            call hdf5_write_data(typegroup, dataname, hdf5_layer, dim_scale(1:rank))

            call hdf5_close_group(typegroup)
            
            ! presursor mapping
            groupname = 'mapping_precursor'
            call hdf5_open_group(casegroup, groupname, typegroup)
            
            do ip = 1, nt%state%dg
                write(unit=tmp, fmt=*) ip
                groupname = 'group_00' // TRIM(ADJUSTL(tmp))
                tmp = ''
                
                call hdf5_open_group(typegroup, groupname, dggroup)
                
                dataname = 'nodal_layer' 
                hdf5_nodal_layer = dist_dnps(ip)%matrix
                rank = 2
                dim_scale(1) = ns%state%nodal
                dim_scale(2) = ns%state%layer
                call hdf5_write_data(dggroup, dataname, hdf5_nodal_layer, dim_scale(1:rank))
                
                dataname = 'nodal_only'
                call dist_dnps(ip)%nodal(mesh, geom, is_function=.FALSE., nodal=hdf5_nodal)
                rank = 1
                dim_scale(1) = ns%state%nodal
                call hdf5_write_data(dggroup, dataname, hdf5_nodal, dim_scale(1:rank))
                
                dataname = 'zone_layer'
                call dist_dnps(ip)%zone_layer(mesh, geom, is_function=.FALSE., zone_layer=hdf5_zone_layer)
                rank = 2
                dim_scale(1) = ns%state%zone
                dim_scale(2) = ns%state%layer
                call hdf5_write_data(dggroup, dataname, hdf5_zone_layer, dim_scale(1:rank))
                
                dataname = 'zone_only'
                call dist_dnps(ip)%zone(mesh, geom, is_function=.FALSE., zone=hdf5_zone)
                rank = 1
                dim_scale(1) = ns%state%zone
                call hdf5_write_data(dggroup, dataname, hdf5_zone, dim_scale(1:rank))
                
                dataname = 'layer'
                call dist_dnps(ip)%layer(mesh, geom, is_function=.FALSE., layer=hdf5_layer)
                rank = 1
                dim_scale(1) = ns%state%layer
                call hdf5_write_data(dggroup, dataname, hdf5_layer, dim_scale(1:rank))
                
                call hdf5_close_group(dggroup)
            end do

            call hdf5_close_group(typegroup)
            
            if (ns%feedback%is_feedback .and. ns%feedback%is_inner)  then
                call Print_thermal_hdf5 (casegroup)
            end if
            
            ! close this case
            call hdf5_close_group(casegroup)
        end if
        
        ! exit hdf5 environment
        call hdf5_file_close(file)
        call h5close_f(hdferr)
        
    end subroutine Print_binary_hdf5
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! write thermal information
    !===============================================================================================
    subroutine Print_thermal_hdf5 (casegroup)
    
        integer(HID_T), intent(in)  :: casegroup
    
        character(len=MAX_WORD_LEN)  :: groupname
        character(len=MAX_WORD_LEN)  :: dataname
        integer(HID_T)  :: typegroup
            
        ! ----------------------------------------------------------------------
        ! thermal: hot-channel
        groupname = 'thermal_hot'
        call hdf5_open_group(casegroup, groupname, typegroup)
            
        dataname = 'flow_velocity' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%flow_velocity, dim=2)
        dim_scale(2) = SIZE(hot_channel%flow_velocity, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%flow_velocity), dim_scale(1:rank))
        
        dataname = 'convection' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%convection, dim=2)
        dim_scale(2) = SIZE(hot_channel%convection, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%convection), dim_scale(1:rank))
        
        dataname = 'rhocoolant' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%rhocoolant, dim=2)
        dim_scale(2) = SIZE(hot_channel%rhocoolant, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%rhocoolant), dim_scale(1:rank))
        
        dataname = 'tcoolant' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tcoolant, dim=2)
        dim_scale(2) = SIZE(hot_channel%tcoolant, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tcoolant), dim_scale(1:rank))
        
        dataname = 'tclad_surf' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tclad_surf, dim=2)
        dim_scale(2) = SIZE(hot_channel%tclad_surf, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tclad_surf), dim_scale(1:rank))
        
        dataname = 'tclad_inner' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tclad_inner, dim=2)
        dim_scale(2) = SIZE(hot_channel%tclad_inner, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tclad_inner), dim_scale(1:rank))
        
        dataname = 'tfuel_surf' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_surf, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_surf, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tfuel_surf), dim_scale(1:rank))
        
        dataname = 'tfuel_center' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_center, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_center, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tfuel_center), dim_scale(1:rank))
        
        dataname = 'tfuel_avg' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_avg, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_avg, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(hot_channel%tfuel_avg), dim_scale(1:rank))
        
        call hdf5_close_group(typegroup)
        
        ! ----------------------------------------------------------------------
        ! thermal: average channel
        groupname = 'thermal_avg'
        call hdf5_open_group(casegroup, groupname, typegroup)
            
        dataname = 'flow_velocity' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%flow_velocity, dim=2)
        dim_scale(2) = SIZE(hot_channel%flow_velocity, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%flow_velocity), dim_scale(1:rank))
        
        dataname = 'convection' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%convection, dim=2)
        dim_scale(2) = SIZE(hot_channel%convection, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%convection), dim_scale(1:rank))
        
        dataname = 'rhocoolant' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%rhocoolant, dim=2)
        dim_scale(2) = SIZE(hot_channel%rhocoolant, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%rhocoolant), dim_scale(1:rank))
        
        dataname = 'tcoolant' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tcoolant, dim=2)
        dim_scale(2) = SIZE(hot_channel%tcoolant, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tcoolant), dim_scale(1:rank))
        
        dataname = 'tclad_surf' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tclad_surf, dim=2)
        dim_scale(2) = SIZE(hot_channel%tclad_surf, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tclad_surf), dim_scale(1:rank))
        
        dataname = 'tclad_inner' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tclad_inner, dim=2)
        dim_scale(2) = SIZE(hot_channel%tclad_inner, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tclad_inner), dim_scale(1:rank))
        
        dataname = 'tfuel_surf' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_surf, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_surf, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tfuel_surf), dim_scale(1:rank))
        
        dataname = 'tfuel_center' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_center, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_center, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tfuel_center), dim_scale(1:rank))
        
        dataname = 'tfuel_avg' 
        rank = 2
        dim_scale(1) = SIZE(hot_channel%tfuel_avg, dim=2)
        dim_scale(2) = SIZE(hot_channel%tfuel_avg, dim=1)
        call hdf5_write_data(typegroup, dataname, TRANSPOSE(avg_channel%tfuel_avg), dim_scale(1:rank))
        
        call hdf5_close_group(typegroup)
        
        ! ----------------------------------------------------------------------
        ! thermal: feedback parameter
        groupname = 'thermal_fdbk'
        call hdf5_open_group(casegroup, groupname, typegroup)
            
        dataname = 'Tf' 
        rank = 2
        dim_scale(1) = SIZE(self_fdbk%Tf%new, dim=1)
        dim_scale(2) = SIZE(self_fdbk%Tf%new, dim=2)
        call hdf5_write_data(typegroup, dataname, self_fdbk%Tf%new, dim_scale(1:rank))
        
        dataname = 'Tm' 
        rank = 2
        dim_scale(1) = SIZE(self_fdbk%Tm%new, dim=1)
        dim_scale(2) = SIZE(self_fdbk%Tm%new, dim=2)
        call hdf5_write_data(typegroup, dataname, self_fdbk%Tm%new, dim_scale(1:rank))
        
        dataname = 'Rho_m' 
        rank = 2
        dim_scale(1) = SIZE(self_fdbk%Rho_m%new, dim=1)
        dim_scale(2) = SIZE(self_fdbk%Rho_m%new, dim=2)
        call hdf5_write_data(typegroup, dataname, self_fdbk%Rho_m%new, dim_scale(1:rank))
        
        call hdf5_close_group(typegroup)
    
    end subroutine Print_thermal_hdf5
    
end module output_hdf5
