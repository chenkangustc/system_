!$
!===================================================================================================
!
!   print distribution information for visualization
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Print_vtk_files
!
!   Public type lists:          No
!
!===================================================================================================
module output_visit
    
    use constants
    use global_state
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use files_dirs,         only : Remove_dir_in_path
    
    use global
    use th_global
	
	use imp_assm_global 
    
    implicit none 
    private
    public  :: Print_vtk_files, Print_vtk_geometry_2D, Print_vtk_distribution_nodal
    
    ! --------------------------------------------------------------------------
    ! sample file name:
    ! visual\forward_nodal.vtk
    ! visual\adjoint_nodal.vtk
    ! visual\_index_xx_nodal.vtk
    ! visual\nodal.visit
    ! --------------------------------------------------------------------------
    
    character(len=MAX_WORD_LEN), parameter  :: VTK_DIR              = 'visual/'                     ! directory for vtk file

    character(len=MAX_WORD_LEN), parameter  :: SUFFIX_NODAL_LAYER   = '_nodal_layer'
    character(len=MAX_WORD_LEN), parameter  :: SUFFIX_NODAL         = '_nodal'
    character(len=MAX_WORD_LEN), parameter  :: SUFFIX_ZONE          = '_zone'
    
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_FORWARD       = 'forward'
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_ADJOINT       = 'adjoint'
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_INDEX         = '_index_'

    character(len=MAX_WORD_LEN), parameter  :: EXTEND_VTK           = '.vtk'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_VISIT         = '.visit'

    ! type for file unit & name mapping
    type, private  :: FileUnit
        integer                      :: unit                                    ! file unit number
        character(len=MAX_WORD_LEN)  :: name                                    ! file name 
        character(len=MAX_WORD_LEN)  :: title                                   ! base title for vtk 
    end type FileUnit
    
    ! type for VisIt output
    type, public  :: VtkOutput
        type(FileUnit)              :: vtk_nodal_layer
        type(FileUnit)              :: vtk_nodal
        type(FileUnit)              :: vtk_zone
        
        type(FileUnit)              :: visit_nodal_layer
        type(FileUnit)              :: visit_nodal
        type(FileUnit)              :: visit_zone
    end type VtkOutput
    
    type(VtkOutput)  :: vtk_files
    logical, save    :: is_visit_exist          = .FALSE.                       ! NOTE: save
    
contains
    !$
    !===============================================================================================
    ! output vtk files
    !===============================================================================================
    subroutine Print_vtk_files (is_adjoint, is_transient, tidx, ctime)
    
        ! intent parameters
        logical, intent(in)  :: is_adjoint
        logical, intent(in)  :: is_transient
        integer, intent(in)  :: tidx
        real(KREAL), intent(in) :: ctime
        
        character(len=MAX_WORD_LEN)  :: tmp_char                                ! transfer integer to character
        character(len=MAX_WORD_LEN)  :: tmp_file                                ! tmp file name without directory in path
        integer  :: io_error
        
        character(len=MAX_WORD_LEN)  :: vtk_file
        character(len=MAX_WORD_LEN)  :: vtk_title
        integer  :: vtk_unit
        integer  :: visit_unit
        
        ! do not output
        if (.NOT. ns%output%is_vtk)  then
            return
        end if
        
        ! ----------------------------------------------------------------------
        ! set file name & title
        
        ! visit file
        vtk_files%visit_nodal_layer%name = TRIM(VTK_DIR) // TRIM(SUFFIX_NODAL_LAYER) // TRIM(EXTEND_VISIT)
        vtk_files%visit_nodal%name       = TRIM(VTK_DIR) // TRIM(SUFFIX_NODAL) // TRIM(EXTEND_VISIT)
        vtk_files%visit_zone%name        = TRIM(VTK_DIR) // TRIM(SUFFIX_ZONE) // TRIM(EXTEND_VISIT)
        
        ! vtk file for steady state forward
        if (.NOT. is_adjoint  .and.  .NOT. is_transient)  then
            vtk_files%vtk_nodal_layer%name = TRIM(VTK_DIR) // TRIM(PREFIX_FORWARD) // TRIM(SUFFIX_NODAL_LAYER) // TRIM(EXTEND_VTK)
            vtk_files%vtk_nodal%name       = TRIM(VTK_DIR) // TRIM(PREFIX_FORWARD) // TRIM(SUFFIX_NODAL) // TRIM(EXTEND_VTK)
            vtk_files%vtk_zone%name        = TRIM(VTK_DIR) // TRIM(PREFIX_FORWARD) // TRIM(SUFFIX_ZONE) // TRIM(EXTEND_VTK)
            
            vtk_files%vtk_nodal_layer%title = TRIM(SUFFIX_NODAL_LAYER)
            vtk_files%vtk_nodal%title       = TRIM(SUFFIX_NODAL)
            vtk_files%vtk_zone%title        = TRIM(SUFFIX_ZONE)
        end if
        
        ! vtk file for steady state adjoint
        if (is_adjoint )  then
            vtk_files%vtk_nodal_layer%name = TRIM(VTK_DIR) // TRIM(PREFIX_ADJOINT) // TRIM(SUFFIX_NODAL_LAYER) // TRIM(EXTEND_VTK)
            vtk_files%vtk_nodal%name       = TRIM(VTK_DIR) // TRIM(PREFIX_ADJOINT) // TRIM(SUFFIX_NODAL) // TRIM(EXTEND_VTK)
            vtk_files%vtk_zone%name        = TRIM(VTK_DIR) // TRIM(PREFIX_ADJOINT) // TRIM(SUFFIX_ZONE) // TRIM(EXTEND_VTK)
            
            vtk_files%vtk_nodal_layer%title = TRIM(SUFFIX_NODAL_LAYER)
            vtk_files%vtk_nodal%title       = TRIM(SUFFIX_NODAL)
            vtk_files%vtk_zone%title        = TRIM(SUFFIX_ZONE)
        end if
        
        ! vtk file for transient state
        if (is_transient )  then
            write(unit=tmp_char, fmt=*)  tidx
            
            vtk_files%vtk_nodal_layer%name = TRIM(VTK_DIR) // TRIM(PREFIX_INDEX) // TRIM(ADJUSTL(tmp_char)) // TRIM(SUFFIX_NODAL_LAYER) // TRIM(EXTEND_VTK)
            vtk_files%vtk_nodal%name       = TRIM(VTK_DIR) // TRIM(PREFIX_INDEX) // TRIM(ADJUSTL(tmp_char)) // TRIM(SUFFIX_NODAL) // TRIM(EXTEND_VTK)
            vtk_files%vtk_zone%name        = TRIM(VTK_DIR) // TRIM(PREFIX_INDEX) // TRIM(ADJUSTL(tmp_char)) // TRIM(SUFFIX_ZONE) // TRIM(EXTEND_VTK)
            
            vtk_files%vtk_nodal_layer%title = TRIM(SUFFIX_NODAL_LAYER)
            vtk_files%vtk_nodal%title       = TRIM(SUFFIX_NODAL)
            vtk_files%vtk_zone%title        = TRIM(SUFFIX_ZONE)
        end if
        
        ! ----------------------------------------------------------------------
        ! open file
        open (status='replace', action='write', iostat=io_error, newunit=vtk_files%vtk_nodal_layer%unit, file=vtk_files%vtk_nodal_layer%name)
        open (status='replace', action='write', iostat=io_error, newunit=vtk_files%vtk_nodal%unit, file=vtk_files%vtk_nodal%name)
        open (status='replace', action='write', iostat=io_error, newunit=vtk_files%vtk_zone%unit, file=vtk_files%vtk_zone%name)
        
        if (is_transient )  then
            if (.NOT. is_visit_exist)  then
                open (status='replace', action='write', iostat=io_error, newunit=vtk_files%visit_nodal_layer%unit, file=vtk_files%visit_nodal_layer%name)
                open (status='replace', action='write', iostat=io_error, newunit=vtk_files%visit_nodal%unit, file=vtk_files%visit_nodal%name)
                open (status='replace', action='write', iostat=io_error, newunit=vtk_files%visit_zone%unit, file=vtk_files%visit_zone%name)
                
                is_visit_exist = .TRUE.
            else 
                open (status='old', action='readwrite', position='append', iostat=io_error, newunit=vtk_files%visit_nodal_layer%unit, file=vtk_files%visit_nodal_layer%name)
                open (status='old', action='readwrite', position='append', iostat=io_error, newunit=vtk_files%visit_nodal%unit, file=vtk_files%visit_nodal%name)
                open (status='old', action='readwrite', position='append', iostat=io_error, newunit=vtk_files%visit_zone%unit, file=vtk_files%visit_zone%name)
            end if
        end if
        
        ! ----------------------------------------------------------------------
        ! write file
        ! for steady state forward
        if (.NOT. is_adjoint  .and.  .NOT. is_transient)  then
            vtk_unit = vtk_files%vtk_nodal_layer%unit
            vtk_title = vtk_files%vtk_nodal_layer%title
            call Print_vtk_geometry_3D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal_layer (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            
            vtk_unit = vtk_files%vtk_nodal%unit
            vtk_title = vtk_files%vtk_nodal%title
            call Print_vtk_geometry_2D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            
            vtk_unit = vtk_files%vtk_zone%unit
            vtk_title = vtk_files%vtk_zone%title
            call Print_vtk_geometry_zone (vtk_unit, vtk_title)
            call Print_vtk_distribution_zone (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
        end if
        
        ! for steady state adjoint
        if (is_adjoint )  then
            vtk_unit = vtk_files%vtk_nodal_layer%unit
            vtk_title = vtk_files%vtk_nodal_layer%title
            call Print_vtk_geometry_3D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal_layer (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            
            vtk_unit = vtk_files%vtk_nodal%unit
            vtk_title = vtk_files%vtk_nodal%title
            call Print_vtk_geometry_2D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            
            vtk_unit = vtk_files%vtk_zone%unit
            vtk_title = vtk_files%vtk_zone%title
            call Print_vtk_geometry_zone (vtk_unit, vtk_title)
            call Print_vtk_distribution_zone (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
        end if
        
        ! for transient state
        if (is_transient )  then
            vtk_unit = vtk_files%vtk_nodal_layer%unit
            vtk_file = vtk_files%vtk_nodal_layer%name
            vtk_title = vtk_files%vtk_nodal_layer%title
            call Print_vtk_geometry_3D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal_layer (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            visit_unit = vtk_files%visit_nodal_layer%unit
            call Remove_dir_in_path (vtk_file, tmp_file)
            write(unit = visit_unit, fmt="(A)")  TRIM(ADJUSTL(tmp_file))
            
            vtk_unit = vtk_files%vtk_nodal%unit
            vtk_file = vtk_files%vtk_nodal%name
            vtk_title = vtk_files%vtk_nodal%title
            call Print_vtk_geometry_2D (vtk_unit, vtk_title)
            call Print_vtk_distribution_nodal (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            visit_unit = vtk_files%visit_nodal%unit
            call Remove_dir_in_path (vtk_file, tmp_file)
            write(unit = visit_unit, fmt="(A)")  TRIM(ADJUSTL(tmp_file))
            
            vtk_unit = vtk_files%vtk_zone%unit
            vtk_file = vtk_files%vtk_zone%name
            vtk_title = vtk_files%vtk_zone%title
            call Print_vtk_geometry_zone (vtk_unit, vtk_title)
            call Print_vtk_distribution_zone (vtk_unit, vtk_title, is_adjoint=is_adjoint, is_transient=is_transient)
            visit_unit = vtk_files%visit_zone%unit
            call Remove_dir_in_path (vtk_file, tmp_file)
            write(unit = visit_unit, fmt="(A)")  TRIM(ADJUSTL(tmp_file))
        end if
        
        ! ----------------------------------------------------------------------
        ! close file
        close (iostat=io_error, status='keep', unit = vtk_files%vtk_nodal_layer%unit)
        close (iostat=io_error, status='keep', unit = vtk_files%vtk_nodal%unit)
        close (iostat=io_error, status='keep', unit = vtk_files%vtk_zone%unit)
        
        if (is_visit_exist )  then
            close (iostat=io_error, status='keep', unit = vtk_files%visit_nodal_layer%unit)
            close (iostat=io_error, status='keep', unit = vtk_files%visit_nodal%unit)
            close (iostat=io_error, status='keep', unit = vtk_files%visit_zone%unit)
        end if
    
    end subroutine Print_vtk_files
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! print 1D geometry information to vtk file per layer
    !===============================================================================================
    subroutine Print_vtk_geometry_zone (unit_, vtk_title)
        
        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        
        ! local
        character(len=MAX_LINE_LEN)  :: string_1, string_2
        integer  :: cell
        integer  :: ir, j
        
        ! ----------------------------------------------------------------------
        ! print title
        write(unit=unit_, fmt="(A)")   '# vtk DataFile Version 3.0'
        write(unit=unit_, fmt="(A)")   '2D distribution for '//TRIM(vtk_title)
        write(unit=unit_, fmt="(A)")   'ASCII'
        write(unit=unit_, fmt="(A)")   'DATASET UNSTRUCTURED_GRID'
        
        ! print point coordinate
        write(unit=unit_, fmt="(A, 3x, I10, 3x, A)")    'POINTS', ns%state%point+1, 'float'
        
        write(unit=unit_, fmt=*) 0.0_KREAL, 0.0_KREAL, 0.0_KREAL
        do ir = 1, ns%state%point
            write(unit=unit_, fmt=*)  geom%coordinate(1,ir), geom%coordinate(2,ir),  0.0_KREAL
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell information
        cell = mesh_vtk%total_zone + mesh_vtk%total_point
        write(unit=unit_, fmt="(A, 3x, I10, 3x, I10)")    'CELLS', mesh_vtk%total_zone, cell
        
        do ir = 1, ns%state%zone
            if (mesh_vtk%mapping(ir) /= 0)  then
                write(unit=unit_, fmt="(1x, I10, TR1, *(I4, TR1))") mesh_vtk%zones(ir)%point, (mesh_vtk%zones(ir)%coordinates(j)%index, j=1, SIZE(mesh_vtk%zones(ir)%coordinates))
            end if
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell type
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_TYPES', mesh_vtk%total_zone
        
        do ir = 1, ns%state%zone
            if (mesh_vtk%mapping(ir) /= 0)  then
                write(unit=unit_, fmt=*)  mesh_vtk%zones(ir)%mesh_type
            end if
        end do
        write(unit=unit_, fmt=*)  '    '
        
    end subroutine Print_vtk_geometry_zone
    
    !$
    !===============================================================================================
    ! print 2D geometry information to vtk file
    !===============================================================================================
    subroutine Print_vtk_geometry_2D (unit_, vtk_title)
        
        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        integer  :: ir
        
        ! ----------------------------------------------------------------------
        ! print title
        write(unit=unit_, fmt="(A)")   '# vtk DataFile Version 3.0'
        write(unit=unit_, fmt="(A)")   '2D distribution for '//TRIM(vtk_title)
        write(unit=unit_, fmt="(A)")   'ASCII'
        write(unit=unit_, fmt="(A)")   'DATASET UNSTRUCTURED_GRID'
        
        ! print point coordinate
        write(unit=unit_, fmt="(A, 3x, I10, 3x, A)")    'POINTS', ns%state%point+1, 'float'
        
        write(unit=unit_, fmt=*) 0.0_KREAL, 0.0_KREAL, 0.0_KREAL
        do ir = 1, ns%state%point
            write(unit=unit_, fmt=*)  geom%coordinate(1,ir), geom%coordinate(2,ir),  0.0_KREAL
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell information
        write(unit=unit_, fmt="(A, 3x, I10, 3x, I10)")    'CELLS', ns%state%nodal, 4*ns%state%nodal
        
        do ir = 1, ns%state%nodal
            write(unit=unit_, fmt=*)  3, mesh%point(1,ir), mesh%point(2,ir), mesh%point(3,ir)
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell type
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_TYPES', ns%state%nodal
        
        do ir = 1, ns%state%nodal
            write(unit=unit_, fmt=*)  5
        end do
        write(unit=unit_, fmt=*)  '    '
        
    end subroutine Print_vtk_geometry_2D
    
    !$
    !===============================================================================================
    ! print 3D geometry information to vtk file
    !===============================================================================================
    subroutine Print_vtk_geometry_3D (unit_, vtk_title)

        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        
        integer  :: bottom
        integer  :: top
        integer  :: ir, ia
        
        ! ----------------------------------------------------------------------
        ! print title
        write(unit=unit_, fmt="(A)")   '# vtk DataFile Version 3.0'
        write(unit=unit_, fmt="(A)")   '3D distribution for '//TRIM(vtk_title)
        write(unit=unit_, fmt="(A)")   'ASCII'
        write(unit=unit_, fmt="(A)")   'DATASET UNSTRUCTURED_GRID'
        
        ! print point coordinate
        write(unit=unit_, fmt="(A, 3x, I10, 3x, A)")    'POINTS', ns%state%point*(ns%state%layer+1)+1, 'float'
        
        write(unit=unit_, fmt=*) 0.0_KREAL, 0.0_KREAL, 0.0_KREAL
        do ia = 0, ns%state%layer
            do ir = 1, ns%state%point
                if (ia == 0)  then
                    write(unit=unit_, fmt=*)  geom%coordinate(1,ir), geom%coordinate(2,ir),  0.0_KREAL
                else
                    write(unit=unit_, fmt=*)  geom%coordinate(1,ir), geom%coordinate(2,ir),  SUM(geom%height(1:ia))
                end if
            end do
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell information
        write(unit=unit_, fmt="(A, 3x, I10, 3x, I10)")    'CELLS', ns%state%nodal*ns%state%layer, 7*ns%state%nodal*ns%state%layer
        
        do ia = 1, ns%state%layer
            bottom = (ia-1)*ns%state%point
            top = (ia)*ns%state%point
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*)  6, bottom+mesh%point(1,ir), bottom+mesh%point(2,ir), bottom+mesh%point(3,ir)  &
                                                &   , top+mesh%point(1,ir), top+mesh%point(2,ir), top+mesh%point(3,ir)
            end do
        end do
        write(unit=unit_, fmt=*)  '    '
        
        ! print cell type
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_TYPES', ns%state%nodal*ns%state%layer
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*)  13
            end do
        end do
        write(unit=unit_, fmt=*)  '    '
        
    end subroutine Print_vtk_geometry_3D
    
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    
    !$
    !===============================================================================================
    ! print 2D distribution information to vtk file per nodal
    !===============================================================================================
    subroutine Print_vtk_distribution_nodal (unit_, vtk_title, is_adjoint, is_transient)

        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        logical, intent(in)  :: is_adjoint
        logical, intent(in)  :: is_transient
        
        real(KREAL)        :: vtk_nodal(ns%state%nodal)
        character(len=MAX_WORD_LEN)  ::  true_title
        character(len=MAX_WORD_LEN)  ::  tmp_char
        integer                  :: ir, ip
        
        ! ----------------------------------------------------------------------
        ! print zone
        true_title = TRIM(vtk_title)//'_zone'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'int'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        ! reversing order
        do ir = 1, ns%state%nodal
            write(unit=unit_, fmt=*)  ns%state%zone + 1 - mesh%zone(ir)
        end do
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print flux
        true_title = TRIM(vtk_title)//'_flux'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        call dist_flux%nodal(mesh, geom, is_function=.FALSE., nodal=vtk_nodal)
        
        do ir = 1, ns%state%nodal
            write(unit=unit_, fmt=*)  vtk_nodal(ir)
        end do
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print power
        if (.NOT. is_adjoint)  then
            true_title = TRIM(vtk_title)//'_power'
            
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            call dist_power%nodal(mesh, geom, is_function=.FALSE., nodal=vtk_nodal)
            
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*)  vtk_nodal(ir)
            end do
            write(unit=unit_, fmt=*)  '    '
            
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
        end if
        
        ! ----------------------------------------------------------------------
        ! print fission rate
        if (.NOT. is_adjoint)  then
            true_title = TRIM(vtk_title)//'_fission_rate'
            
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            call dist_fission_rate%nodal(mesh, geom, is_function=.FALSE., nodal=vtk_nodal)
            
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*)  vtk_nodal(ir)
            end do
            write(unit=unit_, fmt=*)  '    '
            
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
        end if
        
        ! ----------------------------------------------------------------------
        ! print precursor
        if (is_transient )  then
            do ip = 1, nt%state%dg
                write(unit=tmp_char, fmt=*)  ip
                write(unit = true_title, fmt=*) TRIM(vtk_title), '_precursor_group_', TRIM(ADJUSTL(tmp_char))
                
                write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal
                write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
                write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
                
                call dist_dnps(ip)%nodal(mesh, geom, is_function=.FALSE., nodal=vtk_nodal)
                
                do ir = 1, ns%state%nodal
                    write(unit=unit_, fmt=*)  vtk_nodal(ir)
                end do
                write(unit=unit_, fmt=*)  '    '
                
                write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
                write(unit=unit_, fmt=*)  '    '
            end do
        end if
        
    end subroutine Print_vtk_distribution_nodal
    
    !$
    !===============================================================================================
    ! print 3D distribution information to vtk file per nodal per layer
    !===============================================================================================
    subroutine Print_vtk_distribution_nodal_layer (unit_, vtk_title, is_adjoint, is_transient)

        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        logical, intent(in)  :: is_adjoint
        logical, intent(in)  :: is_transient
        
        character(len=MAX_WORD_LEN)  :: true_title
        character(len=MAX_WORD_LEN)  :: tmp_char
        integer                  :: ir, ia, ip, iz
        
        ! ----------------------------------------------------------------------
        ! print material
        true_title = TRIM(vtk_title)//'_material'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'int'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*)  mat_info%loading(mesh%zone(ir), ia)
            end do
        end do
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print flux
        true_title = TRIM(vtk_title)//'_flux'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        do ia = 1, ns%state%layer
            do ir = 1, ns%state%nodal
                write(unit=unit_, fmt=*) dist_flux%matrix(ir,ia)
            end do
        end do
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
        ! ----------------------------------------------------------------------
        ! print power
        if (.NOT. is_adjoint)  then
            true_title = TRIM(vtk_title)//'_power'
            
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    write(unit=unit_, fmt=*) dist_power%matrix(ir,ia)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
        end if
        
        ! ----------------------------------------------------------------------
        ! print fission rate
        if (.NOT. is_adjoint)  then
            true_title = TRIM(vtk_title)//'_fission_rate'
            
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    write(unit=unit_, fmt=*) dist_fission_rate%matrix(ir,ia)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
        end if
        
!        ! ----------------------------------------------------------------------
!        ! print precursor
!        if (is_transient )  then
!            do ip = 1, nt%state%dg
!                write(unit=tmp_char, fmt=*)  ip
!                write(unit = true_title, fmt=*) TRIM(vtk_title), '_precursor_group_', TRIM(ADJUSTL(tmp_char))
!                
!                write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
!                write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
!                write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
!                
!                do ia = 1, ns%state%layer
!                    do ir = 1, ns%state%nodal
!                        write(unit=unit_, fmt=*) dist_dnps(ip)%matrix(ir,ia)
!                    end do
!                end do
!                write(unit=unit_, fmt=*)  '    '
!                
!                write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
!                write(unit=unit_, fmt=*)  '    '
!            
!            end do
!        end if
        
        ! ----------------------------------------------------------------------
        ! thermal information
        if (ns%feedback%is_feedback .and. ns%feedback%is_inner)  then
            ! hot-channel
            true_title = TRIM(vtk_title)//'_hot_'//'flow_velocity'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%flow_velocity(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'convection'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%convection(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'rhocoolant'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%rhocoolant(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tcoolant'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tcoolant(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tclad_surf'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tclad_surf(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tclad_inner'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tclad_inner(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tfuel_surf'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tfuel_surf(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tfuel_center'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tfuel_center(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_hot_'//'tfuel_avg'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) hot_channel%tfuel_avg(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            ! average channel
            true_title = TRIM(vtk_title)//'_avg_'//'flow_velocity'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) avg_channel%flow_velocity(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'convection'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) avg_channel%convection(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'rhocoolant'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) avg_channel%rhocoolant(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tcoolant'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    !write(unit=unit_, fmt=*) avg_channel%tcoolant(ia, iz)
					write(unit=unit_, fmt=*) assm1(iz)%thermal%Tcoolant(ia)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tclad_surf'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    !write(unit=unit_, fmt=*) avg_channel%tclad_surf(ia, iz)
					write(unit=unit_, fmt=*) assm1(iz)%thermal%Tsc(ia)
				end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tclad_inner'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    !write(unit=unit_, fmt=*) avg_channel%tclad_inner(ia, iz)
					write(unit=unit_, fmt=*) assm1(iz)%thermal%Tgs(ia)
				end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tfuel_surf'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    write(unit=unit_, fmt=*) avg_channel%tfuel_surf(ia, iz)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tfuel_center'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    !write(unit=unit_, fmt=*) avg_channel%tfuel_center(ia, iz)
					write(unit=unit_, fmt=*) assm1(iz)%thermal%tfuel_center(ia)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
            
            true_title = TRIM(vtk_title)//'_avg_'//'tfuel_avg'
            write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', ns%state%nodal * ns%state%layer
            write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
            write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
            
            do ia = 1, ns%state%layer
                do ir = 1, ns%state%nodal
                    iz = mesh%zone(ir)
                    !write(unit=unit_, fmt=*) avg_channel%tfuel_avg(ia, iz)
					write(unit=unit_, fmt=*) assm1(iz)%thermal%Tfuel(ia)
                end do
            end do
            write(unit=unit_, fmt=*)  '    '
            write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
            write(unit=unit_, fmt=*)  '    '
        end if
        
    end subroutine Print_vtk_distribution_nodal_layer
    
    !$
    !===============================================================================================
    ! print 1D distribution information to vtk file per layer
    !===============================================================================================
    subroutine Print_vtk_distribution_zone (unit_, vtk_title, is_adjoint, is_transient)

        ! intent parameters
        integer, intent(in)  :: unit_
        character(len=*), intent(in)  :: vtk_title
        logical, intent(in)  :: is_adjoint
        logical, intent(in)  :: is_transient
        
        real(KREAL)        :: vtk_zone(ns%state%zone)
        real(KREAL)        :: active_factor
        character(len=MAX_WORD_LEN)  ::  true_title
        character(len=MAX_WORD_LEN)  ::  tmp_char
        integer                  :: ir, ip, iz
        
        active_factor = REAL(ns%state%layer) / REAL(ns%state%layer - ns%state%layer_top - ns%state%layer_bottom)
        
        ! ----------------------------------------------------------------------
        ! print flux
        true_title = TRIM(vtk_title)//'_flux'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', mesh_vtk%total_zone
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        call dist_flux%zone (mesh, geom, is_function=.FALSE., zone=vtk_zone)
        
        do iz = 1, ns%state%zone
            if (mesh_vtk%mapping(iz) /= 0)  then
                write(unit=unit_, fmt=*)  vtk_zone(iz)
            end if
        end do
        
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
        ! print power
        true_title = TRIM(vtk_title)//'_power'
        
        write(unit=unit_, fmt="(A, 3x, I10)")    'CELL_DATA', mesh_vtk%total_zone
        write(unit=unit_, fmt="(A, 3x, A, 3x, A)") 'SCALARS', TRIM(true_title), 'float'
        write(unit=unit_, fmt="(A, 3x, A)") 'LOOKUP_TABLE', TRIM(true_title)//'_table'
        
        call dist_power%zone (mesh, geom, is_function=.FALSE., zone=vtk_zone)
        vtk_zone = vtk_zone * active_factor
        
        do iz = 1, ns%state%zone
            if (mesh_vtk%mapping(iz) /= 0)  then
                write(unit=unit_, fmt=*)  vtk_zone(iz)
            end if
        end do
        
        write(unit=unit_, fmt=*)  '    '
        
        write(unit=unit_, fmt="(A)")  'POINT_DATA 1'
        write(unit=unit_, fmt=*)  '    '
        
    end subroutine Print_vtk_distribution_zone
        
end module output_visit
