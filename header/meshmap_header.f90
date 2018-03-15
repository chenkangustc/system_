!$
!===================================================================================================
!
!   class for mesh mapping
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          AxialMap
!
!===================================================================================================
module meshmap_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: AxialMap
    
    ! --------------------------------------------------------------------------
    ! type for axial mesh mapping, by volume weight
    type  AxialMap
        integer, public                               :: mesh_nk
        integer, public                               :: mesh_th
        real(KREAL), public, allocatable              :: h_nk(:)
        real(KREAL), public, allocatable              :: h_th(:)
        real(KREAL), public, allocatable              :: wtnk2th(:, :)          ! weight nk belongs th
        real(KREAL), public, allocatable              :: wtth2nk(:, :)          ! weight th belongs nk
        logical, public                               :: is_set  = .FALSE.
    contains
        procedure, public  :: setnk => Set_AxialMap_NK
        procedure, public  :: setth => Set_AxialMap_TH
        procedure, public  :: setmap => Set_AxialMap_mapinfo
        procedure, public  :: print => Print_AxialMap
        procedure, public  :: clean => Free_AxialMap
        procedure, public  :: nk2th => Convert_AxialMap_NK2TH
        procedure, public  :: th2nk => Convert_AxialMap_TH2NK
    end type AxialMap
        
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Set_AxialMap_NK, Set_AxialMap_TH, Set_AxialMap_mapinfo, Print_AxialMap, Free_AxialMap
    private  :: Convert_AxialMap_NK2TH, Convert_AxialMap_TH2NK
        
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_AxialMap_NK (this, mesh_nk, h_nk)
        
        class(AxialMap), intent(in out)  :: this
        integer, intent(in)              :: mesh_nk
        real(KREAL), intent(in)          :: h_nk(:)
        
        integer  :: i_allocate
        
        this%mesh_nk = mesh_nk
        
        if (allocated(this%h_nk))       deallocate(this%h_nk)
        allocate(this%h_nk(this%mesh_nk), stat=i_allocate)
        this%h_nk = h_nk
    
    end subroutine Set_AxialMap_NK
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_AxialMap_TH (this, mesh_th, h_th)
        
        class(AxialMap), intent(in out)  :: this
        integer, intent(in)              :: mesh_th
        real(KREAL), intent(in)          :: h_th(:)
        
        integer  :: i_allocate
        
        this%mesh_th = mesh_th
        
        if (allocated(this%h_th))       deallocate(this%h_th)
        allocate(this%h_th(this%mesh_th), stat=i_allocate)
        this%h_th = h_th
    
    end subroutine Set_AxialMap_TH
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_AxialMap_mapinfo (this)
        
        class(AxialMap), intent(in out)  :: this
        
        real(KREAL)  :: Uth, Unk
        real(KREAL)  :: Lth, Lnk
        real(KREAL)  :: wt
        integer  :: ink, ith
        integer  :: i_allocate
        logical  :: is_seperate
        
        if (.NOT. this%is_set)  then
            allocate(this%wtnk2th(this%mesh_nk, this%mesh_th), stat=i_allocate)
            allocate(this%wtth2nk(this%mesh_th, this%mesh_nk), stat=i_allocate)
            this%wtnk2th = REAL_ZERO
            this%wtth2nk = REAL_ZERO
            
            do ith = 1, this%mesh_th
                Uth = SUM(this%h_th(1: ith))
                Lth = Uth - this%h_th(ith)
                
                do ink = 1, this%mesh_nk
                    Unk = SUM(this%h_nk(1: ink))
                    Lnk = Unk - this%h_nk(ink)
                    
                    is_seperate = (Unk <= Lth) .OR. (Lnk >= Uth)
                    if (.NOT. is_seperate)  then
                        wt = (MIN(Uth, Unk) - MAX(Lth, Lnk)) / (Uth - Lth)
                        this%wtnk2th(ink, ith) = wt
                    end if
                end do
            end do
            
            do ink = 1, this%mesh_nk
                Unk = SUM(this%h_nk(1: ink))
                Lnk = Unk - this%h_nk(ink)
                
                do ith = 1, this%mesh_th
                    Uth = SUM(this%h_th(1: ith))
                    Lth = Uth - this%h_th(ith)
                    
                    is_seperate = (Unk <= Lth) .OR. (Lnk >= Uth)
                    if (.NOT. is_seperate)  then
                        wt = (MIN(Uth, Unk) - MAX(Lth, Lnk)) / (Unk - Lnk)
                        this%wtth2nk(ith, ink) = wt
                    end if
                end do
            end do
            
            this%is_set = .TRUE.
        end if
        
    end subroutine Set_AxialMap_mapinfo
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_AxialMap (this, unit_)
        
        class(AxialMap), intent(in)   :: this
        integer, intent(in)           :: unit_
        
        integer  :: i, j
        
        do j = 1, SIZE(this%wtnk2th, dim=2)
            write(unit=unit_, fmt="(1x, *(ES12.5, TR3))")  this%wtnk2th(:, j)
        end do
        
        do j = 1, SIZE(this%wtth2nk, dim=2)
            write(unit=unit_+1, fmt="(1x, *(ES12.5, TR3))")  this%wtth2nk(:, j)
        end do
    
    end subroutine Print_AxialMap
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_AxialMap (this)
        
        class(AxialMap), intent(in out)   :: this
        
        if (allocated(this%h_nk))           deallocate(this%h_nk)
        if (allocated(this%h_th))           deallocate(this%h_th)
        if (allocated(this%wtnk2th))        deallocate(this%wtnk2th)
        if (allocated(this%wtth2nk))        deallocate(this%wtth2nk)
    
    end subroutine Free_AxialMap
    
    !$
    !===============================================================================================
    ! dim=2 --> Axial
    !===============================================================================================
    subroutine Convert_AxialMap_NK2TH (this, nkin, thout)
        
        class(AxialMap), intent(in)   :: this
        real(KREAL), intent(in out)   :: nkin(:, :)
        real(KREAL), intent(in out)   :: thout(:, :)
        
        integer  :: ink, ith
        
        do ith = 1, this%mesh_th
            thout(:, ith) = REAL_ZERO
            
            do ink = 1, this%mesh_nk
                thout(:, ith) = thout(:, ith) + nkin(:, ink) * this%wtnk2th(ink, ith)
            end do
        end do
    
    end subroutine Convert_AxialMap_NK2TH
    
    !$
    !===============================================================================================
    ! dim=2 --> Axial
    !===============================================================================================
    subroutine Convert_AxialMap_TH2NK (this, thin, nkout)
    
        class(AxialMap), intent(in)   :: this
        real(KREAL), intent(in out)   :: thin(:, :)
        real(KREAL), intent(in out)   :: nkout(:, :)
        
        integer  :: ink, ith
        
        do ink = 1, this%mesh_nk
            nkout(:, ink) = REAL_ZERO
            
            do ith = 1, this%mesh_th
                nkout(:, ink) = nkout(:, ink) + thin(:, ith) * this%wtth2nk(ith, ink)
            end do
        end do
    
    end subroutine Convert_AxialMap_TH2NK
        
end module meshmap_header
