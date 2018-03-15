!$
!===================================================================================================
!
!   class for detector response calculation 
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          Detector  
!
!===================================================================================================
module detector_header 

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use timestep_header,            only : TimeStepInfo
    use state_header,               only : SteadyState
    use geometry_header,            only : Geometry, Meshing
    use contain_header,             only : GroupsFlux
    
    implicit none
    private
    public  :: Detector 

    ! --------------------------------------------------------------------------
    ! type for detector info, detValue = const * product(sigma*phi)
    type, private  :: detector_info_tp
        integer, public                              :: idx      = 1            ! detector index 
        integer, public                              :: ixs      = 1            ! response xsec 
        integer, public                              :: iz       = 1            ! position in radial
        integer, public                              :: ia       = 1            ! position in axial 
        real(KREAL), public                          :: const    = 0.0D0        ! multiply contants
        real(KREAL), public                          :: val      = 0.0D0        ! detector value 
        integer, public                              :: ir       = -1           ! nodal index 
        integer, public                              :: is       = -1           ! angular index 
    end type detector_info_tp                        
                                                     
    type  Detector                                   
        integer, public                              :: nxs      = 0            ! # of xsec defined 
        integer, public                              :: ndet     = 0            ! # of detector defined 
        real(KREAL), public, allocatable             :: detxs(:, :)
        type(detector_info_tp), public, allocatable  :: detpoint(:)
    contains
        procedure, public  :: alloc => Allocate_Detector
        procedure, public  :: clean => Free_Detector
        procedure, public  :: get => Get_Detector_Value
        procedure, public  :: phead => Print_Detector_head
        procedure, public  :: pvalue => Print_Detector_value
    end type Detector
    
contains
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Allocate_Detector (this, ns)
        
        class(Detector), intent(in out)   :: this
        type(SteadyState), intent(in)     :: ns 
        integer  :: i_allocate
        
        ! check status first
        call this%clean ()
        
        allocate(this%detxs(ns%state%ng, this%nxs), stat=i_allocate)
        allocate(this%detpoint(this%ndet), stat=i_allocate)
        
    end subroutine Allocate_Detector
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Free_Detector (this)
        
        class(Detector), intent(in out)  :: this
        
        if (allocated(this%detxs))          deallocate(this%detxs)
        if (allocated(this%detpoint))       deallocate(this%detpoint)
    
    end subroutine Free_Detector  
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Get_Detector_Value (this, ns, geom, mesh, flux)
        
        class(Detector), intent(in out)  :: this
        type(SteadyState), intent(in)    :: ns 
        type(Geometry), intent(in)       :: geom
        type(Meshing), intent(in)        :: mesh 
        type(GroupsFlux), intent(in)     :: flux 
        
        real(KREAL)  :: gflux(ns%state%zone, ns%state%layer, ns%state%ng)
        real(KREAL)  :: val 
        integer  :: iz, ir, ia, ig, idx, is 
        integer  :: i 
        
        gflux = 0.0D0
        do ia = 1, ns%state%layer 
            do ir = 1, ns%state%nodal 
                iz = mesh%zone(ir)
                do ig = 1, ns%state%ng
                    gflux(iz, ia, ig) = gflux(iz, ia, ig) + flux%ngs(ig)%scalar(ir,ia) * geom%area(ir)
                end do 
            end do 
        end do 
        do iz = 1, ns%state%zone
            gflux(iz, :, :) = gflux(iz, :, :) / geom%zone_area(iz)
        end do
        
        do i = 1, this%ndet 
            iz = this%detpoint(i)%iz
            ia = this%detpoint(i)%ia 
            idx = this%detpoint(i)%ixs 
            ir = this%detpoint(i)%ir
            is = this%detpoint(i)%is
            val = 0.0D0 
            
            if (ir <= 0)  then
                do ig = 1, ns%state%ng
                    val = val + gflux(iz, ia, ig) * this%detxs(ig, idx)
                end do 
            else if ((ir > 0) .AND. (is <= 0))  then
                do ig = 1, ns%state%ng
                    val = val + flux%ngs(ig)%scalar(ir,ia) * this%detxs(ig, idx)
                end do 
            else if ((ir > 0) .AND. (is > 0))  then
                do ig = 1, ns%state%ng
                    val = val + flux%ngs(ig)%angular(ir,ia,is) * this%detxs(ig, idx)
                end do 
            end if 
            
            this%detpoint(i)%val = val * this%detpoint(i)%const
        end do 
    
    end subroutine Get_Detector_Value
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_Detector_head (this, unit_, casename)
        
        class(Detector), intent(in)   :: this
        integer, intent(in)           :: unit_
        character(len=*), intent(in)  :: casename 
        
        integer  :: i 
        
        write(unit=unit_, fmt="(1x, A)") '(Begin)'
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
        write(unit=unit_, fmt="(1x, A, 1x)", advance='no') 'Case name is:'
        write(unit=unit_, fmt="('[', A, ']')") TRIM(casename)
        
        write(unit=unit_, fmt="(1x, A)") TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt="(1X, A, /)") 'Detector parameters per time step:'
        
        if (this%nxs <= 0)  return 
        if (this%ndet <= 0)  return 
        
        ! line-1
        write(unit=unit_, fmt='(1x, A)', advance='NO') '  step     time     |'
        do i = 1, this%ndet
            write(unit=unit_, fmt='(I3, A)', advance='NO') this%detpoint(i)%idx, ' #          |'
        end do 
        write(unit=unit_, fmt=*)
        
        ! line-2
        write(unit=unit_, fmt='(1x, A)', advance='NO') '                    |'
        do i = 1, this%ndet
            write(unit=unit_, fmt='(A, I4, A, I3, A)', advance='NO') 'iz=', this%detpoint(i)%iz, ' ia=', this%detpoint(i)%ia, ' |'
        end do 
        write(unit=unit_, fmt=*)
        
        ! line-3
        write(unit=unit_, fmt='(1x, A)', advance='NO') '                    |'
        do i = 1, this%ndet
            write(unit=unit_, fmt='(A, I4, A, I3, A)', advance='NO') 'ir=', this%detpoint(i)%ir, ' is=', this%detpoint(i)%is, ' |'
        end do 
        write(unit=unit_, fmt=*)
        
        ! line-4
        write(unit=unit_, fmt='(2x, A)', advance='NO') '---------------------'
        do i = 1, this%ndet
            write(unit=unit_, fmt='(A)', advance='NO') '----------------'
        end do
        write(unit=unit_, fmt=*)
        
    end subroutine Print_Detector_head
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_Detector_value (this, unit_, cidx, ctime)
        
        class(Detector), intent(in)     :: this
        integer, intent(in)             :: unit_ 
        integer, intent(in)             :: cidx
        real(KREAL), intent(in)         :: ctime 
        
        integer  :: i 
        
        if (this%nxs <= 0)  return 
        if (this%ndet <= 0)  return 
        
        write(unit=unit_, fmt='(1x, I4, TR3, ES12.5)', advance='NO') cidx, ctime
        do i = 1, this%ndet
            write(unit=unit_, fmt='(TR4, ES12.5)', advance='NO') this%detpoint(i)%val
        end do 
        write(unit=unit_, fmt=*)
    
    end subroutine Print_Detector_value
    
end module detector_header 
