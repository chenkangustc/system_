!$
!===================================================================================================
!
!   module for power distribution from transport result
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          LinearPower
!
!===================================================================================================
module gth_power_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use gth_geometry_header,         only : ThermalScale, ThermalGeometry, ThermalAssemblyGeometry
    
    implicit none
    private 
    public  :: LinearPower
    
    ! --------------------------------------------------------------------------
    ! type for linear power density distribution
    type  LinearPower
        real(KREAL), public, allocatable      :: power(:, :)                    ! (Watt)--power level per th-nodal
        real(KREAL), public, allocatable      :: fq_lattice(:)                  ! shape factor from assembly calculation
        real(KREAL), public, allocatable      :: fq_core(:, :)                  ! shape factor from core result
                                              
        real(KREAL), public                   :: normal_old                     ! old value of power amplitude 
        real(KREAL), public                   :: normal_current                 ! new value of power amplitude
                                              
        real(KREAL), public, allocatable      :: avg_linear(:, :)               ! (W/m) --average linear power of fuel pin (absolutely, not normalized)
        real(KREAL), public, allocatable      :: max_linear(:, :)               ! (W/m) --max linear power of fuel pin (absolutely, not normalized)
        real(KREAL), public                   :: gamma                          ! (%) --fraction of coolant direct heating
    contains
        procedure, public  :: alloc => Alloc_LinearPower
        procedure, public  :: clean => Free_LinearPower
        procedure, public  :: print => Print_LinearPower
        procedure, public  :: set_fq_lattice => Set_fq_lattice_info
        procedure, public  :: fq_FA => Get_assembly_fq
        procedure, public  :: fq_pin => Get_fuelpin_fq
        procedure, public  :: fq_point => Get_point_fq
        generic, public    :: set_power => Set_assembly_power_3D, Set_assembly_power_0D
        procedure          :: Set_assembly_power_3D
        procedure          :: Set_assembly_power_0D
    end type LinearPower
    
    ! --------------------------------------------------------------------------
    ! private the real function name
    private  :: Alloc_LinearPower, Free_LinearPower, Print_LinearPower, Set_fq_lattice_info
    private  :: Get_assembly_fq, Get_fuelpin_fq, Get_point_fq
    private  :: Set_assembly_power_3D, Set_assembly_power_0D
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_LinearPower (this, nth)
    
        class(LinearPower), intent(in out)  :: this
        type(ThermalScale), intent(in)      :: nth
        integer  :: i_allocate
        
        ! check allocated first
        call this%clean ()
        
        allocate(this%power(nth%na, nth%nr), stat=i_allocate)
        allocate(this%fq_lattice(nth%nr), stat=i_allocate)
        allocate(this%fq_core(nth%na, nth%nr), stat=i_allocate)

        allocate(this%avg_linear(nth%na, nth%nr), stat=i_allocate)
        allocate(this%max_linear(nth%na, nth%nr), stat=i_allocate)
        
        this%power       = 0.0D0
        this%fq_lattice  = 1.0D0                                              ! this should be further developped
        this%fq_core     = 1.0D0
        
        this%avg_linear  = 0.0D0
        this%max_linear  = 0.0D0
        
        this%normal_old     = 1.0D0
        this%normal_current = 1.0D0
        this%gamma          = 0.0D0
    
    end subroutine Alloc_LinearPower
    
    !$
    !===============================================================================================
    ! finalizer for class of LinearPower
    !===============================================================================================
    subroutine Free_LinearPower (this)
        
        class(LinearPower), intent(in out)  :: this
        
        if (allocated(this%power))              deallocate(this%power)
        if (allocated(this%fq_lattice))         deallocate(this%fq_lattice)
        if (allocated(this%fq_core))            deallocate(this%fq_core)
        
        if (allocated(this%avg_linear))         deallocate(this%avg_linear)
        if (allocated(this%max_linear))         deallocate(this%max_linear)
    
    end subroutine Free_LinearPower
        
    !$
    !===============================================================================================
    ! set current power level per assembly
    !===============================================================================================
    subroutine Print_LinearPower (this, to_avg, to_max)
        
        class(LinearPower), intent(in out)      :: this
        integer, intent(in)                     :: to_avg
        integer, intent(in), optional           :: to_max
        
        ! local variables
        integer  :: unit_
        integer  :: i, j
        
        unit_ = to_avg
        do i = 1, SIZE(this%avg_linear, dim=1)
            write(unit=unit_, fmt="(1x, I4, TR3, *(ES13.6, TR3))")  i, this%avg_linear(i, :)
        end do
        
        if (PRESENT(to_max))  then
            unit_ = to_max
        end if
        do i = 1, SIZE(this%avg_linear, dim=1)
            write(unit=unit_, fmt="(1x, I4, TR3, *(ES13.6, TR3))")  i, this%max_linear(i, :)
        end do
        
    end subroutine Print_LinearPower
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_fq_lattice_info (this, fq)
        
        class(LinearPower), intent(in out)    :: this
        real(KREAL), intent(in)               :: fq(:)
    
        this%fq_lattice = fq
        
    end subroutine Set_fq_lattice_info
    
    !$
    !===============================================================================================
    ! get fq factor by assembly average power density
    !===============================================================================================
    subroutine Get_assembly_fq (this, fq, channel)
        
        class(LinearPower), intent(in)   :: this
        real(KREAL), intent(in out)      :: fq
        integer, intent(in out)          :: channel
        
        ! local
        real(KREAL)  :: max_power    = 0.0D0
        real(KREAL)  :: total_power  = 0.0D0
        real(KREAL)  :: local_avg    = 0.0D0
        integer          :: count  = 0
        integer          :: index  = 0
        integer          :: ir
        
        do ir = 1, SIZE(this%avg_linear, dim=2)
            local_avg = SUM(this%avg_linear(:, ir))
            
            if (local_avg >= EPS_ZERO)  then
                count = count + 1
                total_power = total_power + local_avg
                
                if (local_avg > max_power)  then
                    max_power = local_avg
                    index = ir
                end if
            end if
        end do
        
        fq = max_power / (total_power / count)
        channel = index
    
    end subroutine Get_assembly_fq
    
    !$
    !===============================================================================================
    ! get fq factor by fuel pin power density
    !===============================================================================================
    subroutine Get_fuelpin_fq (this, fq, channel)
        
        class(LinearPower), intent(in)   :: this
        real(KREAL), intent(in out)      :: fq
        integer, intent(in out)          :: channel
        
        ! local
        real(KREAL)  :: max_power    = 0.0D0
        real(KREAL)  :: total_power  = 0.0D0
        real(KREAL)  :: local_avg    = 0.0D0
        real(KREAL)  :: local_max    = 0.0D0
        integer          :: count  = 0
        integer          :: index  = 0
        integer          :: ir
        
        do ir = 1, SIZE(this%avg_linear, dim=2)
            local_avg = SUM(this%avg_linear(:, ir))
            local_max = SUM(this%max_linear(:, ir))
            
            if (local_avg >= EPS_ZERO)  then
                count = count + 1
                total_power = total_power + local_avg
                
                if (local_max > max_power)  then
                    max_power = local_max
                    index = ir
                end if
            end if
        end do
        
        fq = max_power / (total_power / count)
        channel = index
    
    end subroutine Get_fuelpin_fq
    
    !$
    !===============================================================================================
    ! get fq factor by point power density
    !===============================================================================================
    subroutine Get_point_fq (this, fq, channel, axial)
        
        class(LinearPower), intent(in)   :: this
        real(KREAL), intent(in out)      :: fq
        integer, intent(in out)          :: channel
        integer, intent(in out)          :: axial
        
        ! local
        real(KREAL)  :: max_power    = 0.0D0
        real(KREAL)  :: total_power  = 0.0D0
        real(KREAL)  :: local_avg    = 0.0D0
        real(KREAL)  :: local_max    = 0.0D0
        integer          :: index_radial = 0
        integer          :: index_axial  = 0
        integer          :: count  = 0
        integer          :: ir, ia
        
        do ia = 1, SIZE(this%avg_linear, dim=1)
            do ir = 1, SIZE(this%avg_linear, dim=2)
                local_avg = this%avg_linear(ia, ir)
                local_max = this%max_linear(ia, ir)
                
                if (local_avg >= EPS_ZERO)  then
                    count = count + 1
                    total_power = total_power + local_avg
                    
                    if (local_max > max_power)  then
                        max_power = local_max
                        index_radial = ir
                        index_axial  = ia
                    end if
                end if
            end do
        end do

        fq = max_power / (total_power / count)
        channel = index_radial
        axial = index_axial
        
    end subroutine Get_point_fq

    !$
    !===============================================================================================
    ! set current power level per assembly
    !===============================================================================================
    subroutine Set_assembly_power_3D (this, nth, geom_th, geom_assm, power, fq_core)
    
        class(LinearPower), intent(in out)          :: this
        type(ThermalScale), intent(in)              :: nth
        type(ThermalGeometry), intent(in)           :: geom_th
        type(ThermalAssemblyGeometry), intent(in)   :: geom_assm(:)
        real(KREAL), intent(in)                     :: power(:, :)
        real(KREAL), intent(in)                     :: fq_core(:, :)
        
        integer  :: ir, ia
        integer  :: itype
        
        do ia = 1, SIZE(this%power, dim=1)
            do ir = 1, SIZE(this%power, dim=2)
                this%power(ia, ir) = power(ia, ir)
                this%fq_core(ia, ir) = fq_core(ia, ir)
            end do
        end do
        
        do ir = 1, nth%nr
            itype = geom_th%geom_type(ir)
            
            if (itype > 0)  then
                ! for average pin
                do ia = 1, nth%na
                    this%avg_linear(ia, ir) = this%power(ia, ir) / (geom_assm(itype)%n_fuelpin * geom_th%height(ia))
                end do
                ! for max pin
                do ia = 1, nth%na
                    this%max_linear(ia, ir) = this%avg_linear(ia, ir) * this%fq_lattice(ir) * this%fq_core(ia, ir)
                end do
            end if
        end do
        
    end subroutine Set_assembly_power_3D
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_assembly_power_0D (this, power, is_initial)
        
        class(LinearPower), intent(in out)    :: this
        real(KREAL), intent(in)               :: power
        logical, intent(in)                   :: is_initial
        
        ! local 
        real(KREAL)  :: amplitude
        
        if (is_initial)  then
            this%normal_old = 1.0D0
            this%normal_current = 1.0D0
            amplitude = this%normal_current / this%normal_old
            
        else
            this%normal_old = this%normal_current
            this%normal_current = power
            amplitude = this%normal_current / this%normal_old
            
            this%avg_linear = amplitude * this%avg_linear
            this%max_linear = amplitude * this%max_linear
        end if
    
    end subroutine Set_assembly_power_0D
    
end module gth_power_header
