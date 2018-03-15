!$
!===================================================================================================
!
!   module for thermal calculation
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          ThermalDesign
!                               ThermalChannel
!                               ThermalHotPoint
!
!===================================================================================================
module gth_thermal_header

    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use stastics
    
    use abstract_property_header,   only : CoolantProperty
    use gth_geometry_header,        only : ThermalScale, ThermalGeometry, ThermalAssemblyGeometry
    
    implicit none
    private 
    public  :: ThermalDesign, ThermalChannel, ThermalHotPoint
    
    ! --------------------------------------------------------------------------
    ! type for TH design parameter
    ! flow & flowrate --> Kg/s, flow density --> kg/(m^2.s)
    type  ThermalDesign
        integer, public                       :: tftype     = 1         ! type for effective Tf
        real(KREAL), public                   :: tfweight   = 0.7D0     ! weight for effective Tf 
        
        logical         :: is_search    = .FALSE.                       ! perform flowrate search
        real(KREAL)     :: tcoolin                                      ! coolant inlet temperature
        real(KREAL)     :: init_tcoolin                                 ! coolant inlet temperature
        real(KREAL)     :: tcoolout                                     ! coolant final outlet temperature
        real(KREAL)     :: flowrate                                     ! total flowrate of the core
                        
        real(KREAL)     :: tmincoolout                                  ! min temperature of outlet coolant
        real(KREAL)     :: tmaxcoolout                                  ! max temperature of outlet coolant
        real(KREAL)     :: tmaxcladsurf                                 ! max temperature of cladding surface (inner)
        real(KREAL)     :: max_velocity                                 ! max velocity of coolant in channel (m/s)
        
        real(KREAL), allocatable     :: TRU_weight(:)                   ! TRU weight fraction distribution
        real(KREAL), allocatable     :: channel_flowrate(:)             ! flow rate per channel, > 0.0 means active channel ? in kg/(m^2.s)
        real(KREAL), allocatable     :: assembly_flow(:)                ! flow per fuel assembly, in kg/s
        logical, allocatable         :: is_active_channel(:)            ! is this channel active core ?
        logical, allocatable         :: is_active_nodal(:, :)           ! is this nodal active core ?
        
        real(KREAL), allocatable     :: init_flow(:)
        real(KREAL), allocatable     :: init_flowrate(:)
        logical, public              :: is_blockage  = .FALSE.          ! is_blockage in transient ?
        logical, allocatable, public :: block_mask(:)                   ! which channel is blockage ?
        real(KREAL), public          :: block_time   = 0.0D0            ! which time blockage happened ?
        real(KREAL), public          :: percentage   = 0.0D0            ! blockage percentage 
    contains                                                              
        procedure, public  :: alloc => Alloc_ThermalDesign              
        procedure, public  :: clean => Free_ThermalDesign
        procedure, public  :: set => Set_ThermalDesign_flow
        procedure, public  :: print => Print_ThermalDesign
    end type ThermalDesign
                                                                                                                                     
    ! type for a channel data                                           
    type  ThermalChannel                                                
        real(KREAL), public                   :: tcoolant_out           ! final outlet coolant temperature
        real(KREAL), public                   :: rhocoolant_out
        real(KREAL), public, allocatable      :: flow_velocity(:, :)    ! velocity of flow
        real(KREAL), public, allocatable      :: convection(:, :)       ! transfer coefficient between clad surface and coolant
        real(KREAL), public, allocatable      :: rhocoolant(:, :)       ! coolant density distribution
        real(KREAL), public, allocatable      :: tcoolant(:, :)         ! coolant temperature distribution
        real(KREAL), public, allocatable      :: hjunction(:, :)        ! coolant enthalpy
                                                  
        real(KREAL), public, allocatable      :: tclad_surf(:, :)       ! cladding surface temperature distribution
        real(KREAL), public, allocatable      :: tclad_inner(:, :)      ! cladding inner surface temperature distribution
                                                  
        real(KREAL), public, allocatable      :: tfuel_surf(:, :)       ! fuel surface temperature distribution
        real(KREAL), public, allocatable      :: tfuel_center(:, :)     ! fuel center temperature distribution
        real(KREAL), public, allocatable      :: tfuel_avg(:, :)        ! fuel average temperature distribution
        real(KREAL), public, allocatable      :: trod(:, :, :)          ! fuel rod temperature mapping 
    contains                                                                      
        procedure, public  :: alloc => Alloc_ChannelThermal
        procedure, public  :: clean => Free_ChannelThermal
        procedure, public  :: update => Update_ChannelThermal
        procedure, public  :: fix_bottom => Set_ThermalChannel_bottom
        procedure, public  :: fix_top => Set_ThermalChannel_top
        procedure, public  :: homo => Get_Homogeneous_feedback

        procedure, public  :: print => Print_ChannelThermal
        procedure, public  :: print_avg => Print_ThermalChannel_average
        procedure, public  :: print_max => Print_ThermalChennel_max
        procedure, public  :: print_rod => Print_ChannelThermal_rod
    end type ThermalChannel
                                                                                  
    ! type for hot channel data
    type, private  :: hotpoint_record_tp
        integer          :: axial   = 1
        integer          :: radial  = 1
        real(KREAL)      :: value   = 0
    end type hotpoint_record_tp
    
    type  ThermalHotPoint
        type(hotpoint_record_tp)  :: coolant
        type(hotpoint_record_tp)  :: clad_surf
        type(hotpoint_record_tp)  :: clad_inner
        type(hotpoint_record_tp)  :: fuel_center
        
        real(KREAL), public     :: tcoolant                                     ! maximum coolant temperature
        real(KREAL), public     :: tclad_surf                                   ! maximum cladding surface temperature in core
        real(KREAL), public     :: tclad_inner                                  ! maximum cladding inner surface temperature in core
        real(KREAL), public     :: tfuel_center                                 ! maximum fuel core temperature in core
        
        logical, public           :: by_fuel_center = .TRUE.                    ! define hot by fuel center ?
        logical, public           :: by_clad_inner  = .FALSE.                   ! define hot by clad inner ?
        type(hotpoint_record_tp)  :: result                                     ! store hot info
        integer, public           :: hot_FA                                     ! initial hot index ?
    contains
        procedure, public  :: set => Set_HotValue
    end type ThermalHotPoint
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_ThermalDesign (this, nth)
    
        class(ThermalDesign), intent(in out)  :: this
        type(ThermalScale), intent(in)     :: nth
        
        integer  :: ia, i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%TRU_weight(nth%nr), stat=i_allocate)
        allocate(this%channel_flowrate(nth%nr), stat=i_allocate)
        allocate(this%assembly_flow(nth%nr), stat=i_allocate)
        
        allocate(this%is_active_channel(nth%nr), stat=i_allocate)
        allocate(this%is_active_nodal(nth%na, nth%nr), stat=i_allocate)
        
        allocate(this%init_flow(nth%nr), stat=i_allocate)
        allocate(this%init_flowrate(nth%nr), stat=i_allocate)
        allocate(this%block_mask(nth%nr), stat=i_allocate)
        
        this%TRU_weight         = 0.382D0
        this%channel_flowrate   = 0.0D0
        this%assembly_flow      = 0.0D0
        
        this%tcoolin       = 553.0D0                                            ! 280 centigrade
        this%tmincoolout   = 693.0D0                                            ! 420 
        this%tmaxcoolout   = 753.0D0                                            ! 480
        this%tmaxcladsurf  = 793.0D0                                            ! 520
        this%max_velocity  = 10.0
        
        ! default True
        this%is_active_channel = .TRUE.
        this%is_active_nodal   = .TRUE.
        this%init_flow         = this%assembly_flow
        this%init_flowrate     = this%channel_flowrate
        this%block_mask        = .FALSE.
        
    end subroutine Alloc_ThermalDesign
    
    !$
    !===============================================================================================
    ! finalizer for class of ThermalDesign
    !===============================================================================================
    subroutine Free_ThermalDesign (this)
        
        class(ThermalDesign), intent(in out)  :: this
        
        if (allocated(this%TRU_weight))           deallocate(this%TRU_weight)
        if (allocated(this%channel_flowrate))     deallocate(this%channel_flowrate)
        if (allocated(this%assembly_flow))        deallocate(this%assembly_flow)
        
        if (allocated(this%is_active_channel))    deallocate(this%is_active_channel)
        if (allocated(this%is_active_nodal))      deallocate(this%is_active_nodal)
        if (allocated(this%block_mask))           deallocate(this%block_mask)
        if (allocated(this%init_flow))            deallocate(this%init_flow)
        if (allocated(this%init_flowrate))        deallocate(this%init_flowrate)
    
    end subroutine Free_ThermalDesign
    
    !$
    !===============================================================================================
    ! get flow flux and flow density by inlet & outlet enthalpy
    !===============================================================================================
    subroutine Set_ThermalDesign_flow (this, a_coolant, nth, geom_th, geom_assm, lnpower)
        
        class(ThermalDesign), intent(in out)        :: this
        class(CoolantProperty), intent(in), pointer :: a_coolant
        type(ThermalScale), intent(in)              :: nth
        type(ThermalGeometry), intent(in)           :: geom_th
        type(ThermalAssemblyGeometry), intent(in)   :: geom_assm(:)
        real(KREAL), intent(in)                     :: lnpower(:, :)                   
        
        ! local parameter
        real(KREAL)  :: massflow                                            ! coolant mass flow flux, in Kg/s
        real(KREAL)  :: massflowrate                                        ! coolant mass flow flux density, in Kg/(m^2.s)
        real(KREAL)  :: channel_power
        integer          :: ia, ir, ityp
            
        ! get active channel & nodal
        this%is_active_channel = .TRUE.
        this%is_active_nodal   = .TRUE.
        
        do ir = 1, nth%nr
            ityp = geom_th%geom_type(ir)
            
            channel_power = 0.0D0
            do ia = 1, nth%na
                channel_power = channel_power + lnpower(ia, ir)*geom_th%height(ia)
            end do
            
            if (ityp == 0 .or. channel_power < EPS_ZERO)  then
                this%channel_flowrate(ir) = 0.0D0
                this%is_active_channel(ir)  = .FALSE.
                this%is_active_nodal(:, ir) = .FALSE.
            end if
        end do
        
        ! set axial active nodal
        do ir = 1, nth%nr
            if (this%is_active_channel(ir) ) then
                do ia = nth%na_end + 1, nth%na
                    this%is_active_nodal(ia, ir) = .FALSE.
                end do
                
                do ia = 1, nth%na_start - 1
                    this%is_active_nodal(ia, ir) = .FALSE.
                end do
            end if
        end do
            
        ! change input from (assembly, Kg/s) to (channel, Kg/(m^2.s))
        do ir = 1, nth%nr
            ityp = geom_th%geom_type(ir)
            
            if (this%is_active_channel(ir) .AND. (.NOT. this%is_search))  then
                this%channel_flowrate(ir) = this%assembly_flow(ir) / (geom_assm(ityp)%flowarea * geom_assm(ityp)%n_fuelpin)
                this%init_flow(ir) = this%assembly_flow(ir)
                this%init_flowrate(ir) = this%channel_flowrate(ir)
            end if
        end do
    
    end subroutine Set_ThermalDesign_flow
    
    !$
    !===============================================================================================
    !
    !===============================================================================================    
    subroutine Print_ThermalDesign (this, unit_)
        
        class(ThermalDesign), intent(in)  :: this
        integer, intent(in)               :: unit_
        integer  :: i, j
        
        do i = 1, SIZE(this%channel_flowrate)
            write(unit=unit_, fmt="(1x, ES12.5, TR3, L)")  this%channel_flowrate(i), this%is_active_channel(i)
        end do
        
        write(unit=unit_, fmt=*) ' '
        do i = 1, SIZE(this%is_active_nodal, dim=1)
            write(unit=unit_, fmt="(1x, *(L, TR3))")  this%is_active_nodal(i, :)
        end do
    
    end subroutine Print_ThermalDesign
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_ChannelThermal (this, nth)
    
        class(ThermalChannel), intent(in out)  :: this
        type(ThermalScale), intent(in)         :: nth
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%flow_velocity(nth%na, nth%nr), stat=i_allocate)
        allocate(this%convection(nth%na, nth%nr), stat=i_allocate)
        allocate(this%hjunction(0:nth%na, nth%nr), stat=i_allocate)             ! junction = nodal + 1
        allocate(this%tcoolant(0:nth%na, nth%nr), stat=i_allocate)
        allocate(this%rhocoolant(0:nth%na, nth%nr), stat=i_allocate)
        
        allocate(this%tclad_surf(0:nth%na, nth%nr), stat=i_allocate)
        allocate(this%tclad_inner(0:nth%na, nth%nr), stat=i_allocate)
        
        allocate(this%trod(nth%n_mesh, 0:nth%na, nth%nr), stat=i_allocate)
        allocate(this%tfuel_surf(0:nth%na, nth%nr), stat=i_allocate)
        allocate(this%tfuel_center(0:nth%na, nth%nr), stat=i_allocate)
        allocate(this%tfuel_avg(0:nth%na, nth%nr), stat=i_allocate)      
        
        this%tcoolant_out  = REAL_ZERO
        this%rhocoolant_out= REAL_ZERO
        this%flow_velocity = REAL_ZERO
        this%convection    = REAL_ZERO
        this%hjunction     = REAL_ZERO
        this%tcoolant      = REAL_ZERO
        this%rhocoolant    = REAL_ZERO
        
        ! NOTE: this value is used for inital guess in iteration
        this%tclad_surf    = CKELVIN
        this%tclad_inner   = CKELVIN

        this%trod          = REAL_ZERO
        this%tfuel_surf    = REAL_ZERO
        this%tfuel_center  = REAL_ZERO
        this%tfuel_avg     = REAL_ZERO
        
    end subroutine Alloc_ChannelThermal
    
    !$
    !===============================================================================================
    ! finalizer for class of ThermalChannel
    !===============================================================================================
    subroutine Free_ChannelThermal (this)
        
        class(ThermalChannel), intent(in out)  :: this
        
        if (allocated(this%flow_velocity))      deallocate(this%flow_velocity)
        if (allocated(this%convection))         deallocate(this%convection)
        if (allocated(this%hjunction))          deallocate(this%hjunction)
        if (allocated(this%tcoolant))           deallocate(this%tcoolant)
        if (allocated(this%rhocoolant))         deallocate(this%rhocoolant)
        
        if (allocated(this%tclad_surf))         deallocate(this%tclad_surf)
        if (allocated(this%tclad_inner))        deallocate(this%tclad_inner)
        
        if (allocated(this%trod))               deallocate(this%trod)
        if (allocated(this%tfuel_surf))         deallocate(this%tfuel_surf)
        if (allocated(this%tfuel_center))       deallocate(this%tfuel_center)
        if (allocated(this%tfuel_avg))          deallocate(this%tfuel_avg)
        
    end subroutine Free_ChannelThermal
    
    !$
    !===============================================================================================
    ! get Doppler temeprature by 1-emperical, 2-volume weighting, 3-chord weighting, 4-simulate
    !===============================================================================================
    subroutine Update_ChannelThermal (this, a_assembly, nth, design, solve, ia, ir)
        
        class(ThermalChannel), intent(in out)     :: this
        type(ThermalAssemblyGeometry), intent(in) :: a_assembly
        type(ThermalScale), intent(in)            :: nth
        type(ThermalDesign), intent(in)           :: design 
        real(KREAL), intent(in)      :: solve(:)
        integer, intent(in)          :: ia
        integer, intent(in)          :: ir
        
        integer, parameter  :: TF_NEA       = 1 
        integer, parameter  :: TF_VOLUME    = 2
        integer, parameter  :: TF_SIMULATE  = 3 
        integer, parameter  :: TF_CHORD     = 4 
        integer, parameter  :: TF_GDTL      = 5 
        
        real(KREAL)  :: r1, r2 
        real(KREAL)  :: tfsub
        real(KREAL)  :: tavg
        real(KREAL)  :: subwt
        real(KREAL)  :: sumwt
        integer  :: i
        
        this%trod(:, ia, ir) = solve
        this%tclad_surf(ia, ir)   = solve(nth%nf+nth%nc+2)
        this%tclad_inner(ia, ir)  = solve(nth%nf+2) 
        this%tfuel_surf(ia, ir)   = solve(nth%nf+1)
        this%tfuel_center(ia, ir) = solve(1)
        
        select case (design%tftype)
        case (TF_NEA)
            this%tfuel_avg(ia, ir) = design%tfweight * this%tfuel_surf(ia, ir) + (1.0D0-design%tfweight) * this%tfuel_center(ia, ir)
        
        case (TF_VOLUME)
            tavg = 0.0D0
            subwt = 0.0D0
            sumwt = 0.0D0
            do i = 2, nth%nf + 1
                r2 = a_assembly%x_surface(i)
                r1 = a_assembly%x_surface(i-1)
                subwt = PI*r2*r2 - PI*r1*r1
                tfsub = 0.5D0 * (solve(i) + solve(i-1))
                tavg = tavg + subwt*tfsub
                sumwt = sumwt + subwt
            end do
            this%tfuel_avg(ia, ir) = tavg / sumwt
        
        case (TF_SIMULATE)
            tavg = 0.0D0
            subwt = 0.0D0
            sumwt = 0.0D0
            do i = 2, nth%nf + 1
                r2 = a_assembly%x_surface(i)
                r1 = a_assembly%x_surface(i-1)
                subwt = PI*r2*r2 - PI*r1*r1
                tfsub = 0.5D0 * (solve(i) + solve(i-1))
                tavg = tavg + subwt*tfsub
                sumwt = sumwt + subwt
            end do
            this%tfuel_avg(ia, ir) = design%tfweight * (tavg/sumwt) + (1.0D0-design%tfweight) * this%tfuel_surf(ia, ir)
            
        case (TF_CHORD)                                                 ! chord = 4V/S ~ 2*r2*(1-(r1/r2)**2)
            tavg = 0.0D0
            subwt = 0.0D0
            sumwt = 0.0D0
            do i = 2, nth%nf + 1
                r2 = a_assembly%x_surface(i)
                r1 = a_assembly%x_surface(i-1)
                subwt = 2.0D0*r2*(1.0D0-(r1/r2)**2)
                tfsub = 0.5D0 * (solve(i) + solve(i-1))
                tavg = tavg + subwt*tfsub
                sumwt = sumwt + subwt
            end do
            this%tfuel_avg(ia, ir) = tavg / sumwt
            
        case (TF_GDTL)
            tavg = 0.0D0
            subwt = 0.0D0
            sumwt = 0.0D0
            do i = 2, nth%nf + 1
                r2 = a_assembly%x_surface(i)
                r1 = a_assembly%x_surface(i-1)
                subwt = r2 - r1
                tfsub = 0.5D0 * (solve(i) + solve(i-1))
                tfsub = SQRT(tfsub)
                tavg = tavg + subwt*tfsub
                sumwt = sumwt + subwt*(1.0/tfsub)
            end do
            this%tfuel_avg(ia, ir) = tavg / sumwt
            
        case default 
        end select 
        
    end subroutine Update_ChannelThermal
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_ThermalChannel_bottom (this, a_coolant, nth, design)
        
        class(ThermalChannel), intent(in out)        :: this
        class(CoolantProperty), intent(in), pointer  :: a_coolant
        type(ThermalScale), intent(in)               :: nth
        type(ThermalDesign), intent(in)              :: design
        
        ! local variables
        real(KREAL)  :: t_in
        integer  :: ia, ir
        
        t_in = design%tcoolin
        
        do ir = 1, nth%nr
            ! set bottom of active channel
            if (design%is_active_channel(ir) ) then
                
                do ia = 0, nth%na_start - 1
                    this%tcoolant(ia, ir)      = t_in
                    this%rhocoolant(ia, ir)    = a_coolant%get_density(t_in)
                    this%hjunction(ia,ir)      = a_coolant%get_enthalpy(t_in)
                    
                    this%tclad_surf(ia, ir)   = t_in
                    this%tclad_inner(ia, ir)  = t_in
                    
                    this%trod(:, ia, ir)      = t_in
                    this%tfuel_surf(ia, ir)   = t_in
                    this%tfuel_center(ia, ir) = t_in
                    this%tfuel_avg(ia, ir)    = t_in
                end do
                
                do ia = 1, nth%na_start - 1
                    this%flow_velocity(ia, ir) = design%channel_flowrate(ir) / a_coolant%get_density(t_in)
                end do
            
            ! set radial reflector region
            else 
                this%tcoolant(:, ir)      = t_in
                this%rhocoolant(:, ir)    = a_coolant%get_density(t_in)
                this%hjunction(:,ir)      = a_coolant%get_enthalpy(t_in)
                this%tcoolant_out         = t_in
                this%flow_velocity(:, ir) = 0.0D0
                
                this%tclad_surf(:, ir)   = t_in
                this%tclad_inner(:, ir)  = t_in
                
                this%trod(:, :, ir)      = t_in
                this%tfuel_surf(:, ir)   = t_in
                this%tfuel_center(:, ir) = t_in
                this%tfuel_avg(:, ir)    = t_in
            end if
            
            this%hjunction(0,ir) = a_coolant%get_density(t_in)
        end do
    
    end subroutine Set_ThermalChannel_bottom
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_ThermalChannel_top (this, nth, design)
        
        class(ThermalChannel), intent(in out)        :: this
        type(ThermalScale), intent(in)               :: nth
        type(ThermalDesign), intent(in)              :: design
        
        ! local variables
        integer  :: ia, ir
        
        do ir = 1, nth%nr
            ! set top of active channel
            if (design%is_active_channel(ir) ) then
                do ia = nth%na_end + 1, nth%na
                    this%tcoolant(ia, ir)      = this%tcoolant(nth%na_end, ir)
                    this%rhocoolant(ia, ir)    = this%rhocoolant(nth%na_end, ir)
                    this%flow_velocity(ia, ir) = this%flow_velocity(nth%na_end, ir)
                    
                    this%tclad_surf(ia, ir)   = this%tclad_surf(nth%na_end, ir)
                    this%tclad_inner(ia, ir)  = this%tclad_inner(nth%na_end, ir)
                    
                    this%trod(:, ia, ir)      = this%trod(:, nth%na_end, ir)
                    this%tfuel_surf(ia, ir)   = this%tfuel_surf(nth%na_end, ir)
                    this%tfuel_center(ia, ir) = this%tfuel_center(nth%na_end, ir)
                    this%tfuel_avg(ia, ir)    = this%tfuel_avg(nth%na_end, ir)
                end do
            end if
        end do
    
    end subroutine Set_ThermalChannel_top
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Get_Homogeneous_feedback (this, a_coolant, nth, design, geom_th, T_fuel, T_coolant, Rho_coolant)
        
        class(ThermalChannel), intent(in)            :: this
        class(CoolantProperty), intent(in), pointer  :: a_coolant
        type(ThermalScale), intent(in)               :: nth
        type(ThermalDesign), intent(in)              :: design
        type(ThermalGeometry), intent(in)            :: geom_th
        real(KREAL), intent(in out)              :: T_fuel
        real(KREAL), intent(in out)              :: T_coolant
        real(KREAL), intent(in out)              :: Rho_coolant
    
        ! local variables
        real(KREAL)  :: rho
        real(KREAL)  :: tmid
        integer  :: ia, ir, cnt
        
        T_fuel = 0.0D0
        cnt = 0
        do ir = 1, nth%nr
            if (design%is_active_channel(ir))  then
                cnt = cnt + 1
                do ia = nth%na_start, nth%na_end
                    tmid = (this%tfuel_surf(ia, ir) + this%tfuel_surf(ia-1, ir)) / 2.0
                    T_fuel = T_fuel + 0.7D0* geom_th%height(ia) * tmid / SUM(geom_th%height(nth%na_start:nth%na_end))
                    tmid = (this%tfuel_center(ia, ir) + this%tfuel_center(ia-1, ir)) / 2.0
                    T_fuel = T_fuel + 0.3D0* geom_th%height(ia) * tmid / SUM(geom_th%height(nth%na_start:nth%na_end))
                end do
            end if
        end do
        T_fuel = T_fuel / cnt
        
        T_coolant = 0.0D0
        cnt = 0
        do ir = 1, nth%nr
            if (design%is_active_channel(ir))  then
                cnt = cnt + 1
                do ia = nth%na_start, nth%na_end
                    tmid = (this%tcoolant(ia, ir) + this%tcoolant(ia-1, ir)) / 2.0
                    T_coolant = T_coolant + geom_th%height(ia) * tmid / SUM(geom_th%height(nth%na_start:nth%na_end))
                end do
            end if
        end do
        T_coolant = T_coolant / cnt
        
        Rho_coolant = 0.0D0
        cnt = 0
        do ir = 1, nth%nr
            if (design%is_active_channel(ir))  then
                cnt = cnt + 1
                do ia = nth%na_start, nth%na_end
                    tmid = (this%tcoolant(ia, ir) + this%tcoolant(ia-1, ir)) / 2.0
                    rho = a_coolant%get_density (tmid)
                    Rho_coolant = Rho_coolant + geom_th%height(ia) * rho / SUM(geom_th%height(nth%na_start:nth%na_end))
                end do
            end if
        end do
        Rho_coolant = Rho_coolant / cnt
    
    end subroutine Get_Homogeneous_feedback
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_ChannelThermal (this, nth, geom_th, unit_)
        
        class(ThermalChannel), intent(in out)  :: this
        type(ThermalGeometry), intent(in)      :: geom_th
        type(ThermalScale), intent(in)         :: nth
        integer, intent(in)  :: unit_
        
        real(KREAL)  :: distance
        integer  :: ia, ir
        
        ! for coolant
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'flow_velocity:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%flow_velocity(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
    
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'tcoolant:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%tcoolant(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
        
        ! for clad
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'tclad_surf:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%tclad_surf(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
    
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'tclad_inner:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%tclad_inner(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
        
        ! for fuel
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'tfuel_center:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%tfuel_center(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
    
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  'tfuel_avg:'
        distance = 0.0D0
        do ia = 1, nth%na
            distance = distance + 0.0D0*geom_th%height(ia)
            write(unit=unit_, fmt="(1x, I4, TR3, F8.4, TR3, *(ES12.5, TR3))")  ia, distance, this%tfuel_avg(ia, :)
            distance = distance + 1.0D0*geom_th%height(ia)
        end do
    
    end subroutine Print_ChannelThermal
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_ThermalChannel_average (this, unit_, design, nth, tidx, ctime)
        
        class(ThermalChannel), intent(in out)  :: this
        integer, intent(in)                    :: unit_
        type(ThermalDesign), intent(in)        :: design
        type(ThermalScale), intent(in)         :: nth
        integer, intent(in)                    :: tidx
        real(KREAL), intent(in)                :: ctime
        
        real(KREAL)  :: tcoolant
        real(KREAL)  :: rhocoolant
        real(KREAL)  :: tclad_surf
        real(KREAL)  :: tclad_inner
        real(KREAL)  :: tfuel_center
        real(KREAL)  :: tfuel_avg
        real(KREAL)  :: outlet_tcoolant
        real(KREAL)  :: outlet_rhocoolant
        real(KREAL), allocatable  :: tmp(:, :)
        integer                   :: ii, jj
        
        ii = SIZE(this%tcoolant, dim=1) - 1
        jj = SIZE(this%tcoolant, dim=2)
        allocate(tmp(ii, jj))
        
        tmp(1:ii, :) = (this%tcoolant(0:ii-1, :) + this%tcoolant(1:ii, :)) / 2.0
        tcoolant   = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        tmp(1:ii, :) = (this%rhocoolant(0:ii-1, :) + this%rhocoolant(1:ii, :)) / 2.0
        rhocoolant  = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        tmp(1:ii, :) = (this%tclad_surf(0:ii-1, :) + this%tclad_surf(1:ii, :)) / 2.0
        tclad_surf   = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        tmp(1:ii, :) = (this%tclad_inner(0:ii-1, :) + this%tclad_inner(1:ii, :)) / 2.0
        tclad_inner  = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        tmp(1:ii, :) = (this%tfuel_center(0:ii-1, :) + this%tfuel_center(1:ii, :)) / 2.0
        tfuel_center = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        tmp(1:ii, :) = (this%tfuel_avg(0:ii-1, :) + this%tfuel_avg(1:ii, :)) / 2.0
        tfuel_avg    = stastics_average_value (tmp, mask=design%is_active_nodal)
        
        outlet_tcoolant = stastics_average_value (this%tcoolant(nth%na, :), mask=design%is_active_channel)
        outlet_rhocoolant = stastics_average_value (this%rhocoolant(nth%na, :), mask=design%is_active_channel)
        
        write(unit=unit_, fmt="(1x, I5, *(TR3, ES12.5))") tidx, ctime, tcoolant, rhocoolant, tfuel_avg,  &
                                                                    &   outlet_tcoolant, outlet_rhocoolant

    end subroutine Print_ThermalChannel_average
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_ThermalChennel_max (this, unit_, design, nth, tidx, ctime)
        
        class(ThermalChannel), intent(in out)  :: this
        integer, intent(in)                    :: unit_
        type(ThermalDesign), intent(in)        :: design
        type(ThermalScale), intent(in)         :: nth
        integer, intent(in)                    :: tidx
        real(KREAL), intent(in)                :: ctime
        
        real(KREAL)  :: tcoolant
        real(KREAL)  :: rhocoolant
        real(KREAL)  :: tclad_surf
        real(KREAL)  :: tclad_inner
        real(KREAL)  :: tfuel_center
        real(KREAL)  :: tfuel_avg
        
        tcoolant   = stastics_max_value (this%tcoolant(1:, :), design%is_active_nodal)
        rhocoolant  = stastics_max_value (this%rhocoolant(1:, :), design%is_active_nodal)
        tclad_surf   = stastics_max_value (this%tclad_surf(1:, :), design%is_active_nodal)
        tclad_inner  = stastics_max_value (this%tclad_inner(1:, :), design%is_active_nodal)
        tfuel_center = stastics_max_value (this%tfuel_center(1:, :), design%is_active_nodal)
        tfuel_avg    = stastics_max_value (this%tfuel_avg(1:, :), design%is_active_nodal)
        
        write(unit=unit_, fmt="(1x, I5, *(TR3, ES12.5))") tidx, ctime, tcoolant, rhocoolant, tfuel_avg

    end subroutine Print_ThermalChennel_max
    
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Print_ChannelThermal_rod (this, unit_, nth, ir)
    
        class(ThermalChannel), intent(in out)  :: this
        integer, intent(in)                    :: unit_
        type(ThermalScale), intent(in)         :: nth
        integer, intent(in)                    :: ir
        
        integer :: ia
        
        write(unit=unit_, fmt="(1x, A, I3)")  "radial assembly is: ", ir
        do ia = 1, UBOUND(this%trod, dim=2)
            write(unit=unit_, fmt="(1x, I4, *(TR3, ES12.5))")  ia, this%trod(:, ia, ir),          &
                &   this%tfuel_center(ia, ir), this%tfuel_avg(ia, ir), this%tfuel_surf(ia, ir), &
                &   this%tclad_inner(ia,ir), this%tclad_surf(ia,ir)
        end do
        
    end subroutine Print_ChannelThermal_rod
        
    !$
    !===============================================================================================
    ! 
    !===============================================================================================
    subroutine Set_HotValue (this, hot_channel, nth, is_initial)
        
        class(ThermalHotPoint), intent(in out)   :: this
        type(ThermalChannel), intent(in)  :: hot_channel
        type(ThermalScale), intent(in)    :: nth
        logical, intent(in), optional     :: is_initial
        
        integer  :: ir, ia
        
        this%coolant%value     = 0.0D0
        this%clad_surf%value   = 0.0D0
        this%clad_inner%value  = 0.0D0
        this%fuel_center%value = 0.0D0
        
        do ir = 1, nth%nr
            do ia = 1, nth%na
                if (hot_channel%tcoolant(ia,ir) > this%coolant%value)  then
                    this%coolant%value  = hot_channel%tcoolant(ia,ir)
                    this%coolant%axial  = ia
                    this%coolant%radial = ir
                end if
            
                if (hot_channel%tclad_surf(ia,ir) > this%clad_surf%value)  then
                    this%clad_surf%value  = hot_channel%tclad_surf(ia,ir)
                    this%clad_surf%axial  = ia
                    this%clad_surf%radial = ir
                end if
        
                if (hot_channel%tclad_inner(ia,ir) > this%clad_inner%value)  then
                    this%clad_inner%value  = hot_channel%tclad_inner(ia,ir)
                    this%clad_inner%axial  = ia
                    this%clad_inner%radial = ir
                end if
        
                if (hot_channel%tfuel_center(ia,ir) > this%fuel_center%value)  then
                    this%fuel_center%value  = hot_channel%tfuel_center(ia,ir)
                    this%fuel_center%axial  = ia
                    this%fuel_center%radial = ir
                end if
            end do
        end do
        
        this%tcoolant     = this%coolant%value     
        this%tclad_surf   = this%clad_surf%value   
        this%tclad_inner  = this%clad_inner%value  
        this%tfuel_center = this%fuel_center%value 
        
        
        if (this%by_fuel_center)  then
            this%result%axial = this%fuel_center%axial
            this%result%radial = this%fuel_center%radial
            this%result%value = this%fuel_center%value
            
            if (PRESENT(is_initial) .and. is_initial)  then
                this%hot_FA = this%fuel_center%radial
            end if
        end if
        
        if (this%by_clad_inner)  then
            this%result%axial = this%clad_inner%axial
            this%result%radial = this%clad_inner%radial
            this%result%value = this%clad_inner%value
            
            if (PRESENT(is_initial) .and. is_initial)  then
                this%hot_FA = this%clad_inner%radial
            end if
        end if
    
    end subroutine Set_HotValue
    
end module gth_thermal_header
