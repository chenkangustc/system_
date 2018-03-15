!$
!===================================================================================================
!
!   this module is for cross section feedback parameter definition;
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          FeedbackParameter
!
!===================================================================================================
module feedback_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use vector_operation,       only : get_vector_error
    use stastics,               only : stastics_average_value
    use exception_header,       only : WarningCollector
    use state_header,           only : SteadyState
    use geometry_header,        only : Geometry
    use link_header,            only : LinkParameter
    
    implicit none
    private
    public  :: FeedbackParameter
    
    ! --------------------------------------------------------------------------
    type, private  :: feedback_info_tp
        real(KREAL)                   :: out   = 100.0D0                        ! outlet value 
        real(KREAL)                   :: max   = 100.0D0                        ! max vlaue in core
        real(KREAL)                   :: avg   = 1.0D0                          ! average value in core (weighted by volume)
        real(KREAL)                   :: limit = HUGE(1.0)                      ! relative error limit in iteration
        real(KREAL)                   :: error = 0.0D0                          ! relative error
        real(KREAL), allocatable      :: old(:, :)                              ! old value in iteration
        real(KREAL), allocatable      :: new(:, :)                              ! new value in iteration
    contains
        procedure, public  :: alloc => Alloc_feedback_info_tp
        procedure, public  :: clean => Free_feedback_info_tp
        procedure, public  :: check => Check_feedback_info_tp
    end type feedback_info_tp
    
    type  FeedbackParameter
        real(KREAL), allocatable      :: Bu(:, :)                               ! this is an input parameter if used
        real(KREAL), allocatable      :: den_Xe(:, :)                           ! density of Xe
        real(KREAL), allocatable      :: den_Sm(:, :)                           ! density of Sm
        type(feedback_info_tp)        :: CB                                     ! only one value, if the distribution can not be obtained from the TH code
        type(feedback_info_tp)        :: Tf
        type(feedback_info_tp)        :: Tm
        type(feedback_info_tp)        :: Rho_m
        type(feedback_info_tp)        :: power                                  ! normalized power density
        logical                       :: is_CB_search     = .FALSE.             ! is perform critical CB search ?
        logical                       :: is_CB_critical   = .FALSE.             ! is critical CB searched ?
        real(KREAL)                   :: CB_critical      = 1000.0D0
        real(KREAL)                   :: CB_target        = 1.0D0
        real(KREAL)                   :: CB_limit         = 2.0D-5              ! comment CB search
        logical                       :: is_th_external   = .FALSE.
        logical                       :: is_th_inner      = .FALSE.
        logical                       :: is_model         = .FALSE.
        character(len=MAX_WORD_LEN)   :: model_name       = 'LRA'
        real(KREAL)                   :: relax_TH         = 0.7                 ! relaxation factor for NK-TH iteration, especially for full power
        real(KREAL)                   :: relax_CB         = 0.7                 ! relaxation factor for CB search
        integer                       :: relax_idx        = 0
    contains
        procedure, public  :: alloc => Alloc_FeedbackParameter
        procedure, public  :: clean => Free_FeedbackParameter
        procedure, public  :: set => Set_FeedbackParameter
        procedure, public  :: update => Update_FeedbackParameter
        procedure, public  :: check => Check_FeedbackParameter
        procedure, public  :: relax => Relax_FeedbackParameter
        procedure, public  :: print => Print_FeedbackParameter
        procedure, public  :: get_state => Get_FeedbackParameter_Nodal
        procedure, public  :: CB_search => Search_FeedbackParameter_CB
    end type FeedbackParameter
    
    ! --------------------------------------------------------------------------
    type(WarningCollector)  :: a_warning
    integer, parameter      :: ERROR_TYPE_ = 2
    
contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_feedback_info_tp (this, ns)
        
        class(feedback_info_tp), intent(in out)  :: this
        type(SteadyState), intent(in)            :: ns
        
        integer  :: i_allocate
    
        ! check allocated status first
        call this%clean ()
        
        allocate(this%old(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%new(ns%state%zone, ns%state%layer), stat=i_allocate)
        
        this%old = REAL_ONE
        this%new = REAL_ZERO
    
    end subroutine Alloc_feedback_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_feedback_info_tp (this)
    
        class(feedback_info_tp), intent(in out)  :: this
        
        if (allocated(this%old))        deallocate(this%old)
        if (allocated(this%new))        deallocate(this%new)
    
    end subroutine Free_feedback_info_tp
    
    !$
    !===============================================================================================
    ! abs error, infinity normal number
    !===============================================================================================
    subroutine Check_feedback_info_tp (this, is_pass)
        
        class(feedback_info_tp), intent(in out)  :: this
        logical, intent(in out)                  :: is_pass
        
        this%error = get_vector_error (this%old, this%new, ERROR_TYPE_)
        if (this%error <= this%limit)  then
            is_pass = .TRUE.
        else
            is_pass = .FALSE.
        end if
        
    end subroutine Check_feedback_info_tp
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Alloc_FeedbackParameter (this, ns)
        
        class(FeedbackParameter), intent(in out)  :: this
        type(SteadyState), intent(in)             :: ns
        
        integer  :: i_allocate
        
        ! check allocated status first
        call this%clean ()
        
        allocate(this%Bu(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%den_Xe(ns%state%zone, ns%state%layer), stat=i_allocate)
        allocate(this%den_Sm(ns%state%zone, ns%state%layer), stat=i_allocate)
        call this%CB%alloc (ns)
        call this%Tf%alloc (ns)
        call this%Tm%alloc (ns)
        call this%Rho_m%alloc (ns)
        call this%power%alloc (ns)
        
        this%Bu = REAL_ZERO 
    
    end subroutine Alloc_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Free_FeedbackParameter (this)
        
        class(FeedbackParameter), intent(in out)  :: this
        
        if (allocated(this%Bu))         deallocate(this%Bu)
        if (allocated(this%den_Xe))     deallocate(this%den_Xe)
        if (allocated(this%den_Sm))     deallocate(this%den_Sm)
        
        call this%CB%clean ()
        call this%Tf%clean ()
        call this%Tm%clean ()
        call this%Rho_m%clean ()
        call this%power%clean ()
    
    end subroutine Free_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Set_FeedbackParameter (this, self_link)
        
        class(FeedbackParameter), intent(in out)  :: this
        type(LinkParameter), intent(in)           :: self_link
        
        character(len=MAX_WORD_LEN)  :: symbol
        integer  :: i
        
        this%CB%old = 1400.0D0
        this%Tf%old = 900.0D0
        this%Tm%old = 580.0D0
        this%Rho_m%old = 720.0D0
        this%power%old = 1.0D0
        
        do i = 1, self_link%n_parameter
            symbol = TRIM(ADJUSTL(self_link%info(i)%name))
            select case (symbol)
            ! if read from burnupMap, do not overrid
            case ('BU')
!                this%Bu = self_link%info(i)%reference
            case ('CB')
                this%CB%old = self_link%info(i)%reference
            case ('TF')
                this%Tf%old = self_link%info(i)%reference
            case ('STF')
                this%Tf%old = self_link%info(i)%reference*self_link%info(i)%reference
            case ('TM')
                this%Tm%old = self_link%info(i)%reference
            case ('RHO_M')
                this%Rho_m%old = self_link%info(i)%reference
            case ('POWER')
                this%power%old = self_link%info(i)%reference
            end select 
        end do
        
        this%CB%new = this%CB%old 
        this%Tf%new = this%Tf%old
        this%Tm%new = this%Tm%old
        this%Rho_m%new = this%Rho_m%old
        this%power%new = this%power%old
    
    end subroutine Set_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Update_FeedbackParameter (this, geom, Tf, Tm, Rho_m, power, hot_Tf, hot_Tm, out_Tm, hot_Rho_m, mask)
        
        class(FeedbackParameter), intent(in out)  :: this
        type(Geometry), intent(in)        :: geom
        real(KREAL), intent(in)           :: Tf(:, :)
        real(KREAL), intent(in), optional :: Tm(:, :)
        real(KREAL), intent(in), optional :: Rho_m(:, :)
        real(KREAL), intent(in), optional :: power(:, :)
        real(KREAL), intent(in), optional :: hot_Tf
        real(KREAL), intent(in), optional :: hot_Tm
        real(KREAL), intent(in), optional :: out_Tm
        real(KREAL), intent(in), optional :: hot_Rho_m
        logical, intent(in), optional     :: mask(:, :)
        
        real(KREAL), allocatable  :: volume(:, :)
        integer  :: iz, ia
        integer  :: i_allocate
        
        allocate(volume(SIZE(geom%zone_area), SIZE(geom%height)), stat=i_allocate)
        volume = 0.0D0
        do ia = 1, SIZE(volume, dim=2)
            do iz = 1, SIZE(volume, dim=1)
                volume(iz, ia) = geom%zone_area(iz) * geom%height(ia)
            end do
        end do
               
        this%CB%old = this%CB%new
        this%Tf%old = this%Tf%new
        this%Tm%old = this%Tm%new
        this%Rho_m%old = this%Rho_m%new
        this%power%old = this%power%new
        
        this%Tf%new = Tf
        this%Tf%avg = stastics_average_value (this%Tf%new, volume, mask)
        if (PRESENT(Tm))  then
            this%Tm%new = Tm
            this%Tm%avg = stastics_average_value (this%Tm%new, volume, mask)
        end if
        if (PRESENT(Rho_m))  then
            this%Rho_m%new = Rho_m
            this%Rho_m%avg = stastics_average_value (this%Rho_m%new, volume, mask)
        end if
        if (PRESENT(power))     this%power%new = power
        if (PRESENT(hot_Tf))    this%Tf%max = hot_Tf
        if (PRESENT(hot_Tm))    this%Tm%max = hot_Tm
        if (PRESENT(out_Tm))    this%Tm%out = out_Tm
        if (PRESENT(hot_Rho_m)) this%Rho_m%max = hot_Rho_m
        
        if (allocated(volume))          deallocate(volume)
    
    end subroutine Update_FeedbackParameter
        
    !$
    !===============================================================================================
    ! check 'Tf' & 'Tm' only 
    !===============================================================================================
    subroutine Check_FeedbackParameter (this, is_pass, to_write)
        
        class(FeedbackParameter), intent(in out)  :: this
        logical, intent(in out)                   :: is_pass
        logical, intent(in)                       :: to_write
        
        logical  :: pass_(5)
        
        pass_ = .TRUE.
        call this%CB%check (pass_(1))
        call this%power%check (pass_(2))
        call this%Tf%check (pass_(3))                                           ! 
        call this%Tm%check (pass_(4))                                           ! 
        call this%Rho_m%check (pass_(5))
        
        is_pass = .FALSE.
        if (ALL(pass_))  then
            is_pass = .TRUE.
        end if
        
        if (.NOT. to_write) then
            write(FILES%MAIN, fmt="(1x, A)") TRIM(CHAR_SUBMARK)
            write(FILES%MAIN, fmt="(2x, A, TR3, A)") "Normal-2 error of Tf", " Error limit"
            write(FILES%MAIN, fmt="(1x, TR8, *(ES13.6, TR3))") this%Tf%error, this%Tf%limit
            write(FILES%MAIN, fmt="(2x, A, TR3, A)") "Normal-2 error of Tm", " Error limit"
            write(FILES%MAIN, fmt="(1x, TR8, *(ES13.6, TR3))") this%Tm%error, this%Tm%limit
        end if
    
    end subroutine Check_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Relax_FeedbackParameter (this)
        
        class(FeedbackParameter), intent(in out)  :: this
        
        real(KREAL)  :: left, right
        
        this%relax_idx = this%relax_idx + 1
        
        if (this%relax_idx >= 2)  then
            left = 1.0 - this%relax_TH
            right = this%relax_TH
        else
            left = 0.0
            right = 1.0 
        end if 
        
!        this%CB%new = left*this%CB%old + right*this%CB%new
        this%power%new = left*this%power%old + right*this%power%new
        this%Tf%new = left*this%Tf%old + right*this%Tf%new
        this%Tm%new = left*this%Tm%old + right*this%Tm%new
        this%Rho_m%new = left*this%Rho_m%old + right*this%Rho_m%new
    
    end subroutine Relax_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Print_FeedbackParameter (this, unit_)
        
        class(FeedbackParameter), intent(in)  :: this
        integer, intent(in)                   :: unit_
        
        integer  :: i, j
        
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  "CB:"
        do j = 1, SIZE(this%Tf%new, dim=2)
            write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  this%CB%new(:, j)
        end do
        
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  "Tf:"
        do j = 1, SIZE(this%Tf%new, dim=2)
            write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  this%Tf%new(:, j)
        end do
    
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  "Tm:"
        do j = 1, SIZE(this%Tm%new, dim=2)
            write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  this%Tm%new(:, j)
        end do
    
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt=*)  "Rho_m:"
        do j = 1, SIZE(this%Rho_m%new, dim=2)
            write(unit=unit_, fmt="(1x, *(ES13.6, TR3))")  this%Rho_m%new(:, j)
        end do
        
        write(unit=unit_, fmt=*)  TRIM(CHAR_SSUBMARK)
        write(unit=unit_, fmt="(1x, A20, TR3, F8.4)")  "relax_TH:", this%relax_TH
        write(unit=unit_, fmt="(1x, A20, TR3, F8.4)")  "relax_CB:", this%relax_CB
    
    end subroutine Print_FeedbackParameter
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Get_FeedbackParameter_Nodal (this, self_link, iz, ia, values)
        
        class(FeedbackParameter), intent(in)  :: this
        type(LinkParameter), intent(in)       :: self_link
        integer, intent(in)                   :: ia
        integer, intent(in)                   :: iz
        real(8), intent(in out)               :: values(:)                      ! real(8) for consistent with Lilac
        
        character(len=MAX_WORD_LEN)  :: info_line
        character(len=MAX_WORD_LEN)  :: symbol
        real(KREAL)    :: left, right
        integer  :: i
        
        do i = 1, self_link%n_parameter
            symbol = TRIM(ADJUSTL(self_link%info(i)%name))
            select case (symbol)
            case ('BU')
                values(i) = this%Bu(iz, ia)
            case ('CB')
                values(i) = this%CB%new(iz, ia)
            case ('TF')
                values(i) = this%Tf%new(iz, ia)
            case ('STF')
                values(i) = SQRT(this%Tf%new(iz, ia))
            case ('TM')
                values(i) = this%Tm%new(iz, ia)
            case ('RHO_M')
                values(i) = this%Rho_m%new(iz, ia)
            case ('POWER')
                values(i) = this%power%new(iz, ia)
            end select
            
            ! extrapolate the fitting
            left = self_link%info(i)%min - 0.02D0*(self_link%info(i)%reference - self_link%info(i)%min)
            right = self_link%info(i)%max + 0.02D0*(self_link%info(i)%max - self_link%info(i)%reference)
            if (values(i) > right)  then
                values(i) = right
            end if
            if (values(i) < left)  then
                values(i) = left 
            end if
        end do
        
    end subroutine Get_FeedbackParameter_Nodal
    
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Search_FeedbackParameter_CB (this, k_eff, is_pass)
        
        class(FeedbackParameter), intent(in out)  :: this
        real(KREAL), intent(in)      :: k_eff
        logical, intent(in out)      :: is_pass
        
        integer, save          :: cnt = 0
        logical, save          :: is_damp = .FALSE.
        real(KREAL), save      :: CB_old, keff_old
        real(KREAL)            :: CB_new, keff_new
        real(KREAL)            :: error
        real(KREAL)            :: CB_LL, keff_LL
        
        is_pass = .FALSE.
        if (ABS(k_eff - this%CB_target) <= this%CB_limit)  then
            this%is_CB_critical = .TRUE.
            is_pass = .TRUE.
            return
        end if
        
        cnt = cnt + 1
        ! save old value 
        if (cnt == 1)  then
            CB_LL = stastics_average_value (this%CB%new, is_zero=.FALSE.)
            keff_LL = k_eff
        else 
            CB_LL = CB_old
            keff_LL = keff_old
        end if 
        
        if (cnt == 1)  then
            error = (k_eff - this%CB_target) * 2500.0D0
            CB_old = stastics_average_value (this%CB%new, is_zero=.FALSE.)
            keff_old = k_eff
            this%CB%new = this%CB%new + error
        else
            CB_new = stastics_average_value (this%CB%new, is_zero=.FALSE.)
            keff_new = k_eff
            if ((cnt >= 5) .OR. ((ABS(k_eff-this%CB_target)) <= 0.0001) .OR. ((ABS(keff_new - keff_old)) <= 0.0001))  then
                is_damp = .TRUE.
            end if
            
            error = (this%CB_target-keff_new) / (keff_new-keff_old) * (CB_new-CB_old) * this%relax_CB     ! based on new
            this%CB%new = CB_new + error
            CB_old = CB_new
            keff_old = keff_new
        end if
        
        write(FILES%MAIN, fmt="(1x, A)")  TRIM(CHAR_SUBMARK)
        write(FILES%MAIN, fmt="(1x, A)")  "Boron Search Results:"
        write(FILES%MAIN, fmt="(1x, *(A, TR3))")  "  No. ",  "LL-CB ", "keff     ", "L-CB ", "keff    ", "==>", "this CB       delta-CB" 
        write(FILES%MAIN, fmt="(1x, I4, TR3, 2(F9.2, TR1, F9.6, TR1), TR6, F9.2)")  cnt, CB_LL, keff_LL, CB_old, keff_old, stastics_average_value (this%CB%new, is_zero=.FALSE.) 
        
        if (.TRUE.)  then 
            write(OUTPUT_UNIT, fmt="(1x, A)")  TRIM(CHAR_SUBMARK)
            write(OUTPUT_UNIT, fmt="(1x, A)")  "Boron Search Results:"
            write(OUTPUT_UNIT, fmt="(1x, *(A, TR3))")  "  No. ",  "LL-CB ", "keff     ", "L-CB ", "keff    ", "==>", "this CB       delta-CB" 
            write(OUTPUT_UNIT, fmt="(1x, I4, TR3, 2(F9.2, TR1, F9.6, TR1), TR6, F9.2, TR5, F10.4)")  cnt, CB_LL, keff_LL, CB_old, keff_old, stastics_average_value (this%CB%new, is_zero=.FALSE.), error
        end if 
        
    end subroutine Search_FeedbackParameter_CB
    
end module feedback_header
    