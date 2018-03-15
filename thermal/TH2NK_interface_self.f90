!$
!===================================================================================================
!
!   module for thermal hydraulic feedback interface
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Perform_TH_self
!
!   Public type lists:          HotRecord_self
!
!===================================================================================================
module TH2NK_interface_self
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use stastics,           only : stastics_max_value, stastics_min_value
    
    use th_global
    use th_parallel_channel
    
    implicit none
    private
    public  :: Perform_TH_self, HotRecord_self
    
    type  HotRecord_self
        integer                 :: axial   = 1
        integer                 :: radial  = 1
        real(8)                 :: value   = 0.0D0
    end type HotRecord_self
    
    type(HotRecord_self)  :: smax_Tfuel
    type(HotRecord_self)  :: smax_Tcoolant
    type(HotRecord_self)  :: smin_Rhocoolant

contains
    !$
    !===============================================================================================
    !
    !===============================================================================================
    subroutine Perform_TH_self(transient_flag, assembly, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
        
        logical, intent(in)      :: transient_flag                              ! .TRUE. --transient
        real(KREAL), intent(in)  :: assembly(:, :)                              ! (nr, na), in W, 各组件功率;
        real(KREAL), intent(in out)  :: Tfuel(:, :)                             ! (nr, na), in K, 各组件平均燃料温度;
        real(KREAL), intent(in out)  :: Tcoolant(:, :)                          ! (nr, na), in K, 各组件平均冷却剂温度;
        real(KREAL), intent(in out)  :: Rhocoolant(:, :)                        ! (nr, na), in Kg/m^3, 各组件平均冷却剂密度;
        real(KREAL), intent(in out)  :: max_Tfuel                               ! in K, 最热组件最大燃料温度;
        real(KREAL), intent(in out)  :: max_Tcoolant                            ! in K, 最热组件最大冷却剂温度;
        real(KREAL), intent(in out)  :: min_Rhocoolant                          ! in Kg/m^3, 最热组件最大冷却剂密度;
        real(KREAL), intent(in)  :: last                                        ! in s, 上一时间点
        real(KREAL), intent(in)  :: current                                     ! in s, 当前时间点
        real(KREAL), intent(in out)  :: toutlet                                 ! in K, 冷却剂出口平均温度
    
        !
        REAL(KREAL)  :: last_, current_
        real(KREAL), allocatable  :: power(:, :)
        real(KREAL), allocatable  :: fq_core(:, :)
        integer  :: nr, na, npin
        integer  :: ir, ia, ipin, itype
        integer  :: i_allocate
        
        last_ = last
        current_ = current
        nr = SIZE(assembly, dim=1)                                              ! 径向的组件数目
        na = SIZE(assembly, dim=2)                                              ! 轴向的节块数目
        
        allocate(power(na, nr), stat=i_allocate)
        allocate(fq_core(na, nr), stat=i_allocate)
        power = 0.0
        
        do ir = 1, nr
            itype = geom_th%geom_type(ir)
            if (itype > 0)  then
                do ia = 1, na
                    power(ia, ir) = assembly(ir, ia)
                end do
            end if
        end do
        
        fq_core = 1.0D0
        if (transient_flag)  then
            call Driving_ParallelChannel_transient (power, fq_core, 1, last_, current_)
        else
            call Driving_ParallelChannel_steady (power, fq_core)
        end if
        
        if (allocated(power))       deallocate(power)
        if (allocated(fq_core))     deallocate(fq_core)
        
        do ia = 1, na
            do ir = 1, nr
			!use (T(z)+T(z-1))/2.0 because FDM
                Tfuel(ir, ia) = (avg_channel%tfuel_avg(ia, ir) + avg_channel%tfuel_avg(ia-1, ir)) / 2.0
                Tcoolant(ir, ia) = (avg_channel%tcoolant(ia, ir) + avg_channel%tcoolant(ia-1, ir)) / 2.0
                Rhocoolant(ir, ia) = (avg_channel%rhocoolant(ia, ir) + avg_channel%rhocoolant(ia-1, ir)) / 2.0
!                Tfuel(ir, ia) = avg_channel%tfuel_avg(ia, ir) 
!                Tcoolant(ir, ia) = avg_channel%tcoolant(ia, ir) 
!                Rhocoolant(ir, ia) = avg_channel%rhocoolant(ia, ir) 
            end do
        end do
        
        smax_Tfuel%axial  = 1
        smax_Tfuel%radial = 1
        smax_Tfuel%value  = 0.0D0
        smax_Tcoolant%axial  = 1
        smax_Tcoolant%radial = 1
        smax_Tcoolant%value  = 0.0D0
        smin_Rhocoolant%axial  = 1
        smin_Rhocoolant%radial = 1
        smin_Rhocoolant%value  = 0.0D0
        do ia = 1, na
            do ir = 1, nr
                if (avg_channel%tfuel_center(ia, ir) > smax_Tfuel%value)  then
                    smax_Tfuel%axial  = ia
                    smax_Tfuel%radial = ir
                    smax_Tfuel%value  = avg_channel%tfuel_center(ia, ir)
                end if 
                
                if (avg_channel%tcoolant(ia, ir) > smax_Tcoolant%value)  then
                    smax_Tcoolant%axial  = ia
                    smax_Tcoolant%radial = ir
                    smax_Tcoolant%value  = avg_channel%tcoolant(ia, ir)
                end if 
                
                if (avg_channel%rhocoolant(ia, ir) < smin_Rhocoolant%value)  then
                    smin_Rhocoolant%axial  = ia
                    smin_Rhocoolant%radial = ir
                    smin_Rhocoolant%value  = avg_channel%tcoolant(ia, ir)
                end if 
            end do
        end do
        
        max_Tfuel = smax_Tfuel%value
        max_Tcoolant = smax_Tcoolant%value
        min_Rhocoolant = smin_Rhocoolant%value
        
        toutlet = avg_channel%tcoolant_out
    
    end subroutine Perform_TH_self

end module TH2NK_interface_self
