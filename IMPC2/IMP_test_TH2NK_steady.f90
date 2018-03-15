!***************************************************************************************
! test if
!  module: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!  pow(na,nr),fq_core(na,nr)     平均功率密度，功率峰因子
!  nr = SIZE(assembly, dim=1)    径向的组件数目
!  na = SIZE(assembly, dim=2)    轴向的节块数目，原输入变量不需要操作赋值的另外用局部变量表达
! 
!***************************************************************************************
    module test_TH2NK_steady
     use imp_re_input_global
     use imp_assm_global
     use imp_driving_pre_process
     use imp_driving_post_process
     use imp_power_header
     use imp_single_channel
	 use global_state
     use constants
    implicit none           
	                                                            !(zone,layer)
        real(KREAL),allocatable:: assembly(:, :)                          ! (nr, na), in W, 各组件功率;
        real(KREAL),allocatable:: Tfuel(:, :)                             ! (nr, na), in K, 各组件平均燃料温度;
        real(KREAL),allocatable:: Tcoolant(:, :)                          ! (nr, na), in K, 各组件平均冷却剂温度;
        real(KREAL),allocatable:: Rhocoolant(:, :)                        ! (nr, na), in Kg/m^3, 各组件平均冷却剂密度;
        real(KREAL):: max_Tfuel                               ! in K, 最热组件最大燃料温度;
        real(KREAL):: max_Tcoolant                            ! in K, 最热组件最大冷却剂温度;
        real(KREAL):: min_Rhocoolant                          ! in Kg/m^3, 最热组件最大冷却剂密度;
        real(KREAL):: last                                        ! in s, 上一时间点
        real(KREAL):: current                                     ! in s, 当前时间点
        real(KREAL):: toutlet                                 ! in K, 冷却剂出口平均温度
    
        real(KREAL), allocatable  :: power(:, :)
        real(KREAL), allocatable  :: fq_core(:, :)
        integer  :: nr, na, npin
        integer  :: ir, ia, ipin, itype
        integer  :: i_allocate
		
        integer M,N,i,j,n_zone

	!enddo
    contains
    subroutine Perform_TH_imp()
	  print*,'start th_imp calculation...'
     	M=size(assm1(1)%thermal%temperature,dim=1)
        N=size(assm1(1)%thermal%temperature,dim=2)
		allocate(assembly(ns%state%zone,ns%state%layer),Tfuel(ns%state%zone,ns%state%layer),Tcoolant(ns%state%zone,ns%state%layer),Rhocoolant(ns%state%zone,ns%state%layer))
        allocate(power(M,N),fq_core(M,N))
		do i=1,ns%state%zone,1
		  do j=1,ns%state%layer,1
		    assembly(i,j)=2.827*1e7 !W
		  enddo
		enddo
        fq_core=0.0
        power=0.0  
	 	
	 !do k=1,assm1%mesh%n_zone,1
	 !assm zone=1 组件1
	 n_zone=1
	   do i=1,assm1(1)%mesh%ny,1
          do j=1,N,1
           if(j<=assm1(1)%mesh%Nf) power(i,j)=assembly(n_zone,i)
          enddo
       enddo
	  fq_core=0.0
     print*,'finish setting pow,start to calculation...'
     call driving_imp_steady(assm1(1),power,fq_core)!power come from other data,so it should be an interface in place with the data
     print*,'finish the th_imp calculation,start to runout..'
     call Run_output() 

     end subroutine Perform_TH_imp
    end module test_TH2NK_steady

