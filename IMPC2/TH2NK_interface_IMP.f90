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
    module TH2NK_interface_IMP
	 use constants
	 use th_global
	 
     use imp_re_input_global
     use imp_assm_global
     use imp_driving_pre_process
     use imp_driving_post_process
     use imp_power_header
     use imp_single_channel
	 use imp_property
    implicit none

     !real(KREAL),allocatable::power(:,:),fq_core(:,:)
     !integer M,N,i,j
    contains
    subroutine Perform_TH_imp(transient_flag, assembly, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)
		logical, intent(in)      :: transient_flag                              ! .TRUE. --transient
        real(KREAL), intent(in)  :: assembly(:, :)                              ! (zone, layer), in W, 各组件功率;
        real(KREAL), intent(in out)  :: Tfuel(:, :)                             ! (nr, na), in K, 各组件平均燃料温度;
        real(KREAL), intent(in out)  :: Tcoolant(:, :)                          ! (nr, na), in K, 各组件平均冷却剂温度;
        real(KREAL), intent(in out)  :: Rhocoolant(:, :)                        ! (nr, na), in Kg/m^3, 各组件平均冷却剂密度;
        real(KREAL), intent(in out)  :: max_Tfuel                               ! in K, 最热组件最大燃料温度;
        real(KREAL), intent(in out)  :: max_Tcoolant                            ! in K, 最热组件最大冷却剂温度;
        real(KREAL), intent(in out)  :: min_Rhocoolant                          ! in Kg/m^3, 最热组件最大冷却剂密度;
        real(KREAL), intent(in)  :: last                                        ! in s, 上一时间点
        real(KREAL), intent(in)  :: current                                     ! in s, 当前时间点
        real(KREAL), intent(in out)  :: toutlet                                 ! in K, 冷却剂出口平均温度
    
        !local
        REAL(KREAL)  :: last_, current_
		real(KREAL)  :: volumn,TVtotal,dr!used to calculate the average temperature of fuel
		real(KREAL)  :: TAoutlet_total,A_total!used to calculate the core average outlet temperature
		real(KREAL)	 :: xs,xg,xf !used to calculate the Tsurface
		real(KREAL)  :: max_T_inner,max_T_outer
        real(KREAL), allocatable  :: power(:, :)
        real(KREAL), allocatable  :: fq_core(:, :)
		
		integer  :: nf,ng,ns,nRadial,ny
        integer  :: nr, na, npin,M,N
        !integer  :: ir, ia, ipin, itype
        !integer  :: i_allocate
		integer  :: i,j,k
        
        last_ = last
        current_ = current
        nr = SIZE(assembly, dim=1)                                              ! 径向的组件数目zone
        na = SIZE(assembly, dim=2)                                              ! 轴向的节块数目layer
        
    	M=size(assm1(1)%thermal%temperature,dim=1)
        N=size(assm1(1)%thermal%temperature,dim=2)
		!allocate(assembly(nr,na),Tfuel(nr,na),Tcoolant(nr,na),Rhocoolant(nr,na))
        allocate(power(M,N),fq_core(M,N))
		fq_core=1.0D0
        power=0.0 

		do i=1,nr,1
		  do j=1,na,1
		    imp_pow(i,j)=assembly(i,j)!W
		  enddo
		enddo
 
	 !do k=1,assm1%mesh%n_zone,1
	 !assm zone=1 组件1
   !as the transient calculate is after steady calculate,so the initial condition is useless
   !the most important is the boundary condition include inlet temperature/inlet velocity and outlet pressre
   !convert the core boundary to the assembly boundary
   !set the inlet temperature
	 do i=1,nr,1
		assm1(i)%th_boundary%T%inlet=design%tcoolin
	 enddo
   !allocate the flowrate
	 call driving_imp_flowAlloc(assm1,design%assembly_flow)

	 
	 do i=1,nr,1!zone start
	  !if(i==19) then
       print*,'zone=',i
	   do j=1,assm1(i)%mesh%ny,1!dy
          do k=1,N,1
		   !covert the power from zone total W to pin w/m^3 
           !print*,'assembly=',assembly(i,j+assm1(i)%mesh%layer_bottom),'height=',assm1(i)%geom%height(j)
           if(k<=assm1(i)%mesh%Nf) power(j,k)=assembly(i,j+assm1(i)%mesh%layer_bottom)/(assm1(i)%geom%N_fuelpin*assm1(i)%geom%height(j)*3.14159*assm1(i)%geom%pellet**2)
          enddo
       enddo
       
       !if(i==50) then 
       !    print*,power
       !endif
	   
	   if (assm1(i)%th_boundary%u%inlet==0.0) then
		  assm1(i)%thermal%velocity=0.0
		  assm1(i)%th_boundary%u%outlet=0.0
		  assm1(i)%thermal%temperature=assm1(i)%th_boundary%T%inlet
		  assm1(i)%th_boundary%T%outlet=assm1(i)%th_boundary%T%inlet
		  assm1(i)%property%rho=get_density(assm1(i)%th_boundary%T%inlet)
	   else
	     if (transient_flag)  then
              call driving_imp_transient(assm1(i),power, fq_core,last_, current_)
          else
              call driving_imp_steady(assm1(i),power,fq_core)
         end if
	   endif	
	   dr=assm1(i)%geom%pellet/assm1(i)%mesh%Nf
	   do j=1,assm1(i)%mesh%Ny,1
	    volumn=0.0
		TVtotal=0.0
	    do k=1,assm1(i)%mesh%Nf,1 !rod average		
			if (k==1) then
				TVtotal=TVtotal+assm1(i)%thermal%temperature(j,k)*3.14*(k*dr)**2*assm1(i)%geom%height(j)
				volumn=volumn+3.14*(k*dr)**2*assm1(i)%geom%height(j)
			else
				TVtotal=TVtotal+assm1(i)%thermal%temperature(j,k)*3.14*((k*dr)**2-((k-1)*dr)**2)*assm1(i)%geom%height(j)
				volumn=volumn+3.14*((k*dr)**2-((k-1)*dr)**2)*assm1(i)%geom%height(j)
			endif
		enddo
		Tfuel(i,j+assm1(i)%mesh%layer_bottom)=TVtotal/volumn
		Tcoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%thermal%temperature(j,N)
		!print*,'Tcoolant=',Tcoolant(i,j+assm1(i)%mesh%layer_bottom)
		Rhocoolant(i,j+assm1(i)%mesh%layer_bottom)=assm1(i)%property%rho(j,N)
       enddo
	  !endif
	 enddo!zone end
	 !calculate the average toutlet
	 TAoutlet_total=0.0
	 A_total=0.0
     do i=1,nr,1
		!volumn average
		TAoutlet_total=TAoutlet_total+assm1(i)%th_boundary%T%outlet*assm1(i)%hydrau%aflow*assm1(i)%th_boundary%u%outlet
		A_total=A_total+assm1(i)%hydrau%aflow*assm1(i)%th_boundary%u%outlet
	 enddo
	 toutlet=TAoutlet_total/A_total
	 max_Tcoolant=0.0
	 max_Tfuel=0.0
	  !calculate max_Tcoolant max_Tfuel
	 do i=1,nr,1
		do j=1,assm1(i)%mesh%Ny,1
			if(Tfuel(i,j)>max_Tfuel) 		max_Tfuel=Tfuel(i,j)
			if(Tcoolant(i,j)>max_Tcoolant)	max_Tcoolant=Tcoolant(i,j)
		enddo
	 enddo
	 !calculate the surface temperature
	 do i=1,nr,1
	 	Nf=assm1(i)%mesh%Nf
		Ng=assm1(i)%mesh%Ng
		Ns=assm1(i)%mesh%Ns
		Ny=assm1(i)%mesh%Ny
		Nradial=Nf+Ng+Ns+1
		xf=assm1(i)%geom%pellet
		xg=assm1(i)%geom%bond
		xs=assm1(i)%geom%cladth		
		do j=1,assm1(i)%mesh%Ny,1
		  !if (assm1(i)%th_boundary%u%inlet==0.0) then!这种可能性可以排除
			!assm1(i)%thermal%Tfg(j)=assm1(i)%thermal%temperature(j,Nf+1)
			!assm1(i)%thermal%Tgs(j)=assm1(i)%thermal%temperature(j,Nf+Ng+1)
			!assm1(i)%thermal%Tsc(j)=assm1(i)%thermal%temperature(j,Nradial)
		  !else
			assm1(i)%thermal%Tcoolant(j)=Tcoolant(i,j)
			assm1(i)%thermal%Tfuel(j)=Tfuel(i,j)
			assm1(i)%thermal%Tfuel_center(j)=assm1(i)%thermal%temperature(j,1)
			assm1(i)%thermal%Tfg(j)=(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf)*assm1(i)%thermal%temperature(j,Nf+1))/(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf))!芯块外边界
			assm1(i)%thermal%Tgs(j)=(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)*assm1(i)%thermal%temperature(j,Nf+Ng)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf+Ng+1))/(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng))!包壳内边界
			assm1(i)%thermal%Tsc(j)=(assm1(i)%property%htc(j)*assm1(i)%thermal%temperature(j,Nradial)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns)*assm1(i)%thermal%temperature(j,Nradial-1))/(assm1(i)%property%htc(j)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns))!包壳外边界
		  !endif
		enddo
	  enddo!surface zone end
	 !write current_,max_Tcoolant,toutlet,max_Tfuel,max_T_inner,max_T_outer
		max_T_inner=0.0
		max_T_outer=0.0
		do i=1,nr,1
			do j=1,assm1(i)%mesh%Ny,1
				if(max_T_inner<assm1(i)%thermal%Tgs(j)) max_T_inner=assm1(i)%thermal%Tgs(j)
				if(max_T_outer<assm1(i)%thermal%Tsc(j)) max_T_outer=assm1(i)%thermal%Tsc(j)
			enddo
        enddo
       
       !open(5,file='.\output\selftimelist.txt',STATUS='new')
       !close(5)
	   !open(6,file='.\output\selftimelist.txt',STATUS='OLD',ACCESS='append')
       !write(6,100) current_,max_Tcoolant,toutlet,max_Tfuel,max_T_inner,max_T_outer
	   !100 Format(1x,F8.1,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)
       ! !write(2,*) assm1%pow%power
       !close(6) 
	  ! open(6,file='.\output\Tfuel.txt')
      ! write(6,*) Tfuel
      ! !write(2,*) assm1%pow%power
      ! close(6)  
	  ! open(7,file='.\output\Tcoolant.txt')
      ! write(7,*) Tcoolant
      ! !write(2,*) assm1%pow%power
      ! close(7) 
	  
	  if (allocated(power))       deallocate(power)
      if (allocated(fq_core))     deallocate(fq_core)
      
     !*********************************************
     !power come from other data,so it should be an interface in place with the data

     !do while(timer1%ctime<timer1%ttotal) 
     ! timer1%ctime=timer1%ctime+timer1%dt
     ! call update_power(power,fq_core,timer1%ltime,timer1%ctime)
     ! call driving_imp_transient(assm1,power, fq_core,timer1%ltime,timer1%ctime)
     ! call timer1%record(assm1%th_boundary%T%outlet,assm1%th_boundary%u%inlet,power(1,1))
     !print*,'ctime=',timer1%ctime
     ! timer1%ltime=timer1%ctime
     !enddo
   
     !call Run_output() 
     
	 !read(*,*)
     end subroutine Perform_TH_imp
    end module TH2NK_interface_IMP

