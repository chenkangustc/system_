module imp_driving_pre_process
    use constants
	use global 
    use global_state
	use imp_assm_global
    use imp_re_input_global
	
	
    implicit none
    private
    public::Sys_pre_process
    public::alloc_assembly
    public::free_assembly
    public::set_assembly
    public::init_assembly
    public::cal_grid
contains
    subroutine Sys_pre_process()
     implicit none
	 integer i,layer_core
     !real(KREAL)::ttotal
     !real(KREAL)::ltime
     !real(KREAL)::ctime
     !integer::Nt
     write(*,*)'start the sys pre process:'
     !alloc
	  layer_core=ns%state%layer-ns%state%layer_bottom-ns%state%layer_top
		if(ns%feedback%is_feedback .and. ns%feedback%is_inner) then
			allocate(assm1(ns%state%zone))
			allocate(imp_pow(ns%state%zone,ns%state%layer))
			do i=1,ns%state%zone,1
				!if(allocated(assm1(i)%geom%height))  deallocate(assm1(i)%geom%height)
				allocate(assm1(i)%geom%height(layer_core))
			enddo
		endif
     !读取参数
     call reInputdata%set()
     !call reInputdata%publish()
     !参数赋值
	 do i=1,ns%state%zone,1
         !print*,assm1(i)
		call set_assembly(assm1(i),reInputdata,ns%state%zone,ns%state%layer,ns%state%layer_bottom,ns%state%layer_top,geom%height)
		print*,'set nf ng ns...'
		print*,'nf=',assm1(i)%mesh%nf,'ng=',assm1(i)%mesh%ng,'ns=',assm1(i)%mesh%ns
		!ttotal=150.0
		!Nt=150
		!call timer1%set(ttotal,Nt)!(ttotal,Nt)
		!分配空间
		call alloc_assembly(assm1(i),ns%state%layer)
		!call timer1%alloc()
		!初始化
		call init_assembly(assm1(i))
	 enddo
	  print*,assm1(1)%geom%rod,assm1(1)%geom%pellet,assm1(1)%geom%bond,assm1(1)%geom%cladth,assm1(1)%geom%pitch,assm1(1)%geom%n_pin,assm1(1)%geom%n_fuelpin

     !    timer1%init(ttotal,Nt,ctime,ltime)
     !ltime=0.0
     !ctime=0.0
     !call timer1%init(ctime,ltime)
    end subroutine Sys_pre_process
    
    subroutine init_assembly(assm)
      implicit none
      type(sys_assembly),intent(in out)::assm
      write(*,*)'init assembly'
      !热物性初始化
      call assm%property%init(assm%mesh%Nf,assm%mesh%Ng,assm%mesh%Ns,assm%mesh%Ny)
      !热工参数初始化
      call assm%Thermal%init(assm%initdata%Ti,assm%initdata%Pi,assm%initdata%ui)
      !边界条件初始化
      call assm%th_boundary%init(assm%initdata%Tin,assm%initdata%uin,assm%initdata%pout)
      assm%th_boundary%p%inlet=assm%thermal%pressure(1)+2500.0
      !热源初始化
      !write(*,*)'init power'
      !assm%pow%power=0.0
      !assm%pow%fq_core=0.0
      !网格
	  assm%geom%pellet=assm%geom%rod-assm%geom%cladth-assm%geom%bond
      call cal_grid(assm)
      call assm%hydrau%cal(assm%geom%pellet,assm%geom%pd)
     endsubroutine init_assembly
    
     subroutine alloc_assembly(assm,layer)
      implicit none
      type(sys_assembly),intent(in out)::assm
	  integer,intent(in)::layer
      !local
      integer::i_allocate
      integer::M,N,Ny
      Ny=assm%mesh%Ny
      M=Ny+1
      N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
      !check allocated first
      call Free_assembly(assm)
      
	  
      allocate(assm%property%rho(0:M,0:N))!(1:M,1:N)
      allocate(assm%property%shc(0:M,0:N))
      allocate(assm%property%ctc(0:M,0:N))
	  allocate(assm%property%dvs(0:M,0:N))
      allocate(assm%property%htc(0:M))
      
      allocate(assm%thermal%Temperature(M-1,N))
      allocate(assm%thermal%Pressure(M-1))
      allocate(assm%thermal%Velocity(Ny-1))
	  allocate(assm%thermal%Tcoolant(Ny))
	  allocate(assm%thermal%Tfuel(Ny))
	  allocate(assm%thermal%Tfuel_center(Ny))
	  allocate(assm%thermal%Tfg(Ny))
	  allocate(assm%thermal%Tgs(Ny))
	  allocate(assm%thermal%Tsc(Ny))
      
      allocate(assm%mesh%r(0:M,0:N))
      allocate(assm%mesh%z(0:M,0:N))
      
      allocate(assm%pow%power(M-1,N))
      allocate(assm%pow%fq_core(M-1,N))
      write(*,*)'alloc data array'
     end subroutine alloc_assembly
     
     subroutine Free_assembly(assm)
      implicit none
      type(sys_assembly),intent(in out)::assm
	  !if(allocated(assm%geom%height))  deallocate(assm%geom%height)
	  
      if(allocated(assm%property%rho))  deallocate(assm%property%rho)
      if(allocated(assm%property%shc))  deallocate(assm%property%shc)
      if(allocated(assm%property%ctc))  deallocate(assm%property%ctc)
      if(allocated(assm%property%htc))  deallocate(assm%property%htc)
      
      if(allocated(assm%thermal%temperature))  deallocate(assm%thermal%temperature)
      if(allocated(assm%thermal%pressure))  deallocate(assm%thermal%pressure)
      if(allocated(assm%thermal%Velocity))  deallocate(assm%thermal%Velocity)
	  if(allocated(assm%thermal%Tcoolant))  deallocate(assm%thermal%Tcoolant)
	  if(allocated(assm%thermal%Tfuel))  deallocate(assm%thermal%Tfuel)
	  if(allocated(assm%thermal%Tfuel_center))  deallocate(assm%thermal%Tfuel_center)
	  if(allocated(assm%thermal%Tfg))  deallocate(assm%thermal%Tfg)
	  if(allocated(assm%thermal%Tgs))  deallocate(assm%thermal%Tgs)
	  if(allocated(assm%thermal%Tsc))  deallocate(assm%thermal%Tsc)
      
      if(allocated(assm%mesh%r))  deallocate(assm%mesh%r)
      if(allocated(assm%mesh%z))  deallocate(assm%mesh%z)
      
      if(allocated(assm%pow%power))  deallocate(assm%pow%power)
      if(allocated(assm%pow%fq_core))  deallocate(assm%pow%fq_core)
     end subroutine Free_assembly
     
     subroutine set_assembly(assm,reInputdata,zone,layer,layer_bottom,layer_top,height)
      implicit none
      type(sys_assembly),intent(in out)::assm
      type(sys_re_input),intent(in)::reInputdata
	  integer,intent(in)::zone
	  integer,intent(in)::layer
	  integer,intent(in)::layer_bottom
	  integer,intent(in)::layer_top
	  real(KREAL),intent(in)::height(:)
	  !local
	  integer layer_core
	  layer_core=layer-layer_bottom-layer_top
	  !integer n_start
	  !integer n_end
	  !n_start=layer_bottom+1
	  !n_end=layer-layer_top
	  !real,intent(in)::height(:)
      write(*,*)'set assmebly as below:'
      !设置几何参数
	  !print*,geom%height
      call assm%geom%set(reInputdata%xf,reInputdata%xg,reInputdata%xs,reInputdata%acf,Height,reInputdata%pd,reInputdata%nFuelPin,reInputdata%npin)
      print*,'assm height(1)=',assm%geom%height(1)
	  !设置网格参数
      call assm%mesh%set(reInputdata%nf,reInputdata%ng,reInputdata%ns,zone,layer_core,layer_bottom,layer_top)
      !设置初始值
      call assm%initdata%set(reInputdata%Ti,reInputdata%Pi,reInputdata%Ui,reInputdata%Tin,reInputdata%Pin,reInputdata%Uin)
      !设置收敛因子
      call assm%confactor_%set(reInputdata%alpha,reInputdata%sigma)
      call assm%hydrau%set(reInputdata%f)
     end subroutine set_assembly
     
     subroutine cal_grid(assm)
       implicit none
       type(sys_assembly),intent(in out)::assm
       !local
       real(KREAL):: Df,Dg,Ds,Dy 
       integer  M,N,i,j,k
     write(*,*)'calculate the grid value...'
     Df=assm%geom%pellet/assm%mesh%Nf
     Dg=assm%geom%Bond/assm%mesh%Ng
     Ds=assm%geom%Cladth/assm%mesh%Ns
     !Dy=assm%geom%Height/assm%mesh%Ny
     M=assm%mesh%Ny+1
     N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
     !print*,'after M, N ,ok'
     Do i=0,M,1
        !print*,i !Dy=1.0
	   if(i>0.and.i<M-1) then 
	     !print*,'1>0 1<M-1'
	     k=assm%mesh%layer_bottom+i
		 !dy=0.01
		 !print*,'k=',k
		  !Dy=assm%geom%height(k)
          !print*,dy
		 !Dy=assm%geom%rod
		 Dy=geom%height(k)/100 !convert the unit from cm to m
		 !print*,'dy=',dy
       elseif(i==0)then 
         Dy=0.0
       elseif(i==M)then
         Dy=0.0
       endif		 
         do j=0,N,1
            if (j==0)then
               assm%mesh%r(i,j)=0.0
            elseif(j==1)then
               assm%mesh%r(i,j)=Df/2.0
            elseif(j>1.and.j<=assm%mesh%Nf) then
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+Df
            elseif(j==assm%mesh%Nf+1)then
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+Df/2.0+Dg/2.0
            elseif (j>assm%mesh%Nf+1.and.j<=assm%mesh%Nf+assm%mesh%Ng) then
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+Dg
            elseif (j==assm%mesh%Nf+assm%mesh%Ng+1)then
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+ Dg/2.0+Ds/2.0
            elseif (j>assm%mesh%Nf+assm%mesh%Ng+1.and.j<=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns)then
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+Ds
            else!流体的径向坐标，没有实际意义
               assm%mesh%r(i,j)=assm%mesh%r(i,j-1)+Ds
            endif
            
            if(i==0)then
              assm%mesh%z(i,j)=0.0
            elseif (i==1)then
              assm%mesh%z(i,j)=Dy/2.0
            elseif (i>1.and.i<M)then
              assm%mesh%z(i,j)=assm%mesh%z(i-1,j)+Dy
            elseif (i==M)then
              assm%mesh%z(i,j)=assm%mesh%z(i-1,j)+Dy/2.0
            endif
			!print*,'r z ok '	 
         enddo
		 
      enddo  
       
     end subroutine cal_grid   

end module imp_driving_pre_process