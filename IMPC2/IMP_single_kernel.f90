module imp_single_kernel
    use constants
    use imp_mathkerel
    use imp_assembly_header
    implicit none
     private
     public::solve_momentum
     public::cal_momentumA
     public::solve_momentumA
     public::solve_pressureCorrection
     public::cal_pmodifyA
     public::solve_pmodifyA
     public::modify_PV
     public::cal_th_convection
     public::cal_th_temperature
     public::solve_temperature
     public::update_property
    contains
    subroutine solve_momentum(assm,flag,rhoi,ui,dt,ap)!(几何网格水力学)（p）（初始条件）（边界条件）（瞬态计算选项）
        type(sys_assembly),intent(in out)::assm
        real(KREAL),intent(in)::flag
        !real(KREAL),intent(in)::pguess(:)
        real(KREAL),intent(in)::rhoi(:)
        real(KREAL),intent(in)::ui(:)
        real(KREAL),intent(in)::dt
        real(KREAL),intent(in out)::ap(:)
        !local
        integer M!M=Ny-1,矢量控制体的个数
        real(KREAL),allocatable::A(:,:),b(:)
        M=assm%mesh%ny-1
        allocate(A(1:M,1:M),b(1:M))
        !pguess=assm%thermal%pressure
        call cal_momentumA(assm,flag,rhoi,ui,dt,A,b,ap)
        call solve_momentumA(M,A,b,assm%thermal%Velocity)

    end subroutine solve_momentum

    !subroutine cal_momentumA(N,f,De,rhoi,rho,uin,ui,ulast,pin,pout,pguess,dx,dt,A,b)!rho是当前迭代步的物性
subroutine cal_momentumA(assm,flag,rhoi,ui,dt,A,b,ap)
     implicit none
     type(sys_assembly),intent(in out)::assm
     real(KREAL),intent(in)::flag
     !real(KREAL),intent(in)::pguess(:)
     real(KREAL),intent(in)::rhoi(:),ui(:)
     real(KREAL),intent(in)::dt
     real(KREAL),intent(in out)::A(:,:)
     real(KREAL),intent(in out)::b(:)
     real(KREAL),intent(in out)::ap(:)
      !local
       integer N,i,k!N是矢量控制体的个数
       integer nr,nz!nr是径向的控制体个数，nz是轴向的标量控制体个数
       real(KREAL):: dx
       real(KREAL):: f,De
       real(KREAL):: uin,pin,pout
       real(KREAL),allocatable::ulast(:),pguess(:)
       real(KREAL),allocatable::rho(:)
       real(KREAL),allocatable::aw(:),ae(:),api(:),bs(:)
       N=assm%mesh%ny-1
       nr=assm%mesh%nf+assm%mesh%ng+assm%mesh%ns
       nz=assm%mesh%ny
       allocate(aw(1:N),ae(1:N),api(1:N),bs(1:N))
       allocate(rho(0:nz+1),ulast(1:N),pguess(1:N+1))
        !dx=assm%geom%height/assm%mesh%Ny
        f=assm%hydrau%fric
        de=assm%hydrau%de
        rho=assm%property%rho(:,nr+1)
        uin=assm%th_boundary%u%inlet
        
        pout=assm%th_boundary%p%outlet
        ulast=assm%thermal%velocity
        pguess=assm%thermal%Pressure
		!pin=assm%th_boundary%p%inlet
		pin=(3*pguess(1)-pguess(2))/2.0
       !dx rhoi,ulast,f,De,rho,uin,pin,api,pout
       !计算各个控制体的常系数和源项
       api=0.0
       do i=1,N,1	     
	     if(i==1)then
           k=i+assm%mesh%layer_bottom		 
		   dx=assm%geom%height(k)+assm%geom%height(k+1)/2
		 elseif(i>1.and.i<N) then	
           k=i+assm%mesh%layer_bottom		 
		   dx=assm%geom%height(k)
		 elseif(i==N)then
		   k=assm%mesh%Ny+assm%mesh%layer_top
		   dx=assm%geom%height(k-1)+assm%geom%height(k)
         endif		 
		 
         if(i==1)then
             aw(i)=0.0
             ae(i)=0.0
             if (flag==1.0) api(i)=1.50*dx/dt*(2*RHOI(i)+RHOI(i+1))/3.0
             ap(i)=api(i)+f/(2*De)*1.50*dx*ulast(i)*(2*RHO(i)+RHO(i+1))/3.0+RHO(0)*uin!replace rho(0)*uin with rho(1)*uin
             bs(i)=pin-pguess(i+1)-(2*RHO(i)+RHO(i+1))/3.0*9.8*1.5*dx+RHO(0)*uin*uin  !replace rho(0)*uin with rho(1)*uin           
         elseif(i>1.and.i<N)then
             aw(i)=ulast(i-1)*RHO(i)
             ae(i)=0.0
             if(flag==1.0) api(i)=dx/dt*(RHOI(i)+RHOI(i+1))/2.0
             ap(i)=aw(i)+api(i)+f/(2*De)*dx*ulast(i)*(RHO(i)+RHO(i+1))/2.0
             bs(i)=pguess(i)-pguess(i+1)-(RHO(i)+RHO(i+1))/2.0*9.8*dx
         elseif(i==N)then
             aw(i)=RHO(N)*ulast(N-1)
             ae(i)=0.0
             if(flag==1.0) api(i)=1.50*dx/dt*(RHOI(N)+2*RHOI(N+1))/3.0
             ap(i)=aw(i)+api(i)+f/(2*De)*1.50*dx*ulast(N)*(RHO(N)+2*RHO(N+1))/3.0
             bs(i)=pguess(i)-pout-(RHO(i)+2*RHO(i+1))/3.0*9.8*1.5*dx            
         endif
       enddo
       
      !用直接发求矩阵，第一步，写出系数矩阵
      A=0.0
       do i=1,N,1
           if(i==1)then
               A(i,i)=ap(i)
           elseif(i>1.and.i<N)then
               A(i,i)=ap(i)
               A(i,i-1)=-aw(i) 
           elseif(i==N)then
               A(i,i)=ap(i)
               A(i,i-1)=-aw(i)
           endif
           b(i)=bs(i)+api(i)*ui(i)
       enddo
    end subroutine cal_momentumA
    

    
subroutine solve_pressureCorrection(assm,flag,ap,rhoi,dt,pmodify,btotal)
    implicit none
    type(sys_assembly),intent(in out)::assm
    real(KREAL),intent(in)::flag
    real(KREAL),intent(in)::ap(:)
    real(KREAL),intent(in)::rhoi(:)
    real(KREAL),intent(in)::dt
    real(KREAL),intent(in out)::pmodify(:)
    real(KREAL),intent(in out)::btotal
    !local
    integer N
    real(KREAL),allocatable::A(:,:),b(:)
    N=size(pmodify)
    allocate(A(N-1,N-1),b(N-1))
    call cal_pmodifyA(assm,flag,ap,rhoi,dt,A,b,btotal)
    call solve_pmodifyA(A,b,pmodify)
end subroutine solve_pressureCorrection
    
subroutine cal_pmodifyA(assm,flag,ap,rhoi,dt,A,b,btotal)
 implicit none
 type(sys_assembly),intent(in out)::assm
 real(KREAL),intent(in)::flag
 real(KREAL),intent(in)::ap(:)
 real(KREAL),intent(in)::rhoi(:)
 real(KREAL),intent(in)::dt
 real(KREAL),intent(in out)::A(:,:)
 real(KREAL),intent(in out)::b(:)
 real(KREAL),intent(in out)::btotal
 !local
 integer Ny,i,nr,k
 real(KREAL):: uin,pout,dx
 real(KREAL),allocatable::rho(:),ulast(:)
 real(KREAL),allocatable::bp(:),be(:),bw(:),bb(:)

 Ny=assm%mesh%Ny
 nr=assm%mesh%nf+assm%mesh%ng+assm%mesh%ns
 allocate(RHO(0:Ny+1),ulast(1:Ny-1))
 allocate(bp(1:Ny),be(1:Ny),bw(1:Ny),bb(1:Ny))

 uin=assm%th_boundary%u%inlet
 pout=assm%th_boundary%p%outlet
 !dx=assm%geom%height/assm%mesh%Ny
 rho=assm%property%rho(:,nr+1)
 ulast=assm%thermal%velocity
        if(flag==0.0) then !steady
            dx=0.0
        endif
 
        do i=1,Ny,1	     
          if(flag==0.0)then
		     dx=0.0
		  elseif(flag==1.0)then
           k=i+assm%mesh%layer_bottom		 
		   dx=assm%geom%height(k)
		  endif

           if(i==1)then
             be(i)=1.5*(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=0.0
             bp(i)=be(i)
             bb(i)=RHO(0)*uin-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt
           elseif(i==2)then
             be(i)=(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=1.5*(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=be(i)+bw(i)
             bb(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt               
           elseif(i>2.and.i<Ny-1)then
             be(i)=(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=be(i)+bw(i)
             bb(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt    
           elseif(i==Ny-1)then
             be(i)=0.0
             bw(i)=(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=bw(i)+(RHO(i)+RHO(i+1))/2.0/ap(i)
             bb(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt  
           elseif(i==Ny)then
             be(i)=0.0
             bw(i)=0.0
             bp(i)=0.0
             bb(i)=0.0
           endif
       enddo
       
       A=0.0
       do i=1,Ny-1,1
           if(i==1)then
               A(i,i)=bp(i)
               A(i,i+1)=-be(i)
               b(i)=bb(i)
           elseif(i>1.and.i<Ny-1)then
               A(i,i)=bp(i)
               A(i,i+1)=-be(i)
               A(i,i-1)=-bw(i)
               b(i)=bb(i)
           elseif(i==Ny-1)then
               A(i,i)=bp(i)
               A(i,i-1)=-bw(i)
               b(i)=bb(i)
           endif
       enddo
       
      btotal=0.0
      do i=1,Ny,1
        btotal=btotal+abs(bb(i))
      enddo
end subroutine cal_pmodifyA



subroutine modify_PV(assm,ap,pmodify)
  implicit none
  type(sys_assembly),intent(in out)::assm
  real(KREAL),intent(in)::ap(:)
  real(KREAL),intent(in out)::pmodify(:)
  !local
  integer i,Ny
  real(KREAL):: alpha
  Ny=assm%mesh%Ny
  alpha=0.8
            do i=1,Ny-1,1
                 assm%thermal%pressure(i)=assm%thermal%pressure(i)+alpha*pmodify(i)
               if(i==1)then
                 assm%thermal%velocity(i)=assm%thermal%velocity(i)+1.5*(pmodify(i)-pmodify(i+1))/ap(i)
               elseif(i>1.and.i<Ny-1)then
                 assm%thermal%velocity(i)=assm%thermal%velocity(i)+(pmodify(i)-pmodify(i+1))/ap(i)
               elseif(i==Ny-1)then
               !assm%thermal%velocity(i)=assm%thermal%velocity(i)+(pmodify(i)-pout)/ap(i)!pmout=0.0
                 assm%thermal%velocity(i)=assm%thermal%velocity(i)+pmodify(i)/ap(i)
               endif
            enddo
            assm%th_boundary%u%outlet=assm%thermal%velocity(Ny-1)
            assm%th_boundary%p%inlet=1.50*assm%thermal%pressure(1)-0.5*assm%thermal%pressure(2)
            assm%thermal%pressure(Ny)=(2*assm%th_boundary%p%outlet+assm%thermal%pressure(Ny-1))/3.0
end subroutine modify_PV

subroutine cal_th_convection(assm)
 implicit none
 type(sys_assembly),intent(in out)::assm
 !lcoal
 integer i,Ny,nr
 real(KREAL):: De,Area,wet,velocity!单纯输入的变量可以用局部变量来替换
 real(KREAL),allocatable::rho(:),dvs(:),shc(:),ctc(:)
 Ny=assm%mesh%Ny
 nr=assm%mesh%nf+assm%mesh%ng+assm%mesh%ns
 allocate(rho(0:Ny+1),dvs(0:Ny+1),shc(0:Ny+1),ctc(0:Ny+1))
 de=assm%hydrau%de
 Area=assm%hydrau%aflow
 wet=assm%hydrau%wet
 rho=assm%property%rho(:,nr+1)
 dvs=assm%property%dvs(:,nr+1)
 shc=assm%property%shc(:,nr+1)
 ctc=assm%property%ctc(:,nr+1)
 
     do i=1,Ny,1
      if (i==1)then
          velocity=assm%thermal%velocity(1)
      elseif(i>1.and.i<Ny)then
          velocity=(assm%thermal%velocity(i-1)+assm%thermal%velocity(i))/2.0
      else
          velocity=assm%thermal%velocity(i-1)
      endif
      call get_convection(De,Area,wet,RHO(i),velocity,DVS(i),SHC(i),CTC(i),assm%property%htc(i))!DVS(i,N)动力粘度 Pa*s
    enddo
    assm%property%htc(0)=assm%property%htc(1)!边界上的对流换热系数
    assm%property%htc(Ny+1)=assm%property%htc(Ny)
end subroutine cal_th_convection
    
subroutine cal_th_temperature(assm,flag,Ti,rhoi,dt)
 implicit none
 type(sys_assembly),intent(in out)::assm
 real(KREAL),intent(in)::flag
 real(KREAL),intent(in)::Ti(:,:),rhoi(:,:)
 real(KREAL),intent(in)::dt
 !local
    real(KREAL):: Area,Xt,xf,xg,xs,Df,Dg,Ds,Dy,uin,Tin
    integer  M,N,i,j,k,Nf,Ng,Ns,Ny
    real(KREAL),allocatable::Tj(:,:)
    real(KREAL),allocatable::RHO(:,:),SHC(:,:),CTC(:,:),DVS(:,:)
    real(KREAL),allocatable::aw(:,:),ae(:,:),ap(:,:),as(:,:),an(:,:),api(:,:),bs(:,:),q(:,:)
     Area=assm%hydrau%aflow
     uin=assm%th_boundary%u%inlet
     Tin=assm%th_boundary%T%inlet
     xf=assm%geom%pellet
     xg=assm%geom%Bond
     xs=assm%geom%Cladth
     Xt=Xf+Xg+Xs !包壳外径
     Nf=assm%mesh%Nf
     Ng=assm%mesh%Ng
     Ns=assm%mesh%Ns
     Ny=assm%mesh%Ny
     M=assm%mesh%Ny+1
     N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
     Df=Xf/Nf
     Dg=Xg/Ng
     Ds=Xs/Ns
     !Dy=assm%geom%height/Ny
     
    allocate(Tj(1:M-1,1:N))
    allocate(RHO(0:M,0:N),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N))
    allocate(aw(1:M-1,1:N),ae(1:M-1,1:N),ap(1:M-1,1:N),as(1:M-1,1:N),an(1:M-1,1:N),api(1:M-1,1:N),bs(1:M-1,1:N),q(1:M-1,1:N))
    rho=assm%property%rho
    shc=assm%property%shc
    ctc=assm%property%ctc
    dvs=assm%property%dvs
    !set the pow
    q=assm%pow%power
    
    api=0.0
    Do i=1,M-1,1
	    Dy=assm%geom%height(i+assm%mesh%layer_bottom)
        Do j=1,N,1
         if (j==1)then!轴对称边界的控制体
          aw(i,j)=0.0
          ae(i,j)=(assm%mesh%r(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Df/dt
          ap(i,j)=ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Df*q(i,j)
         elseif (j>1.and.j<Nf)then!fuel内部控制体
          aw(i,j)=(assm%mesh%r(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=(assm%mesh%r(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Df*q(i,j)
         elseif(j==Nf)then!f-g边界左侧控制体
          aw(i,j)=(assm%mesh%r(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=2*(assm%mesh%r(i,j)+Df/2.0)*(Df+Dg)/(Df/CTC(i,j)+Dg/CTC(i,j+1))/(Df+Dg)
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Df*q(i,j)  
         elseif (j==Nf+1)then!f-g边界右侧控制体
          aw(i,j)=2*(assm%mesh%r(i,j)-Dg/2.0)*(df+dg)/(df/CTC(i,j-1)+dg/CTC(i,j))/(df+dg)
          ae(i,j)=(assm%mesh%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Dg*q(i,j)
         elseif(j>Nf+1.and.j<Nf+Ng)then!g气隙内部控制体
          aw(i,j)=(assm%mesh%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=(assm%mesh%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng)then!g-c边界左侧控制体
          aw(i,j)=(assm%mesh%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=2*(assm%mesh%r(i,j)+Dg/2.0)*(Dg+Ds)/(Dg/CTC(i,j)+Ds/CTC(i,j+1))/(Dg+Ds)
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng+1)then!g-c边界右侧控制体
          aw(i,j)=2*(assm%mesh%r(i,j)-Ds/2.0)*(dg+ds)/(dg/CTC(i,j-1)+ds/CTC(i,j))/(dg+ds)
          ae(i,j)=(assm%mesh%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Ds*q(i,j)
         elseif(j>Nf+Ng+1.and.j<Nf+Ng+Ns)then!c包壳内部控制体
          aw(i,j)=(assm%mesh%r(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(assm%mesh%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns)then!s-fluid边界左侧控制体
          aw(i,j)=(assm%mesh%r(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(assm%mesh%r(i,j)+Ds/2.0)/(1.0/assm%property%htc(i)+ds/(2.0*CTC(i,j)))
          as(i,j)=0.0
          an(i,j)=0.0
          if(flag==1.0) api(i,j)=RHO(i,j)*SHC(i,j)*assm%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=assm%mesh%r(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns+1)then!fluid控制体
          
          if(i==1)then!流体入口的控制体
           aw(i,j)=Dy/SHC(i,j)*2.0*3.14*Xt/Area*1.0/(1.0/assm%property%htc(i)+Ds/(2*CTC(i,j-1)))
           ae(i,j)=0.0
           as(i,j)=0.0
           an(i,j)=0.0
           if(flag==1.0) api(i,j)=RHOI(i,j)*Dy/dt
           ap(i,j)=aw(i,j)+api(i,j)+RHO(i-1,j)*uin
           bs(i,j)=Dy/SHC(i,j)*q(i,j)+RHO(i-1,j)*uin*assm%th_boundary%T%inlet          
          else!流体内部以及出口控制体
           aw(i,j)=Dy/SHC(i,j)*2.0*3.14*Xt/Area*1.0/(1.0/assm%property%htc(i)+Ds/(2*CTC(i,j-1)))
           ae(i,j)=0.0
           as(i,j)=0.0
           an(i,j)=0.5*(RHO(i,j)+RHO(i-1,j))*assm%thermal%velocity(i-1)
           if(flag==1.0) api(i,j)=RHOI(i,j)*Dy/dt
           ap(i,j)=an(i,j)+aw(i,j)+api(i,j)
           bs(i,j)=Dy/SHC(i,j)*q(i,j)        
          endif
         endif              
      enddo
    enddo

     
     do k=1,5000,1
       do i=1,M-1,1
           do j=1,N,1
               if(k==1) then
                Tj(i,j)=Ti(i,j)
               else
                if(j==1)then
                 Tj(i,j)=(ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j>1.and.j<N)then
                 Tj(i,j)=(aw(i,j)*Tj(i,j-1)+ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j==N) then
                  if(i==1)then
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                  else
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+an(i,j)*Tj(i-1,j)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)  
                 endif
                endif
              endif
           enddo
       enddo
     enddo
     
     do i=1,M-1,1!as for the solid,no need to know the inlet and outlet temperature
         do j=1,N,1
            !if(i>0.and.i<M)then
            !    if(j==0)then
            !      assm%thermal%Temperature(i,j)=Tj(i,j+1)
            !    else
                  assm%thermal%Temperature(i,j)=Tj(i,j)
            !    endif
            !elseif(i==0)then
            !    if(j==0)then
            !      assm%thermal%Temperature(i,j)=Tj(i+1,j+1)
            !    elseif(j>0.and.j<N)then
            !      assm%thermal%Temperature(i,j)=Tj(i+1,j)
            !    elseif(j==N)then!入口
            !      assm%thermal%Temperature(i,j)=Tin
            !    endif
           !elseif(i==M)then
           !     if(j==0)then
           !       assm%thermal%Temperature(i,j)=Tj(i-1,j+1)
           !     else
           !       assm%thermal%Temperature(i,j)=Tj(i-1,j)
           !     endif
            !endif           
         enddo
     enddo
	 assm%th_boundary%T%outlet=assm%thermal%Temperature(M-1,N)
end subroutine cal_th_temperature

subroutine solve_temperature(assm,flag,Ti,rhoi,dt)
 implicit none
 type(sys_assembly),intent(in out)::assm
 real(KREAL),intent(in)::flag
 real(KREAL),intent(in)::Ti(:,:)
 real(KREAL),intent(in)::rhoi(:,:)
 real(KREAL):: dt
 !local
 call cal_th_convection(assm)
 call cal_th_temperature(assm,flag,Ti,rhoi,dt)
end subroutine solve_temperature


subroutine update_property(assm,drho)
  implicit none
  type(sys_assembly),intent(in out)::assm
  real(KREAL):: drho
  !local
  integer i,M,N,Ny
  real(KREAL),allocatable::rhof(:)!用于存放上一迭代步的密度
  Ny=assm%mesh%Ny
  M=Ny+1
  N=assm%mesh%Nf+assm%mesh%Ng+assm%mesh%Ns+1
  allocate(rhof(0:M))
      do i=0,M,1
        RHOF(i)=assm%property%rho(i,N)
      enddo
      !print*,RHOF
      do i=1,M-1,1
         assm%property%rho(i,N)=get_density(assm%thermal%Temperature(i,N))
      enddo
      !print*,'temperature=',assm%thermal%Temperature(:,N)
      assm%property%rho(0,N)=assm%property%rho(1,N)
      assm%property%rho(M,N)=assm%property%rho(M-1,N)
     ! print*,'rho=',assm%property%rho(:,N)
      drho=0.0
      do i=0,M,1
        drho=drho+abs((assm%property%rho(i,N)-RHOF(i))/RHOF(i))
        !print*,'drho=',drho
      enddo
      !print*,'drho=',drho
end subroutine update_property
    subroutine solve_momentumA(N,A,b,u)
     implicit none
     integer,intent(in)::N
     real(KREAL),intent(in)::A(:,:)
     real(KREAL),intent(in)::b(:)
     real(KREAL),intent(in out)::u(:)
     !local
     integer i
     real(KREAL),allocatable::uc(:)
     allocate(uc(1:N))
     
     call tdma(N,A,b,uc)
     
     do i=1,N,1
        if(uc(i)>0.0)then
        u(i)=uc(i)
        else
        u(i)=0.0
        endif
      enddo
    end subroutine solve_momentumA
    
subroutine solve_pmodifyA(A,b,pmodify)
  implicit none
  real(KREAL),intent(in)::A(:,:),b(:)
  real(KREAL),intent(in out)::pmodify(:)
  !local
  integer M,N,i
  real(KREAL),allocatable::pmm(:)
  M=size(b)
  N=size(pmodify)
  allocate(pmm(M))
  call tdma(M,A,b,pmm)
  do i=1,N,1
     if(i<N)then
        pmodify(i)=pmm(i)
      else
        pmodify(i)=0.0
      endif
   enddo
end subroutine solve_pmodifyA
end module imp_single_kernel