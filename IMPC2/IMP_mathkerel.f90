module imp_mathkerel
    use constants
    implicit none
    private
    public::tdma
    public::get_convection
    public::get_Nusselt
    public::get_hyconstant
    
contains
    subroutine tdma(N,A,B,u)
      real(KREAL):: A(N,N)
      real(KREAL):: B(N)
      real(KREAL):: u(N)
      integer i,N
      real(KREAL),dimension(N)::aa,bb,cc,dd,x,y

      do i=1,N,1
          if(i==1)then
          aa(1)=0.0
          bb(1)=A(1,1)
          cc(1)=A(1,2)
          dd(1)=B(1)
          elseif (i==N) then
           aa(N)=A(N,N-1)
           bb(N)=A(N,N)
           cc(N)=0.0
           dd(N)=B(N)
          else
           aa(i)=A(i,i-1)
           bb(i)=A(i,i)
           cc(i)=A(i,i+1)
           dd(i)=B(i)
          endif
      enddo
      
      x(1)=cc(1)/bb(1)
      y(1)=dd(1)/bb(1)
      
      do i=2,N,1
          x(i)=cc(i)/(bb(i)-x(i-1)*aa(i))
          y(i)=(dd(i)-y(i-1)*aa(i))/(bb(i)-x(i-1)*aa(i))
      enddo
      
      u(N)=y(N)
      
      do i=N-1,1,-1
          u(i)=y(i)-x(i)*u(i+1)
      enddo
    endsubroutine tdma    
    
    subroutine get_convection(lenth,flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,convection) !lenth是特征长度,取水力学直径  Nu=(h*x)/k
     real(KREAL):: lenth
     real(KREAL):: flow_area
     real(KREAL):: wetted_perimeter
     real(KREAL):: density
     real(KREAL):: velocity
     real(KREAL):: De
     real(KREAL):: viscosity,capacity,convection,conductivity
     real(KREAL):: Pr,Re,Pe,Nu
     
     call get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)  
     convection=Nu*conductivity/lenth
    end subroutine get_convection
    
    subroutine get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)!液态重金属的努赛尔数
     real(KREAL):: flow_area
     real(KREAL):: wetted_perimeter    
     real(KREAL):: density
     real(KREAL):: velocity
     real(KREAL):: De
     real(KREAL):: Re 
     real(KREAL):: viscosity
     real(KREAL):: capacity
     real(KREAL):: conductivity
     real(KREAL):: Pr
     real(KREAL):: Pe
     real(KREAL):: Nu
     
     De=4*flow_area/wetted_perimeter
     Re=4*density*velocity*De/viscosity
     Pr=viscosity*capacity/conductivity
     Pe=Re*Pr
     Nu=5.0+2.5D-2*Pe**0.8
    end subroutine get_Nusselt
    !
    subroutine get_hyconstant(rc,pd,Aflow,wet,de)
       real(KREAL):: rc,p,pd !r是包壳外半径 p是对边距
       real(KREAL):: Aflow,Ashell,Atotal,wet,de
       
           p=pd*2*rc
           Ashell=3.14*rc*rc
           Atotal=0.5*sqrt(3.0)*p*p
           Aflow=Atotal-Ashell
           wet=2*3.14*rc+p/sqrt(3.0)*6
           de=4*Aflow/wet
    endsubroutine get_hyconstant
end module imp_mathkerel
