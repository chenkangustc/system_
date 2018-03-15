module imp_assm_header
    use imp_mathkerel
    use imp_property
	use constants
    implicit none
    type,public::AssmGeom
      
	  !real(KREAL):: rod            !元件半径
      real(KREAL):: pellet         !芯块半径
      real(KREAL):: Bond           !元件气隙厚度
      real(KREAL):: Cladth         !元件外壳厚度
      !real(KREAL):: AssmCladth    !组件外壳厚度
      real(KREAL):: pitch          !组件外对边距（包含包壳厚度）
      real(KREAL):: pd             !燃料元件PD比        
      real(KREAL):: rod            !元件半径
      real(KREAL):: area           !芯块横截面积
	  real(KREAL),allocatable:: Height(:)       !组件高度（活性区）
	  integer N_fuelpin            !燃料pin的个数
      integer N_pin                !pin总数  
    contains
      procedure,public::set=>set_assmgeom
	  !procedure,public::alloc=>alloc_assmgeom
	  !procedure,public::free=>free_assmgeom
      procedure,public::print=>print_assmgeom
    end type AssmGeom
    
    type,public::Assmmesh
        integer nf
        integer ng
        integer ns
		integer n_zone
        integer ny
        integer layer_bottom
        integer layer_top
        real(KREAL),allocatable::r(:,:)
        real(KREAL),allocatable::z(:,:)
      contains
      procedure,public::set=>set_assmmesh
      procedure,public::print=>print_assmmesh
      end type Assmmesh
    
    type,public::boundary
        real(KREAL):: inlet
        real(KREAL):: outlet
    end type boundary
    
    type,public::iteration
        real(KREAL),allocatable::Temperature(:,:) !pvt
        real(KREAL),allocatable::pressure(:)
        real(KREAL),allocatable::velocity(:)
    end type iteration
      
    type,public::th_boundary
       type(boundary)::p
       type(boundary)::u
       type(boundary)::T
       !type(boundary)::rho
       contains
       procedure,public::init=>init_th_boundary
       procedure,public::update=>update_th_boundary
       !procedure,public::init=>init_th_boundary !设置出口的边界条件
    end type th_boundary
    
    type,public::hydraulic
        real(KREAL):: fric
        real(KREAL):: aflow
        real(KREAL):: wet
        real(KREAL):: de
    contains
        procedure,public::set=>set_hydraulic
        procedure,public::print=>print_hydraulic
        procedure,public::cal=>cal_hydraulic
    end type hydraulic
    
        
    type,public::material!热物性和水力学参数
        real(KREAL),allocatable::rho(:,:)!热物性
        real(KREAL),allocatable::shc(:,:)
        real(KREAL),allocatable::ctc(:,:)
        real(KREAL),allocatable::dvs(:,:)
        real(KREAL),allocatable::htc(:)
     contains
       procedure,public::init=>init_material
    end type material
     
    type,public::thermal!                
        real(KREAL),allocatable::Temperature(:,:) !pvt
        real(KREAL),allocatable::pressure(:)
        real(KREAL),allocatable::velocity(:)
		real(KREAL),allocatable::Tcoolant(:)
		real(KREAL),allocatable::Tfuel(:)
		real(KREAL),allocatable::Tfuel_center(:)
		real(KREAL),allocatable::Tfg(:)
		real(KREAL),allocatable::Tgs(:)
		real(KREAL),allocatable::Tsc(:)
      contains
       procedure,public::init=>init_thermal
      end type thermal
      
    type,public::AssmInit
      real(KREAL):: Ti!初始温度
      real(KREAL):: Pi!初始压力
      real(KREAL):: Ui!初始速度
      real(KREAL):: Tin
      !real(KREAL):: Pin
      real(KREAL):: Pout
      real(KREAL):: Uin
    contains
      procedure,public::set=>set_assminit
      procedure,public::print=>print_assminit
    end type AssmInit
    
    type,public::confactor
      real(KREAL):: alpha!初始温度
      real(KREAL):: sigma!初始压力
    contains
      procedure,public::set=>set_confactor
      procedure,public::print=>print_confactor
    end type confactor
    
    type,public::assmpow
        real(KREAL),allocatable::power(:,:)
        real(KREAL),allocatable::fq_core(:,:)
    contains
        procedure,public::set=>set_power
        procedure,public::print=>print_power
    end type
        
    type,public::sys_time
        type(th_boundary)::th_boundary
        type(thermal)::thermal
        type(material)::material
        type(assmpow)::power
    end type sys_time
    
     private::set_assmgeom
     private::print_assmgeom
     private::set_assmmesh
     private::print_assmmesh
     private::set_assminit
     private::print_assminit
     private::set_confactor
     private::print_confactor
     private::set_hydraulic
     private::print_hydraulic
     private::set_power
     private::print_power
     private::init_th_boundary!会随时间变化的量用init
     private::update_th_boundary
     private::init_material
     private::init_thermal
     private::cal_hydraulic
     
     !private::cal_grid
    contains
    subroutine set_assmgeom(this,rod,Bond,Cladth,pitch,Height,pd,N_fuelpin,N_pin)
       implicit none
       class(assmgeom),intent(in out)::this
       real(KREAL),intent(in)::rod      !元件半径
       real(KREAL),intent(in)::Bond         !元件气隙厚度
       real(KREAL),intent(in)::Cladth     	!元件外壳厚度
       !real(KREAL)::,intent(in)::AssmCladth !组件外壳厚度
       real(KREAL),intent(in)::pitch     !组件外对边距（包含包壳厚度）
       real(KREAL),intent(in)::Height(:)         !组件高度（活性区）
       real(KREAL),intent(in)::pd
       integer,intent(in)::N_fuelpin       !元件的个数
	   integer,intent(in)::N_pin
       this%rod=rod
       this%Bond=Bond
       this%Cladth=Cladth
       !this%AssmCladth=AssmCladth
       this%pitch=pitch
       this%Height=Height*0.01
       !this%Height=Height
       this%pd=pd
       this%N_fuelpin=N_fuelpin
	   this%N_pin=N_pin
       call this%print()
     end subroutine set_assmgeom
     
     subroutine print_assmgeom(this)
      implicit none
      class(assmgeom),intent(in out)::this
      write(*,*)'set geom as below:'
      write(*,100) this%rod,this%Bond,this%Cladth,this%pitch,this%Height(1),this%pd,this%N_fuelpin
      100 format(1x,'rod=',F8.4,3x,'Bond=',F8.4,3x,'Cladth=',F8.4,3x,1x,'pitch=',F7.3,3x,'Height=',F7.3,3x,'pd=',F7.3,3x,'N_fuelpin=',I3/)
     end subroutine print_assmgeom
     
     subroutine print_assmmesh(this)
       implicit none
       class(assmmesh),intent(in out)::this
       write(*,*)'set mesh as below:'
       write(*,100) this%ny,this%layer_bottom,this%layer_top
       100 format(1x,'Ny=',I3,3x,'layer_bottom=',I3,3x,'layer_top=',I3/)
     end subroutine print_assmmesh
     
     subroutine print_hydraulic(this)
       implicit none
       class(hydraulic),intent(in out)::this
       write(*,*)'set fric as below:'
       write(*,100) this%fric
       100 format(1x,'fric=',F7.3/)
     end subroutine print_hydraulic
     
     subroutine print_assminit(this)
       implicit none
       class(assminit),intent(in out)::this
       write(*,*)'set assminit as below:'
       write(*,100) this%Pi,this%Ui,this%Ti,this%Pout,this%uin,this%Tin
       100 format(1x,'Pi=',F8.1,3x,'Ui=',F5.2,3x,'Ti=',F5.1,3x,'Pout=',F8.1,3x,'Uin=',F5.2,3x,'Tin=',F5.1/)             
       !write(*,100) this%Pi,this%Ui,this%Ti,this%Pin,this%uin,this%Tin
       !100 format(1x,'Pi=',F8.1,3x,'Ui=',F5.2,3x,'Ti=',F5.1,3x,'Pin=',F8.1,3x,'Uin=',F5.2,3x,'Tin=',F5.1/)
     end subroutine print_assminit
     
     subroutine print_confactor(this)
       implicit none
       class(confactor),intent(in out)::this
       write(*,*)'set confactor as below:'
       write(*,100) this%alpha,this%sigma
       100 format(1x,'alpha=',F7.4,3x,'sigma=',F7.4/)
     end subroutine print_confactor
     
     subroutine set_assmmesh(this,nf,ng,ns,n_zone,ny,layer_bottom,layer_top)
        implicit none
        class(assmmesh),intent(in out)::this
		integer,intent(in)::nf
		integer,intent(in)::ng
		integer,intent(in)::ns
		integer,intent(in)::n_zone
        integer,intent(in)::ny
        integer,intent(in)::layer_bottom
        integer,intent(in)::layer_top
		this%nf=nf
		this%ng=ng
		this%ns=ns
		this%n_zone=n_zone
        this%ny=ny
        this%layer_bottom=layer_bottom
        this%layer_top=layer_top
        call this%print()
     end subroutine set_assmmesh
     
     subroutine init_th_boundary(this,Tin,uin,pout)
       implicit none
       class(th_boundary),intent(in out)::this
       real(KREAL),intent(in)::Tin
       real(KREAL),intent(in)::uin
       real(KREAL),intent(in)::pout
       write(*,*)'init the PVT inlet boundary value'
       this%T%inlet=Tin
       this%u%inlet=uin
       this%p%outlet=pout
       !this%rho%inlet=get_density(Tin)
     end subroutine init_th_boundary
     
     subroutine update_th_boundary(this,ctime)
       implicit none
       class(th_boundary),intent(in out)::this
       real(KREAL),intent(in)::ctime
       if(this%u%inlet>0.6)then
        this%u%inlet=2.0*1.0/(1.0+(ctime/5.5))
       elseif(this%u%inlet<=0.6)then
        this%u%inlet=0.6
       endif
       print*,'update the uin,uin=',this%u%inlet
     end subroutine update_th_boundary
     !subroutine init_material(this,LBE,he,T91)  
     subroutine init_material(this,Nf,Ng,Ns,Ny)
       implicit none
       class(material),intent(in out)::this
       integer,intent(in)::Nf
       integer,intent(in)::Ng
       integer,intent(in)::Ns
       integer,intent(in)::Ny
       integer M,N!local
       integer i,j!local
       write(*,*)'init the material value'
       M=Ny+1
       N=Nf+Ng+Ns+1
          Do i=0,M,1
             Do j=0,N,1
                 if (j>=0.and.j<=Nf)then !芯块热物性 UO2
                  !RHOI(i,j)=10980
                  this%RHO(i,j)=10980
                  this%SHC(i,j)=300.0
                  this%CTC(i,j)=4.33
                  this%DVS(i,j)=0.0
                 elseif (j>Nf.and.j<=Nf+Ng) then!气隙热物性 He
                  !RHOI(i,j)=1.785
                  this%RHO(i,j)=1.785
                  this%SHC(i,j)=1.260
                  this%CTC(i,j)=0.124
                  this%DVS(i,j)=0.0
                 elseif (j>Nf+Ng.and.j<=Nf+Ng+Ns)then!包壳热物性 Ti
                  !RHOI(i,j)=7900
                  this%RHO(i,j)=7900
                  this%SHC(i,j)=502.42
                  this%CTC(i,j)=18.84
                  this%DVS(i,j)=0.0
                 else!流体物性
                  !RHOI(i,j)=10470  
                  this%RHO(i,j)=10470
                  this%SHC(i,j)=159
                  this%CTC(i,j)=3.61
                  this%DVS(i,j)=5.0 !动力粘度，而非运动粘度
                 endif                 
             enddo
         enddo    
     end subroutine init_material
     
     subroutine init_thermal(this,Temperature,Pressure,Velocity)
        implicit none
        class(thermal),intent(in out)::this
        real(KREAL),intent(in)::Temperature
        real(KREAL),intent(in)::Pressure
        real(KREAL),intent(in)::Velocity
		!local
		integer i,N
        write(*,*)'init the thermal value'
        this%Temperature=Temperature
        !this%Pressure=Pressure
        this%Velocity=Velocity
		N=size(this%Pressure)
		do i=1,N,1
		   !this%Pressure(i)=Pressure-5000.0*(i-1)
		   this%Pressure(i)=Pressure+500.0*(N-i)
		enddo
		this%Tfuel_center=0.0
		this%Tcoolant=0.0
		this%Tfuel=0.0
		this%Tfg=0.0
		this%Tgs=0.0
		this%Tsc=0.0
     end subroutine init_thermal
     
     subroutine set_AssmInit(this,Ti,Pi,Ui,Tin,Pout,Uin)
        implicit none
        class(AssmInit),intent(in out)::this
        real(KREAL),intent(in)::Ti
        real(KREAL),intent(in)::Pi
        real(KREAL),intent(in)::Ui
        real(KREAL),intent(in)::Tin
        real(KREAL),intent(in)::Pout
        real(KREAL),intent(in)::Uin
        this%Ti=Ti
        this%Pi=Pi
        this%Ui=Ui
        this%Tin=Tin
        this%Pout=Pout
        this%Uin=Uin
        call this%print()
     end subroutine set_AssmInit
     
     subroutine set_confactor(this,alpha,sigma)
        implicit none
        class(confactor),intent(in out)::this
        real(KREAL),intent(in)::alpha
        real(KREAL),intent(in)::sigma
        this%alpha=alpha
        this%sigma=sigma
        call this%print
     end subroutine set_confactor
 
     subroutine set_hydraulic(this,fric)
        implicit none
        class(hydraulic),intent(in out)::this
        real(KREAL),intent(in)::fric
        this%fric=fric
        call this%print
     end subroutine set_hydraulic
     
     subroutine cal_hydraulic(this,rc,pd)
        implicit none
        class(hydraulic),intent(in out)::this
        real(KREAL),intent(in)::rc,pd
        print*,'calculate the hydraulic constants...'
        call get_hyconstant(rc,pd,this%Aflow,this%wet,this%de)
     end subroutine cal_hydraulic
     
     subroutine set_power(this,power,fq_core)
        implicit none
        class(assmpow),intent(in out)::this
        real(KREAL),intent(in)::power(:,:)
        real(KREAL),intent(in)::fq_core(:,:)
        write(*,*)'set power...'
        this%power=power
        this%fq_core=fq_core
        !call this%print
     end subroutine set_power
     
    subroutine print_power(this)
        implicit none
        class(assmpow),intent(in out)::this
        write(*,*)'set power as below:'
        write(*,100) this%power
        100     format(1x,'power=',I3)
        
        write(*,*)'set fq_core as below:'
        write(*,101) this%fq_core
        101 format(1x,'fq_core=',I3)
     end subroutine print_power
     
end module imp_assm_header