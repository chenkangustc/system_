!
!value
!type::IHX
!
!type subroutine
! 			private::init_IHX
! 			private::alloc_IHX
! 			private::free_IHX
! 			private::cal_thermal_steady
! 			private::cal_thermal_transient
!*******************************************************************************
module Imp_IHX_header
	use constants
	use Imp_property
	implicit none
	
	type IHX
		!geom
		real(KREAL),allocatable::Length(:)
		real(KREAL)::Lsingle
		real(KREAL)::AreaOuter
		integer::Ntube
		real(KREAL)::Rtube
		real(KREAL)::AreaTubeSingle
		real(KREAL)::AreaTubeTotal
		!mesh
		integer::N
		!hydraulic
		real(KREAL)::DeOuter
		real(KREAL)::FricOuter
		real(KREAL)::Kouter
		real(KREAL)::WetOuter
		real(KREAL)::Qouter
		real(KREAL)::betaOuter
		!material
		real(KREAL),allocatable::rho(:)
		!thermal
		real(KREAL),allocatable::T(:)
		!initial
		real(KREAL)::Ti

	  contains
		procedure,public::init=>init_IHX
        procedure,public::alloc=>alloc_IHX
        procedure,public::free=>free_IHX
		procedure,public::thcals=>cal_thermal_steady
		procedure,public::thcalt=>cal_thermal_transient
	end type IHX
	
	private::init_IHX
    private::alloc_IHX
    private::free_IHX
	private::cal_thermal_steady
	private::cal_thermal_transient
  contains
	subroutine init_IHX(this)
		implicit none
		class(IHX),intent(in out)::this
		!local
		integer::i,N,Ntube
		real(KREAL)::Lsingle,rtube,wetOuter,Qouter,Kiouter
		real(KREAL),allocatable::temperature(:)
		
		Lsingle=this%Lsingle
		Qouter=this%Qouter
		Kiouter=this%Kouter/this%N
		rtube=this%Rtube
		Ntube=this%Ntube
		N=size(this%T)
		allocate(temperature(1:N))
		
		this%T(:)=this%Ti
		this%AreaTubeSingle=PI*rtube*rtube
		this%AreaTubeTotal=Ntube*this%AreaTubeSingle
	    this%DeOuter=4*Qouter/WetOuter	
		this%betaOuter=0.0	
		do i=1,N,1
			temperature(i)=this%T(i)
			this%rho(i)=get_density(temperature(i))
			this%Length(i)=Lsingle/N
			!betaOuter
			this%betaOuter=this%betaOuter+0.5*(this%fricouter*this%length(i)/this%DeOuter+Kiouter)*1/(this%rho(i)*this%areaouter**2)
		end do
		if(allocated(temperature)) deallocate(temperature)
	end subroutine init_IHX
	
	subroutine alloc_IHX(this)
      implicit none
      class(IHX),intent(in out)::this
      !local
      integer::N
      N=this%N
      !integer,intent(in)::N
      !check allocated first
      call Free_IHX(this)
      allocate(this%length(1:N))
      allocate(this%rho(1:N))
      allocate(this%T(1:N))
    end subroutine alloc_IHX
    
    subroutine free_IHX(this)
      implicit none
      class(IHX),intent(in out)::this
      if(allocated(this%length)) deallocate(this%length)
      if(allocated(this%rho)) deallocate(this%rho)
      if(allocated(this%T)) deallocate(this%T)
    end subroutine free_IHX
	
	subroutine cal_thermal_steady(this)
		implicit none
		class(IHX),intent(in out)::this
		
	end subroutine cal_thermal_steady
	
	subroutine cal_thermal_transient(this)
		implicit none
		class(IHX),intent(in out)::this
		
	end subroutine cal_thermal_transient
end module Imp_IHX_header