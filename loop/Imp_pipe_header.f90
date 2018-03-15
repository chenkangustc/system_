module Imp_pipe_header
    use constants
    use imp_property
  
    type pipe
	    !geom 
        real(KREAL)::Length
        real(KREAL)::r
        real(KREAL)::theta
		real(KREAL)::area
		!mesh 
		integer::N
		!hydraulic
		real(KREAL)::De
        real(KREAL)::fric
        real(KREAL)::K
		real(KREAL)::Q
		real(KREAL)::beta
		!thermal
		real(KREAL)::T
		!material
		real(KREAL)::rho
    contains
        !procedure,public::set=>set_pipe
        procedure,public::init=>init_pipe
        ! procedure,public::alloc=>alloc_pipe
        ! procedure,public::free=>free_pipe    
    end type pipe
    !private::set_pipe
    private::init_pipe
    ! private::alloc_pipe
    ! private::free_pipe
  contains
      subroutine init_pipe(this)
      implicit none
      class(pipe),intent(in out)::this
      !local
      integer::i,N 
      real(KREAL)::Length,Aflow,wet,rpipe
      N=this%N
      Length=this%Length
	  Aflow=this%Q
	  rpipe=this%r
      do i=1,N,1
          ! this%geom%length(i)=Length/N
          this%rho=get_density(this%T)
      end do
	  !area
	  this%area=PI*rpipe*rpipe
      !de
	  wet=2*PI*rpipe
      this%de=4.0*Aflow/wet
	  !beta
	  this%beta=0.5*(this%fric*this%length/this%de+this%K)*1.0/(this%rho*this%area**2)
    end subroutine init_pipe
    ! subroutine alloc_pipe(this,N)
      ! implicit none
      ! class(pipe),intent(in out)::this
      ! !local
      ! integer::N
      ! N=this%mesh%N
      ! !integer,intent(in)::N
      ! !check allocated first
      ! call Free_pipe(this)
      ! allocate(this%geom%length(1:N))
      ! allocate(this%property%rho(1:N))
      ! allocate(this%thermal%T(1:N))
    ! end subroutine alloc_pipe
    
    ! subroutine free_pipe(this)
      ! implicit none
      ! class(pipe),intent(in out)::this
      ! if(allocated(this%geom%length)) deallocate(this%geom%length)
      ! if(allocated(this%property%rho)) deallocate(this%property%rho)
      ! if(allocated(this%thermal%T)) deallocate(this%thermal%T)
    ! end subroutine free_pipe
    

end module Imp_pipe_header