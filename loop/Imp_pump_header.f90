module  Imp_pump_header
	use constants
	use imp_property
	implicit none
    
    type pump
		real(KREAL)::I
		!OperatingCondition
		real(KREAL)::yita
		real(KREAL)::omega
		real(KREAL)::Q
		real(KREAL)::omegae
		real(KREAL)::Qe
		real(KREAL)::He!额定扬程
		!material
		real(KREAL)::rho
		!thermal
		real(KREAL)::T
	  contains
		procedure,public::init=>init_pump
	end type pump
	private::init_pump
  contains
    subroutine init_pump(this)
		implicit none
		class(pump),intent(in out)::this
		!local
        real(KREAL)::temperature
		this%Q=this%Qe
		this%omega=this%omegae
		temperature=this%T
		this%rho=get_density(temperature)
	end subroutine init_pump
end module  Imp_pump_header