module imp_timer_header
    use constants
    implicit none
    private
    public::sys_timer
    
    type::sys_timer
        !private
        real(KREAL)::ltime
        real(KREAL)::ctime
        integer::Nt
        real(KREAL)::dt
        real(KREAL)::ttotal
        real(KREAL),allocatable::Tout(:)
        real(KREAL),allocatable::Pow(:)
        real(KREAL),allocatable::Uin(:)
    contains
        procedure,public::set=>set_timer
        procedure,public::init=>init_timer
        procedure,public::record=>record_timer
        procedure,public::alloc=>alloc_timer
        procedure,public::clc=>clc_timer
    end type sys_timer
    private::set_timer
    private::init_timer
    private::record_timer
    private::alloc_timer
    private::clc_timer
    contains
    subroutine set_timer(this,ttotal,Nt)
     class(sys_timer),intent(inout)::this
     real(KREAL),intent(in)::ttotal
     integer,intent(in)::Nt
     this%ttotal=ttotal
     this%Nt=Nt
     this%dt=this%ttotal/this%Nt
    end subroutine set_timer
    
    subroutine init_timer(this,ctime,ltime)
     class(sys_timer),intent(inout)::this
     real(KREAL),intent(in)::ctime
     real(KREAL),intent(in)::ltime    

     this%ctime=ctime
     this%ltime=ltime
     this%tout=0.0
     this%uin=0.0
     this%pow=0.0
    end subroutine init_timer
    
    subroutine record_timer(this,tout,uin,pow)
     class(sys_timer),intent(inout)::this
     real(KREAL),intent(in)::tout
     real(KREAL),intent(in)::uin
     real(KREAL),intent(in)::pow
     !local
     integer i
     i=this%ctime/this%dt
     this%Tout(i)=tout
     this%uin(i)=uin
     this%Pow(i)=pow
    end subroutine record_timer
    
    subroutine alloc_timer(this)
     class(sys_timer),intent(inout)::this
     integer Nt
     Nt=this%Nt
     call this%clc
     allocate(this%tout(Nt))
     allocate(this%uin(Nt))
     allocate(this%pow(Nt))
    end subroutine alloc_timer
    
    subroutine clc_timer(this)
     class(sys_timer),intent(inout)::this
     if(allocated(this%tout))  deallocate(this%tout)
     if(allocated(this%uin))  deallocate(this%uin)
     if(allocated(this%pow))  deallocate(this%pow)

    end subroutine clc_timer
    


end module imp_timer_header