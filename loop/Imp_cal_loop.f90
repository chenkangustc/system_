module Imp_cal_loop
    use Imp_loop_global
	use constants
	implicit none
	contains
	subroutine driving_loop_transient(current)
        implicit none
        real(KREAL),intent(in)::current
        
		call cal_loop_hydraulic(current)
	end subroutine driving_loop_transient
	
	subroutine cal_loop_hydraulic(current)
		implicit none
		real(KREAL),INTENT(in)::current!time
		!local
		real(KREAL)::alpha,beta,LAsum
		real(KREAL)::flowrate
		integer i!i is the num of the device
		LAsum=PipePR%Length/PipePR%Area+core%Length/core%Area+PipeRI%Length/PipeRI%Area+IHX1%Lsingle/IHX1%AreaOuter+PipeIP%Length/PipeIP%Area
		alpha=LAsum+pump1%rho*pump1%yita*pump1%I*pump1%omegae**2/pump1%Qe**2
		call cal_beta(beta,0)
		flowrate=alpha/(beta*current+alpha/pump1%Qe)
        call set_flowrate(flowrate)
	end subroutine cal_loop_hydraulic
	
	subroutine cal_beta(beta,formula)
		implicit none
		real(KREAL),intent(in out)::beta
		integer,intent(in)::formula!if formula==1 then K/fric else he/we
		select case(formula)
        case(0)!rho keep constant
			beta=pump1%rho*9.80*pump1%He/pump1%Qe**2
		case(1)!rho changes	
            beta=PipePR%beta+core%beta+PipeRI%beta+IHX1%betaOuter+PipeIP%beta
		end select
    end subroutine cal_beta
    
	subroutine test_loop_hydraulic()
		implicit none
		!local
		real(KREAL)::alpha,beta,LAsum
		real(KREAL)::flowrate,ffe
		integer i,current!i is the num of the device
		!LAsum=PipePR%Length/PipePR%Area+core%Length/core%Area+PipeRI%Length/PipeRI%Area+IHX1%Lsingle/IHX1%AreaOuter+PipeIP%Length/PipeIP%Area
		LAsum=959.0!m-1
		pump1%Qe=125.0
        pump1%He=41.9
		pump1%yita=0.756
		pump1%omegae=49*PI
		pump1%I=3.7
		pump1%rho=1000.0
		alpha=LAsum+pump1%rho*pump1%yita*pump1%I*pump1%omegae**2/pump1%Qe**2
		call cal_beta(beta,0)
		
		open(1,file='.\flowrate.txt')
		do i=0,25,1
			current=i
			flowrate=alpha/(beta*current+alpha/pump1%Qe)
            call set_flowrate(flowrate)
			print*, flowrate
            ffe=flowrate/pump1%Qe
			write(1,100) current,ffe
			100 Format(I6,' ',F10.2)
		enddo
		close(1)
    end subroutine test_loop_hydraulic
    
    subroutine set_flowrate(flowrate)
        implicit none
        real(KREAL),intent(in)::flowrate
            pump1%Q=flowrate
            PipePR%Q=flowrate
            core%Q=flowrate
            PipeRI%Q=flowrate
            IHX1%Qouter=flowrate
            PipeIP%Q=flowrate
    end subroutine set_flowrate
end module Imp_cal_loop
