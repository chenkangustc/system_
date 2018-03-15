module testNK2TH
	use global_state
    use imp_assm_global
	use TH2NK_interface_imp,        only : Perform_TH_imp
	implicit none
    contains
    subroutine driving_testNK2TH()

        real(KREAL)  :: power(ns%state%zone, ns%state%layer)
        real(KREAL)  :: fq(ns%state%zone, ns%state%layer)
        
        logical  :: transient_flag = .FALSE.
        real(KREAL)  :: Tfuel(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Tcoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: Rhocoolant(ns%state%zone, ns%state%layer) 
        real(KREAL)  :: toutlet
        real(KREAL)  :: max_Tfuel 
        real(KREAL)  :: max_Tcoolant 
        real(KREAL)  :: min_Rhocoolant 
        real(KREAL)  :: last 
        real(KREAL)  :: current 
		!local
		integer Nradial,i_zone,nTime,i
        real(KREAL) ::tTotal,dtime
        
        Tfuel = 0.0; Tcoolant = 0.0; Rhocoolant = 0.0; 
        max_Tfuel = 0.0; max_Tcoolant = 0.0; min_Rhocoolant = 0.0; 
        last = 0.0; current = 0.0;
	    
	    transient_flag=.FALSE.
	   open(unit=1,file='.\output\powDistribution.txt')
       !read(1,*) xf,xg,xs,height,nf,ng,ns,ny,f,Tin,pout,Tic,uic,tmax,nt,sigma,sigmab,alpha
       !read(1,*) this%xf,this%xg,this%xs,this%acf,this%height,this%pd,this%npin,this%nf,this%ng,this%ns,this%f,this%Tin,this%pout,this%uin,this%pin,this%Ti,this%ui,this%pi,this%alpha,this%sigma
       read(1,100) power
	   100 Format(F15.5)
	   close(1)
	   
	   !print*,'power=',power
	   if(transient_flag==.FALSE.) then
            if (ns%feedback%is_inner)  then
	          call Perform_TH_imp(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
			endif  
	   else       
			  transient_flag=.FALSE.
			  call Perform_TH_imp(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  	   
			  power=0.8*power
			  transient_flag=.TRUE.
       		  tTotal=1.0
			  nTime=10
			  dtime=tTotal/nTime
			  do i=1,nTime,1
				  current=current+dtime
				  call Perform_TH_imp(transient_flag, power, Tfuel, Tcoolant, Rhocoolant, max_Tfuel, max_Tcoolant, min_Rhocoolant, last, current, toutlet)  
				  print*,'time=',current,'maxTfuel=',max_Tfuel,'maxTcoolant=',max_Tcoolant
				  last=current
			  enddo
	   endif
	      
       i_zone=19
	   Nradial=assm1(i_zone)%mesh%Nf+assm1(i_zone)%mesh%Ng+assm1(i_zone)%mesh%Ns+1
	   print*,'maxTfuel=',max_Tfuel,'maxTcoolant=',max_Tcoolant
	   print*,'Tcoolant=',Tfuel(i_zone,:)
	   print*,'zone=',i_zone,'Tinlet=',assm1(i_zone)%th_boundary%T%inlet
	   print*,'zone=',i_zone,'Temperature=',assm1(i_zone)%thermal%temperature(:,Nradial)
       print*,'zone=',i_zone,'uinlet=',assm1(i_zone)%th_boundary%u%inlet
	   print*,'zone=',i_zone,'velocity=',assm1(i_zone)%thermal%velocity
	   print*,'zone=',i_zone,'pressure=',assm1(i_zone)%thermal%pressure
       read(*,*)
    endsubroutine driving_testNK2TH
endmodule testNK2TH