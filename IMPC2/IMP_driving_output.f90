module imp_driving_post_process
    use imp_assm_global
	use constants
    implicit none
    private
    public::Run_output
	public::Free_imp_thermal
contains

	subroutine Free_imp_thermal()
		if(allocated(assm1))   deallocate(assm1)
		if(allocated(imp_pow)) deallocate(imp_pow)
	endsubroutine Free_imp_thermal
	 
    subroutine Run_output()
	   !local
	   integer i,j,k
	   integer zone_,Nf,Ng,Ns,Ny,Nradial
	   real(KREAL)::temperature,pressure,velocity,density,htc
	   real(KREAL)::coordinate_r,coordinate_z
	   real(KREAL)::Tfg,Tgs,Tsc,xf,xg,xs
	   !integer,allocatable::Nf_(:),Ng_(:),Ns_(:),Ny_(:),Nradial(:)
	   !real,allocatable::temperature(:,:,:),pressure(:),velocity(:)
	   !real,allocatable::pressure(:),velocity(:)
	   zone_=size(assm1)
	   ! allocate(Nf_(zone_),Ng_(zone_),Ns_(zone_),Ny_(zone_),Nradial(zone_))
	   ! do i=1,zone_,1
	    ! Nf_(i)=assm1(i)%mesh%nf
		! ng_(i)=assm1(i)%mesh%ng
		! ns_(i)=assm1(i)%mesh%ns
		! ny_(i)=assm1(i)%mesh%ny	
		! Nradial(i)=Nf_(i)+ng_(i)+ns_(i)+1
	   ! enddo
	   !allocate(temperature(zone_,ny_(i),Nradial(i)),pressure(ny_(i)),velocity(ny_(i)))
	   ! allocate(pressure(ny_(i)),velocity(ny_(i)))
        ! do i=1,zone_,1
		  ! do j=1,ny_(i),1
		    ! do k=1,Nradial(i),1
			 ! !temperature(i,j,k)=assm1(i)%thermal%temperature(j,k)
			 ! if (k==Nradial(i)) then
			    ! pressure(j)=assm1(i)%thermal%pressure(j)
				! if (j==1) then
				  ! velocity(1)=assm1(i)%thermal%velocity(1)
				! elseif(j>1.and.j<=Ny_(i)-1) then
				  ! velocity(j)=(assm1(i)%thermal%velocity(j-1)+assm1(i)%thermal%velocity(j))/2.0
				! elseif(j==Ny_(i)) then
				  ! velocity(j)=assm1(i)%thermal%velocity(j-1)
				! endif
			 ! endif			 
			! enddo
		  ! enddo
		! enddo
		
        ! open(1,file='.\output\mesh.txt')
        ! write(1,*) zone_,assm1(1)%mesh%Nf,assm1(1)%mesh%Ng,assm1(1)%mesh%Ns,assm1(1)%mesh%Ny
        ! close(1) 
        
        open(1,file='.\output\solid.txt')
        do i=1,zone_,1
		    Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns
            do j=0,Ny+1,1
			  do k=0,Nradial,1
			    if(k==0)then
				  if(j==0)then
                     temperature=assm1(i)%thermal%temperature(1,1)		
				  elseif(j==Ny+1)then
				     temperature=assm1(i)%thermal%temperature(Ny,1)
				  else
					 temperature=assm1(i)%thermal%temperature(j,1)
				  endif
				else
				  if(j==0)then
				     temperature=assm1(i)%thermal%temperature(1,k)
				  elseif(j==Ny+1)then
				     temperature=assm1(i)%thermal%temperature(Ny,k)
			      else
				     temperature=assm1(i)%thermal%temperature(j,k)
				  endif
				endif
				 coordinate_r=assm1(i)%mesh%r(j,k)
				 coordinate_z=assm1(i)%mesh%z(j,k)
                 write(1,100) i,j,k,temperature,coordinate_r,coordinate_z,Nf,Ng,Ns
				 100 Format(I6,' ',I6,' ',I6,' ',F10.2,' ',F10.4,' ',F10.4,' ',I6,' ',I6,' ',I6)
			 enddo
		 enddo
	    enddo
        close(1) 
		
		open(2,file='.\output\fluid.txt')
		do i=1,zone_,1
		    Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns+1
		  do j=0,Ny+1,1		   
		   if(j==0)then
			temperature=assm1(i)%th_boundary%T%inlet
			pressure=assm1(i)%th_boundary%p%inlet
			velocity=assm1(i)%th_boundary%u%inlet
		   elseif(j==Ny+1)then
		    temperature=assm1(i)%th_boundary%T%outlet
			pressure=assm1(i)%th_boundary%p%outlet
			velocity=assm1(i)%th_boundary%u%outlet
		   else
		    temperature=assm1(i)%thermal%temperature(j,Nradial)
			pressure=assm1(i)%thermal%pressure(j)
			if (j==1) then
			 velocity=assm1(i)%thermal%velocity(1)
			elseif(j>1.and.j<=Ny-1) then
			 velocity=(assm1(i)%thermal%velocity(j-1)+assm1(i)%thermal%velocity(j))/2.0
			elseif(j==Ny) then
			 velocity=assm1(i)%thermal%velocity(j-1)
			endif
		   endif
		   density=assm1(i)%property%rho(j,Nradial)
		   htc=assm1(i)%property%htc(j)
		   coordinate_z=assm1(i)%mesh%z(j,Nradial)
		   write(2,101) i,j,temperature,velocity,pressure,density,htc,coordinate_z
		   101 Format(I6,' ',I6,' ',F10.2,' ',F10.2,' ',F10.1,' ',F10.1,' ',F10.2,' ',F10.4)
		  enddo
		enddo
		close(2)
		
		open(3,file='.\output\mesh.txt')
		do i=1,zone_,1
		    Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns+1
			write(3,102) zone_,Ny,Nradial,Nf,Ng,Ns
			102 Format(I6,' ',I6,' ',I6,' ',I6,' ',I6,' ',I6)
		enddo
		close(3)
		
		open(4,file='.\output\Tsurface.txt')
		do i=1,zone_,1
			Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns+1
			xf=assm1(i)%geom%pellet
			xg=assm1(i)%geom%bond
			xs=assm1(i)%geom%cladth
		 do j=1,Ny,1
		    Tfg=(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf)*assm1(i)%thermal%temperature(j,Nf+1))/(assm1(i)%property%ctc(j,Nf)*(Xg/Ng)+assm1(i)%property%ctc(j,Nf+1)*(Xf/Nf))!芯块外边界
			Tgs=(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)*assm1(i)%thermal%temperature(j,Nf+Ng)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng)*assm1(i)%thermal%temperature(j,Nf+Ng+1))/(assm1(i)%property%ctc(j,Nf+Ng)*(Xs/Ns)+assm1(i)%property%ctc(j,Nf+Ng+1)*(Xg/Ng))!包壳内边界
			Tsc=(assm1(i)%property%htc(j)*assm1(i)%thermal%temperature(j,Nradial)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns)*assm1(i)%thermal%temperature(j,Nradial-1))/(assm1(i)%property%htc(j)+2*assm1(i)%property%ctc(j,Nradial-1)/(Xs/Ns))!包壳外边界
			coordinate_z=assm1(i)%mesh%z(j,Nradial)
			write(4,103) i,j,Tfg,Tgs,Tsc,coordinate_z
			103 Format(I6,' ',I6,' ',F10.2,' ',F10.2,' ',F10.2,' ',F10.4)
		 enddo
		enddo
		close(4)
		
		open(5,file='.\output\power.txt')
		do i=1,zone_,1
			Nf=assm1(i)%mesh%Nf
			Ng=assm1(i)%mesh%Ng
			Ns=assm1(i)%mesh%Ns
			Ny=assm1(i)%mesh%Ny
			Nradial=Nf+Ng+Ns+1
			do j=1,Ny,1
				coordinate_z=assm1(i)%mesh%z(j,Nradial)
				write(5,104)i,j,imp_pow(i,j),coordinate_z
				104 Format(I6,' ',I6,' ',F15.4,' ',F10.4)
			enddo
		enddo
		close(5)
	! do i=0,M,1
        ! Tfg(i)=(CTC(i,Nf)*(Xg/Ng)*T(i,Nf)+CTC(i,Nf+1)*(Xf/Nf)*T(i,Nf+1))/(CTC(i,Nf)*(Xg/Ng)+CTC(i,Nf+1)*(Xf/Nf))!芯块外边界
        ! Tgs(i)=(CTC(i,Nf+Ng)*(Xs/Ns)*T(i,Nf+Ng)+CTC(i,Nf+Ng+1)*(Xg/Ng)*T(i,Nf+Ng+1))/(CTC(i,Nf+Ng)*(Xs/Ns)+CTC(i,Nf+Ng+1)*(Xg/Ng))!包壳内边界
        ! Tsf(i)=(htc(i)*T(i,N)+2*CTC(i,N-1)/(Xs/Ns)*T(i,N-1))/(htc(i)+2*CTC(i,N-1)/(Xs/Ns))!包壳外边界
    ! enddo
					  ! if(k==Nradial) then
			    ! write(2,*) i,j,k,nf_(i),ng_(i),ns_(i),assm1(i)%thermal%temperature(j,k),velocity(j),pressure(j)
			  ! else 
			    ! write(2,*) i,j,k,nf_(i),ng_(i),ns_(i),assm1(i)%thermal%temperature(j,k),0.0,0.0
			  ! endif
		
		! if(allocated(Nf_))  deallocate(Nf_)
		! if(allocated(Ng_))  deallocate(Ng_)
		! if(allocated(Ns_))  deallocate(Ns_)
		! if(allocated(Ny_))  deallocate(Ny_)
		! if(allocated(Nradial))  deallocate(Nradial)
		! !if(allocated(temperature))  deallocate(temperature)
		! if(allocated(pressure))  deallocate(pressure)
		! if(allocated(velocity))  deallocate(velocity)

	  ! open(3,file='.\output\pressure2.txt')
      ! write(3,*) assm1(1)%thermal%pressure
 
      ! close(3) 
	  ! open(4,file='.\output\velocity2.txt')
      ! write(4,*) assm1(1)%thermal%velocity
  
      ! close(4) 
	  ! open(5,file='.\output\pow.txt')
      ! write(5,*) assm1(1)%pow%power

      ! close(5)       
      !open(3,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tout.txt')
      !write(3,*) timer1%tout
      !close(3) 
      !open(4,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tpow.txt')
      !write(4,*) timer1%pow
      !close(4) 
      !open(5,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tuin.txt')
      !write(5,*) timer1%uin
      !close(5)

      !print*,assm1%th_boundary%p%outlet
      !liquid PVT、rho distribution
      !solid temperature distribution
      
	  ! print*,assm1(1)%th_boundary%u%inlet,assm1(1)%Thermal%velocity,assm1(1)%th_boundary%u%outlet
      ! print*,assm1(1)%th_boundary%p%inlet,assm1(1)%Thermal%pressure,assm1(1)%th_boundary%p%outlet
      ! print*,assm1(1)%th_boundary%T%inlet,assm1(1)%Thermal%temperature(:,assm1(1)%mesh%Nf+assm1(1)%mesh%Ng+assm1(1)%mesh%Ns+1),assm1(1)%th_boundary%T%outlet
      ! print*,assm1(1)%geom%height
	end subroutine Run_output
end module imp_driving_post_process
