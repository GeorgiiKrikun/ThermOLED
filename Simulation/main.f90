program main
	
	use constants
	implicit none
	
	
	include 'mpif.h'
	! Control over iterations
	integer :: number_of_main_iterations = 10000, number_of_phi_iterations=1,number_of_sg_iterations=100,&
	smalliterationinbig=10000
	
	integer :: i,j,k,cut, lbu=249,Nx,Ny,Nz,Ncutx,Ncuty,Ncutz,Ntx,Nty,Ntz,status(MPI_STATUS_SIZE)&
		,ierr,rc,n_tasks,len,compositioniterator,superiterator
	
	
		
	integer, dimension(:,:), allocatable :: cutproperties
	integer :: backup=0
	real :: starttime,endtime
	
	double precision :: lengthx=100d-9,lengthy=100d-9,lengthz=100d-9,V = 6.0d0,errorphi=1.0d0,errorSG=1.0d0,&
	maxheating,errorT=1.0d0,totalheating,totaloutflow,timestep = 10.0d-13,ptimestep = 1.0d-20,nnewmax,receive,&
	iterationTmax=0.0
	logical :: heatstep = .true.
	double precision, dimension(:,:,:), allocatable :: heatfluxx,heatfluxy,heatfluxz,&
		phi,n, p,EFx,EFy,EFz,ephi,hphi,T,currentx,currenty,currentz,ntempcurrentx,ntempcurrenty,&
		ntempcurrentz,ptempcurrentx,ptempcurrenty,ptempcurrentz,heating, changescurrent,changesT&
		,changesphi,changesn,changesp,changespcurrent,changesncurrent,recombination,relativerecombination
	double precision, dimension(:), allocatable :: maxerrorphi,maxerrorSG,maxerrorT,maxT
	
	double precision, dimension(:,:,:), allocatable :: phi1,n1, p1,EFx1,EFy1,EFz1,T1,currentx1,currenty1,currentz1,heating1
	
	CHARACTER(len=32) :: arg
	character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
        integer :: skipheatiteration = 0
	
	call MPI_INIT(ierr)
	if(ierr.ne.MPI_SUCCESS) then
		write(*,'(a)') 'ERROR in start: MPI initialisation failed'
		call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
		stop
	end if
	
	call MPI_COMM_RANK(MPI_COMM_WORLD,task_ID,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,n_tasks,ierr)
	call MPI_GET_PROCESSOR_NAME(hostname,len,ierr)

	if (task_ID == 0) then
		call CPU_TIME(starttime)
	end if
	
	call obtainmaterialcomposition
	if (task_ID == 0) then
		call GETARG(1,arg)
		read(arg,*) Nx
		call GETARG(2,arg)
		read(arg,*) Ny
		call GETARG(3,arg)
		read(arg,*) Nz
		call GETARG(4,arg)
		read(arg,*) Ncutx
		call GETARG(5,arg)
		read(arg,*) Ncuty
		call GETARG(6,arg)
		read(arg,*) Ncutz
		call GETARG(7,arg)
		read(arg,*) V

		call writecalculationinformation(Nx,Ny,Nz)

		allocate(heatfluxx(Nx,Ny,Nz),heatfluxy(Nx,Ny,Nz),heatfluxz(Nx,Ny,Nz),phi(Nx,Ny,Nz),n(Nx,Ny,Nz),&
		p(Nx,Ny,Nz),EFx(Nx,Ny,Nz),EFy(Nx,Ny,Nz),EFz(Nx,Ny,Nz),ephi(Nx,Ny,Nz),hphi(Nx,Ny,Nz),T(Nx,Ny,Nz),&
		currentx(Nx,Ny,Nz),currenty(Nx,Ny,Nz),currentz(Nx,Ny,Nz),ntempcurrentx(Nx,Ny,Nz),ntempcurrenty(Nx,Ny,Nz),&
		ntempcurrentz(Nx,Ny,Nz),ptempcurrentx(Nx,Ny,Nz),ptempcurrenty(Nx,Ny,Nz),ptempcurrentz(Nx,Ny,Nz),&
		heating(Nx,Ny,Nz),changescurrent(Nx,Ny,Nz),changesT(Nx,Ny,Nz),changesphi(Nx,Ny,Nz),changesn(Nx,Ny,Nz),&
		changesp(Nx,Ny,Nz),changespcurrent(Nx,Ny,Nz),changesncurrent(Nx,Ny,Nz),recombination(Nx,Ny,Nz), &
		relativerecombination(Nx,Ny,Nz))
		
		allocate(phi1(Nx,Ny,Nz),n1(Nx,Ny,Nz),p1(Nx,Ny,Nz),EFx1(Nx,Ny,Nz),EFy1(Nx,Ny,Nz),EFz1(Nx,Ny,Nz),T1(Nx,Ny,Nz),&
		currentx1(Nx,Ny,Nz),currenty1(Nx,Ny,Nz),currentz1(Nx,Ny,Nz), heating1(Nx,Ny,Nz))
		
		allocate(maxerrorphi(0:n_tasks-1),maxerrorSG(0:n_tasks-1),maxerrorT(0:n_tasks-1),maxT(0:n_tasks-1))
! 		allocate(cutproperties(n_tasks-1,6))
				
		if (backup>0) then
			call quantitiesbackup(phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T,backup)
		else
			do j = 1,Ny
				do i =1,Nx				
					!phi(i,j,1) = V*(1-2*0.25*( dsqrt((REAL(i,8)-(Nx+1)/2d0)**2 + (REAL(j,8)-(Ny+1)/2d0)**2)/((Nx+1)*(Ny+1)/4.0d0)))
! 					phi(:,j,1) = V*(1-2*0.013*dabs( (Real(Ny,8)+1d0)/2d0-Real(j,8) )/dabs((Real(Ny,8)-1d0)/2d0))
					phi(:,j,1)=V
				end do
			end do
			!print*, MINVAL(phi(1,:,1))
			!call SLEEP(10)
			phi(:,:,Nz)=0
			n=1d21
			p=1d21
			T=300.0d0
		end if

		
		allocate(subgridx(0:Nx+1), subgridy(0:Ny+1), subgridz(Nz))
		do i=1,Nx
			subgridx(i) = lengthx*(i-1)/(Nx-1)
		end do
		subgridx(0) = subgridx(1) - subgridx(Nx)+subgridx(Nx-1)
		subgridx(Nx+1) = subgridx(Nx) + subgridx(2)-subgridx(1)
		
		do j=1,Ny
			subgridy(j) = lengthy*(j-1)/(Ny-1)
		end do
		
		subgridy(0) = subgridy(1) - subgridy(Ny) + subgridy(Ny-1)
		subgridy(Ny+1) = subgridy(Ny) + subgridy(2) - subgridy(1)

		do k=1,Nz
			subgridz(k) = lengthz*(k-1)/(Nz-1)
		end do
		
		print*, subgridx
		print*, subgridy
		print*, subgridz
		
		allocate(epsR(0:Nx+1,0:Ny+1,Nz),rho(0:Nx+1,0:Ny+1,Nz),cp(0:Nx+1,0:Ny+1,Nz),kappa(0:Nx+1,0:Ny+1,Nz),HOMO(0:Nx+1,0:Ny+1,Nz)&
		,LUMO(0:Nx+1,0:Ny+1,Nz),psigma(0:Nx+1,0:Ny+1,Nz),nsigma(0:Nx+1,0:Ny+1,Nz),materials(0:Nx+1,0:Ny+1,Nz))
		do k = 1,Nz
			do compositioniterator = 1,composition%matnum
				if ((subgridz(k) >= composition%matstart(compositioniterator)) .and. &
				(subgridz(k) < composition%matend(compositioniterator))) then
					epsR(:,:,k) = composition%epsR(compositioniterator)
					rho(:,:,k) = composition%rho(compositioniterator)
					cp(:,:,k) = composition%cp(compositioniterator)
					kappa(:,:,k) = composition%kappa(compositioniterator)
					HOMO(:,:,k) = composition%penergy(compositioniterator)
					LUMO(:,:,k) = composition%nenergy(compositioniterator)
					psigma(:,:,k) = composition%psigma(compositioniterator)
					nsigma(:,:,k) = composition%nsigma(compositioniterator)
					materials(:,:,k) = compositioniterator
				end if
			end do
		end do
				
! 		materials = composition%matnum
! 		DETERMINATION OF CUT PROPERTIES
		
		allocate(cutproperties(n_tasks-1,6))
		call cuteverything(Nx,Ny,Nz,Ncutx,Ncuty,Ncutz,cutproperties)
		
		do i = 1,n_tasks-1
			print*, cutproperties(i,:)
		end do
		
		call MPI_BCAST(Nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(Ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(Nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		
		do i = 1, Ncutx+1
			do j = 1, Ncuty+1
				do k = 1, Ncutz+1
					cut = k+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					call MPI_SEND(cutproperties(cut,:),6,MPI_INTEGER,cut,1,MPI_COMM_WORLD,ierr)
				end do
			end do
		end do
		
		allocate(threadpoints(3,3,3))
		
! 		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		do i = 1, Ncutx+1
			do j = 1, Ncuty+1
				do k = 1, Ncutz+1
					cut = k+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					threadpoints = 0
					if (i .ne. 1) then
						threadpoints(1,2,2) = k+(j-1)*(Ncutz+1)+(i-2)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(1,2,2) = k+(j-1)*(Ncutz+1)+(Ncutx)*(Ncutz+1)*(Ncuty+1)
					end if
					if (j .ne. 1) then
						threadpoints(2,1,2) = k+(j-2)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(2,1,2) = k+(Ncuty)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					end if
					if (k .ne. 1) then
						threadpoints(2,2,1) = k-1+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(2,2,1) = -1
					end if
					if (i .ne. (Ncutx+1)) then
						threadpoints(3,2,2) = k+(j-1)*(Ncutz+1)+(i)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(3,2,2) = k+(j-1)*(Ncutz+1)
					end if
					if (j .ne. (Ncuty+1)) then
						threadpoints(2,3,2) = k+(j)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(2,3,2) = k+(i-1)*(Ncutz+1)*(Ncuty+1)
					end if
					if (k .ne. (Ncutz+1)) then
						threadpoints(2,2,3) = k+1+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					else
						threadpoints(2,2,3) = -2
					end if
					
					Ntx = cutproperties(cut,4)-cutproperties(cut,1)+1
					Nty = cutproperties(cut,5)-cutproperties(cut,2)+1
					Ntz = cutproperties(cut,6)-cutproperties(cut,3)+1
					
					call MPI_SEND3D_INT(threadpoints,cut,2)
					
					call MPI_SEND(subgridx(cutproperties(cut,1):cutproperties(cut,4)),&
					Ntx,MPI_DOUBLE,cut,3,MPI_COMM_WORLD,ierr)
					call MPI_SEND(subgridy(cutproperties(cut,2):cutproperties(cut,5)),&
					Nty,MPI_DOUBLE,cut,4,MPI_COMM_WORLD,ierr)
					call MPI_SEND(subgridz(cutproperties(cut,3):cutproperties(cut,6)),&
					Ntz,MPI_DOUBLE,cut,5,MPI_COMM_WORLD,ierr)
					call MPI_SEND3D_DP(epsR(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,6)
					call MPI_SEND3D_INT(materials(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,7)
					call MPI_SEND3D_DP(psigma(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,8)
					call MPI_SEND3D_DP(nsigma(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,9)
					call MPI_SEND3D_DP(rho(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,10)
					call MPI_SEND3D_DP(cp(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,11)
					call MPI_SEND3D_DP(kappa(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,12)
					call MPI_SEND3D_DP(HOMO(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,13)
					call MPI_SEND3D_DP(LUMO(cutproperties(cut,1):cutproperties(cut,4),&
					cutproperties(cut,2):cutproperties(cut,5),cutproperties(cut,3):cutproperties(cut,6)),&
					cut,14)
					call MPI_SEND3D_DP(phi(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,15)
					call MPI_SEND3D_DP(EFx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,16)
					call MPI_SEND3D_DP(EFy(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,17)
					call MPI_SEND3D_DP(EFz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,18)
					call MPI_SEND3D_DP(n(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,19)
					call MPI_SEND3D_DP(p(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,20)
					call MPI_SEND3D_DP(currentx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,21)
					call MPI_SEND3D_DP(currenty(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,22)
					call MPI_SEND3D_DP(currentz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,23)
					call MPI_SEND3D_DP(heating(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,24)
					call MPI_SEND3D_DP(T(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,25)
					call MPI_SEND3D_DP(relativerecombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,26)
					call MPI_SEND3D_DP(recombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
					cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)),&
					cut,27)
					
				end do
			end do
		end do
	end if
	
	if (task_ID .ne. 0) then
		allocate(cutproperties(1,6))
		call MPI_RECV(cutproperties(1,:),6,MPI_INTEGER,0,1,MPI_COMM_WORLD,status,ierr)
		Print*, task_ID, cutproperties, Ntx,Nty,Ntz
		Ntx = cutproperties(1,4)-cutproperties(1,1)+1
		Nty = cutproperties(1,5)-cutproperties(1,2)+1
		Ntz = cutproperties(1,6)-cutproperties(1,3)+1
		allocate(heatfluxx(Ntx,Nty,Ntz),heatfluxy(Ntx,Nty,Ntz),heatfluxz(Ntx,Nty,Ntz),phi(Ntx,Nty,Ntz),n(Ntx,Nty,Ntz),&
		p(Ntx,Nty,Ntz),EFx(Ntx,Nty,Ntz),EFy(Ntx,Nty,Ntz),EFz(Ntx,Nty,Ntz),ephi(Ntx,Nty,Ntz),hphi(Ntx,Nty,Ntz),T(Ntx,Nty,Ntz),&
		currentx(Ntx,Nty,Ntz),currenty(Ntx,Nty,Ntz),currentz(Ntx,Nty,Ntz),ntempcurrentx(Ntx,Nty,Ntz),ntempcurrenty(Ntx,Nty,Ntz),&
		ntempcurrentz(Ntx,Nty,Ntz),ptempcurrentx(Ntx,Nty,Ntz),ptempcurrenty(Ntx,Nty,Ntz),ptempcurrentz(Ntx,Nty,Ntz),&
		heating(Ntx,Nty,Ntz),changescurrent(Ntx,Nty,Ntz),changesT(Ntx,Nty,Ntz),changesphi(Ntx,Nty,Ntz),changesn(Ntx,Nty,Ntz),&
		changesp(Ntx,Nty,Ntz),changespcurrent(Ntx,Nty,Ntz),changesncurrent(Ntx,Nty,Ntz),relativerecombination(Ntx,Nty,Ntz),&
		recombination(Ntx,Nty,Ntz))

		allocate(threadpoints(3,3,3))
		call MPI_RECV3D_INT(threadpoints,0,2)

		allocate(subgridx(Ntx),subgridy(Nty),subgridz(Ntz))
		call MPI_RECV(subgridx,Ntx,MPI_DOUBLE,0,3,MPI_COMM_WORLD,status,ierr)
		call MPI_RECV(subgridy,Nty,MPI_DOUBLE,0,4,MPI_COMM_WORLD,status,ierr)
		call MPI_RECV(subgridz,Ntz,MPI_DOUBLE,0,5,MPI_COMM_WORLD,status,ierr)

		allocate(epsR(Ntx,Nty,Ntz),materials(Ntx,Nty,Ntz),psigma(Ntx,Nty,Ntz),nsigma(Ntx,Nty,Ntz),&
		rho(Ntx,Nty,Ntz),cp(Ntx,Nty,Ntz),kappa(Ntx,Nty,Ntz),HOMO(Ntx,Nty,Ntz),LUMO(Ntx,Nty,Ntz))
		
		call MPI_RECV3D_DP(epsR,0,6)
		call MPI_RECV3D_INT(materials,0,7)
		call MPI_RECV3D_DP(psigma,0,8)
		call MPI_RECV3D_DP(nsigma,0,9)
		call MPI_RECV3D_DP(rho,0,10)
		call MPI_RECV3D_DP(cp,0,11)
		call MPI_RECV3D_DP(kappa,0,12)
		call MPI_RECV3D_DP(HOMO,0,13)
		call MPI_RECV3D_DP(LUMO,0,14)
		call MPI_RECV3D_DP(phi(2:Ntx-1,2:Nty-1,:),0,15)
		call MPI_RECV3D_DP(EFx(2:Ntx-1,2:Nty-1,:),0,16)
		call MPI_RECV3D_DP(EFy(2:Ntx-1,2:Nty-1,:),0,17)
		call MPI_RECV3D_DP(EFz(2:Ntx-1,2:Nty-1,:),0,18)
		call MPI_RECV3D_DP(n(2:Ntx-1,2:Nty-1,:),0,19)
		call MPI_RECV3D_DP(p(2:Ntx-1,2:Nty-1,:),0,20)
		call MPI_RECV3D_DP(currentx(2:Ntx-1,2:Nty-1,:),0,21)
		call MPI_RECV3D_DP(currenty(2:Ntx-1,2:Nty-1,:),0,22)
		call MPI_RECV3D_DP(currentz(2:Ntx-1,2:Nty-1,:),0,23)
		call MPI_RECV3D_DP(heating(2:Ntx-1,2:Nty-1,:),0,24)
		call MPI_RECV3D_DP(T(2:Ntx-1,2:Nty-1,:),0,25)
		call MPI_RECV3D_DP(relativerecombination(2:Ntx-1,2:Nty-1,:),0,26)
		call MPI_RECV3D_DP(recombination(2:Ntx-1,2:Nty-1,:),0,27)
		
		
		
		call MPI_TRIPLE_ROTATE(phi,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(EFx,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(EFy,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(EFz,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(n,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(p,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(currentx,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(currenty,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(currentz,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(heating,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(T,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(relativerecombination,threadpoints,Ntx,Nty,Ntz)
		call MPI_TRIPLE_ROTATE(recombination,threadpoints,Ntx,Nty,Ntz)
		
	end if
		
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 	!ITERATION STARTS
! 	if (task_ID .ne. 0) then
! 		if (backup == 0) then
! 			do i=1,10000
! 				call PoisonIteration(Ntx,Nty,Ntz,phi,n,p,ptimestep,errorPhi,V)
! 				call MPI_TRIPLE_ROTATE(phi,threadpoints,Ntx,Nty,Ntz)
! 				print*, "phi", i, errorPhi
! 			end do
! 						
! 			do i=1,1000000
! 				call PoisonIteration(Ntx,Nty,Ntz,phi,n,p,ptimestep,errorPhi,V)
! 				call MPI_TRIPLE_ROTATE(phi,threadpoints,Ntx,Nty,Ntz)
! 
! 				call getElectricField(Ntx,Nty,Ntz,phi,EFx,EFy,EFz)
! 				call MPI_TRIPLE_ROTATE(EFx,threadpoints,Ntx,Nty,Ntz)
! 				call MPI_TRIPLE_ROTATE(EFy,threadpoints,Ntx,Nty,Ntz)
! 				call MPI_TRIPLE_ROTATE(EFz,threadpoints,Ntx,Nty,Ntz)
! 
! 				ephi = phi + LUMO
! 				hphi = phi + HOMO
! 
! 				call sgiteration(Ntx,Nty,Ntz,n,p,phi,phi+LUMO,phi+HOMO,T,currentx,currenty,currentz,timestep,errorSG,&
! 	ntempcurrentx,ntempcurrenty,ntempcurrentz,ptempcurrentx,ptempcurrenty,ptempcurrentz,EFx,EFy,EFz,relativerecombination,&
! 	recombination)
! 				call MPI_TRIPLE_ROTATE(n,threadpoints,Ntx,Nty,Ntz)
! 				call MPI_TRIPLE_ROTATE(p,threadpoints,Ntx,Nty,Ntz)
! 				call MPI_TRIPLE_ROTATE(relativerecombination,threadpoints,Ntx,Nty,Ntz)
! 				call MPI_TRIPLE_ROTATE(recombination,threadpoints,Ntx,Nty,Ntz)
! 				if (mod(i,1000) == 0) then
!                                     print*, "Scharfetter-Gummel", i/1000, errorSG
!                                 end if
! 			end do
! 		end if
! 	end if
! 	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	call SLEEP(2)
	if (task_ID .ne. 0) then
		do i=1,number_of_main_iterations
			call MPI_SEND(errorphi,1,MPI_DOUBLE,0,i,MPI_COMM_WORLD,ierr)
			call MPI_SEND(errorSG,1,MPI_DOUBLE,0,2*i,MPI_COMM_WORLD,ierr)
			call MPI_SEND(errorT,1,MPI_DOUBLE,0,3*i,MPI_COMM_WORLD,ierr)
			call MPI_SEND(maxval(reshape(T, (/ Ntx*Nty*Ntz /) )),1,MPI_DOUBLE,0,4*i,MPI_COMM_WORLD,ierr)	

			if (threadpoints(2,2,1) == -1) then
				call MPI_SEND3D_DP(phi(2:Ntx-1,2:Nty-1,1:Ntz-1),0,1)
				call MPI_SEND3D_DP(n(2:Ntx-1,2:Nty-1,1:Ntz-1),0,2)
				call MPI_SEND3D_DP(p(2:Ntx-1,2:Nty-1,1:Ntz-1),0,3)
				call MPI_SEND3D_DP(EFx(2:Ntx-1,2:Nty-1,1:Ntz-1),0,4)
				call MPI_SEND3D_DP(EFy(2:Ntx-1,2:Nty-1,1:Ntz-1),0,5)
				call MPI_SEND3D_DP(EFz(2:Ntx-1,2:Nty-1,1:Ntz-1),0,6)
				call MPI_SEND3D_DP(currentx(2:Ntx-1,2:Nty-1,1:Ntz-1),0,7)
				call MPI_SEND3D_DP(currenty(2:Ntx-1,2:Nty-1,1:Ntz-1),0,8)
				call MPI_SEND3D_DP(currentz(2:Ntx-1,2:Nty-1,1:Ntz-1),0,9)
				call MPI_SEND3D_DP(heating(2:Ntx-1,2:Nty-1,1:Ntz-1),0,10)
				call MPI_SEND3D_DP(T(2:Ntx-1,2:Nty-1,1:Ntz-1),0,11)
				call MPI_SEND3D_DP(relativerecombination(2:Ntx-1,2:Nty-1,1:Ntz-1),0,12)
				call MPI_SEND3D_DP(recombination(2:Ntx-1,2:Nty-1,1:Ntz-1),0,13)
				
			elseif (threadpoints(2,2,3) == -2) then
				call MPI_SEND3D_DP(phi(2:Ntx-1,2:Nty-1,2:Ntz),0,1)
				call MPI_SEND3D_DP(n(2:Ntx-1,2:Nty-1,2:Ntz),0,2)
				call MPI_SEND3D_DP(p(2:Ntx-1,2:Nty-1,2:Ntz),0,3)
				call MPI_SEND3D_DP(EFx(2:Ntx-1,2:Nty-1,2:Ntz),0,4)
				call MPI_SEND3D_DP(EFy(2:Ntx-1,2:Nty-1,2:Ntz),0,5)
				call MPI_SEND3D_DP(EFz(2:Ntx-1,2:Nty-1,2:Ntz),0,6)
				call MPI_SEND3D_DP(currentx(2:Ntx-1,2:Nty-1,2:Ntz),0,7)
				call MPI_SEND3D_DP(currenty(2:Ntx-1,2:Nty-1,2:Ntz),0,8)
				call MPI_SEND3D_DP(currentz(2:Ntx-1,2:Nty-1,2:Ntz),0,9)
				call MPI_SEND3D_DP(heating(2:Ntx-1,2:Nty-1,2:Ntz),0,10)
				call MPI_SEND3D_DP(T(2:Ntx-1,2:Nty-1,2:Ntz),0,11)
				call MPI_SEND3D_DP(relativerecombination(2:Ntx-1,2:Nty-1,2:Ntz),0,12)
				call MPI_SEND3D_DP(recombination(2:Ntx-1,2:Nty-1,2:Ntz),0,13)

			else
				call MPI_SEND3D_DP(phi(2:Ntx-1,2:Nty-1,2:Ntz-1),0,1)
				call MPI_SEND3D_DP(n(2:Ntx-1,2:Nty-1,2:Ntz-1),0,2)
				call MPI_SEND3D_DP(p(2:Ntx-1,2:Nty-1,2:Ntz-1),0,3)
				call MPI_SEND3D_DP(EFx(2:Ntx-1,2:Nty-1,2:Ntz-1),0,4)
				call MPI_SEND3D_DP(EFy(2:Ntx-1,2:Nty-1,2:Ntz-1),0,5)
				call MPI_SEND3D_DP(EFz(2:Ntx-1,2:Nty-1,2:Ntz-1),0,6)
				call MPI_SEND3D_DP(currentx(2:Ntx-1,2:Nty-1,2:Ntz-1),0,7)
				call MPI_SEND3D_DP(currenty(2:Ntx-1,2:Nty-1,2:Ntz-1),0,8)
				call MPI_SEND3D_DP(currentz(2:Ntx-1,2:Nty-1,2:Ntz-1),0,9)
				call MPI_SEND3D_DP(heating(2:Ntx-1,2:Nty-1,2:Ntz-1),0,10)
				call MPI_SEND3D_DP(T(2:Ntx-1,2:Nty-1,2:Ntz-1),0,11)
				call MPI_SEND3D_DP(relativerecombination(2:Ntx-1,2:Nty-1,2:Ntz-1),0,12)
				call MPI_SEND3D_DP(recombination(2:Ntx-1,2:Nty-1,2:Ntz-1),0,13)
				
			end if
				
			do k = 1,smalliterationinbig
                                
                                call PoisonIteration(Ntx,Nty,Ntz,phi,n,p,ptimestep,errorPhi,V)
				call MPI_TRIPLE_ROTATE(phi,threadpoints,Ntx,Nty,Ntz)
						
				call getElectricField(Ntx,Nty,Ntz,phi,EFx,EFy,EFz)		
				call MPI_TRIPLE_ROTATE(EFx,threadpoints,Ntx,Nty,Ntz)
				call MPI_TRIPLE_ROTATE(EFy,threadpoints,Ntx,Nty,Ntz)
				call MPI_TRIPLE_ROTATE(EFz,threadpoints,Ntx,Nty,Ntz)
				
				
				ephi = phi + LUMO
				hphi = phi + HOMO
				
				if (i > number_of_phi_iterations) then
                                    call sgiteration(Ntx,Nty,Ntz,n,p,phi,phi+LUMO,phi+HOMO,T,&
                                        currentx,currenty,currentz,timestep,errorSG,&
                                        ntempcurrentx,ntempcurrenty,ntempcurrentz,&
                                        ptempcurrentx,ptempcurrenty,ptempcurrentz,&
                                        EFx,EFy,EFz,relativerecombination,recombination)
                                        
                                    call MPI_TRIPLE_ROTATE(n,threadpoints,Ntx,Nty,Ntz)
                                    call MPI_TRIPLE_ROTATE(p,threadpoints,Ntx,Nty,Ntz)
                                    call MPI_TRIPLE_ROTATE(relativerecombination,threadpoints,Ntx,Nty,Ntz)
                                    call MPI_TRIPLE_ROTATE(recombination,threadpoints,Ntx,Nty,Ntz)
                                end if
                                
                                if (i>number_of_sg_iterations+number_of_phi_iterations) then
                                    call getHeating(Ntx,Nty,Ntz,currentx,currenty,currentz,n,p,phi,T,heating,maxheating,&
                                        EFx,EFy,EFz,recombination)
                                    call MPI_TRIPLE_ROTATE(heating,threadpoints,Ntx,Nty,Ntz)
                                    
                                    if  (maxval(kappa) > 0.1) then
                                    
                                        call HeatIteration(Ntx,Nty,Ntz,T,heating,timestep,300.0d0, errorT,totalheating,&
                                                totaloutflow,&
                                                heatfluxx, heatfluxy, heatfluxz)
                                    else  
                                            call HeatIteration(Ntx,Nty,Ntz,T,heating,100*timestep,300.0d0, errorT,totalheating,&
                                                totaloutflow,&
                                                heatfluxx, heatfluxy, heatfluxz)
                                    end if
                                    
                                    call MPI_TRIPLE_ROTATE(T,threadpoints,Ntx,Nty,Ntz)
                                    
                                end if    
                        end do
                end do        
	end if
		
	if (task_ID == 0) then
		do superiterator = 1,(number_of_main_iterations)
						
			do i = 1, Ncutx+1
				do j = 1, Ncuty+1
					do k = 1, Ncutz+1
						cut = k+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
						call MPI_RECV(maxerrorphi(cut),1,MPI_DOUBLE,cut,superiterator,MPI_COMM_WORLD,status,ierr)
						call MPI_RECV(maxerrorSG(cut),1,MPI_DOUBLE,cut,2*superiterator,MPI_COMM_WORLD,status,ierr)
						call MPI_RECV(maxerrorT(cut),1,MPI_DOUBLE,cut,3*superiterator,MPI_COMM_WORLD,status,ierr)
						call MPI_RECV(maxT(cut),1,MPI_DOUBLE,cut,4*superiterator,MPI_COMM_WORLD,status,ierr)
					end do
				end do
			end do	
		
			
 			
 			do i = 1, Ncutx+1
				do j = 1, Ncuty+1
					do k = 1, Ncutz+1
						cut = k+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
						if (cutproperties(cut,3) == 1) then
							call MPI_RECV3D_DP(phi(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,1)
							call MPI_RECV3D_DP(n(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,2)
							call MPI_RECV3D_DP(p(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,3)
							call MPI_RECV3D_DP(EFx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,4)
							call MPI_RECV3D_DP(EFy(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,5)
							call MPI_RECV3D_DP(EFz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,6)
							call MPI_RECV3D_DP(currentx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,7)
							call MPI_RECV3D_DP(currenty(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,8)
							call MPI_RECV3D_DP(currentz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,9)
							call MPI_RECV3D_DP(heating(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,10)
							call MPI_RECV3D_DP(T(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,11)
							call MPI_RECV3D_DP(relativerecombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,12)
							call MPI_RECV3D_DP(recombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3):cutproperties(cut,6)-1),&
							cut,13)
							
						elseif (cutproperties(cut,6) == Nz) then
							call MPI_RECV3D_DP(phi(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,1)
							call MPI_RECV3D_DP(n(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,2)
							call MPI_RECV3D_DP(p(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,3)							
							call MPI_RECV3D_DP(EFx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,4)
							call MPI_RECV3D_DP(EFy(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,5)
							call MPI_RECV3D_DP(EFz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,6)
							call MPI_RECV3D_DP(currentx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,7)
							call MPI_RECV3D_DP(currenty(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,8)
							call MPI_RECV3D_DP(currentz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,9)
							call MPI_RECV3D_DP(heating(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,10)
							call MPI_RECV3D_DP(T(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,11)
							call MPI_RECV3D_DP(relativerecombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,12)
							call MPI_RECV3D_DP(recombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)),&
							cut,13)
							
						else
							call MPI_RECV3D_DP(phi(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,1)
							call MPI_RECV3D_DP(n(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,2)
							call MPI_RECV3D_DP(p(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,3)
							call MPI_RECV3D_DP(EFx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,4)
							call MPI_RECV3D_DP(EFy(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,5)
							call MPI_RECV3D_DP(EFz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,6)
							call MPI_RECV3D_DP(currentx(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,7)
							call MPI_RECV3D_DP(currenty(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,8)
							call MPI_RECV3D_DP(currentz(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,9)
							call MPI_RECV3D_DP(heating(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,10)
							call MPI_RECV3D_DP(T(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,11)
							call MPI_RECV3D_DP(relativerecombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,12)
							call MPI_RECV3D_DP(recombination(cutproperties(cut,1)+1:cutproperties(cut,4)-1,&
							cutproperties(cut,2)+1:cutproperties(cut,5)-1,cutproperties(cut,3)+1:cutproperties(cut,6)-1),&
							cut,13)
							
						end if
				
					end do
				end do
			end do
			call makebackup(phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T,superiterator,&
			relativerecombination,recombination)

			maxerrorphi(0) = 0.0d0
			maxerrorSG(0) = 0.0d0
			maxerrorT(0) = 0.0d0
			maxT(0) = 0.0d0
			maxerrorphi(0) = maxval(dabs(maxerrorphi))
			maxerrorSG(0) = maxval(dabs(maxerrorSG))
			maxerrorT(0) = maxval(dabs(maxerrorT))
			maxT(0) = maxval(dabs(maxT))
 			print*,superiterator, maxerrorphi(0),maxerrorSG(0), maxerrorT(0),maxT(0)
			
			if ((maxerrorphi(0)<1d-10) .and. (maxerrorSG(0)<1d-10) .and. (maxerrorT(0)<1d-10)) then
				call CPU_TIME(endtime)				
				print*, "CPU_TIME",endtime-starttime
				call MPI_Abort(MPI_COMM_WORLD, 25,ierr)
			end if
		end do
	end if
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!	if (task_ID == 0) then
!		call CPU_TIME(endtime)
!		totalheating = 0.0d0
!		totaloutflow = 0.0d0
!		do i = 1,n_tasks-1
!			call MPI_RECV(receive,1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,status,ierr)
!!			call MPI_RECV(receive,1,MPI_DOUBLE,i,i+n_tasks,MPI_COMM_WORLD,status,ierr)
	!		totaloutflow = totaloutflow + receive
	!	end do
! 		print*, "heating/cooling relative difference:", (totalheating/totaloutflow-1.0d0)
!	else
!		call MPI_SEND(totalheating,1,MPI_DOUBLE,0,task_ID,MPI_COMM_WORLD,ierr)
!		call MPI_SEND(totaloutflow,1,MPI_DOUBLE,0,task_ID+n_tasks,MPI_COMM_WORLD,ierr)
!	end if
	
	contains
	
	subroutine cuteverything(Nx,Ny,Nz,Ncutx,Ncuty,Ncutz,cutproperties)
		integer, intent(in) :: Nx,Ny,Nz,Ncutx,Ncuty,Ncutz
		integer, dimension(:,:),intent(inout) :: cutproperties
		integer :: subx(Ncutx+2),suby(Ncuty+2),subz(Ncutz+2)
		integer :: i,j,k,cut
		
		subx(1)=0
		cut = floor((Nx+2.0)/(Ncutx+1.0))
		do i = 2,size(subx)-1
			subx(i)=subx(i-1)+cut
		end do
		
		subx(size(subx))=Nx+1
		suby(1)=0
		cut = floor((Ny+2.0)/(Ncuty+1.0))
		do i = 2,size(suby)-1
			suby(i)=suby(i-1)+cut
		end do
		suby(size(suby))=Ny+1
		
		subz(1)=1
		cut = floor((Nz)/(Ncutz+1.0))
		do i = 2,size(subz)-1
			subz(i)=subz(i-1)+cut
		end do
		subz(size(subz))=Nz
		
		do i = 1, Ncutx+1
			do j = 1, Ncuty+1
				do k = 1, Ncutz+1
					cut = k+(j-1)*(Ncutz+1)+(i-1)*(Ncutz+1)*(Ncuty+1)
					if (i .ne. 1) then
						cutproperties(cut,1) = subx(i)-1
					else
						cutproperties(cut,1) = subx(i)
					end if
					if (j .ne. 1) then
						cutproperties(cut,2) = suby(j)-1
					else
						cutproperties(cut,2) = suby(j)
					end if
					if (k .ne. 1) then
						cutproperties(cut,3) = subz(k)-1
					else
						cutproperties(cut,3) = subz(k)
					end if			
					cutproperties(cut,4) = subx(i+1)
					cutproperties(cut,5) = suby(j+1)
					cutproperties(cut,6) = subz(k+1)	
					end do
			end do
		end do
		

	end subroutine cuteverything	
	
	subroutine makebackup(phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T,Iiteration,relativerecombination,recombination)
		double precision, dimension(:,:,:),intent(in) :: phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T,&
		relativerecombination,recombination
		integer :: Nx,Ny,Nz
		integer,intent(in) :: Iiteration
		character(len=9) :: iteration
		
		Nx = SIZE(T(:,1,1))
		Ny = SIZE(T(1,:,1))
		Nz = SIZE(T(1,1,:))
		
		write(iteration,"(I9)") Iiteration
			
		open(unit = 555, file = "backup/phi"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) phi
		close(555)
		
		open(unit = 555, file = "backup/EFx"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) EFx
		close(555)
		
		open(unit = 555, file = "backup/EFy"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) EFy
		close(555)
		
		open(unit = 555, file = "backup/EFz"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) EFz
		close(555)
		
		open(unit = 555, file = "backup/n"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) n
		close(555)
		
		open(unit = 555, file = "backup/p"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) p
		close(555)
		
		open(unit = 555, file = "backup/currentx"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) currentx
		close(555)
		
		open(unit = 555, file = "backup/currenty"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) currenty
		close(555)
		
		open(unit = 555, file = "backup/currentz"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) currentz
		close(555)
		
		open(unit = 555, file = "backup/heating"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) heating
		close(555)
		
		open(unit = 555, file = "backup/T"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) T
		close(555)
		
		open(unit = 555, file = "backup/relativerecombination"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) relativerecombination
		close(555)
		
		open(unit = 555, file = "backup/recombination"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) recombination
		close(555)

		open(unit = 555, file = "backup/jouleheating"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*) currentz*EFz+currenty*EFy+currentx*EFx
		close(555)	

		open(unit = 555, file = "backup/recombinationheating"//trim(ADJUSTL(iteration))//".bkp",action = "write", status = "replace")
		write(555,*)  0.75 * (Erecombination)*recombination
		close(555)	
	end subroutine makebackup
		
		
	subroutine quantitiesbackup(phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T,Iiteration)
		double precision, dimension(:,:,:) :: phi,EFx,EFy,EFz,n,p,currentx,currenty,currentz,heating,T
		integer :: Nx,Ny,Nz
		integer :: Iiteration
		character(len=9) :: iteration
		
		Nx = SIZE(T(:,1,1))
		Ny = SIZE(T(1,:,1))
		Nz = SIZE(T(1,1,:))
		
		write(iteration,"(I9)") Iiteration
			
		open(unit = 555, file = "backup/phi"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) phi
		close(555)
		
		open(unit = 555, file = "backup/EFx"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) EFx
		close(555)
		
		open(unit = 555, file = "backup/EFy"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) EFy
		close(555)
		
		open(unit = 555, file = "backup/EFz"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) EFz
		close(555)
		
		open(unit = 555, file = "backup/n"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) n
		close(555)
		
		open(unit = 555, file = "backup/p"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) p
		close(555)
		
		open(unit = 555, file = "backup/currentx"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) currentx
		close(555)
		
		open(unit = 555, file = "backup/currenty"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) currenty
		close(555)
		
		open(unit = 555, file = "backup/currentz"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) currentz
		close(555)
		
		open(unit = 555, file = "backup/heating"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) heating
		close(555)
		
		open(unit = 555, file = "backup/T"//trim(ADJUSTL(iteration))//".bkp",action = "read")
		read(555,*) T
		close(555)
	
	end subroutine quantitiesbackup	
	
	subroutine MPI_SEND2D_INT(array,destination,tag)
		integer,intent(in) :: array(:,:),destination,tag
		integer, allocatable :: buffer(:)
		integer :: Nx,Ny,i,j		
		
		Nx = SIZE(array(:,1))
		Ny = SIZE(array(1,:))
		
		allocate(buffer(Nx*Ny+2))

		do i = 1,Nx
			do j = 1,Ny
				buffer(i+(j-1)*Nx) = array(i,j)
			end do
		end do
		buffer(Nx*Ny+1) = Nx
		buffer(Nx*Ny+2) = Ny
		
		call MPI_SEND(buffer,Nx*Ny+2,MPI_INTEGER,destination,tag,MPI_COMM_WORLD,ierr)		

	end subroutine MPI_SEND2D_INT

	subroutine MPI_SEND3D_INT(array,destination,tag)
		integer,intent(in) :: array(:,:,:),destination,tag
		integer, allocatable :: buffer(:)
		integer :: Nx,Ny,Nz,i,j,k
		
		Nx = SIZE(array(:,1,1))
		Ny = SIZE(array(1,:,1))
		Nz = SIZE(array(1,1,:))
		
		allocate(buffer(Nx*Ny*Nz+3))

		do i = 1,Nx
			do j = 1,Ny
				do k = 1,Nz
					buffer(i+(j-1)*Nx+(k-1)*Nx*Ny) = array(i,j,k)
				end do
			end do
		end do
		buffer(Nx*Ny*Nz+1) = Nx
		buffer(Nx*Ny*Nz+2) = Ny
		buffer(Nx*Ny*Nz+3) = Nz

		call MPI_SEND(buffer,Nx*Ny*Nz+3,MPI_INTEGER,destination,tag,MPI_COMM_WORLD,ierr)		

	end subroutine MPI_SEND3D_INT

	subroutine MPI_RECV3D_INT(array,source,tag)
		integer,intent(in) :: source,tag
		integer,dimension(:,:,:) :: array
		integer, allocatable :: buffer(:)
		integer :: Nx,Ny,Nz,i,j,k,buffersize,ierr
		INTEGER :: status(MPI_STATUS_SIZE)
		print*, "receiving",source,tag		
		Nx = SIZE(array(:,1,1))
		Ny = SIZE(array(1,:,1))
		Nz = SIZE(array(1,1,:))

		call MPI_PROBE(source,tag,MPI_COMM_WORLD,status,ierr)
		call MPI_GET_COUNT(status,MPI_INTEGER,buffersize,ierr)
		allocate(buffer(buffersize))

		call MPI_RECV(buffer,buffersize,MPI_INTEGER,source,tag,MPI_COMM_WORLD,status,ierr)
		if ((buffer(SIZE(buffer)-2) .ne. Nx) .or. (buffer(SIZE(buffer)-1) .ne. Ny) &
					.or. (buffer(SIZE(buffer)) .ne. Nz)) then
! 			print*, "ERROR IN MPI_RECV3D_INT, sent array size differs from recieved array"
			return
		else
			do i = 1,Nx
				do j = 1,Ny
					do k = 1,Nz
						array(i,j,k) = buffer(i+(j-1)*Nx+(k-1)*Nx*Ny)
					end do
				end do
			end do
		end if


	end subroutine MPI_RECV3D_INT

	subroutine MPI_SEND2D_DP(array,destination,tag)
		integer,intent(in) :: destination,tag
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,i,j		
		double precision,intent(in) :: array(:,:)
		
		Nx = SIZE(array(:,1))
		Ny = SIZE(array(1,:))
		
		allocate(buffer(Nx*Ny+2))

		do i = 1,Nx
			do j = 1,Ny
				buffer(i+(j-1)*Nx) = array(i,j)
			end do
		end do
		buffer(Nx*Ny+1) = real(Nx,8)
		buffer(Nx*Ny+2) = real(Ny,8)
		
		call MPI_SEND(buffer,Nx*Ny+2,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD,ierr)		

	end subroutine MPI_SEND2D_DP

	subroutine MPI_RECV2D_DP(array,source,tag)
		integer,intent(in) :: source,tag
		double precision,dimension(:,:) :: array
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,i,j,buffersize,ierr
		INTEGER :: status(MPI_STATUS_SIZE)
		
		Nx = SIZE(array(:,1))
		Ny = SIZE(array(1,:))

		call MPI_PROBE(source,tag,MPI_COMM_WORLD,status,ierr)
		call MPI_GET_COUNT(status,MPI_DOUBLE,buffersize,ierr)
		allocate(buffer(buffersize))
		call MPI_RECV(buffer,buffersize,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,status,ierr)
		if ((buffer(SIZE(buffer)-1) .ne. real(Nx,8)) .or. (buffer(SIZE(buffer)) .ne. real(Ny,8)))then
			print*, "ERROR IN MPI_RECV2D_DP, sent array size differs from recieved array"
			return
		else
			do i = 1,Nx
				do j = 1,Ny
					array(i,j) = buffer(i+(j-1)*Nx)
				end do
			end do
		end if

	end subroutine MPI_RECV2D_DP

	subroutine MPI_SEND2D_DP_NB(array,destination,tag)
		integer,intent(in) :: destination,tag
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,i,j		
		double precision,intent(in) :: array(:,:)
		
		Nx = SIZE(array(:,1))
		Ny = SIZE(array(1,:))
		
		allocate(buffer(Nx*Ny+2))

		do i = 1,Nx
			do j = 1,Ny
				buffer(i+(j-1)*Nx) = array(i,j)
			end do
		end do
		buffer(Nx*Ny+1) = real(Nx,8)
		buffer(Nx*Ny+2) = real(Ny,8)
		
		call MPI_ISEND(buffer,Nx*Ny+2,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD,ierr)		

	end subroutine MPI_SEND2D_DP_NB

	subroutine MPI_RECV2D_DP_NB(array,source,tag)
		integer,intent(in) :: source,tag
		double precision,dimension(:,:) :: array
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,i,j,buffersize,ierr
		INTEGER :: status(MPI_STATUS_SIZE)
		
		Nx = SIZE(array(:,1))
		Ny = SIZE(array(1,:))

		call MPI_PROBE(source,tag,MPI_COMM_WORLD,status,ierr)
		call MPI_GET_COUNT(status,MPI_DOUBLE,buffersize,ierr)
		allocate(buffer(buffersize))
		call MPI_IRECV(buffer,buffersize,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,status,ierr)
		if ((buffer(SIZE(buffer)-1) .ne. real(Nx,8)) .or. (buffer(SIZE(buffer)) .ne. real(Ny,8)))then
			print*, "ERROR IN MPI_RECV2D_DP, sent array size differs from recieved array"
			return
		else
			do i = 1,Nx
				do j = 1,Ny
					array(i,j) = buffer(i+(j-1)*Nx)
				end do
			end do
		end if

	end subroutine MPI_RECV2D_DP_NB
	
	subroutine MPI_SEND3D_DP(array,destination,tag)
		integer,intent(in) :: destination,tag
		double precision,dimension(:,:,:),intent(in) :: array
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,Nz,i,j,k

		Nx = SIZE(array(:,1,1))
		Ny = SIZE(array(1,:,1))
		Nz = SIZE(array(1,1,:))

		allocate(buffer(Nx*Ny*Nz+3))

		do i = 1,Nx
			do j = 1,Ny
				do k = 1,Nz
					buffer(i+(j-1)*Nx+(k-1)*Nx*Ny) = array(i,j,k)
				end do
			end do
		end do

		buffer(Nx*Ny*Nz+1) = real(Nx,8)
		buffer(Nx*Ny*Nz+2) = real(Ny,8)
		buffer(Nx*Ny*Nz+3) = real(Nz,8)
		call MPI_SEND(buffer,Nx*Ny*Nz+3,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD,ierr)		
	end subroutine MPI_SEND3D_DP

	subroutine MPI_RECV3D_DP(array,source,tag)
		integer,intent(in) :: source,tag
		double precision,dimension(:,:,:) :: array
		double precision, allocatable :: buffer(:)
		integer :: Nx,Ny,Nz,i,j,k,buffersize,ierr
		INTEGER :: status(MPI_STATUS_SIZE)
		
		Nx = SIZE(array(:,1,1))
		Ny = SIZE(array(1,:,1))
		Nz = SIZE(array(1,1,:))

		call MPI_PROBE(source,tag,MPI_COMM_WORLD,status,ierr)
		call MPI_GET_COUNT(status,MPI_DOUBLE,buffersize,ierr)
		allocate(buffer(buffersize))
		call MPI_RECV(buffer,buffersize,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,status,ierr)

		if ((buffer(SIZE(buffer)-2) .ne. real(Nx,8)) .or. (buffer(SIZE(buffer)-1) .ne. real(Ny,8)) &
					.or. (buffer(SIZE(buffer)) .ne. real(Nz,8))) then
			print*, "ERROR IN MPI_RECV3D_DP, sent array size differs from recieved array"
			print*, buffer(SIZE(buffer)-2), Nx
			print*, buffer(SIZE(buffer)-1), Ny
			print*, buffer(SIZE(buffer)), Nz
			return
		else
			do i = 1,Nx
				do j = 1,Ny
					do k = 1,Nz
						array(i,j,k) = buffer(i+(j-1)*Nx+(k-1)*Nx*Ny)
					end do
				end do
			end do
		end if


	end subroutine MPI_RECV3D_DP

	subroutine MPI_ROTATE(arraytosend,arraytoget,wheretosend,wheretoget,tag)
		double precision, dimension(:,:) :: arraytosend,arraytoget
		integer,intent(in) :: wheretosend,wheretoget,tag
		if (wheretosend	> 0) then
			call MPI_SEND2D_DP(arraytosend,wheretosend,tag)
		end if
		if (wheretoget	> 0) then
			call MPI_RECV2D_DP(arraytoget,wheretoget,tag)
		end if
	end subroutine MPI_ROTATE
	
	subroutine MPI_TRIPLE_ROTATE(array,threadpoints,Ntx,Nty,Ntz)
		double precision, dimension(:,:,:) :: array
		integer ,intent(in) :: Ntx, Nty, Ntz,threadpoints(3,3,3)
		
		call MPI_ROTATE(array(:,2,:),array(:,Nty,:),threadpoints(2,1,2),threadpoints(2,3,2),1)
		call MPI_ROTATE(array(:,Nty-1,:),array(:,1,:),threadpoints(2,3,2),threadpoints(2,1,2),2)
		
		call MPI_ROTATE(array(2,:,:),array(Ntx,:,:),threadpoints(1,2,2),threadpoints(3,2,2),1)
		call MPI_ROTATE(array(Ntx-1,:,:),array(1,:,:),threadpoints(3,2,2),threadpoints(1,2,2),2)
		
		call MPI_ROTATE(array(:,:,2),array(:,:,Ntz),threadpoints(2,2,1),threadpoints(2,2,3),1)
		call MPI_ROTATE(array(:,:,Ntz-1),array(:,:,1),threadpoints(2,2,3),threadpoints(2,2,1),2)

	end subroutine MPI_TRIPLE_ROTATE

	subroutine writecalculationinformation(Nx,Ny,Nz)
		integer, intent(in) :: Nx,Ny,Nz
		open(unit = 555, file = "scripts/calculationinformation.info",action = "write", status = "replace")
		write(555,*) Nx,Ny,Nz
		close(555)		
	end subroutine writecalculationinformation
	
end program main

