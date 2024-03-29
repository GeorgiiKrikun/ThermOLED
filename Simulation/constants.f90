module constants
	implicit none
	
	double precision, parameter :: q=1.6d-19,eps=8.85d-12, kt = 1.38d-23,langrec = 0.0d0, &
	Erecombination = 1.6d-19, pi=3.1415926,hplank = 6.62d-34, Rydberg = 1.097d7, &
	BorRadius = 5.29d-11,hbar = 1.054d-34,newtonh_bot = 5.0d3,newtonh_top = 5.0d3
	double precision, dimension(:,:,:), allocatable :: epsR,rho,cp,kappa,HOMO,LUMO,psigma,nsigma
	double precision, dimension(:), allocatable :: subgridx, subgridy, subgridz
	
	
	integer, dimension(:,:,:), allocatable :: threadpoints,materials
	integer :: task_ID
	
	type materialproperties
		integer :: matnum
		double precision, allocatable, dimension(:) :: nmu0, nsigma, nenergy
		double precision, allocatable, dimension(:) :: pmu0, psigma, penergy
		double precision, allocatable, dimension(:) :: lattice, epsR, rho, cp, kappa,matstart,matend,layer
	end type materialproperties
	
	type(materialproperties) :: composition
	
	contains
	
! 	loading material composition from input/input.dat

	subroutine callbackup(Nx,Ny,Nz,backupdata,filename)
		integer, intent(in) :: Nx,Ny,Nz
		double precision, intent(in), dimension(Nx,Ny,Nz) :: backupdata
		character(len=*), intent(in) :: filename
		integer :: i,j,k
		open(unit = 555, file = filename,action = "write", status = "replace")
		
		write(555,*) Nx,Ny,Nz
		
		do i = 1,Nx
			do j = 1,Ny
				do k = 1,Nz
					write(555,*) i,j,k,backupdata(i,j,k)
				end do
			end do
		end do
			
		close(555)
		
	end subroutine callbackup

	subroutine loadbackup(Nx,Ny,Nz,backupdata,filename)
		integer, intent(in) :: Nx,Ny,Nz
		double precision, intent(inout), dimension(Nx,Ny,Nz) :: backupdata
		character(len=*), intent(in) :: filename
		integer :: i,j,k,Nxcheck,Nycheck,Nzcheck,IOSTATUS
		double precision :: readdata
		open(unit = 555, file = filename,action = "read")
		
		read(555,*) Nxcheck,Nycheck,Nzcheck
		if (.not.((Nx == Nxcheck) .and. (Ny == Nycheck) .and. (Nz == Nzcheck))) then
			print*, "ERROR WHILE LOADING DATA. PROBLEMS WITH DIMENSIONALITY"
		end if
		
		do 
			read(555,*,IOSTAT = IOSTATUS) i,j,k,readdata
			if (IOstatus == 0) then
				backupdata(i,j,k) = readdata
			else
				exit
			end if
		end do
			
		close(555)
		
	end subroutine loadbackup

	subroutine obtainmaterialcomposition()
		integer :: i,n
		
		open(unit = 101, file = "input/input.dat")
		read(101,*)
		read(101,*) composition%matnum
		n = composition%matnum
				
		allocate(composition%pmu0(n),composition%psigma(n),&
			composition%penergy(n),composition%nmu0(n),composition%nsigma(n)&
			,composition%nenergy(n),composition%lattice(n),composition%epsR(n),&
			composition%rho(n),composition%cp(n),composition%kappa(n),&
			composition%matstart(n),composition%matend(n))
		print*, n
		do i=1,n
			read(101,*) composition%pmu0(i),composition%psigma(i),composition%penergy(i),&
				composition%nmu0(i),composition%nsigma(i),composition%nenergy(i),&
				composition%lattice(i),composition%epsR(i),&
				composition%rho(i),composition%cp(i),composition%kappa(i),&
				composition%matstart(i),composition%matend(i)
			print*, composition%pmu0(i),composition%psigma(i),composition%penergy(i),&
				composition%nmu0(i),composition%nsigma(i),composition%nenergy(i),&
				composition%lattice(i),composition%epsR(i),&
				composition%rho(i),composition%cp(i),composition%kappa(i),&
				composition%matstart(i),composition%matend(i)
		end do
		
		call Sleep(1)
		
		close(101)
	end subroutine obtainmaterialcomposition
	
	double precision function mu(n,T,E,charge,material) !charge = 1 - holes, -1 - electrons 
		double precision, intent(in) :: n,T,E
		integer, intent(in) :: charge,material
		
		double precision :: c1=1.8d-9,c2=0.42,delta,sigma,shat,mu0,lattice
		
		if (charge == 1) then
			mu0 = composition%pmu0(material)
			sigma = composition%psigma(material)*kt*300
		elseif (charge == -1) then
			mu0 = composition%nmu0(material)
			sigma = composition%nsigma(material)*kt*300
		end if
		
! 		sigma = 2.0d0*kt*300
		lattice = composition%lattice(material)
		
		shat = sigma/(kt*T)
		
		delta = 2*(dlog(shat**2-shat)-0.32663425997828094)/shat**2
		mu0 = mu0*c1*dexp(-c2*shat**2) !T
		mu0 = mu0*dexp(0.5*(shat**2-shat)*(2*n*lattice**3)**delta) !n,p
		mu0 = mu0*dexp(0.44*(shat**1.5 - 2.2)*(dsqrt(1.0d0+0.8d0*(E*q*lattice/sigma)**2.0d0)-1.0d0))
		mu = mu0
! 		mu = 3.67d-7
!  		print*, mu,material
! 		print*, "here2"
! 		if (abs(T - 300.0d0)<1.0d-6) then
! 			mu = 3.67d-6
! 		else
! 			mu = 3.67d-7
! 		end if
				
	end function mu

	!Diffusion constant
	double precision function Dc(n,T,E,charge,material) !charge = 1 - holes, -1 - electrons 
		double precision, intent(in) :: n,T,E	
		integer, intent(in) :: charge,material
		
	 	Dc = mu(n,T,E,charge,material)*kt*T/q
	end function Dc

	double precision function S(n,T,material)
		double precision, intent(in) :: n, T
		integer, intent(in) :: material
		
! 		S = 1.0d-3
		S = 0
	end function S
	
	double precision function dSdT(n,T,material)
		double precision, intent(in) :: n, T
		integer, intent(in) :: material
		dSdT = 0.0d0
	end function dSdT
end module constants
