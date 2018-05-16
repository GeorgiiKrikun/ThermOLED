subroutine PoisonIteration(Nx,Ny,Nz,phi,n,p,timestep,error,V)
	use constants
	implicit none
	
	integer :: i,j,k
	integer, intent(in) :: Nx,Ny,Nz
	double precision, intent(in),dimension(Nx,Ny,Nz) :: n,p
	double precision, intent(inout),dimension(Nx,Ny,Nz) :: phi
	double precision, intent(inout) :: error
	double precision, intent(in) :: V
	double precision, dimension(Nx,Ny,Nz) :: newphi
	double precision :: uphi,dphi,lphi,rphi,ophi,iphi,pphi,timestep
	double precision :: udist,ddist,ldist,rdist,odist,idist, pp, nn
	newphi = phi
	
	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 2,Nz-1
				uphi = phi(i,j,k-1)
				dphi = phi(i,j,k+1)
				lphi = phi(i-1,j,k)
				rphi = phi(i+1,j,k)
				ophi = phi(i,j-1,k)
				iphi = phi(i,j+1,k)
				pphi = phi(i,j,k)
				
				udist = dabs(subgridz(k)-subgridz(k-1))
				ddist = dabs(subgridz(k)-subgridz(k+1))
				ldist = dabs(subgridx(i)-subgridx(i-1))
				rdist = dabs(subgridx(i)-subgridx(i+1))
				odist = dabs(subgridy(j)-subgridy(j-1))
				idist = dabs(subgridy(j)-subgridy(j+1))
				
				nn = n(i,j,k)
				pp = p(i,j,k)
				newphi(i,j,k) = -q*(nn-pp)/(eps*epsR(i,j,k))
				newphi(i,j,k) = newphi(i,j,k) + 2*(rphi/rdist+lphi/ldist-pphi*(1d0/rdist+1d0/ldist))/(rdist+ldist)
				newphi(i,j,k) = newphi(i,j,k) + 2*(dphi/ddist+uphi/udist-pphi*(1d0/ddist+1d0/udist))/(ddist+udist)
				newphi(i,j,k) = newphi(i,j,k) + 2*(iphi/idist+ophi/odist-pphi*(1d0/idist+1d0/odist))/(idist+odist)
				newphi(i,j,k) = pphi+newphi(i,j,k)*timestep
			end do
		end do
	end do
	
! 	if (threadpoints(2,2,1) == -1) then
! 		do i = 2,Nx-1
! 			do j = 2,Ny-1
! 				newphi(i,j,1) = phi(i,j,1)
! 			end do
! 		end do
! 	end if
	
! 	if (threadpoints(2,2,3) == -2) then
! 		do i = 2,Nx-1
! 			do j = 2,Ny-1
! 				newphi(i,j,Nz) = phi(i,j,Nz)
! 			end do
! 		end do
! 	end if
	
	error = 0.0d0
	
	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 2,Nz-1
				if (phi(i,j,k) .ne. 0d0) then
					if (((newphi(i,j,k)/phi(i,j,k)-1.0d0) > error) .and. ((i-1)*(j-1)*(k-1) .ne. 0d0) &
					.and. ((i-Nx)*(j-Ny)*(k-Nz) .ne. 0d0)) then
						error = abs(newphi(i,j,k)/phi(i,j,k)-1.0d0)
					end if
				end if
			end do
		end do
	end do
	phi = newphi
	
end subroutine PoisonIteration

! subroutine PoisonIteration(Nx,Ny,Nz,phi,n,p,error,V)
! 	use constants
! 	implicit none
! 	
! 	integer :: i,j,k
! 	integer, intent(in) :: Nx,Ny,Nz
! 	double precision, intent(in),dimension(Nx,Ny,Nz) :: n,p
! 	double precision, intent(inout),dimension(Nx,Ny,Nz) :: phi
! 	double precision, intent(inout) :: error
! 	double precision, intent(in) :: V
! 	double precision, dimension(Nx,Ny,Nz) :: newphi
! 	double precision :: uphi,dphi,lphi,rphi,ophi,iphi,pphi
! 	double precision :: udist,ddist,ldist,rdist,odist,idist, pp, nn
! 	newphi = phi
! 	
! 	do i = 2,Nx-1
! 		do j = 2,Ny-1
! 			do k = 2,Nz-1
! 				uphi = phi(i,j,k-1)
! 				dphi = phi(i,j,k+1)
! 				lphi = phi(i-1,j,k)
! 				rphi = phi(i+1,j,k)
! 				ophi = phi(i,j-1,k)
! 				iphi = phi(i,j+1,k)
! 				
! 				udist = dabs(subgridz(k)-subgridz(k-1))
! 				ddist = dabs(subgridz(k)-subgridz(k+1))
! 				ldist = dabs(subgridx(i)-subgridx(i-1))
! 				rdist = dabs(subgridx(i)-subgridx(i+1))
! 				odist = dabs(subgridy(j)-subgridy(j-1))
! 				idist = dabs(subgridy(j)-subgridy(j+1))
! 				
! 				nn = n(i,j,k)
! 				pp = p(i,j,k)
! 				
! 				newphi(i,j,k) = q*(nn-pp)/(eps*epsR(i,j,k))
! 				newphi(i,j,k) = newphi(i,j,k) - 2*((rphi)/(rdist)+(lphi)/(ldist))/(rdist+ldist) - &
! 				2*((uphi)/(udist)+(dphi)/(ddist))/(udist+ddist) - &
! 				2*((ophi)/(odist)+(iphi)/(idist))/(odist+idist)
! 				if ((i == 3) .and. (j == 2) .and. (k == 3)) then
!						write(*,'(4D10.2)') 2*((rphi)/(rdist)+(lphi)/(ldist))/(rdist+ldist),&
!						2*((ophi)/(odist)+(iphi)/(idist))/(odist+idist),&
!						2*((uphi)/(udist)+(dphi)/(ddist))/(udist+ddist),&
!						- 2*((rphi)/(rdist)+(lphi)/(ldist))/(rdist+ldist) - &
!						2*((uphi)/(udist)+(dphi)/(ddist))/(udist+ddist) - &
!						2*((ophi)/(odist)+(iphi)/(idist))/(odist+idist)
! 					
! 				end if
! 				
! 				newphi(i,j,k) = newphi(i,j,k)/&
! 				(-2*(1/ldist + 1/rdist)/(ldist+rdist)-2*(1/udist + 1/ddist)/(udist+ddist)-2*(1/idist + 1/odist)/(idist+odist))
! 				if ((i == 3) .and. (j == 2) .and. (k == 3)) then
!						write(*,'(4D10.2)') 2*(1/ldist + 1/rdist)/(ldist+rdist),2*(1/idist + 1/odist)/(idist+odist),&
!						2*(1/udist + 1/ddist)/(udist+ddist),&
!	 					(-2*(1/ldist + 1/rdist)/(ldist+rdist)-2*(1/udist + 1/ddist)/(udist+ddist)-2*(1/idist + 1/odist)/(idist+odist))
! 				end if
! 				
! 			end do
! 		end do
! 	end do
! 	
! 	if (threadpoints(2,2,1) == -1) then
! 		do i = 2,Nx-1
! 			do j = 2,Ny-1
! 				newphi(i,j,1) = V	
! 			end do
! 		end do
! 	end if
! 	
! 	if (threadpoints(2,2,3) == -2) then
! 		do i = 2,Nx-1
! 			do j = 2,Ny-1
! 				newphi(i,j,Nz) = 0.0d0
! 			end do
! 		end do
! 	end if
! 		
! 	error = 0.0d0
! 	do i = 1,Nx
! 		do j = 1,Ny
! 			do k = 1,Nz
! 				if (((newphi(i,j,k)/phi(i,j,k)-1.0d0) > error) .and. ((i-1)*(j-1)*(k-1) .ne. 0) &
! 				.and. ((i-Nx)*(j-Ny)*(k-Nz) .ne. 0)) then
! 					error = abs(newphi(i,j,k)/phi(i,j,k)-1.0d0)
! 				end if
! 			end do
! 		end do
! 	end do
! 		
! 	phi = newphi
! 		
! end subroutine PoisonIteration

subroutine GetElectricField(Nx,Ny,Nz,phi,EFx,EFy,EFz)
	use constants
	implicit none
	
	integer :: i,j,k
	integer, intent(in) :: Nx,Ny,Nz
	double precision, intent(in), dimension(Nx,Ny,Nz) :: phi
	double precision,  dimension(Nx,Ny,Nz) :: EFx,EFy,EFz
	double precision :: uphi,dphi,lphi,rphi,ophi,iphi,pphi
	double precision :: udist,ddist,ldist,rdist,odist,idist
	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 1,Nz
				if ((threadpoints(2,2,1) .ne. -1) .and. (k == 1)) then
					cycle
				end if
				if ((threadpoints(2,2,3) .ne. -2) .and. (k == Nz)) then
					cycle
				end if
				
				lphi = phi(i-1,j,k)
				rphi = phi(i+1,j,k)
				ophi = phi(i,j-1,k)
				iphi = phi(i,j+1,k)
				pphi = phi(i,j,k)
				
				ldist = dabs(subgridx(i)-subgridx(i-1))
				rdist = dabs(subgridx(i)-subgridx(i+1))
				odist = dabs(subgridy(j)-subgridy(j-1))
				idist = dabs(subgridy(j)-subgridy(j+1))
				
				if ((k .ne. 1 ) .and. (k .ne. Nz)) then
					uphi = phi(i,j,k-1)
					dphi = phi(i,j,k+1)
					udist = dabs(subgridz(k)-subgridz(k-1))
					ddist = dabs(subgridz(k)-subgridz(k+1))
					EFz(i,j,k) = -(dphi-uphi)/((ddist+udist))
				elseif (k == 1 ) then
					dphi = phi(i,j,k+1)
					ddist = dabs(subgridz(k)-subgridz(k+1))
					EFz(i,j,k) = -(dphi-pphi)/(ddist)
				elseif (k == Nz ) then
					uphi = phi(i,j,k-1)
					udist = dabs(subgridz(k)-subgridz(k-1))
					EFz(i,j,k) = -(pphi-uphi)/(udist)
				end if
				EFx(i,j,k) = -(rphi-lphi)/((ldist+rdist))
				EFy(i,j,k) = -(iphi-ophi)/((idist+odist))
				
			end do
		end do
	end do
	
	
end subroutine GetElectricField
	
subroutine GetHeating(Nx,Ny,Nz,currentx,currenty,currentz,n,p,phi,T,heating,maxheating,EFx,EFy,EFz,recombination)
	use constants
	implicit none
	
	integer, intent(in) :: Nx, Ny, Nz
	double precision, dimension(Nx,Ny,Nz), intent(in) :: currentx, currenty, currentz, n, p, phi, T, EFx, EFy,EFz,&
		recombination
	double precision, dimension(Nx,Ny,Nz) :: heating
	
	integer :: i,j,k,layer,ilayer,olayer,dlayer,ulayer,llayer,rlayer
	
	double precision :: ldist = 0, rdist = 0, udist = 0, ddist = 0, idist = 0, odist = 0, maxheating
	
	double precision :: hh,EExx,jjxx,EEyy,jjyy,EEzz,jjzz,pphi,&
						lT,rT,dT,uT,TT,iT,oT,QQxx,QQyy,QQzz, nn, pp 
	maxheating = 0.0d0
	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 1,Nz
	
				if ((threadpoints(2,2,1) .ne. -1) .and. (k == 1)) then
					cycle
				end if
				if ((threadpoints(2,2,3) .ne. -2) .and. (k == Nz)) then
					cycle
				end if
				pphi = phi(i,j,k)
				TT = T(i,j,k)
				nn = n(i,j,k)
				pp = p(i,j,k)
				hh = 0.0d0
				
				ldist = dabs(subgridx(i)-subgridx(i-1))
				rdist = dabs(subgridx(i)-subgridx(i+1))
				odist = dabs(subgridy(j)-subgridy(j-1))
				idist = dabs(subgridy(j)-subgridy(j+1))
				if (k .ne. 1) then
					ulayer = materials(i,j,k-1)
					uT = T(i,j,k-1)
					udist = dabs(subgridz(k)-subgridz(k-1))
				end if
				if (k .ne. Nz) then
					dlayer = materials(i,j,k+1)
					dT = T(i,j,k+1)
					ddist = dabs(subgridz(k)-subgridz(k+1))
				end if
				llayer = materials(i-1,j,k)
				lT = T(i-1,j,k)
				rlayer = materials(i+1,j,k)
				rT = T(i+1,j,k)
				olayer = materials(i,j-1,k)
				oT = T(i,j-1,k)
				ilayer = materials(i,j+1,k)
				iT = T(i,j+1,k)
				
				EExx = EFx(i,j,k)
				EEyy = EFy(i,j,k)
				EEzz = EFz(i,j,k)
				
				!Joule heating
				hh = hh + EExx*currentx(i,j,k) + EEyy*currenty(i,j,k) + EEzz*currentz(i,j,k)
				!Thomson heating
				QQxx = - (rT-lT)/(rdist+ldist)
				QQyy = - (iT-oT)/(idist+odist)
				if ((k > 1) .and. (k < Nz)) then
					QQzz = -(dT-uT)/(ddist+udist)
				elseif (k == 1) then
					QQzz = -(dT - TT)/(ddist)
				elseif (k == Nz) then
					QQzz = -(TT - uT)/(udist)
				end if
				hh = hh -q*dSdT(nn,TT,layer)*TT*(currentx(i,j,k)*QQxx+currenty(i,j,k)*QQyy+currentz(i,j,k)*QQzz)
				!RECOMBINATION HEATING
				hh = hh + 0.75 * (Erecombination)*recombination(i,j,k)
				heating(i,j,k) = hh
				if (hh > maxheating) then
					maxheating = hh
				end if
				
			end do
		end do
	end do
	
	
end subroutine GetHeating

subroutine HeatIteration(Nx,Ny,Nz,T,heating, deltat,enviromentalT, error,totalheating,totaloutflow,heatfluxx, heatfluxy, heatfluxz)
	use constants
	implicit none
	
	
	integer, intent(in) :: Nx, Ny, Nz
	double precision, dimension(Nx,Ny,Nz), intent(in) :: heating
	double precision, dimension(Nx,Ny,Nz), intent(inout) :: T
	double precision, dimension(Nx,Ny,Nz) :: newT,heatfluxx,heatfluxy,heatfluxz
	double precision :: deltat,rrho,ccp,enviromentalT,error
	
	integer :: k,i,j
	double precision :: nominator, denominator
	double precision :: uT,dT,rT,lT,iT,oT
	double precision :: udist,ddist,rdist,ldist,idist,odist
	double precision :: ukappa,dkappa,rkappa,lkappa,kkappa,ikappa,okappa
	double precision :: hheating,TT,nnewT,qx0,qx1,qy0,qy1,qz0,qz1,divqx,divqy,divqz
	double precision :: temp1,temp2,maxmax,maxk,maxi,maxj
	double precision :: totalheating,totaloutflow
	
	totalheating = 0.0d0
	totaloutflow = 0.0d0
	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 1,Nz
				
				if ((threadpoints(2,2,1) .ne. -1) .and. (k == 1)) then
					cycle
				end if
				
				if ((threadpoints(2,2,3) .ne. -2) .and. (k == Nz)) then
					cycle
				end if
				
				if (k .ne. 1) then
					ukappa = kappa(i,j,k-1)
					uT = T(i,j,k-1)
					udist = dabs(subgridz(k)-subgridz(k-1))
				end if
				
				if (k .ne. Nz) then
					dkappa = kappa(i,j,k+1)
					dT = T(i,j,k+1)
					ddist = dabs(subgridz(k)-subgridz(k+1))
				end if
				
				lkappa = kappa(i-1,j,k)
				lT = T(i-1,j,k)
				ldist = dabs(subgridx(i)-subgridx(i-1))
			
				rkappa = kappa(i+1,j,k)
				rT = T(i+1,j,k)
				rdist = dabs(subgridx(i)-subgridx(i+1))
							
				okappa = kappa(i,j-1,k)
				oT = T(i,j-1,k)
				odist = dabs(subgridy(j)-subgridy(j-1))
			
				ikappa = kappa(i,j+1,k)
				iT = T(i,j+1,k)
				idist = dabs(subgridy(j)-subgridy(j+1))
				
				hheating = heating(i,j,k)
				TT = T(i,j,k)
				kkappa = kappa(i,j,k)
				ccp = cp(i,j,k)
				rrho = rho(i,j,k)
				
				qx1 = -(kkappa+rkappa)*(rT-TT)/(2*rdist)
				qx0 = -(kkappa+lkappa)*(TT-lT)/(2*ldist)
				divqx = 2*(qx1-qx0)/(ldist+rdist)
				heatfluxx(i,j,k) = (qx0+qx1)/2
				
				qy1 = -(kkappa+ikappa)*(iT-TT)/(2*idist)
				qy0 = -(kkappa+okappa)*(TT-oT)/(2*odist)
				divqy = 2*(qy1-qy0)/(idist+odist)
				heatfluxy(i,j,k) = (qy0+qy1)/2
								
				if ((k .ne. 1) .and. (k .ne. Nz)) then
				
					qz1 = -(kkappa+dkappa)*(dT-TT)/(2*ddist)
					qz0 = -(kkappa+ukappa)*(TT-uT)/(2*udist)
					divqz = (qz1-qz0)/ddist
					heatfluxz(i,j,k) = (qy0+qy1)/2
					nnewT = TT + deltat*(-divqx - divqy - divqz + hheating)/(rrho*ccp)
! 					print*,deltat*(-divqx - divqy - divqz + hheating)/(rrho*ccp)
					
					heatfluxz(i,j,k) = (qz0+qz1)/2
					totalheating = totalheating + hheating*(idist+odist)*(rdist+ldist)*(ddist+udist)/8
					
				elseif (k == 1) then
					
					qz1 = -(kkappa+dkappa)*(dT-TT)/(2*ddist)
					qz0 = newtonh*(enviromentalT-TT)
					divqz = (qz1-qz0)/ddist
					nnewT = TT + deltat*(-divqx - divqy - divqz + hheating)/(rrho*ccp)
					heatfluxz(i,j,k) = (qz0+qz1)/2
					totalheating = totalheating + hheating*(idist+odist)*(rdist+ldist)*(ddist)/4
					totaloutflow = totaloutflow - qz0*(idist + odist)*(rdist + ldist)/4

				elseif (k == Nz) then
					
					qz1 = newtonh*(TT-enviromentalT)
					qz0 = -(kkappa+ukappa)*(TT-uT)/(2*udist)
					divqz = (qz1-qz0)/udist
					nnewT = TT + deltat*(-divqx - divqy - divqz + hheating)/(rrho*ccp)
					heatfluxz(i,j,k) = (qz0+qz1)/2
					totalheating = totalheating + hheating*(idist+odist)*(rdist+ldist)*(ddist)/4
					totaloutflow = totaloutflow + qz1*(idist + odist)*(rdist + ldist)/4
					
				end if
				newT(i,j,k) = nnewT
			end do
		end do
	end do
	
	error = 0.0d0
	
	do i = 2, Nx-1
		do j = 2, Ny-1
			do k = 1, Nz
				
				if ((threadpoints(2,2,1) .ne. -1) .and. (k == 1)) then
					cycle
				end if
				
				if ((threadpoints(2,2,3) .ne. -2) .and. (k == Nz)) then
					cycle
				end if
				
				TT = T(i,j,k)
				
				nnewT = newT(i,j,k)
				
				if (dabs(nnewT/TT-1.0d0) > 0.0d0) then
					error = dabs(nnewT/TT-1.0d0)
				end if
			
			end do
		end do
	end do
		
	T = newT
end subroutine HeatIteration
