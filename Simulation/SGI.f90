subroutine sgiteration(Nx,Ny,Nz,n,p,phi,ephi,hphi,T,currentx,currenty,currentz,deltatin,nnewmax,&
ntempcurrentx,ntempcurrenty,ntempcurrentz,ptempcurrentx,ptempcurrenty,ptempcurrentz,&
EFx,EFy,EFz,relativerecombination,recombinationperm3)
	
	use constants
	implicit none
	
	integer, intent(in) :: Nx,Ny,Nz
	double precision, dimension(Nx,Ny,Nz), intent(in) :: ephi, hphi, T,EFx,EFy,EFz,phi
	double precision, dimension(Nx,Ny,Nz) :: n,p,currentx,currenty,currentz,ntempcurrentx,&
	ntempcurrenty,ntempcurrentz,ptempcurrentx,ptempcurrenty,ptempcurrentz,nnew, pnew,&
	recombinationperm3,relativerecombination
	double precision, intent(in) :: deltatin
	double precision :: deltat
	
	integer :: i,j,k, ni,nj,nk,ni0,nj0,kmax,imax,jmax
	double precision :: iteration, DD,mmu, nm, Emf, Ef, Seebeck, nnewmax,vt, SS,Ebottom,Etop
		
	double precision :: uphi,dphi,lphi,rphi,ophi,iphi,pphi
	double precision :: udist,ddist,ldist,rdist,odist,idist
	double precision :: uephi,dephi,lephi,rephi,oephi,iephi,pephi
	double precision :: uhphi,dhphi,lhphi,rhphi,ohphi,ihphi,phphi
	double precision :: un,dn,ln,rn,on,in,nn
	double precision :: up,dp,lp,rp,op,ip,pp
	double precision :: recombination
	double precision :: uE,dE,rE,lE,oE,iE,EE
	double precision :: uT,dT,rT,lT,oT,iT,TT, tbc=1.0d21
	double precision :: PhiTop=0.5, PhiBottom=0.5, ymaxtop, ymaxbottom, umaxtop,umaxbottom
	double precision :: xcurrent0, xcurrent1, ycurrent0, ycurrent1, zcurrent0,zcurrent1
	double precision :: ppsigma=0, lpsigma=0, rpsigma=0, upsigma=0, dpsigma=0, opsigma=0,ipsigma=0
	double precision :: nnsigma=0, lnsigma=0, rnsigma=0, unsigma=0, dnsigma=0, onsigma=0,insigma=0
	integer :: layer,dlayer,ulayer,rlayer,llayer,olayer,ilayer
	double precision :: pycurrent0,pycurrent1,nycurrent0,nycurrent1
	
	open(333,FILE = "check")

	do i = 2,Nx-1
		do j = 2,Ny-1
			do k = 1,Nz
				
				if ((threadpoints(2,2,1) .ne. -1) .and. (k == 1)) then
					cycle
				end if
				
				if ((threadpoints(2,2,3) .ne. -2) .and. (k == Nz)) then
					cycle
				end if
				
				layer = materials(i,j,k)
				nn = n(i,j,k)
				pp = p(i,j,k)
				EE = dsqrt(EFx(i,j,k)**2+EFy(i,j,k)**2+EFz(i,j,k)**2)
				pephi = ephi(i,j,k)
				phphi = hphi(i,j,k)
				ppsigma = psigma(i,j,k)
				nnsigma = nsigma(i,j,k)
				TT = T(i,j,k)
				pphi = phi(i,j,k)
				
				if (k .ne. 1) then
					uhphi = hphi(i,j,k-1)
					uephi = ephi(i,j,k-1)
					ulayer = materials(i,j,k-1)
					uphi = phi(i,j,k-1)
					un = n(i,j,k-1)
					up = p(i,j,k-1)
					unsigma = nsigma(i,j,k-1)
					upsigma = psigma(i,j,k-1)
					uE = dsqrt(EFx(i,j,k-1)**2+EFy(i,j,k-1)**2+EFz(i,j,k-1)**2)
					uT = T(i,j,k-1)
					udist = dabs(subgridz(k)-subgridz(k-1))
				end if
					
				if (k .ne. Nz) then
					dhphi = hphi(i,j,k+1)
					dephi = ephi(i,j,k+1)
					dlayer = materials(i,j,k+1)
					dphi = phi(i,j,k+1)
					dn = n(i,j,k+1)
					dp = p(i,j,k+1)
					dnsigma = nsigma(i,j,k+1)
					dpsigma = psigma(i,j,k+1)
					dE = dsqrt(EFx(i,j,k+1)**2+EFy(i,j,k+1)**2+EFz(i,j,k+1)**2)
					dT = T(i,j,k+1)
					ddist = dabs(subgridz(k)-subgridz(k+1))
				end if
				
				if ( k==1 ) then
					udist = ddist
				elseif (k == Nz) then
					ddist = udist
				end if
				
				lhphi = hphi(i-1,j,k)
				lephi = ephi(i-1,j,k)
				llayer = materials(i-1,j,k)
				lphi = phi(i-1,j,k)
				ln = n(i-1,j,k)
				lp = p(i-1,j,k)
				lnsigma = nsigma(i-1,j,k)
				lpsigma = psigma(i-1,j,k)
				lE = dsqrt(EFx(i-1,j,k)**2+EFy(i-1,j,k)**2+EFz(i-1,j,k)**2)
				lT = T(i-1,j,k)
				ldist = dabs(subgridx(i)-subgridx(i-1))
					
				rhphi = hphi(i+1,j,k)
				rephi = ephi(i+1,j,k)
				rlayer = materials(i+1,j,k)
				rphi = phi(i+1,j,k)
				rn = n(i+1,j,k)
				rp = p(i+1,j,k)
				rnsigma = nsigma(i+1,j,k)
				rpsigma = psigma(i+1,j,k)
				rE = dsqrt(EFx(i+1,j,k)**2+EFy(i+1,j,k)**2+EFz(i+1,j,k)**2)
				rT = T(i+1,j,k)
				rdist = dabs(subgridx(i)-subgridx(i+1))
				
				ohphi = hphi(i,j-1,k)
				oephi = ephi(i,j-1,k)
				olayer = materials(i,j-1,k)
				ophi = phi(i,j-1,k)
				on = n(i,j-1,k)
				op = p(i,j-1,k)
				onsigma = nsigma(i,j-1,k)
				opsigma = psigma(i,j-1,k)
				oE = dsqrt(EFx(i,j-1,k)**2+EFy(i,j-1,k)**2+EFz(i,j-1,k)**2)
				oT = T(i,j-1,k)
				odist = dabs(subgridy(j)-subgridy(j-1))
				
				ihphi = hphi(i,j+1,k)
				iephi = ephi(i,j+1,k)
				ilayer = materials(i,j+1,k)
				iphi = phi(i,j+1,k)
				in = n(i,j+1,k)
				ip = p(i,j+1,k)
				insigma = nsigma(i,j+1,k)
				ipsigma = psigma(i,j+1,k)
				iE = dsqrt(EFx(i,j+1,k)**2+EFy(i,j+1,k)**2+EFz(i,j+1,k)**2)
				iT = T(i,j+1,k)
				idist = dabs(subgridy(j)-subgridy(j+1))
				Etop = EFz(i,j,k)
				Ebottom = EFz(i,j,k)
!_______________UPDATING SIGMA TO CREATE NEW POTENTIAL____________________________________________
! 				print*, "nnsigma",nnsigma
				
				nnsigma = -(2*nnsigma**2*kt*TT)/q
				ppsigma = -(2*ppsigma**2*kt*TT)/q
				if (lT .ne. 0) then
					lnsigma = -(2*lnsigma**2*kt*lT)/q
					lpsigma = -(2*lpsigma**2*kt*lT)/q
				else
					lnsigma = nnsigma
					lpsigma = ppsigma
				end if
				if (rT .ne. 0) then
					rnsigma = -(2*rnsigma**2*kt*rT)/q
					rpsigma = -(2*rpsigma**2*kt*rT)/q
				else
					rnsigma = nnsigma
					rpsigma = ppsigma
				end if
				if (dT .ne. 0) then
					dnsigma = -(2*dnsigma**2*kt*dT)/q
					dpsigma = -(2*dpsigma**2*kt*dT)/q
				else
					dnsigma = nnsigma
					dpsigma = ppsigma
				end if
				if (uT .ne. 0) then
					unsigma = -(2*unsigma**2*kt*uT)/q
					upsigma = -(2*upsigma**2*kt*uT)/q
				else
					unsigma = nnsigma
					upsigma = ppsigma
				end if
				if (oT .ne. 0) then
					onsigma = -(2*onsigma**2*kt*oT)/q
					opsigma = -(2*opsigma**2*kt*oT)/q
				else
					onsigma = nnsigma
					opsigma = ppsigma
				end if
				
				if (iT .ne. 0) then
					insigma = -(2*insigma**2*kt*iT)/q
					ipsigma = -(2*ipsigma**2*kt*iT)/q
				else
					insigma = nnsigma
					ipsigma = ppsigma
				end if

!_________________________________________________________________________________________________

				deltat = deltatin
				recombination = GR(nn,pp,TT,EE,layer,epsR(i,j,k))
				recombinationperm3(i,j,k) = recombination
				relativerecombination(i,j,k) = recombination/(nn+pp)
! 				recombinationpernm3(i,j,k) = recombination
				if ((k .ne. 1) .and. (k .ne. Nz)) then
					SS = (S(nn,TT,layer) + S(rn,rT,rlayer))/2
					DD = (Dc(nn,TT,EFx(i,j,k),-1,layer)+Dc(rn,rT,EFx(i+1,j,k),-1,rlayer))/2
					xcurrent1 = -(DD/rdist)*(rn*B(q*(rephi-pephi + SS*(rT-TT)+(rnsigma-nnsigma))/(kt*(TT+rT)/2))-&
					nn*B(-q*(rephi-pephi+SS*(rT-TT)+(rnsigma-nnsigma))/(kt*(TT+rT)/2)))
					
					SS = (S(nn,TT,layer) + S(ln,lT,llayer))/2
					DD = (Dc(nn,TT,EFx(i,j,k),-1,layer)+Dc(ln,lT,EFx(i-1,j,k),-1,llayer))/2
					xcurrent0 = -(DD/ldist)*(nn*B(q*(pephi-lephi+SS*(TT-lT)+(nnsigma-lnsigma))/(kt*(TT+lT)/2))-&
					ln*B(-q*(pephi-lephi+SS*(TT-lT)+(nnsigma-lnsigma))/(kt*(TT+lT)/2)))
					
					SS = (S(nn,TT,layer) + S(dn,dT,dlayer))/2
					DD = (Dc(nn,TT,EFz(i,j,k),-1,layer)+Dc(dn,dT,EFz(i,j,k+1),-1,dlayer))/2
					zcurrent1 = -(DD/ddist)*(dn*B(q*(dephi-pephi+SS*(dT-TT)+(dnsigma-nnsigma))/(kt*(TT+dT)/2))-&
					nn*B(-q*(dephi-pephi+SS*(dT-TT)+(dnsigma-nnsigma))/(kt*(TT+dT)/2)))
					
					SS = (S(nn,TT,layer) + S(un,uT,ulayer))/2
					DD = (Dc(nn,TT,EFz(i,j,k),-1,layer)+Dc(un,uT,EFz(i,j,k-1),-1,ulayer))/2
					zcurrent0 = -(DD/udist)*(nn*B(q*(pephi-uephi+SS*(TT-uT)+(nnsigma-unsigma))/(kt*(TT+uT)/2))-&
					un*B(-q*(pephi-uephi+SS*(TT-uT)+(nnsigma-unsigma))/(kt*(TT+uT)/2)))
					
					SS = (S(nn,TT,layer) + S(in,iT,ilayer))/2
					DD = (Dc(nn,TT,EFy(i,j,k),-1,layer)+Dc(in,iT,EFy(i,j+1,k),-1,ilayer))/2
					ycurrent1 = -(DD/idist)*(in*B(q*(iephi-pephi + SS*(iT-TT)+(insigma-nnsigma))/(kt*(TT+iT)/2))-&
					nn*B(-q*(iephi-pephi+SS*(iT-TT)+(insigma-nnsigma))/(kt*(TT+iT)/2)))
					
					SS = (S(nn,TT,layer) + S(on,oT,olayer))/2
					DD = (Dc(nn,TT,EFy(i,j,k),-1,layer)+Dc(on,oT,EFy(i,j-1,k),-1,olayer))/2
					ycurrent0 = -(DD/odist)*(nn*B(q*(pephi-oephi+SS*(TT-oT)+(nnsigma-onsigma))/(kt*(TT+oT)/2))-&
					on*B(-q*(pephi-oephi+SS*(TT-oT)+(nnsigma-onsigma))/(kt*(TT+oT)/2)))
					
					nnew(i,j,k) = nn - deltat*(recombination + 2*(xcurrent1-xcurrent0)/(ldist+rdist)+&
					2*(ycurrent1-ycurrent0)/(idist+odist)+2*(zcurrent1-zcurrent0)/(udist+ddist))
					ntempcurrentx(i,j,k) = (xcurrent0+xcurrent1)/2
					ntempcurrenty(i,j,k) = (ycurrent0+ycurrent1)/2
					ntempcurrentz(i,j,k) = (zcurrent0+zcurrent1)/2
					
				elseif (k == 1) then	
					nnew(i,j,k) = tbc
					
				elseif (k == Nz) then
					
					DD = (Dc(nn,TT,EFx(i,j,k),-1,layer)+Dc(rn,rT,EFx(i+1,j,k),-1,rlayer))/2
					SS = (S(nn,TT,layer) + S(rn,rT,rlayer))/2
					xcurrent1 = -(DD/rdist)*(rn*B(q*(rephi-pephi + SS*(rT-TT)+(rnsigma-nnsigma))/(kt*(TT+rT)/2))-&
					nn*B(-q*(rephi-pephi+SS*(rT-TT)+(rnsigma-nnsigma))/(kt*(TT+rT)/2)))
					
					SS = (S(nn,TT,layer) + S(ln,lT,llayer))/2
					DD = (Dc(nn,TT,EFx(i,j,k),-1,layer)+Dc(ln,lT,EFx(i-1,j,k),-1,llayer))/2
					xcurrent0 = -(DD/ldist)*(nn*B(q*(pephi-lephi+SS*(TT-lT)+(nnsigma-lnsigma))/(kt*(TT+lT)/2))-&
					ln*B(-q*(pephi-lephi+SS*(TT-lT)+(nnsigma-lnsigma))/(kt*(TT+lT)/2)))

					DD = (Dc(nn,TT,EFy(i,j,k),-1,layer)+Dc(in,iT,EFy(i,j+1,k),-1,ilayer))/2
					SS = (S(nn,TT,layer) + S(in,iT,ilayer))/2
					ycurrent1 = -(DD/idist)*(in*B(q*(iephi-pephi + SS*(iT-TT)+(insigma-nnsigma))/(kt*(TT+iT)/2))-&
					nn*B(-q*(iephi-pephi+SS*(iT-TT)+(insigma-nnsigma))/(kt*(TT+iT)/2)))
					
					SS = (S(nn,TT,layer) + S(on,oT,olayer))/2
					DD = (Dc(nn,TT,EFy(i,j,k),-1,layer)+Dc(on,oT,EFy(i,j-1,k),-1,olayer))/2
					ycurrent0 = -(DD/odist)*(nn*B(q*(pephi-oephi+SS*(TT-oT)+(nnsigma-onsigma))/(kt*(TT+oT)/2))-&
					on*B(-q*(pephi-oephi+SS*(TT-oT)+(nnsigma-onsigma))/(kt*(TT+oT)/2)))
					
					SS = (S(nn,TT,layer) + S(un,uT,ulayer))/2
					DD = (Dc(nn,TT,EFz(i,j,k),-1,layer)+Dc(un,uT,EFz(i,j,k-1),-1,ulayer))/2
					zcurrent0 = -(DD/udist)*(nn*B(q*(pephi-uephi+SS*(TT-uT)+(nnsigma-unsigma))/(kt*(TT+uT)/2))-&
					un*B(-q*(pephi-uephi+SS*(TT-uT)+(nnsigma-unsigma))/(kt*(TT+uT)/2)))
					
					Ebottom = EFz(i,j,k)
					ymaxbottom = (q/(Ebottom*16*pi*eps*epsR(i,j,k)))
					if (Ebottom .ge. 0) then
						umaxbottom = q*phiBottom - sqrt(q**3*Ebottom/(4*pi*eps*epsR(i,j,k)))
					else
						umaxbottom = q*phiBottom
					end if
					zcurrent1 = (-jforward(TT,Umaxbottom,q*phibottom,Ebottom,eps*epsR(i,j,k))+&
					jbackflow(TT,nn,q*phibottom,-1) +&
					jbackdrift(Ebottom,TT,nn,-1,1))/q
					
					nnew(i,j,k) = nn - deltat*(recombination + 2*(xcurrent1-xcurrent0)/(ldist+rdist)+&
					2*(ycurrent1-ycurrent0)/(idist+odist)+2*(zcurrent1-zcurrent0)/(udist+ddist))
					ntempcurrentx(i,j,k) = (xcurrent0+xcurrent1)/2
					ntempcurrenty(i,j,k) = (ycurrent0+ycurrent1)/2
					ntempcurrentz(i,j,k) = (zcurrent0+zcurrent1)/2
					
				end if
					
				if ((k .ne. 1) .and. (k .ne. Nz)) then
					SS = (S(pp,TT,layer) + S(rp,rT,rlayer))/2
					DD = (Dc(pp,TT,EFx(i,j,k),1,layer)+Dc(rp,rT,EFx(i+1,j,k),1,rlayer))/2
					xcurrent1 = -(DD/rdist)*(rp*B(-q*(rhphi-phphi+SS*(rT-TT)-rpsigma+ppsigma)/(kt*(TT+rT)/2))-&
					pp*B(q*(rhphi-phphi+ SS*(rT-TT)-rpsigma+ppsigma)/(kt*(TT+rT)/2)))
					
					SS = (S(pp,TT,layer) + S(lp,lT,llayer))/2
					DD = (Dc(pp,TT,EFx(i,j,k),1,layer)+Dc(lp,lT,EFx(i-1,j,k),1,llayer))/2
					xcurrent0 = -(DD/ldist)*(pp*B(-q*(phphi-lhphi+SS*(TT-lT)-ppsigma+lpsigma)/(kt*(TT+lT)/2))-&	
					lp*B(q*(phphi-lhphi+SS*(TT-lT)-ppsigma+lpsigma)/(kt*(TT+lT)/2)))
					
					SS = (S(pp,TT,layer) + S(ip,iT,ilayer))/2
					DD = (Dc(pp,TT,EFy(i,j,k),1,layer)+Dc(ip,iT,EFy(i,j+1,k),1,ilayer))/2
					ycurrent1 = -(DD/idist)*(ip*B(-q*(ihphi-phphi+SS*(iT-TT)-ipsigma+ppsigma)/(kt*(TT+iT)/2))-&
					pp*B(q*(ihphi-phphi+ SS*(iT-TT)-ipsigma+ppsigma)/(kt*(TT+iT)/2)))
					
					SS = (S(pp,TT,layer) + S(op,oT,olayer))/2
					DD = (Dc(pp,TT,EFy(i,j,k),1,layer)+Dc(op,oT,EFy(i,j-1,k),1,olayer))/2
					ycurrent0 = -(DD/odist)*(pp*B(-q*(phphi-ohphi+SS*(TT-oT)-ppsigma+opsigma)/(kt*(TT+oT)/2))-&	
					op*B(q*(phphi-ohphi+SS*(TT-oT)-ppsigma+opsigma)/(kt*(TT+oT)/2)))
					
					SS = (S(pp,TT,layer) + S(dp,dT,dlayer))/2
					DD = (Dc(pp,TT,EFz(i,j,k),1,layer)+Dc(dp,dT,EFz(i,j,k+1),1,dlayer))/2
					zcurrent1 = -(DD/ddist)*(dp*B(-q*(dhphi-phphi+SS*(dT-TT)-dpsigma+ppsigma)/(kt*(TT+dT)/2))-&
					pp*B(q*(dhphi-phphi+SS*(dT-TT)-dpsigma+ppsigma)/(kt*(TT+dT)/2)))
					
					SS = (S(pp,TT,layer) + S(up,uT,ulayer))/2
					DD = (Dc(pp,TT,EFz(i,j,k),1,layer)+Dc(up,uT,EFz(i,j,k-1),1,ulayer))/2
					zcurrent0 = -(DD/udist)*(pp*B(-q*(phphi-uhphi+SS*(TT-uT)-ppsigma+upsigma)/(kt*(TT+uT)/2))-&
					up*B(q*(phphi-uhphi+SS*(TT-uT)-ppsigma+upsigma)/(kt*(TT+uT)/2)))
					
					pnew(i,j,k) = pp - deltat*(recombination + 2*(xcurrent1-xcurrent0)/(ldist+rdist)+&
					2*(ycurrent1-ycurrent0)/(odist+idist) + 2*(zcurrent1-zcurrent0)/(udist+ddist))
					ptempcurrentx(i,j,k) = (xcurrent0+xcurrent1)/2
					ptempcurrenty(i,j,k) = (ycurrent0+ycurrent1)/2
					ptempcurrentz(i,j,k) = (zcurrent0+zcurrent1)/2
					
				elseif (k == Nz) then
					
					pnew(i,j,k) = tbc
					
				elseif (k == 1) then	
					SS = (S(pp,TT,layer) + S(rp,rT,rlayer))/2
					DD = (Dc(pp,TT,EFx(i,j,k),1,layer)+Dc(rp,rT,EFx(i+1,j,k),1,rlayer))/2
					xcurrent1 = -(DD/rdist)*(rp*B(-q*(rhphi-phphi+SS*(rT-TT)-rpsigma+ppsigma)/(kt*(TT+rT)/2))-&
					pp*B(q*(rhphi-phphi+ SS*(rT-TT)-rpsigma+ppsigma)/(kt*(TT+rT)/2)))
					
					SS = (S(pp,TT,layer) + S(lp,lT,llayer))/2
					DD = (Dc(pp,TT,EFx(i,j,k),1,layer)+Dc(lp,lT,EFx(i-1,j,k),1,llayer))/2
					xcurrent0 = -(DD/ldist)*(pp*B(-q*(phphi-lhphi+SS*(TT-lT)-ppsigma+lpsigma)/(kt*(TT+lT)/2))-&
					lp*B(q*(phphi-lhphi+SS*(TT-lT)-ppsigma+lpsigma)/(kt*(TT+lT)/2)))
					
					SS = (S(pp,TT,layer) + S(ip,iT,ilayer))/2
					DD = (Dc(pp,TT,EFy(i,j,k),1,layer)+Dc(ip,iT,EFy(i,j+1,k),1,ilayer))/2
					ycurrent1 = -(DD/idist)*(ip*B(-q*(ihphi-phphi+SS*(iT-TT)-ipsigma+ppsigma)/(kt*(TT+iT)/2))-&
					pp*B(q*(ihphi-phphi+ SS*(iT-TT)-ipsigma+ppsigma)/(kt*(TT+iT)/2)))
					
					SS = (S(pp,TT,layer) + S(op,oT,olayer))/2
					DD = (Dc(pp,TT,EFy(i,j,k),1,layer)+Dc(op,oT,EFy(i,j-1,k),1,olayer))/2
					ycurrent0 = -(DD/odist)*(pp*B(-q*(phphi-ohphi+SS*(TT-oT)-ppsigma+opsigma)/(kt*(TT+oT)/2))-&	
					op*B(q*(phphi-ohphi+SS*(TT-oT)-ppsigma+opsigma)/(kt*(TT+oT)/2)))
					
					SS = (S(pp,TT,layer) + S(dp,dT,dlayer))/2
					DD = (Dc(pp,TT,EFz(i,j,k),1,layer)+Dc(dp,dT,EFz(i,j,k+1),1,dlayer))/2
					zcurrent1 = -(DD/ddist)*(dp*B(-q*(dhphi-phphi+SS*(dT-TT)-dpsigma+ppsigma)/(kt*(TT+dT)/2))-&
					pp*B(q*(dhphi-phphi+SS*(dT-TT)-dpsigma+ppsigma)/(kt*(TT+dT)/2)))
					
					Etop = EFz(i,j,k)
					ymaxtop = (q/(Etop*16*pi*eps*epsR(i,j,k)))
					
					if (Etop .ge. 0) then
						umaxtop = q*phiTop - sqrt(q**3*Etop/(4*pi*eps*epsR(i,j,k)))
					else
						umaxtop = q*phiTop
					end if
					
					zcurrent0 = (jforward(TT,Umaxtop,q*phitop,Etop,eps*epsR(i,j,k))& 
					-jbackflow(TT,pp,q*phitop,+1)&
					-jbackdrift(Etop,TT,pp,+1,1))/q
					
					pnew(i,j,k) = pp - deltat*(recombination + 2*(xcurrent1-xcurrent0)/(ldist+rdist)+&
					2*(ycurrent1-ycurrent0)/(odist+idist) + 2*(zcurrent1-zcurrent0)/(udist+ddist))
					ptempcurrentx(i,j,k) = (xcurrent0+xcurrent1)/2
					ptempcurrenty(i,j,k) = (ycurrent0+ycurrent1)/2
					ptempcurrentz(i,j,k) = (zcurrent0+zcurrent1)/2
				end if
					
			end do	
		end do		
	end do			
	
	nnewmax = 0.0d0
	
	do i = 2, Nx-1
		do j = 2, Ny-1
			do k = 2, Nz-1
				if (dabs(nnew(i,j,k)/n(i,j,k)-1.0d0) > nnewmax) then
					nnewmax = dabs(nnew(i,j,k)/n(i,j,k)-1.0d0)
				end if
			end do
		end do
	end do
	
	n = nnew
	p = pnew
	
	currentx = (ptempcurrentx-ntempcurrentx)*q
	currenty = (ptempcurrenty-ntempcurrenty)*q
	currentz = (ptempcurrentz-ntempcurrentz)*q
	if (threadpoints(2,2,1) == -1) then
		currentx(:,:,1) = currentx(:,:,2)
		currenty(:,:,1) = currenty(:,:,2)
		currentz(:,:,1) = currentz(:,:,2)
	end if
	
	if (threadpoints(2,2,3) == -2) then
		currentx(:,:,Nz) = currentx(:,:,Nz-1)
		currenty(:,:,Nz) = currenty(:,:,Nz-1)
		currentz(:,:,Nz) = currentz(:,:,Nz-1)
	end if
	
	contains
	
	double precision function GR(n,p,T,E,material,dielconst)
		use constants
		double precision, intent(in) :: n,p,T,E,dielconst
		double precision :: gamma, n0 = 0.0d10, p0 = 0.0d10
		integer, intent(in) :: material
		gamma = q*max(mu(n,T,E,-1,material),mu(p,T,E,+1,material))/(eps*dielconst)
		
		GR = gamma * (n * p - n0*p0)

	end function GR

	double precision function jbackdrift(E,T,p,charge,material)
		use constants
		double precision, intent(in) :: E,p,T
		integer, intent(in) :: charge,material
		
		if (E > 0) then
			jbackdrift = 0.0d0
			return
		end if
		
		jbackdrift = q*dabs(E)*p*mu(p,T,E,charge,material)
		
	end function jbackdrift

	double precision function jbackflow(T,p,Ain,charge)
		use constants
		double precision, intent(in) :: T,p,Ain
		integer, intent(in) :: charge
		double precision :: A,n
		
		n=1.0d20
		A=1.20173d6
		
		jbackflow = A*T**2*dexp(-Ain/(kt*T))*p/n
! 		print*, "BACK",jbackflow
	end function jbackflow

	double precision function jforward(T,Umax,Ain,E,epstot)
		double precision, intent(in) :: T,Umax,Ain,E,epstot
		if (E .le. 0) then
			jforward = 0
		else
			jforward = jthermal(T,Umax,E)+jtunneling(Ain,E,epstot)
		end if
			
	end function jforward

	double precision function jthermal(T,Umax,E)
		use constants
		double precision, intent(in) :: T,Umax,E
		double precision :: A
		A=1.20173d6
		if (E .ge. 0) then
			jthermal = A*T**2*dexp(-Umax/(kt*T))
		else
			jthermal = 0.0d0
		end if
		
		
! 		print*,"THERMAL INJECTION:", jthermal
	end function jthermal

	double precision function jtunneling(Ain,E,epstot)
		use constants
		double precision, intent(in) :: Ain,E,epstot
		double precision :: h,hh,v,t,Jt,Et

		if (E .le. 0) then
			jtunneling = 0
			return
		end if
		
		if (task_ID == 1) then
! 			print*, "here",E,Ain
		end if
		
		h = q*dsqrt(q*E/(4*pi*epstot))/Ain

		if (task_ID == 1) then
! 			print*, "here",h
		end if
		
		hh = dsqrt((1.0d0-h)/(1.0d0+h))
		if (hh < 0.0d0 .or. hh > 1.0d0) then
			print*, "something wrong with h in tunneling current, please check again"
			jtunneling = 0
			read*,
			return
		end if
		if (task_ID == 1) then
! 			print*, "here",1+h,hh
		end if
		v = dsqrt(1+h)*(EllipticE(hh) - h*EllipticK(hh))
		if (task_ID == 1) then
! 			print*, "here"
		end if
		t = ((1+h)*EllipticE(hh) - h * EllipticK(hh))/dsqrt(1+h)
		if (task_ID == 1) then
! 			print*, "here"
		end if
		
		
		Et = (4 * Ain**(1.5d0))/(3*q*BorRadius*dsqrt(Rydberg))
		Jt = q*Ain**2/(9*pi**2*BorRadius**2*hbar*Rydberg)
		jtunneling = Jt*(E/(Et*t))**2*dexp(-(Et*v)/E)
		
	end function jtunneling

	double precision function nbottomconcentration()
		nbottomconcentration = 1.0d10
	end function nbottomconcentration

	double precision function ptopconcentration()
		ptopconcentration = 1.0d10
	end function ptopconcentration

	double precision function ntopcurrent()
		ntopcurrent = 0.0d0
	end function ntopcurrent

	double precision function pbottomcurrent()
		pbottomcurrent = 0.0d0
	end function pbottomcurrent

pure	double precision function B(x)
		double precision, intent(in) :: x
		if ((abs(x)>=1.0d-3) .and. (abs(x)<=41)) then
			B = -x/(1-exp(x))	
		elseif(abs(x)<1.0d-3) then
			B = 1 - (x/2)
		elseif(x>41) then
			B=0.0d0
		elseif(x<-41) then
			B=-x
		end if
	end function B

double precision function ellipticE(x)
		double precision, intent(in) :: x
		double precision :: a, g,o, anew, gnew,onew, e
		a = (1+dsqrt(1-x**2))/2
		g = dsqrt(dsqrt(1-x**2))
		e = 1.0d-15
		if (x == 1.0d0) then
			EllipticE = 1.0d0
			return
		end if
		if (x < 0.0d0 .or. x > 1.0d0) then
			print*, "x < 0.0d0 or x > 1.0d0, PRESS ANY KEY TO EXIT, x:", x
			read*,
			return
		end if
		
		do
			
			anew = (a+g)/2
			gnew = dsqrt(a*g)
			if (anew-gnew < e) then
				exit
			end if
			a = anew
			g = gnew
		end do
		
		ellipticE = 3.1415926/(2*anew)
		
		
		a = 1-x**2
		g = 1
		o = 0
		do
			anew = (a + g)/2
			gnew = o + dsqrt(((a-o)*(g-o)))
			onew = o - dsqrt(((a-o)*(g-o)))
			if (anew-gnew < e) then
				exit
			end if
			a = anew
			g = gnew
			o = onew
		end do
		ellipticE = ellipticE*anew

	end function ellipticE


	double precision function ellipticK(x)
		double precision, intent(in) :: x
		double precision :: a, g, anew, gnew, e
		if (x < 0.0d0 .or. x .ge. 1.0d0) then
			print*, "x < 0.0d0 or x >= 1.0d0, PRESS ANY KEY TO EXIT,x :", x
			read*,
			return
		end if
		
		a = (1+dsqrt(1-x**2))/2
		g = dsqrt(dsqrt(1-x**2))
		e = 1.0d-15
		
		do
			anew = (a+g)/2
			gnew = dsqrt(a*g)
			if (anew-gnew < e) then
				exit
			end if
			a = anew
			g = gnew
		end do
		
		ellipticK = 3.1415926/(2*anew)
	end function ellipticK
	
end subroutine sgiteration
	
