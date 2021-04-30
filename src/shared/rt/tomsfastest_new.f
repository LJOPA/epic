**********A subroutine to calculate D flux / D Pressure in an **********
***atmosphere.  It uses the Exponential integrals to do the ************
***angle averaging prior to the radiative transfer.  It is *************
***approximately 5 times faster than DISORT when the number*************
***of levels in the model is 75.****************************************
******formulation taken from Goukenleuque et al. 2000.*****************


	subroutine tomsfastdfdtau(nlay,pi2,
     &taut,extab3,extab4,wn,t,p,wnumlow,wnumhigh,
     &fluxtotdif,fluxbase,fluxtot,flux)

	Implicit None

	integer nlay,j,i,start,incs,ab1,ab2,abbase,z,zz
	real*8 taut(nlay),wn,fluxbase,
     &t(nlay),p(nlay),pi2,tauuse(0:nlay),ause(nlay),
     &fluxtot(nlay),flux(nlay),tau(0:nlay),
     &B,fluxtotdif(nlay-1),wnumlow,wnumhigh,
     &extab3(400000),extab4(400000)


	tau(0)=0.
	do j=1,nlay
	  tau(j)=taut(j)+tau(j-1)
        end do
        j=1
	start=j
	do while (tau(j).le.5.e-4.and.j.lt.nlay) 
          tau(j)=0.
          j=j+1
          start=j
        end do
	do j=start,nlay
          z=j-start+1
          zz=nlay+start-j
	  tauuse(z)=tau(zz)
	  ause(z)=B(wn,t(zz))   !a(nlay+start-j)
	end do
	incs=nlay-start+1
	do j=1,incs
	  flux(j)=0.
	  do i=1,incs-1
	    ab1=min(1+nint(abs(tauuse(i)-tauuse(j))*10000.),400000)
	    ab2=min(1+nint(abs(tauuse(i+1)-tauuse(j))*10000.),400000)
	    flux(j)=(ause(i+1)*extab3(ab2)-ause(i)*  !note
     &extab3(ab1)+(ause(i)-ause(i+1))*dabs((extab4(ab1)-
     &extab4(ab2))/(tauuse(i)-tauuse(i+1))))+flux(j)
	  end do
	  abbase=min(1+nint(abs(tauuse(1)-tauuse(j))*10000.),400000)
	  flux(j)=(flux(j)+ause(1)*extab3(abbase))*
     &(wnumhigh-wnumlow)*0.001*pi2 !flux (W/m^2)
	  fluxtot(incs+start-j)=flux(j)
	end do
	do i=1,start-1
	   fluxtot(i)=flux(incs)
	end do
!could pull out this division by delta p and do it after summing up all the flux totals
	do i=1,nlay-1
	fluxtotdif(i)=(fluxtot(i+1)-fluxtot(i))/(p(i+1)-p(i)) !DFlux/DP @ play layer
	end do
	fluxbase=fluxtot(nlay) 
	return
	end 
