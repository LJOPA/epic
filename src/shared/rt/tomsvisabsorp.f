*************This subroutine assumes that the solar flux is striking
*************the planet as a direct columnated beam. 	

        Subroutine tomsvisabsorp(nlay,tau,cosz,flux,dfdt,
     &   fluxbase,fluxinc2,s,delg)
	Implicit None
	integer nlay,i,s,j,dum
	real tau(10,nlay-1),cosz,flux,dfdt(nlay)
	real*8 etau,fluxinc(nlay),fluxinc2(nlay),airmass,
     &        fluxbase,delg(10)

	airmass=1./cosz
****Assume that 0 opacity above first pressure level

        dfdt(1)=0.
	fluxbase=0.

	  fluxinc(1)=delg(1)*flux/airmass
          fluxinc2(1)=fluxinc(1)
	  do i=1,nlay-1
	    etau=exp(-tau(1,i)*airmass)
            dfdt(i+1)=fluxinc(i)*(1.-etau)
	    fluxinc(i+1)=fluxinc(i)*etau
            fluxinc2(i+1)=fluxinc(i)
	  end do
	  fluxbase=fluxinc(nlay)+fluxbase

        do j=2,s
	  fluxinc(1)=delg(j)*flux/airmass
	  fluxinc2(1)=fluxinc(1)+fluxinc2(1)
	  do i=1,nlay-1
	    etau=exp(-tau(j,i)*airmass)
            dfdt(i+1)=dfdt(i+1)+fluxinc(i)*(1.-etau)
	    fluxinc(i+1)=fluxinc(i)*etau
            fluxinc2(i+1)=fluxinc(i)+fluxinc2(i+1)
	  end do
	  fluxbase=fluxinc(nlay)+fluxbase
        enddo

	return
	end 
