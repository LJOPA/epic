       subroutine brunodiurn(lat,re,rp,plandec,sdis,cos_psi,
     &   lsubs,avecos_day,avecos)
       implicit none
       integer l,i
       real*8 lat,lsubs,pi,pi180,re,rp,plandec,theta_s,
     &  avecos_day,avecos,
     &   sdis,theta,phi,phi_s,phi_e,cos_psi,d,q,test
       real fbeam
cc 'avecos_day' is average cosine of the zenith angle integrated
cc        from sunrise to sunset
cc 'Q' is the mean daily solar insolation in units o W/m^2/day

! -26.7 .le. theta_s .le. 26.7
      pi = acos(-1.)
      pi180 = pi/180.
          theta_s=plandec*pi180  !rad
	  phi_s=lsubs*pi180  !rad
      d=sdis   !AU
      theta=lat*pi180  !rad
      phi=0.  !let longitude be 0.
      test=-((re/rp)**2.)*tan(theta)*tan(theta_s)
       if(test .ge. 1.)then
           phi_e=0.
       elseif(test .gt. -1. .and.
     &        test .lt. 1.) then
           phi_e=acos(-((re/rp)**2.)*tan(theta)*tan(theta_s))
       elseif(test .le. -1.) then
           phi_e=pi
       endif
      cos_psi= (cos(theta)*cos(theta_s)*cos(phi-phi_s)+
     &       ((re/rp)**2.)*sin(theta)*sin(theta_s))/
     &       ((cos(theta))**2.+((re/rp)**4.)*(sin(theta))**2.)**(0.5)
      avecos_day= (cos(theta)*cos(theta_s)*sin(phi_e)+
     &      (re/rp)**2.*phi_e*sin(theta)*sin(theta_s))/
     & (pi*((cos(theta))**2.+((re/rp)**4.)*(sin(theta))**2.)**(0.5))
      avecos=avecos_day*pi/phi_e
c      Q=(fbeam)*(avecos_day)
      end
