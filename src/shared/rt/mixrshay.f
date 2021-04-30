*---------------------------------------------------------------------------

	subroutine mixr(nlay,ngas,ig,nr,ntp,pr,xmr,p,mix) 
	
** Inputs:  nlay, ngas:  number of layers, number of gases
*	    nr:     	 number of points defining a distribution
*	    pr: 	 pressure levels defining inflexion pts in distribution
*	    xmr:  	 mixing ratio at inflexion pts
*	    p:		 presures at nlay levels
** Output:  xmix: 	 mixing ratio at each pressure p

	Implicit none
	
	Integer nlay,nr,ig,ntp,j,i,ngas
	Real*8 pr(ntp),xmr(ntp),mix(ngas,nlay),p(nlay),b,sl

	
	j = 1 
	do i=1,nlay
 300	   If (p(i).lt.pr(1)) then
                mix(ig,i)=dexp(((log(xmr(1))-log(xmr(2)))/
     &(log(pr(1))-
     &log(pr(2))))*log(p(i))+log(xmr(1))-(log(pr(1))*
     &((log(xmr(1))-
     &log(xmr(2)))/(log(pr(1))-log(pr(2))))))
	   else if (p(i).gt.pr(nr)) then
                mix(ig,i)=dexp(((log(xmr(nr-1))-log(xmr(nr)))/
     &(log(pr(nr-1))-
     &log(pr(nr))))*log(p(i))+log(xmr(nr-1))-
     &(log(pr(nr-1))*(
     &(log(xmr(nr-1))-
     &log(xmr(nr)))/(log(pr(nr-1))-log(pr(nr))))))
	   else if (p(i).ge.pr(j).and.p(i).le.pr(j+1)) then
	     sl=(log(xmr(j+1))-log(xmr(j)))/
     .          (log(pr(j+1))-log(pr(j)))
             b=exp(log(xmr(j))-sl*log(pr(j)))
	     mix(ig,i)=exp((log(p(i)))*sl+log(b)) 
	   else
	     j = j+1
	     go to 300
	   End if 
	end do
	

	Return 
	End 	  
	  	
