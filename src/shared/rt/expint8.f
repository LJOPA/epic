


	function expint(n,x)
	integer n,maxit
	real*8 expint,x,eps,fpmin,euler
	parameter (maxit=100,eps=1.e-7,fpmin=1.e-30,
     &euler=0.5772156649)
	integer i,ii,nm1
	real*8 a,b,c,d,del,fact,h,psi
	nm1=n-1
	if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.
     &or.n.eq.1))) then
		print*,'bad arguments in expint'
	else if (n.eq.0) then
		expint=exp(-x)/x
	else if (x.eq.0.)then
		expint=1./nm1
	else if (x.gt.1.) then
		b=x+n
		c=1./fpmin
		d=1./b
		h=d
		do i=1,maxit
			a=-i*(nm1+i)
			b=b+2
			d=1./(a*d+b)
			c=b+a/c
			del=c*d
			h=h*del
			if (abs(del-1.).lt.eps) then
				expint=h*exp(-x)
				return
			endif
		end do
		print*,'continued fracion failed in expint'
	else
		if (nm1.ne.0) then
			expint=1./nm1
		else
			expint=-log(x)-Euler
		end if
		fact=1.
		do i=1,maxit
			fact=-fact*x/i
			if (i.ne.nm1) then
				del=-fact/(i-nm1)
			else
				psi=-euler
				do ii=1,nm1
					psi=psi+1./ii
				enddo
				del=fact*(-log(x)+psi)
			end if
			expint=expint+del
			if (abs(del).lt.abs(expint)*eps) return
		end do
		print*,'series failed in expint'
	end if
	return
	end 

