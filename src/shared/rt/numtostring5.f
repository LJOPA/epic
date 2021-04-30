

	subroutine numtostring5(x,xx)

	Implicit None
	integer x
	character*1 z1,xs1
	character*2 z2,xs2
	character*3 z3,xs3
	character*4 z4,xs4
	character*5 xx

	z1='0'
	z2='00'
	z3='000'
	z4='0000'
	if (x.le.9) then
	  write(xs1,'(I1)'),x
	  xx=z4//xs1
	elseif (x.lt.100) then
	  write(xs2,'(I2)'),x
	  xx=z3//xs2
	elseif (x.lt.1000) then
	  write(xs3,'(I3)'),x
	  xx=z2//xs3
	elseif (x.lt.10000) then
	  write(xs4,'(I4)'),x
	  xx=z1//xs4
	else
	  write(xx,'(I5)'),x
	endif 
	
	return
	end
