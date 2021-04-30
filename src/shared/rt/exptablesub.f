

	subroutine exptable(extab3,extab4)

	Implicit None

	integer n,i
        real*8 x,expint,extab3(400000),extab4(400000)

	open(10,file='exp3.out',status='unknown')
        n=3
	do i=1,400000
        x=0.0001*(i-1)
        extab3(i)=expint(n,x)
	write(10,*) extab3(i),n,x
	end do
	close(10)

	open(10,file='exp4.out',status='unknown')
        n=4
	do i=1,400000
        x=(0.0001*(i-1))
        extab4(i)=expint(n,x)
	write(10,*) extab4(i),n,x
	end do
	close(10)

	return
	end
