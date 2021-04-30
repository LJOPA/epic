        subroutine tominterp(x,y,nx,newx,newy,nnew) 

	implicit none

        integer i,j,nnew,nx
        real*8 x(nx),y(nx),newx(nnew),newy(nnew),mp,b

        do i=1,nnew
		if (newx(i).lt.x(1)) then
                        mp=(y(1)-y(1+1))/(log(x(1))-log(x(1+1)))
                        b=y(1)-log(x(1))*mp
                        newy(i)=log(newx(i))*mp+b
                else
                do j=1,nx-1
		if (newx(i).eq.x(j)) then
			newy(i)=y(j)
                else if (newx(i).gt.x(j).and.newx(i).lt.x(j+1)) then
                        mp=(y(j)-y(j+1))/(log(x(j))-log(x(j+1)))
                        b=y(j)-log(x(j))*mp
                        newy(i)=log(newx(i))*mp+b
                end if
                end do
		if (newx(i).gt.x(nx)) then
                        mp=(y(nx-1)-y(nx))/(log(x(nx-1))-log(x(nx)))
                        b=y(nx-1)-log(x(nx-1))*mp
                        newy(i)=log(newx(i))*mp+b
                end if
		end if
        end do
	return
	end

