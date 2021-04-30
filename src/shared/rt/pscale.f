        subroutine pscale(p1,scale,p,nlay) 

        integer i
        real*8 p(nlay),scale,p1

!  Set up the pressure grid 

        p(1)  = p1
        do i=2,nlay
          p(i) = p(i-1)*exp(scale) 
        end do 
        return
        end 
