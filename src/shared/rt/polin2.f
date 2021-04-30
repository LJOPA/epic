      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      REAL*8 dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      Parameter (NMAX=570000,MMAX=570000)
C     Uses polint
      INTEGER j,k
      REAL*8 ymtmp(MMAX),yntmp(NMAX)
!print*,'hello'
      do j=1,m
         do k=1,n
            yntmp(k)=ya(j,k)
         enddo
         call polint(x2a,yntmp,n,x2,ymtmp(j),dy)  
      enddo
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END
