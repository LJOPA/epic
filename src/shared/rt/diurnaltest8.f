       subroutine diurnaltest(qtr,stgkji,stgdayi,
     &           daily_f,qtr_use,maxcly,nlay,latstr,workdir,DF,
     &           fluxbaseday,fluxbase)

       implicit none

       integer guess,nlay,qtr,maxcly
       Parameter (guess=88)

       real*8 DF(maxcly,qtr),slope,qtr_use(qtr),pi,pi180,
     &   b,y(maxcly,guess),daily_f(maxcly),
     &   qtr_guess(guess),deltax,sumy(maxcly),dum,
     &   fluxbaseday(qtr),fluxbase,yy(guess)

       integer m,n,lk,l,ll,k
       character*34 workdir
       character*3 stgkji,latstr
       character*5 stgdayi

      pi = acos(-1.)
      pi180 = pi/180.

      do ll=1,guess
        qtr_guess(ll)=0.
      enddo

      qtr_guess(1)=qtr_use(1)

      do ll=2,guess
        qtr_guess(ll)=qtr_guess(ll-1)-(qtr_use(1)-qtr_use(qtr))/(guess)
      enddo

      m=guess-1
      n=m-1
      open(12,file=workdir//'yfile'//stgkji//'.out',
     &         form='unformatted',
     &         access='direct',recl=4*guess)
      do l=1,maxcly
        do lk=1,guess
  	  dum=qtr_guess(lk)
          call locate(qtr_use,qtr,dum,k)
          if(k.ge.10.or.k.le.0)then
            print*,'BIG WARNING: problem with cosine values'
          endif
          slope=(DF(l,k)-DF(l,k+1))/
     &    (qtr_use(k)-qtr_use(k+1))
          b=DF(l,k)-slope*qtr_use(k)
          y(l,lk)=slope*qtr_guess(lk)+b !unit DF
        end do
        write(12,rec=l) (real(y(l,lk)),lk=1,guess)
      enddo
      close(12)
          
      deltax=(qtr_guess(1)-qtr_guess(guess))/m       !h=(b-a)/n
      deltax=deltax*106.6          !106.6=10.66*3600./360. Converts angle into time in seconds.

      print*,'deltax',deltax

      open(15,file=workdir//'daily_f2.out',form='unformatted',
     &     access='direct',recl=4*maxcly)
            
      do l=1,maxcly
        sumy(l)=0.
        do lk=1,guess-1
          sumy(l)=(y(l,lk)+y(l,lk+1))/2.+sumy(l)
        enddo
        sumy(l)=sumy(l)*deltax
        daily_f(l)=2.*sumy(l)
      enddo     
      write(15,rec=1) (real(daily_f(l)),l=1,maxcly)
      close(15)

      do lk=1,guess
        dum=qtr_guess(lk)
        call locate(qtr_use,qtr,dum,k)
        if (k.ge.10.or.k.le.0) then
          print*,'BIG WARNING: problem with cosine values'
        endif
        slope=(fluxbaseday(k)-fluxbaseday(k+1))/
     &  (qtr_use(k)-qtr_use(k+1))
        b=fluxbaseday(k)-slope*qtr_use(k)
        yy(lk)=slope*qtr_guess(lk)+b !unit DF
      end do
      sumy(1)=0.
      do lk=1,guess-1
        sumy(1)=(yy(lk)+yy(lk+1))/2.+sumy(1)
      enddo
      sumy(1)=sumy(1)*deltax
      fluxbase=2.*sumy(1)

      return
      END

