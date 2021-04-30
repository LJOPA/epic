        Subroutine irwin_rewrite(pnum,tnum,info,delg,pbar,t,wnirw,
     &wnirwl,wnirwh)

        implicit none

        integer pnum,tnum,il,n,m,j,i,l,ll
        character*70 dum
        REAL*8 err,delg(10),pbar(pnum),t(tnum),
     &info(10,pnum,tnum,1501),wnirw(1501),wnirwh(1501),wnirwl(1501)
c  This subroutine opens the Irwin data file and extracts the information.
c The Irwin file consists of blocks of temperatures and for each temp(19),
c there is 20 pressures each with 10 k-coeff values. The statistical 
c weights(10) remain constant for each and every set of 10 k-coeff and are 
c used in conjunction with the k-coeff. They weight the "contribution" of
c each of the 10 k-coeff relative to the total transmission at each pressure.

        open (unit=7,file='other/irwinH_iii.par',status='old')
        do l=1,15
           read(7,*) dum 
        enddo
        do ll=1,10
           read(7,'(F8.7)') delg(ll)   !g_ordinates
        enddo
        do ll=1,1
           read(7,*) dum    
        enddo 
        do n=1,pnum
           read(7,'(F15.13)') err
           pbar(n)=err*1.01325
        enddo
        do n=1,2
           read(7,*) dum    
        enddo
        do m=1,tnum
           read(7,'(F7.4)') t(m)
        enddo
        do ll=1,9
           read(7,*) dum
        enddo

         do i=1,1501
           do il=1,tnum !temperature   (50,410,20)
             do ll=1,7
               read(7,fmt='(A70)') dum
             enddo
             do j=1,pnum               !number of P
	       read(7,*) err,info(1,j,il,i),info(2,j,il,i),
     &info(3,j,il,i),info(4,j,il,i),info(5,j,il,i)
               read(7,*) err,info(6,j,il,i),info(7,j,il,i),
     &info(8,j,il,i),info(9,j,il,i),info(10,j,il,i)
             enddo
           enddo
           wnirw(i)=(i-1)*5.+2000.
	   wnirwl(i)=(i-1)*5.+1997.5
	   wnirwh(i)=(i-1)*5.+2002.5
*	   print*,wnirw(i)
        enddo
        
      close(7) 
      return  
      end

