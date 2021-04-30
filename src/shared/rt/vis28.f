         SUBROUTINE vis2(nlay,nwn,t,wnf,wvnmhi_vhm,
     &       wvnmlo_vhm,wnot,dtaucc,ff_ch4,cd,nunwn,workdir)

         Implicit none 
         integer nunwn,nlay,nwn,k,i,l,po
         REAL*8 wnf,dtaucc(nlay,nwn),wnfvlmx,wn1vlmx,
     &       wnfvhm,wn1vhm,wn1vlm,wnfvlm,ef,wn(nunwn),wnot(nwn),
     &       wvnmhi_vhm(nwn),wvnmlo_vhm(nwn),t(nlay),cd(nlay),
     &       ff_ch4(1,nlay),vac,air,methvis(4750)
         character*34 workdir
         wnfvlmx=3.3289E4           !300.4nm
         wn1vlmx=1.9231E4           !520.0
         wnfvhm=1.9227E4            !520.1
         wn1vhm=1.0050E4            !995.0
         wnfvlm=1.0048E4            !995.2
         wn1vlm=9523.81                   !1050.0

         wn(1)=wnf
         wnot(1)=wn(1)
         ef=0.1
         do i=2,nunwn
            wn(i)=1/((995.0-ef)*0.0000001)
            wnot(i)=wn(i)
            ef=ef+0.1
         enddo
         do i=2,nunwn-1
           wvnmhi_vhm(i)=(wn(i)+wn(i+1))/2.
           wvnmlo_vhm(i)=(wn(i)+wn(i-1))/2.
         enddo
         
!-------------------------------------------------------------
          open(43,file='nsubs/1995hi_karko.tab',status='old')
*          open(40,file=workdir//'vis2_tau.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
          do po=1,4750
            read(43,'(1x,f5.1,2x,f6.2,f8.4)') vac,air,
     &        methvis(po)
          enddo
          close(43)
          do l=1,nlay
             Do k=0,4749
                dtaucc(l,k+1) = cd(l)*methvis(4750-k)*ff_ch4(1,l)
             End do
*             write(40,rec=l) (real(dtaucc(l,k)),k=1,nunwn)
          END DO
*          close(40)
          return
          END

