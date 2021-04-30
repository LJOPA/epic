       SUBROUTINE vis3(nlay,nwn,t,wnf,wvnmhi_vlmx,wvnmlo_vlmx,
     &            wnot,dtaucc,ff_ch4,cd,nunwn,workdir)

         implicit none
         integer nunwn,nlay,nwn,k,po,i,l
         REAL*8 wnf,wnot(nwn),wnfvlmx,wn1vlmx,
     &   wnfvhm,ff_ch4(1,nlay),
     &   wn1vlm,wnfvlm,ef,wn(nunwn),wn1vhm,wvnmhi_vlmx(nwn), 
     &   wvnmlo_vlmx(nwn),t(nlay),cd(nlay),
     &   dtaucc(nlay,nwn),methvisl(1875),vac,air
         character*34 workdir

         wnfvlmx=3.3289E4           !300.4nm
         wn1vlmx=1.9231E4           !520.0
         wnfvhm=1.9227E4            !520.1
         wn1vhm=1.0050E4            !995.0
         wnfvlm=1.0048E4            !995.2
         wn1vlm=9523.81                   !1050.0


          wn(1)=wnf
          wnot(1)=wn(1)
          ef=0.4
          do i=2,nunwn
             wn(i)=1/((520.0-ef)*0.0000001)
             wnot(i)=wn(i)
             ef=ef+0.4
          enddo
          do i=2,nunwn-1
             wvnmhi_vlmx(i)=(wn(i)+wn(i+1))/2.
             wvnmlo_vlmx(i)=(wn(i)+wn(i-1))/2.             
          enddo

      open(45,file='nsubs/1995lo_karko.tab',status='old') 
*          open(40,file=workdir//'vis3_tau.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
         do po=1,1875
           read(45,'(f6.1,1x,f7.2,1x,f7.4)') vac,air,methvisl(po)
         enddo
      close(45)
       
          do l=1,nlay
             Do k=0,549
                dtaucc(l,k+1) = (cd(l))*methvisl(550-k)*ff_ch4(1,l)
             End do
*             write(40,rec=l) (real(dtaucc(l,k)),k=1,nunwn)
          ENDDO
*          close(40)
        return

      END

