        SUBROUTINE nirh(xc,nlay,t,kj,wnf,wvnmhi_n,wvnmlo_n,wnot,
     &             dtaucc,tlay_o,lat_ratio,lat_mixch4,
     &             nwn,nunwn,hit_geisa,hitran,R)
        implicit none
        integer j,nunwn,nlay,nwn
       REAL*8 t(nlay),wvnmhi_n(nwn),tlay_o(nlay),
     & lat_ratio(nlay),lat_mixch4(nlay),
     & R,wn(nunwn),wnf(nwn),h2,
     & t1(nlay),t2(nlay),wvnmlo_n(nwn),
     & newc,Ec,wnot(nwn),dtaucc(nlay,nwn),dum
       Real tauc2(nunwn),tauc(nunwn)
       Logical ct,ch4,hitran,hit_geisa
          INTEGER rec1,rec2,kj,strm,i,ij,l,xc(nlay)
          CHARACTER*3 str3,str
          CHARACTER*1 str1
          CHARACTER*2 str2
!          ct=.true.
         ch4=.true.
           do l=1,nlay
            do ij=1,nunwn
              dtaucc(l,ij)=0.
            enddo
           enddo
          if (nwn.ne.nunwn) then
           print*,'error in nir8_hitran.f: nwn and nunwn are different'
          endif
          h2 = 1.43879
          strm=kj 
          print*,'kj',kj
          wn(1)=wnf(1)       !setup wno scale
          wnot(1)=wn(1)
          do i=2,nunwn
            wn(i)=(wn(i-1)/R)+wn(i-1)
            wnot(i)=wn(i)
          enddo
          do i=2,nunwn-1  
            wvnmlo_n(i)=(wn(i)+wn(i-1))/2.
            wvnmhi_n(i)=(wn(i)+wn(i+1))/2.
          enddo
             if(strm.lt.10) then
               write (str1,'(I1)'),strm
               str='00'//str1
             elseif(strm.ge.10 .and. strm.lt.100) then
               write (str2,'(I2)'),strm
               str='0'//str2
             elseif(strm.ge.100)then
               write (str3,'(I3)'),strm
               str=str3
             endif
!          if (ct)then
!         open(46,
!     &file=
!     & 'SAT/nsubs/LAT/p40_opacityHR/opacity_mirct'//str//'.out',
!     &     form='unformatted',
!     &     access='direct',recl=4*nunwn)
!          do l=1,nlay
!                rec1=xc(l)*(nlay)
!                rec2=(xc(l)+1.)*(nlay)
!                t1(l)=tlay_o(l)+(xc(l)*10.)
!                t2(l)=tlay_o(l)+(10.+(xc(l)*10.))
!             read(46,rec=l+rec1) (tauc(ij),ij=1,nunwn)
!             read(46,rec=l+rec2) (tauc2(ij),ij=1,nunwn)
!             do ij=1,nunwn
!                 if (tauc2(ij).eq.0 .or. tauc(ij).eq.0) then
!                    newc=0.
!                 else
!                 Ec=(log(tauc(ij)/tauc2(ij)))*
!     &              1./((1/t2(l)-1/t1(l))*h2)
!                 newc=(tauc(ij)*(exp(h2*Ec*
!     &              ((1./t1(l))-(1./t(l))))))*lat_ratio(l)
!                 endif
!                 dtaucc(l,ij)=newc+dtaucc(l,ij)
!                 if (dtaucc(l,ij).lt.0.)then
!                   print*,dtaucc(l,ij),l,ij,'ct'
!                   print*,xc(l),rec1,rec2,t1(l),t2(l)
!                  pause
!                 endif
!             enddo
!          enddo
!          close(46)
!         else
!           do l=1,nlay
!            do ij=1,nunwn
!              dtaucc(l,ij)=0.
!            enddo
!           enddo
!        endif !ct
!        if(ct) then
!         open(46,
!     &file=
!     &  'SAT/nsubs/LAT/p40_opacityHR/opacity_mirct'//str//'.out',
!     &      form='unformatted',
!     &      access='direct',recl=4*nunwn)
!          do l=1,nlay
!                rec1=xc(l)*(nlay)
!                rec2=(xc(l)+1.)*(nlay)
!                t1(l)=tlay_o(l)+(xc(l)*10.)
!                t2(l)=tlay_o(l)+(10.+(xc(l)*10.))
!             read(46,rec=l+rec1) (tauc(ij),ij=1,nunwn)
!             read(46,rec=l+rec2) (tauc2(ij),ij=1,nunwn)
!             do ij=1,nunwn
!                 if (tauc2(ij).eq.0 .or. tauc(ij).eq.0) then
!                    newc=0.
!                 else
!                 Ec=(log(tauc(ij)/tauc2(ij)))*
!     &              1./((1/t2(l)-1/t1(l))*h2)
!                 newc=(tauc(ij)*(exp(h2*Ec*
!     &              ((1./t1(l))-(1./t(l))))))*lat_ratio(l)
!                 endif
!                 dtaucc(l,ij)=newc
!             enddo
!          enddo
!          close(46)
!         else
!           do l=1,nlay
!            do ij=1,nunwn
!              dtaucc(l,ij)=0. 
!            enddo
!           enddo
!         endif !ct
         newc=0.
         if(ch4)then
         if(HITRAN)then
            open(44,
     &      file='SAT/nsubs/LAT/p40_opacityHR/HITRAN/opacity_nirch4_'
     &        //str//'.out',
     &      form='unformatted',
     &      access='direct',recl=4*nunwn)
         endif
         if(HIT_GEISA)then
           print*,'HITRAN-GEISA NIR',kj
           if (kj.lt.256) then
             open(44,
     &         file='SAT/nsubs/LAT/p40_opacityHR/HITRAN/opacity_nirch4_'
     &           //str//'.out',form='unformatted',access='direct',
     &           recl=4*nunwn)
            elseif(kj.ge.256) then
          open(44,
     &    file='SAT/nsubs/LAT/p40_opacityHR/HITRAN/opacity_nirch4_geisa'
     &           //str//'.out',form='unformatted',access='direct',
     &           recl=4*nunwn)
             endif
         endif

          do l=1,nlay
                rec1=xc(l)*(nlay)
                rec2=(xc(l)+1.)*(nlay)
                t1(l)=tlay_o(l)+(xc(l)*10.)
                t2(l)=tlay_o(l)+(10.+(xc(l)*10.))
             read(44,rec=l+rec1) (tauc(ij),ij=1,nunwn)
             read(44,rec=l+rec2) (tauc2(ij),ij=1,nunwn)
             do ij=1,nunwn
                 if (tauc2(ij).eq.0 .or. tauc(ij).eq.0) then
                    newc=0.
                 else
                 Ec=(log(tauc(ij)/tauc2(ij)))*
     &              1./((1/t2(l)-1/t1(l))*h2)
                 newc=(tauc(ij)*(exp(h2*Ec*
     &              ((1./t1(l))-(1./t(l))))))*lat_mixch4(l)
                 endif
                 dtaucc(l,ij)=newc+dtaucc(l,ij)
             enddo
          enddo
          close(44)
          endif   
          do ij=1,nunwn
            do l=1,nlay
          dtaucc(l,ij)=dtaucc(l,ij)
          if (dtaucc(l,ij) .lt. 0.) then
            print*,dtaucc(l,ij),l,ij
          endif
           
              enddo
          enddo
       open(13,file='nir_tau_hitran.out',form='unformatted',
     &      access='direct',recl=4*nwn)
           do j=1,nlay
              write(13,rec=j) (real(dtaucc(j,l)),l=1,nwn)
            enddo
            close(13)
         print*,'done nirh'
         return
         END
