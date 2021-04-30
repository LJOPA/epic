        Subroutine irwin(l,ff,nunwn,nlay,cd,tirw_arr,pirw_arr,tirw,
     &       pirw,trans_arr,tnum,pnum,info,play_o,t,irwintest,
     &       wnout,nwn,absfact,trans_arr3,kcheck,delg)
        implicit none
        integer pnum,tnum,nlay,nunwn
        integer i,ii,l,jj,jk,j,nwn,m,s

       Real*8 cd(nlay),tirw_arr(2,nlay),trans_arr(nunwn,2,2),
     &     pirw_arr(2,nlay),info(10,pnum,tnum,1501),
     &     pirw_lo(nlay),tirw_hi(nlay),pirw_hi(nlay),tirw(tnum), 
     &     pirw(pnum),ff(1,nlay),slope,b,tirw_lo(nlay),dum,
     &     trans(10),k1(10),k1h(10),k1hp(10),wnout(nwn),col(0:nlay),
     &     k1ht(10),play_o(nlay),t(nlay),delg(10),wavel(119),
     &     k2(10),pix(119),wno(119),afact(119),knew(119),wavell(119),
     &     ka(2),pixl(119),wnol(119),afactl(119),knewl(119),wnolo,
     &     pa(2),afactnew,knew2,wnohi,k1h1(10),k1hp1(10),k1ht1(10),
     &     dy,k11(10),afacthi,afactlo,klo,khi,trans_arr3(nunwn,10,2,2),
     &     trans1(10),trans2(10),trans3(10)
      Logical irwintest,absfact,kcheck

c  This subroutine opens the Irwin data file and extracts the information.
c The Irwin file consists of blocks of temperatures and for each temp(19),
c there is 20 pressures each with 10 k-coeff values. The statistical 
c weights(10) remain constant for each and every set of 10 k-coeff and are 
c used in conjunction with the k-coeff. They weight the "contribution" of
c each of the 10 k-coeff relative to the total transmission at each pressure.
c  the info_arr(#,p,t,wn). #'s below:
c  info(1): k1
c  info(2): k2
c  info(3): k3
c  info(4): k4
c  info(5): k5
c  info(6): k6
c  info(7): k7
c  info(8): k8
c  info(9): k9
c  info(10): k10
cccc 08/30/07  Adding Bruno's mods to the k. There are 'abs_factors' and
cccc           'k_abs'. 'k_abs' replace weak-line absorption coeff.
cccc           'abs_factors' multiply existing coeff. Need to interpolate.
         if(absfact) then
            open(65,file='SAT/nsubs/kmod.csv',status='old')
 300        format(e9.2,1x,e9.2,1x,e9.2,1x,e9.2)
            do m=1,119
               read(65,300) wavel(m),wno(m),afact(m),knew(m)
            enddo
            close(65)
         endif           
*         delg(1) = 3.333600E-02
*         delg(2) = 7.472600E-02
*         delg(3) = 0.1095430
*         delg(4) = 0.1346330
*         delg(5) = 0.1477620
*         delg(6) = 0.1477620
*         delg(7) = 0.1346330
*         delg(8) = 0.1095430
*         delg(9) = 7.472600E-02
*         delg(10) =3.333600E-02

        
*        do l=1,nlay
         call locate(pirw,pnum,play_o(l),jk)
         call locate(tirw,tnum,t(l),jj)
         tirw_lo(l)=tirw(jj)      !lo & hi bracket value of T & P
         pirw_lo(l)=pirw(jk)
         tirw_hi(l)=tirw(jj+1)
         pirw_hi(l)=pirw(jk+1)
         if (play_o(l).lt.pirw(1)) then
           jk=1
           pirw_lo(l)=play_o(l)
           pirw_hi(l)=pirw(1)
         endif  
         if (t(l).lt.50.)then
           jj=1
           tirw_lo(l)=t(l)
           tirw_hi(l)=tirw(1) 
         endif  
          tirw_arr(1,l)=tirw_lo(l)
          tirw_arr(2,l)=tirw_hi(l)
          pirw_arr(1,l)=pirw_lo(l)
          pirw_arr(2,l)=pirw_hi(l)
          do i=559,1501 !,nunwn  !wavelength
             ii=i-558
            if(irwintest)then         !A
!             IF(wnout(i).lt.4800.)then   !B
                do j=1,10
                   k1h(j)=info(j+3,jk+1,jj+1,i)  !hi temp, hi P
                   k1hp(j)=info(j+3,jk+1,jj,i)  !lo temp,hi P
                   k1ht(j)=info(j+3,jk,jj+1,i)   !hi temp, lo P
                   k1(j)=info(j+3,jk,jj,i)      !lo temp, lo P
                 end do
!              IF(wnout(i).ge.4800.)then   !B
                 do j=1,10         !test to see if Bruno is write about
                   k1h(j)=k1h(j)/2. !Irwin K's being too big
                   k1hp(j)=k1hp(j)/2.
                   k1ht(j)=k1ht(j)/2.
                   k1(j)=k1(j)/2.
                 enddo
!              ENDIF   
                   !B
            elseif(absfact)then    !A     Bruno's changes to k's
          if(wnout(ii).ge.6170. .and. wnout(ii).le.12000.) then!D
                call locate(wno,119,wnout(ii),s)
            !     print*,wnout(i),'wnout'
                 wnolo=wno(s)
                 wnohi=wno(s+1)
                 afactlo=afact(s)
                 afacthi=afact(s+1)
                 klo=knew(s)
                 khi=knew(s+1)
                 slope=(afacthi-afactlo)/(wnohi-wnolo)
                 b=(afacthi)-slope*wnohi
                 afactnew=slope*wnout(ii)+b
                 do j=1,10
                   k1h1(j)=(info(j+3,jk+1,jj+1,i))*afactnew  !hi t, hi P
                   k1hp1(j)=(info(j+3,jk+1,jj,i))*afactnew  !lo t,hi P
                   k1ht1(j)=(info(j+3,jk,jj+1,i))*afactnew   !hi temp, lo P
                   k11(j)=(info(j+3,jk,jj,i))*afactnew      !lo temp, lo P
                 end do
                 slope=(khi-klo)/(wnohi-wnolo)
                 b=(khi)-slope*wnohi
                 knew2=slope*wnout(ii)+b
                 do j=1,10               !weak line updates
                   k1h(j)=k1h1(j)+knew2
                   k1hp(j)=k1hp1(j)+knew2
                   k1ht(j)=k1ht1(j)+knew2
                   k1(j)=k11(j)+knew2
                 enddo
              else
                 do j=1,10
                   k1h(j)=info(j+3,jk+1,jj+1,i)  !hi temp, hi P
                   k1hp(j)=info(j+3,jk+1,jj,i)  !lo temp,hi P
                   k1ht(j)=info(j+3,jk,jj+1,i)   !hi temp, lo P
                   k1(j)=info(j+3,jk,jj,i)      !lo temp, lo P
                 end do
              endif  !D
            else    !A
               do j=1,10
                   k1h(j)=info(j+3,jk+1,jj+1,i)  !hi temp, hi P
                   k1hp(j)=info(j+3,jk+1,jj,i)  !lo temp,hi P
                   k1ht(j)=info(j+3,jk,jj+1,i)   !hi temp, lo P
                   k1(j)=info(j+3,jk,jj,i)      !lo temp, lo P
                 end do
            endif  !A
            if (play_o(l).le.pirw(1)) then
                do j=1,10
                   slope=((k1(j))-(k1hp(j)))/((pirw(1))-
     &                      (pirw(2)))
                   b=(k1(j))-slope*(pirw(1))
                   k1hp(j)=k1(j)
                   k1(j)=(slope*(play_o(l))+b)
                end do
                do j=1,10
                   slope=((k1ht(j))-(k1h(j)))/((pirw(1))-
     &                      (pirw(2)))
                   b=(k1ht(j))-slope*(pirw(1))
                   k1h(j)=k1ht(j)
                   k1ht(j)=(slope*(play_o(l))+b)
                    !messed around with linear/log interp/extrap
                    !doesn't make a difference in the heating rates
                enddo
                do j=1,10
                    if(k1(j).lt.0.) then
                      k1(j)=1.e-25
                    elseif(k1ht(j).lt.0.) then
                      k1ht(j)=1.e-25
                    elseif(k1hp(j).lt.0.) then
                      k1hp(j)=1.e-25
                    elseif(k1h(j).lt.0.) then
                      k1h(j)=1.e-25
                    endif
                 enddo
                 pirw_arr(1,l)=play_o(l)
                 pirw_arr(2,l)=pirw(1)
                  endif
                  if (t(l).lt.50.) then
                    do j=1,10
                     slope=(k1(j)-k1ht(j))/((tirw(jj))-
     &                        (tirw(jj+1)))
                     b=k1(j)-slope*(tirw(jj))
                     k1ht(j)=k1(j)
                     k1(j)=slope*t(l)+b
                    end do
                    do j=1,10
                     slope=(k1hp(j)-k1h(j))/((tirw(jj))-
     &                        (tirw(jj+1)))
                     b=k1hp(j)-slope*(tirw(jj))
                     k1h(j)=k1hp(j)
                     k1hp(j)=slope*t(l)+b
                    end do
                   end if
                 trans_arr(ii,1,1)=0.
                 trans_arr(ii,1,2)=0.
                 trans_arr(ii,2,1)=0.
                 trans_arr(ii,2,2)=0.
               do j=1,10
                 trans_arr3(ii,j,1,1)=0.
                 trans_arr3(ii,j,1,2)=0.
                 trans_arr3(ii,j,2,1)=0.
                 trans_arr3(ii,j,2,2)=0.
              enddo
             if(kcheck)then
               do j=1,10
                 trans(j)=k1(j)*cd(l)*ff(1,l)
                 trans_arr3(ii,j,1,1)=trans(j)
                 trans1(j)=(k1ht(j)*cd(l)*ff(1,l))
                 trans_arr3(ii,j,1,2)=trans1(j)
                 trans2(j)=(k1h(j)*cd(l)*ff(1,l))
                 trans_arr3(ii,j,2,2)=trans2(j)
                 trans3(j)=(k1hp(j)*cd(l)*ff(1,l))
                 trans_arr3(ii,j,2,1)=trans3(j)

               !   print*,trans_arr3(ii,j,2,1),k1hp(j)
               enddo 
             else
                 do j=1,10
                   trans(j)=exp(-k1(j)*cd(l)*ff(1,l))
                   trans_arr(ii,1,1)=trans_arr(i,1,1)+trans(j)
!                   print*,trans_arr(i,1,1),'1'
                   trans1(j)=exp(-k1ht(j)*cd(l)*ff(1,l))
                   trans_arr(ii,1,2)=trans_arr(i,1,2)+trans1(j)
!                   print*,trans_arr(i,1,2),'2'
                   trans2(j)=exp(-k1h(j)*cd(l)*ff(1,l))
                   trans_arr(ii,2,2)=trans_arr(i,2,2)+trans2(j)
!                  print*,trans_arr(i,2,2),'3'
                   trans3(j)=exp(-k1hp(j)*cd(l)*ff(1,l))
                   trans_arr(ii,2,1)=trans_arr(i,2,1)+trans3(j)
!                   print*,trans_arr(i,2,1),'4'
                 end do
            endif
         IF(kcheck)then
!          do j=1,10
!           if(trans_arr3(ii,j,1,1).gt.1.0.and.
!     &        trans_arr3(ii,j,1,1).lt.1.003)then
!              trans_arr3(ii,j,1,1)=1.
!           endif
!           if(trans_arr3(ii,j,1,2).gt.1.0.and.
!     &        trans_arr3(ii,j,1,2).lt.1.003)then
!              trans_arr3(ii,j,1,2)=1.
!           endif
!           if(trans_arr3(ii,j,2,1).gt.1.0.and.
!     &        trans_arr3(ii,j,2,1).lt.1.003)then
!              trans_arr3(ii,j,2,1)=1.
!           endif
!           if(trans_arr3(ii,j,2,2).gt.1.0.and.
!     &        trans_arr3(ii,j,2,2).lt.1.003)then
!              trans_arr3(ii,j,2,2)=1.
!           endif
!             if (trans_arr3(ii,j,1,1).le.0.)then
!                 trans_arr3(ii,j,1,1)=1.e-22
!             endif
!             if (trans_arr3(ii,j,1,2).le.0.)then
!                 trans_arr3(ii,j,1,2)=1.e-22
!             endif
!             if (trans_arr3(ii,j,2,1).le.0.)then
!                 trans_arr3(ii,j,2,1)=1.e-22
!             endif
!             if (trans_arr3(ii,j,2,2).le.0.)then
!                 trans_arr3(ii,j,2,2)=1.e-22
!             endif
!          enddo
         ELSE

         if(trans_arr(ii,1,1).gt.1.0.and.trans_arr(ii,1,1).lt.1.003)then
              trans_arr(ii,1,1)=1.
         endif
         if(trans_arr(ii,1,2).gt.1.0.and.trans_arr(ii,1,2).lt.1.003)then
              trans_arr(ii,1,2)=1.
         endif
         if(trans_arr(ii,2,1).gt.1.0.and.trans_arr(ii,2,1).lt.1.003)then
              trans_arr(ii,2,1)=1.
         endif
         if(trans_arr(ii,2,2).gt.1.0.and.trans_arr(ii,2,2).lt.1.003)then
              trans_arr(ii,2,2)=1.
         endif
                   if (trans_arr(ii,1,1).le.0.)then
                      trans_arr(ii,1,1)=1.e-22
                   endif
                   if (trans_arr(ii,1,2).le.0.)then
                      trans_arr(ii,1,2)=1.e-22
                   endif
                   if (trans_arr(ii,2,1).le.0.)then
                      trans_arr(ii,2,1)=1.e-22
                   endif
                   if (trans_arr(ii,2,2).le.0.)then
                      trans_arr(ii,2,2)=1.e-22
                   endif
         endif              
         ENDDO
*           enddo
        return 
        end

