          SUBROUTINE nir(wnout,nlay,nwn,t,wnf,wvnmhi_n,wvnmlo_n,
     &         dtaucc,workdir,play_o,pnum,tnum,info_arr,nunwn,
     &         ff_ch4,dlay,cd,irwinhalf,absfact,kcheck,dtaucc5,delg,
     &tirw,pirw)

         implicit none
         integer jj,nunwn,nlay,nwn,pnum,tnum,s,k,kl,h,ii,i,l,j
         REAL*8 kll,wnf(nwn),wvnmhi_n(nwn),
     &   wvnmlo_n(nwn),t(nlay),dtaucc(nlay,nwn),
     &   play_o(nlay),cd(nlay),dum,ef,dtaucc5(10,nlay,1501),
     &   tirw(tnum),pirw(pnum),pirw_arr(2,nlay),trans_arr5(2,2),
     &   dlay(nlay),pirw_arr2(2),trans_arr2(2,2),tra,dy,col(0:nlay),
     &   delg(10),trans_arr4(10,2,2),
     &   tirw_arr(2,nlay),tirw_arr2(2),trans_arr(nunwn,2,2),
     &   info_arr(10,pnum,tnum,1501),trans_arr3(nunwn,10,2,2),
     &   trans(nlay,nwn),wnout(nwn),ff_ch4(1,nlay),
     &   slope,b,dtauc5(10,nlay,nunwn,2),trans5(10,nlay,nwn)
            character*34 workdir
            Logical irwinhalf,absfact,kcheck
!t(nlay) IS tlay(nlay), just so you know
!------------------------------------------------------------------------
!             Read in Irwin parameters 
!               k-coeff [km/amagat], p[atm],t[K]
!------------------------------------------------------------------------
           do jj=1,10
             do j=1,nlay
               do l=1,nunwn
                  dtaucc5(jj,j,l)=0.
               enddo
             enddo
            enddo


            h=h+nlay
            wnout(1)=wnf(1)
            ef=5.
            do i=2,nunwn
               wnout(i)=4790.+ef
               ef=ef+5.
            enddo
            do i=2,nunwn-1
               wvnmhi_n(i)=(wnout(i)+wnout(i+1))/2.
               wvnmlo_n(i)=(wnout(i)+wnout(i-1))/2.
            enddo

            do j=1,nlay
            call irwin(j,ff_ch4,nunwn,nlay,cd,tirw_arr,pirw_arr,
     &tirw,pirw,trans_arr,tnum,pnum,info_arr,play_o,t,irwinhalf,wnout,
     &nwn,absfact,trans_arr3,kcheck,delg)
            if (kcheck) then
              do l=1,nunwn
                do kl=1,2
                  do ii=1,2
                    do jj=1,10
                      trans_arr4(jj,kl,ii)=trans_arr3(l,jj,kl,ii)
                      pirw_arr2(kl)=pirw_arr(kl,j)
                      tirw_arr2(ii)=tirw_arr(ii,j)
                    enddo
                   enddo
                 enddo
                   do jj=1,10
                     do ii=1,2
                       trans_arr5(2,ii)=trans_arr4(jj,2,ii)
                       trans_arr5(1,ii)=trans_arr4(jj,1,ii)
                       slope=((trans_arr5(2,ii))-
     &                        (trans_arr5(1,ii)))/
     &                       (log(pirw_arr2(2))-log(pirw_arr2(1)))
                       b=(trans_arr5(2,ii))-
     &                       slope*log(pirw_arr2(2))
                       dtauc5(jj,j,l,ii)=(slope*log(play_o(j))+b)
                     enddo
                   slope=((dtauc5(jj,j,l,2))-(dtauc5(jj,j,l,1)))/
     &                      (tirw_arr2(2)-tirw_arr2(1))
                     b=(dtauc5(jj,j,l,2))-slope*tirw_arr2(2)
                     dtaucc5(jj,j,l)=(slope*t(j)+b)
                   enddo
               enddo
             ELSE
               do l=1,nunwn
                  do kl=1,2
                    do ii=1,2
                       trans_arr2(kl,ii)=trans_arr(l,kl,ii)
                       pirw_arr2(kl)=pirw_arr(kl,j)
                       tirw_arr2(ii)=tirw_arr(ii,j)
                    enddo
                   enddo
          call polin2(pirw_arr2,tirw_arr2,trans_arr2,2,2,play_o(j),t(j),
     &               tra,dy)
                   dtaucc(j,l)=(tra)
                
               enddo
           ENDIF
                      
            enddo
           return
           END
