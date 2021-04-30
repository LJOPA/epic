        Subroutine irwin_only(l,ff,nunwn,nlay,cd,tirw,
     &       pirw,trans_arr,tnum,pnum,info,play_o,t,irwintest,
     &       wnout,nwn,absfact,trans_arr3,kcheck,delg,col)
        implicit none
        integer pnum,tnum,nlay,nunwn
        integer i,ii,l,jj,jk,j,nwn,m,s

       Real*8 cd(nlay),trans_arr(nunwn,2,2),
     &     info(10,pnum,tnum,nunwn),
     &     pirw_lo(nlay),tirw_hi(nlay),pirw_hi(nlay),tirw(tnum), 
     &     pirw(pnum),ff(1,nlay),slope,b,dum,cdff,
     &     k1(10),k1h(10),k1hp(10),wnout(nwn),col(0:nlay),
     &     k1ht(10),play_o(nlay),t(nlay),delg(10),wavel(119),
     &     k2(10),pix(119),wno(119),afact(119),knew(119),wavell(119),
     &     ka(2),pixl(119),wnol(119),afactl(119),knewl(119),wnolo,
     &     pa(2),afactnew,knew2,wnohi,k1h1(10),k1hp1(10),k1ht1(10),
     &     dy,k11(10),afacthi,afactlo,klo,khi,trans_arr3(nunwn,10,2,2)
      Logical irwintest,absfact,kcheck

c  This subroutine opens the Irwin data file and extracts the information.
c The Irwin file consists of blocks of temperatures and for each temp(19),
c there is 20 pressures each with 10 k-coeff values. The statistical 
c weights(10) remain constant for each and every set of 10 k-coeff and are 
c used in conjunction with the k-coeff. They weight the "contribution" of
c each of the 10 k-coeff relative to the total transmission at each pressure.
c  the info_arr(#,p,t,wn). #'s below:
c  info(1): k1   info(6): k6
c  info(2): k2   info(7): k7
c  info(3): k3   info(8): k8
c  info(4): k4   info(9): k9
c  info(5): k5   info(10): k10
cccc 08/30/07  Adding Bruno's mods to the k. There are 'abs_factors' and
cccc           'k_abs'. 'k_abs' replace weak-line absorption coeff.
cccc           'abs_factors' multiply existing coeff. Need to interpolate.
        
 	 cdff=cd(l)*ff(1,l)
         call locate(pirw,pnum,play_o(l),jk)
         call locate(tirw,tnum,t(l),jj)
         if (play_o(l).lt.pirw(1)) then
           jk=1
         endif  
         if (t(l).lt.50.)then
           jj=1
         endif  
          do i=1,nunwn  !wavelength
             ii=i
            do j=1,10
              k1h(j)=info(j,jk+1,jj+1,i)  !hi temp, hi P
              k1hp(j)=info(j,jk+1,jj,i)  !lo temp,hi P
              k1ht(j)=info(j,jk,jj+1,i)   !hi temp, lo P
              k1(j)=info(j,jk,jj,i)      !lo temp, lo P
            end do
            do j=1,10
              trans_arr3(i,j,1,1)=0.
              trans_arr3(i,j,1,2)=0.
              trans_arr3(i,j,2,1)=0.
              trans_arr3(i,j,2,2)=0.
            enddo
            do j=1,10
              trans_arr3(i,j,1,1)=k1(j)*cdff
              trans_arr3(i,j,1,2)=k1ht(j)*cdff
              trans_arr3(i,j,2,2)=k1h(j)*cdff
              trans_arr3(i,j,2,1)=k1hp(j)*cdff
            enddo 
        return 
        end

