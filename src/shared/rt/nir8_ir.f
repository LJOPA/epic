          SUBROUTINE nir_ir(nlay,nwn,t,play_o,pnum,tnum,
     &info_arr,nnwn,ff_ch4,cd,dtaucc5,tirw,pirw)

         Implicit None

         integer jj,nnwn,nlay,nwn,pnum,tnum,s,k,kl,h,
     &ii,i,l,j,jk,ll
         REAL*8 t(nlay),play_o(nlay),cd(nlay),
     &dtaucc5(10,nlay,nnwn),tirw(tnum),pirw(pnum),
     &info_arr(10,pnum,tnum,nnwn),ff_ch4(1,nlay),m,b,
     &dtauc5(10,nlay,nnwn,2)

!t(nlay) IS tlay(nlay), just so you know
!------------------------------------------------------------------------
!             Read in Irwin parameters 
!               k-coeff [km/amagat], p[atm],t[K]
!------------------------------------------------------------------------

            h=h+nlay

*            open(14,file=workdir//'nir_play.out',form='unformatted', 
*     &          access='direct',recl=4*nlay)
!       open(13,file=workdir//'nir_tau_eq270n.out',form='unformatted', 
!     &          access='direct',recl=4*nnwn)
*            open(15,file=workdir//'nir_trans.out',form='unformatted',
*     &          access='direct',recl=4*nnwn)
            do l=1,nlay
              call locate(pirw,pnum,play_o(l),jk)
              call locate(tirw,tnum,t(l),jj)
              if (play_o(l).lt.pirw(1)) then
                jk=1
              endif  
              if (t(l).lt.50.)then
                jj=1
              endif  
*info_arr(j,jk,jj,ll)     !lo temp, lo P
*info_arr(j,jk,jj+1,ll)   !hi temp, lo P
*info_arr(j,jk+1,jj+1,ll) !hi temp, hi P
*info_arr(j,jk+1,jj,ll)  !lo temp,hi P
              do ll=1,nnwn
                do j=1,10
                  do ii=1,2
                    m=(info_arr(j,jk,jj+(ii-1),ll)-
     &info_arr(j,jk+1,jj+(ii-1),ll))/(pirw(jk)-pirw(jk+1))
                    b=info_arr(j,jk,jj+(ii-1),ll)-m*pirw(jk)
                    dtauc5(j,l,ll,ii)=m*play_o(l)+b
                  end do 
                  m=(dtauc5(j,l,ll,2)-dtauc5(j,l,ll,1))/
     &(tirw(jj+1)-tirw(jj))
                  b=dtauc5(j,l,ll,2)-m*tirw(jj+1)
                  dtaucc5(j,l,ll)=(m*t(j)+b)*cd(l)*ff_ch4(1,l)
                enddo
              enddo
	    enddo
!            close(13)

           return
           END
