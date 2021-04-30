       SUBROUTINE midtexpb(nlay,t,play,ptau,npres,kj,wnhi,wnlo,wn,
     &   dtauc,nwn,nunwn,cd,ff_ch4,ff_c2h2,ff_c2h6,qh2,he,dlay)

       Implicit none

       integer j,nunwn,nlay,nwn,npres,k,intt(nlay)
       REAL*8 t(nlay),wnhi(nwn),wn(nwn),qh2,
     & wnlo(nwn),ptau(npres),h22,he,dlay(nlay),m,b,cdlay,
     & dtauc(nlay,nwn),mid_c3h8,ff_ch4(1,nlay),ff_c2h2(1,nlay),
     & ff_c2h6(1,nlay),play(nlay),cd(nlay),d25,h2,
     & invt,invt1,invt2,Ec,test,qh2he,ffh2(nlay),temp

       Real tauc(5,nunwn),dum1(nunwn)
       INTEGER rec1,rec2,kj,i,ij,l,tj
       CHARACTER*1 str1,str2*2,str*3,file1*30

       h22=qh2*qh2
       h2=1.43879
       qh2he=qh2*he
       do i=1,nlay
          ffh2(i)=h2*ff_ch4(1,i)
       enddo
       file1='/data/seasonaltom/code/opacity'

       tj=kj+1
       if (tj.lt.10) then
         write (str1,'(I1)'),tj
         str='00'//str1
       elseif (tj.lt.100) then
         write (str2,'(I2)'),tj
         str='0'//str2
       else
         write (str,'(I3)'),tj
       endif

       d25=1./25.
       do i=1,nlay
         intt(i)=int(t(i)/25)
       end do

       open(45,file=file1//'/data/opacity_mir0_'//str//'.dat',
     &form='unformatted',access='direct',recl=4*nunwn)
       open(46,file=file1//'/data/opacity_mir2_'//str//'.dat'
     &,form='unformatted',access='direct',recl=4*nunwn)
       open(47,file=file1//'/data/opacity_mir3_'//str//'.dat',
     &form='unformatted',access='direct',recl=4*nunwn)
       open(48,file=file1//'/data/opacity_mir4_'//str//'.dat'
     &,form='unformatted',access='direct',recl=4*nunwn)



       read(45,rec=1) (dum1(i),i=1,nunwn)
       do i=1,nunwn
         wnlo(i)=dum1(i)
       end do
       read(45,rec=2) (dum1(i),i=1,nunwn)
       do i=1,nunwn
         wn(i)=dum1(i)
       end do
       read(45,rec=3) (dum1(i),i=1,nunwn)
       do i=1,nunwn
         wnhi(i)=dum1(i)
       end do
************************************************
       do l=1,nlay
***********************************************
         temp=25.*intt(l)
!I think I could come up with a way to speed this up and
!it would significantly increase the speed of the code
         do i=1,3
           rec1=3+i+(intt(l)-3)*3.
           read(45,rec=rec1) (tauc(4,ij),ij=1,nunwn)
           read(45,rec=rec1+3) (tauc(5,ij),ij=1,nunwn)
           do ij=1,nunwn                !remove
             m=(tauc(5,ij)-tauc(4,ij))*d25
	     b=tauc(5,ij)-m*temp
             tauc(i,ij)=m*t(l)+b
           enddo                
	 enddo
         cdlay=cd(l)*dlay(l)
         do i=1,nunwn
           dtauc(l,i)=cdlay*(tauc(1,i)*h22+tauc(2,i)*qh2he+
     &tauc(3,i)*ffh2(l))
         enddo
********************************************************************
         rec1=3+l+npres*(intt(l)-3)
         invt1=1./(25.*intt(l))
         invt2=1./(25.*intt(l)+25.)
         invt=1./(t(l))
	if (tj.ge.78.and.tj.lt.175) then
         read(46,rec=rec1) (tauc(1,ij),ij=1,nunwn)
         read(46,rec=rec1+npres) (tauc(2,ij),ij=1,nunwn)
         do ij=1,nunwn
           Ec=log(tauc(1,ij)/tauc(2,ij))/(h2*(invt2-invt1))
           dtauc(l,ij)=tauc(1,ij)*exp(h2*Ec*
     &(invt1-invt))*ff_ch4(1,l)*cd(l)+dtauc(l,ij)
         enddo
        endif
*********************************************************************
	if ((tj.le.75).or.((tj.ge.127).and.tj.le.172)) then
         read(47,rec=rec1) (tauc(1,ij),ij=1,nunwn)
         read(47,rec=rec1+npres) (tauc(2,ij),ij=1,nunwn)
         do ij=1,nunwn
           Ec=log(tauc(1,ij)/tauc(2,ij))/(h2*(invt2-invt1))
           dtauc(l,ij)=tauc(1,ij)*exp(h2*Ec*
     &(invt1-invt))*ff_c2h2(1,l)*cd(l)+dtauc(l,ij)
         enddo
        endif
***********************************************************************
        if (tj.ge.42.and.tj.le.82) then
         read(48,rec=rec1) (tauc(1,ij),ij=1,nunwn)
         read(48,rec=rec1+npres) (tauc(2,ij),ij=1,nunwn)
         do ij=1,nunwn
           Ec=log(tauc(1,ij)/tauc(2,ij))/(h2*(invt2-invt1))
           dtauc(l,ij)=tauc(1,ij)*exp(h2*Ec*
     &(invt1-invt))*ff_c2h6(1,l)*cd(l)+dtauc(l,ij)
         enddo
	endif
       enddo 

       close(45)
       close(46)
       close(47)
       close(48)

         return
         END
