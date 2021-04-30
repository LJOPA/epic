          SUBROUTINE fir(nlay,nwn,t,wnhi,wnlo,wn,
     &          dtauc,nunwn,qh2,he,ff_ch4,cd,dlay)

          implicit none
          integer nunwn,nlay,nwn,l,i,ij,rec1
          real tauf(nunwn),tauf2(nunwn),dum1(nunwn)
          real*8 t(nlay),wnhi(nwn),wnlo(nwn),cd(nlay),dlay(nlay),
     &         wn(nwn),dtauc(nlay,nwn),dtau(3,nunwn),
     &         m,b,cdlay,h22,qh2,he,ff_ch4(1,nlay)
          character*30 file1

	  file1='/data/seasonaltom/code/opacity'
	  h22=qh2**2
          open(45,file=file1//'/data/opacity_fir0_1.dat'
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
            do l=1,nlay
	      do i=1,3
                rec1=3+i+(int(t(l)/25)-3)*3.
                read(45,rec=rec1) (tauf(ij),ij=1,nunwn)
                read(45,rec=rec1+3) (tauf2(ij),ij=1,nunwn)
              do ij=1,nunwn
                m=(tauf2(ij)-tauf(ij))/(25.)
		b=tauf2(ij)-m*25.*int(t(l)/25)
                dtau(i,ij)=m*t(l)+b
              enddo
	      enddo
	      cdlay=cd(l)*dlay(l)
	      do i=1,nunwn
		dtauc(l,i)=cdlay*(dtau(1,i)*h22+dtau(2,i)*qh2*he+
     &dtau(3,i)*qh2*ff_ch4(1,l))
	      enddo
            enddo
          close(45)


          return
         end

