       implicit none
       integer nr,ntp,nspec,nwn,nlay,nl,tpc,ngas,nh2
       Parameter (ntp=80,nr=ntp,nspec=1,nwn=16000,nlay=70,nl=100000,
     &    tpc=20,ngas=1,nh2=150)
       REAL*8 p(nlay),r,d(nlay),xgeo(ngas),mw(ngas),nexp(ngas),
     &   wn1(nspec),wnf(nspec),delwn(nspec),wn(nspec,nwn),cof,
     &   f(nspec,nl,ngas),s(nspec,nl,ngas),a(nspec,nl,ngas),
     &   e(nspec,nl,ngas),th2(nh2),xmr2(nr),xmr3(nr),
     &   pr(ntp),tr(ntp),xm,he,h2,dratio,t(nlay),g,col(0:nlay),
     &   cd(nlay),mixh6fact,mixh8fact,mixh2fact,mix10fact,
     &   mix11fact,mixch4fact,mix13fact,amfact,lat,tlay(nlay),
     &   play(nlay),dlay(nlay),mtau(ngas,nspec,nlay,nwn),ff(ngas,nlay),
     &   scale,xmr1(nr),xmr4(nr),xmr5(nr),xmr6(nr),xmr7(nr),
     &   xmr8(nr),pxmr1(nr),pxmr2(nr),pxmr3(nr),pxmr4(nr),ab2(nh2),
     &   pxmr5(nr),pxmr6(nr),pxmr7(nr),mixnh3fact,ab1(nh2),
     &   pxmr8(nr),abshe(nh2,nwn),absh(nh2,nwn),absch4(nh2,nwn),
     &   bs1(nlay),bs2(nlay),tauc(0:tpc,nlay,nwn),tmin
       character*5 planet,oname
       Integer lnl,i,ii,j,nppt,npt(nspec,ngas),h,k,gas,ij,jk,
     &   il,iw,gasnum,nfile,hh,nunwn,kj,xc,lambda,nfile1
       Character*40 opacityfile_ch4,ptfile_ch4,opacityfile,ptfile
	character*80  gasf(ngas)
       Character*1 str1
       Character*3 str3,str
       Character*2  str2
       
*---------------------------------------------------------------------
 
        print*,'Pick region: (3)MIDIR, (4)FIR, (5)NIR'
        read*,lambda 
!        lambda=5
        if (lambda.eq.3) then
           print*,'Choose molecule: (1) CH4, (2) C2H2, (3) C2H6, (0)CT'
          read*,gas
!        gas=5
        elseif(lambda.eq.4) then
        gas=0
        elseif(lambda.eq.5) then
        gas=1   !CH4
        endif
	gasnum=gas+1   !This only works so far with ch4,c2h2,c2h6
        print*,'molecule:',gas,'lambda:',lambda
*---------------------------------------------------------------------
        open(33,file='sat_info.dat', status='old')
        read(33,*) planet
        read(33,*) he             !fraction of He  (plan Sci comp)
        read(33,*) h2             ! fraction of H2  (plan sci comp)
        read(33,*) xm             ! mean molec weight (kg) (PlanSciComp)
        read(33,*) dratio         !=for ch4  D/H value/1.5e-4 (for earth)
        close(33)

        open(31,file='var_script_shay',status='old')
        read(31,*) mixh2fact
        read(31,*) mixch4fact
        read(31,*) mixh6fact
        read(31,*) mixh8fact
        read(31,*) amfact
        read(31,*) mix10fact
        read(31,*) mix11fact
        read(31,*) mix13fact
        read(31,*) mixnh3fact
        read(31,*) lat
        close(31)
        call gravsub(planet,lat,g)
        g=g/100.
        cof=25.0
        xc=0
        wnf(1)=0.
        if(gas.eq.0)then
          oname='ct'
          if(lambda.eq.3)then
            wnf(1) =600.
            nfile=172         !111
            nfile1=1
          elseif(lambda.eq.4)then
            wnf(1) = 0.
            nfile=1
            nfile1=1
          endif
        elseif(gas.eq.1)then
          oname='ch4'
          if(lambda.eq.5)then
            wnf(1) =2000. 
            nfile1=173   !173
            nfile=321   !321
            gasf(1) = 'model/06_hit04.par'
          else              
            nfile1=82   !53
            nfile=172    !111
            wnf(1)=953.1361774178714
             gasf(1) = 'model/ch4_345t'
          endif
          nexp(1) = .45
          mw(1) = 16.e-3
          xgeo(1) = 1.5
        elseif (gas.eq.2) then                 !c2h2
          nexp(1) = .62
          oname='c2h2'
          nfile=81      !52
          nfile1=1
          mw(1) = 26.E-3        ! molecular weight (kg) of C2H2
          wnf(1) = 600.
          gasf(1) = 'model/c2h2_345t'
          xgeo(1) = 1.
        elseif (gas.eq.3) then                 !c2h6
          nexp(1) = .94
          nfile=81  !52
          nfile1=1
          oname='c2h6'
          mw(1) = 30.E-3
          wnf(1) = 600.
          gasf(1) = 'model/c2h6_345t'
          xgeo(1) = 1.5
        elseif (gas.eq.4) then                 !c2h4
          nfile=1
          oname='c2h4'
          nexp(1) = 1.
          mw(1) = 28.e-3
          wnf(1) =wnf(1)
      gasf(1) = 'molecules/c2h4_345t'
          xgeo(1) = 1.5
        elseif (gas.eq.5) then                 !c3h8
          nfile1=1
          nfile=51
          mw(1)=44.E-3
          oname='c3h8'
          wnf(1) = 600.
          gasf(1) = 'molecules/c3h8_345t'
          xgeo(1)=1.5
          nexp(1)=1.0
        elseif (gas.eq.6) then                 !ch3d
          nfile1=1  !53
          nfile=172    !111
          wnf(1)=600.      !952.53119878
          mw(1)=17.E-3
          xgeo(1)=1.5
          nexp(1)=.45
          gasf(1)='/home/grad38/sholmes/SAT/line_p/molecules/ch3d_345t'
        endif
        DO KJ=nfile1,nfile
         print*,kj,'nfile'
          h=0
          hh=0
          wn1(1) = wnf(1)
          DO jk=0,110,10
             xc=int(jk/10)
!	print*,xc
             if (kj.lt.10) then
                write (str1,'(I1)'),kj
                str='00'//str1
             elseif (kj.ge.10 .and. kj.lt.100) then
                write (str2,'(I2)'),kj
                str='0'//str2
             elseif (kj.ge.100) then
                write (str3,'(I3)'),kj
                str=str3
             endif
             do lnl=1,1
               wn(lnl,1)=wn1(lnl)
               if(lambda.ne.4)then
                 nunwn=nwn
                if(lambda.eq.5)then    !nir
                   R=1.52e6
                else
                   R=2800000.    !oops, overkill- only need 1.6e6 !mid
                endif
                write(6,201) wn1(1)
 201    format(' wno start: ',f7.1,' cm-1')
                 do i=2,nunwn
                   wn(lnl,i)=(wn(lnl,i-1)/R)+wn(lnl,i-1)
                 end do
                 delwn(lnl)=wn(lnl,2)-wn(lnl,1)
                 wnf(lnl)=wn(lnl,nunwn)
               elseif(lambda.eq.4)then
                 nunwn=600
                 do i=2,nunwn
                   wn(lnl,i)=(wn(lnl,i-1))+1.
                   delwn(lnl)=wn(lnl,2)-wn(lnl,1)
                   wnf(lnl)=wn(lnl,nunwn)
                 enddo 
               endif
              end do
! Pressure Scale:
        p(1)    =       1.02E-7
        scale   =       0.2273          ! pressure units [bars]
        call pscale(p(1),scale,p,nlay)
*------------------------------------------------------------------------
*  MidIR: Read in Geisa parameters                                      |
*------------------------------------------------------------------------
        if(lambda.ne.4 .and. gas.ne.0)then
             call hitranever(1,gasf,1,1,nl,ngas,cof,wn1,wnf,f,s,a,
     &               e,nppt)
             npt(1,1)=nppt
        endif
*------------------------------------------------------------------------
          open(15,file='mixsatmoses2000sep.dat',status='old')
 181      format(10e10.3)
          do i=1,ntp        ! read T-P profile
            read(15,181) pr(ntp+1-i),tr(ntp+1-i),xmr1(ntp+1-i),
     &       xmr2(ntp+1-i),xmr3(ntp+1-i),xmr4(ntp+1-i),
     &       xmr5(ntp+1-i),xmr6(ntp+1-i),xmr7(ntp+1-i),
     &       xmr8(ntp+1-i)

            pr(ntp+1-i)=pr(ntp+1-i)/1000.
            tr(ntp+1-i)=tr(ntp+1-i)-50.+jk
            pxmr1(ntp+1-i)=pr(ntp+1-i)
            pxmr2(ntp+1-i)=pr(ntp+1-i)
            pxmr3(ntp+1-i)=pr(ntp+1-i)
            pxmr4(ntp+1-i)=pr(ntp+1-i)
            pxmr5(ntp+1-i)=pr(ntp+1-i)
            pxmr6(ntp+1-i)=pr(ntp+1-i)
            pxmr7(ntp+1-i)=pr(ntp+1-i)
            pxmr8(ntp+1-i)=pr(ntp+1-i)
            xmr1(ntp+1-i)=xmr1(ntp+1-i)*mixch4fact
            xmr2(ntp+1-i)=xmr2(ntp+1-i)*mixh2fact
            xmr3(ntp+1-i)=xmr3(ntp+1-i)*mixh6fact
            xmr4(ntp+1-i)=xmr4(ntp+1-i)*mixh8fact
            xmr5(ntp+1-i)=xmr5(ntp+1-i)*mix10fact
            xmr6(ntp+1-i)=xmr6(ntp+1-i)*mix11fact
            xmr7(ntp+1-i)=xmr7(ntp+1-i)*mix13fact
            xmr8(ntp+1-i)=xmr8(ntp+1-i)*mixnh3fact
          enddo
          close(15)
          call tominterp(pr,tr,ntp,p,t,nlay)
          col(0) = 0.
          do j=1,nlay
            call column(p(j),t(j),xm,g,col(j),d(j))
            cd(j) = col(j) - col(j-1)
          end do
          play(1) = 0.5*p(1)                              ! LAYER VALUES
          do i=2,nlay 
            play(i) = (p(i-1) + p(i))/2.
          end do
          call tominterp(p,t,nlay,play,tlay,nlay)
          call tominterp(p,d,nlay,play,dlay,nlay)
          if (gas.eq.1) then                 !ch4
           do i=1,1
            call mixr(nlay,ngas,i,nr,ntp,pxmr1,xmr1,play,ff)
           enddo
          endif
          if (gas.eq.2) then                 !c2h2
           do i=1,ngas
            call mixr(nlay,ngas,i,nr,ntp,pxmr2,xmr2,play,ff)
           enddo
          endif
          if (gas.eq.3) then                 !c2h6
           do i=1,1
            call mixr(nlay,ngas,i,nr,ntp,pxmr3,xmr3,play,ff)
           enddo
          endif
          if (gas.eq.4) then                 !c2h4
           do i=1,1
            call mixr(nlay,ngas,i,nr,ntp,pxmr6,xmr6,play,ff)
           enddo
          endif
          if (gas.eq.5) then                 !c3h8
           do i=1,1
            call mixr(nlay,ngas,i,nr,ntp,pxmr4,xmr4,play,ff)
           enddo
          endif

* --------------------------------------------------------------------
*   Line by Line calculation: optical depth per layer per freq.
* --------------------------------------------------------------------
       IF (lambda.ne.3 .and. gas.ne.0)then
         call lbylshay(gasnum,1,1,tlay,play,cd,nlay,1,nlay,wn,nunwn,cof,
     .           ngas,xgeo,mw,f,s,a,e,nexp,nl,npt,ff,mtau)
       ENDIF
* --------------------------------------------------------------------
*   Continuum Absorption from Collision-Induced Dipole Moments part1
* --------------------------------------------------------------------
       if(gas.eq.0)then
         tmin = 40.
         do ii=1,nh2 
           th2(ii) = tmin + 2.*(ii)
         end do
         call h2hesubver(1,1,wn,nunwn,th2,nh2,absh,abshe,absch4)

         DO iw=1,nunwn
            do il=1,nh2
               ab1(il) = absh(il,iw)
               ab2(il) = abshe(il,iw)
            end do
            call tominterp(th2,ab1,nh2,t,bs1,nlay)
            call tominterp(th2,ab2,nh2,t,bs2,nlay)
            Do k=1,nlay
                tauc(xc,k,iw) = (col(k)-col(k-1)) * dlay(k) *
     &            (bs1(k)*(h2**2) + bs2(k)*(he*h2))
            Enddo
         ENDDO
        endif
       IF(lambda.eq.3 .and. gas.eq.0)then
          open(46,file='output/opacity_mirct'//str//'.out',
     &       form='unformatted',
     &       access='direct',recl=4*nunwn)
          do i=1,nlay
            write(46,rec=i+h) (real(tauc(xc,i,ij)),ij=1,nunwn)
          enddo
           close(46)
           h=h+nlay
       ELSEIF(lambda.ne.4 .and. gas.ne.0) then
          if(gas.eq.1)then
            if(lambda.eq.3)then
       opacityfile_ch4='output/opacity_mirch4_'//str//'.out'
         elseif(lambda.eq.5)then
       opacityfile_ch4='output/opacity_nirch4_geisa'//str//'.out'
            endif  
            ptfile_ch4='output/pt_mirch4.out'
          elseif(gas.eq.2)then
       opacityfile='output/opacity_mirc2h2'//str//'.out'
          elseif(gas.eq.3)then
       opacityfile='output/opacity_mirc2h6'//str//'.out'
          elseif(gas.eq.4)then
       opacityfile='output/opacity_mirc2h4'//str//'.out'
          elseif(gas.eq.5)then
       opacityfile='output/opacity_mirc3h8'//str//'.out'
          endif
          if(gas.eq.1)then
             open(45,file=opacityfile_ch4,form='unformatted',
     &          access='direct',recl=4*nunwn)
             if(xc.eq.0)then
       open(47,file='output/pt_ch4.out',form='unformatted',
     &             access='direct',recl=4*nlay)
                write(47,rec=1+hh) (real(p(ii)),ii=1,nlay)
                write(47,rec=2+hh) (real(t(ii)),ii=1,nlay)
                close(47)
                hh=hh+2
             endif
             do i=1,nlay
                write(45,rec=i+h) (real(mtau(1,1,i,ij)),ij=1,nunwn)
             enddo
             close(45)
             h=h+nlay
          elseif(gas.ne.0 .or. gas.ne.1)then
             open(45,file=opacityfile,form='unformatted',
     &          access='direct',recl=4*nunwn)
             do i=1,nlay
                write(45,rec=i+h) (real(mtau(1,1,i,ij)),ij=1,nunwn)
             enddo
             close(45)
             h=h+nlay
          endif
        ELSEIF(lambda.eq.4)then
             opacityfile='output/opacity_firct.out'
             ptfile='output/pt_firct.out'
             open(45,file=opacityfile,form='unformatted',
     &          access='direct',recl=4*nunwn)
             if(xc.eq.0)then
               open(47,file=ptfile,form='unformatted',
     &            access='direct',recl=4*nlay)
               write(47,rec=1+hh) (real(p(ii)),ii=1,nlay)
               write(47,rec=2+hh) (real(t(ii)),ii=1,nlay)
               close(47)
               hh=hh+2
             endif
             do i=1,nlay
               write(45,rec=i+h) (real(tauc(xc,i,ij)),ij=1,nunwn)
             enddo
             close(45)
             h=h+nlay
        ENDIF
 45    enddo
       ENDDO
       END

