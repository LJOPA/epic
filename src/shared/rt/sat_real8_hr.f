*-------------------------------------------------------------
*	Program started 03/01/2005 by SBH, using T's         | 
*		program as a starting point		     |
*-------------------------------------------------------------
*    [NOTES:] 
*-------------------------------------------------------------
* started: 05/30/05					     |
!--DISORT PARAMETERS/VARIABLES (using DISOTEST.f as template)    


************go to /tom   Need to teach code how to search pressure
***********and temperature at the same time, need to read in pres
***********and temp

	Implicit None
   
       INTEGER  nlyr,maxcly,maxmom,maxphi,maxulv,maxumu,ibcnd,
     &  nmom,numu,nstr,nphi,ntau,nnwn,nlay,pnum,tnum,lae,wae,
     &  nlt,lsbs,npres,ntemp
       Parameter (nlay=70,nnwn=943,maxphi=3,maxmom=299,npres=70,
     &  maxulv=nlay,maxumu=48,maxcly=nlay-1,pnum=20,tnum=19,lae=25,
     &  wae=191,nlt=10,lsbs=37,ntemp=10)
       integer totyear,last,nmol,nr,ntp,nspec,nwn,ngas,startpix
       Parameter (ntp=70,nr=ntp,nspec=1,nwn=16000,ngas=1)
       real day(24148)
       REAL*8 p_real(nlay),t_o_real(nlay),extab3(400000),
     &extab4(400000)
       REAL*8 r,d(nlay),wnf,lsubs(24148),fluxbasetot,
     &  psave2(10,nwn),cd(nlay),lat,value,decr,fluxbase,
     &  cmu1,cmu2,cmu3,angler,g,plandec(24148),sdis(24148),
     &  dtaucc5(10,nlay,1501),lsub,wn1n,angle,qh2,ryd,CpR,Cp,
     &  dtdt(maxcly),p(nlay),t_o(nlay),t(nlay),fluxbaseday(10),qhe,
     &  DF_big(maxcly,10),sumdfdp(maxcly),delp(maxcly),
     &  wnout(nwn),qtr_use(10),wnfn,wmidmax,tnews(nlay),pi2,
     &  dtsum2(10,nwn),sdfdp(maxcly),dec,req,rpole,pi,pi180,obl,scale, 
     &  wn1vlm,wnfvlm,wn1vhm,wnfvhm,wn1vlmx,wnfvlmx,pr(ntp),tr(ntp),
     &  fluxtotsum,dtau_sum(nwn),cosz(10),tfluxdiff(maxcly),
     &  fluxtottheat(nlay),psave(nwn),dtsum(nwn),fluxtot_vis3(nwn),
     &  wmidmin,wfirmin,wfirmax,daily_fnvu(maxcly),
     &  wn1uv,wnfuv,daily_fnvutot(maxcly),clouddiv,fluxtotcheck,
     &  sumdf2(nlay-1,nwn),satsec,ringfrac,daily_fm(maxcly),
     &  Qday,latr,glatr,sumdfDT(nlay),dtdt_lay(nlay),xm,
     &  play(nlay),playm1(maxcly),tlay(nlay),cloudgrey2(nlay),
     &  xmr1(nr),xmr2(nr),xmr3(nr),xmr4(nr),xmr5(nr),xmr6(nr),
     &  fluxtottcool(nlay),xmr7(nr),xmr8(nr),pjt(nlay),
     &  tjt(nlay),starttau,play_o(nlay),tlay_o(nlay),pm1(maxcly),
     &  sdfdpnuv(maxcly),sdfdpnuv_lay(maxcly),cloudgrey(nlay),
     &  dta(nlay),dlay(nlay),amfact,info_arr(10,pnum,tnum,1501),
     &  cos_psi,ff_ch4(1,nlay),ptaulay(npres),
     &  ff_c2h2(1,nlay),ff_c2h4(1,nlay),ff_c2h6(1,nlay),wnirwl(1501),
     &  einit(1816),einitc(1815),dtauc(nlay,nwn),totalfbeam,wnirw(1501),
     &  wnsol(1816),wnsolc(1816),binsize(1815),sumflux,dum,wnirwh(1501),
     &  cpsi,fluxtot(nwn),fluxtot_nir(nwn),fluxtottom(nlay),
     &  sfflux(nlay,nwn),fluxtot_uv(nwn),fluxtot_vis1(nwn),
     &  fflux(nlay),ff_c3h8(1,nlay),ttau(ntemp),ptau(npres),
     &  fluxtott(nlay,nwn),p_aero(wae,lae,nlay),lat_aero(wae,lae,nlay),
     &  lambda_aero(wae,lae,nlay),tau_aeroabs(wae,lae,nlay),
     &  lat_chk(25),pbar_aero(wae,lae,nlay),wno_aero(wae,lae,nlay),
     &  wno_aerol(wae),pbar_aerol(nlay),tau_aeroabsl(nlay,wae),
     &  err,dtau_aerol(nlay,nwn),dtau_aero,tau_aerol(nlay,nwn),
     &  tau_aeronew(nlay,nwn),wno_aerolr(wae),latxj(10),pxx(81),
     &  dfwnc(maxcly),dfwn(maxcly,nwn),dislim,bb,
     &  dfb(nlay),Fbinc(nlay),fluxtottv(nlay),den(81,nlt),
     &  delg(10),mix(81,7,nlt,lsbs),tau_aeroabslr(nlay,wae),
     &  dtau_sum2(10,nwn),fluxtottcoolwn(maxcly,nwn),
     &  resnir,res,press(81),temp(81,nlt),mixj(81,7,nlt,lsbs),
     &  fluxtot_vis2(nwn),tirw(tnum),pirw(pnum),
     &  wn2n,nn,ljsbs(lsbs),fflatc2h2(nlay,lsbs),fflatc2h6(nlay,lsbs),
     &  avecos_day,ffxc2h2(nlay,10,lsbs),fflatch4(nlay,lsbs),
     &  etauFb(10),Fb(nlay),avecos,ffxc2h6(nlay,10,lsbs),
     &  wmid,nir1,nir2,wmid1,wmid2,ffxch4(nlay,10,lsbs),
     &  dtauc3(nlay,nwn),dtauc_tom(nlay),fluxtot_niri(nwn),
     &  fluxtottheatwn(nlay,nwn),fluxtot_nirhi(nwn),fluxtot_nirh(nwn)
!--DISORT PARAMETERS/VARIABLES (using DISOTEST.f as template)    
       Real accur,alb,btemp,fbeam,FISOT,phi(maxphi),
     &   pmom(0:maxmom,maxcly),phi0,ssalb(maxcly),temper(0:maxcly),
     &   temis,ttemp,wvnmh,wvnml,umu(maxumu),umu0,utau(maxulv),
     &   rfldir(maxulv),rfldn(maxulv),flup(maxulv),dfdt(nlay),
     &   uavg(maxulv),uu(maxumu,maxulv,maxphi),albmed(maxumu),
     &   trnmed(maxumu),dtauc_dis(10,maxcly),wnout_lo(nwn),
     &   dtauc_t(nlay),wnout_hi(nwn),sumdf(maxcly),tinit(nlay)
       Character header*127
       Logical cloud,lamber,plank,onlyfl,prnt(5),usrang,usrtau,
     &tempread,juliehot,scattering,aerosol,bruno,irwinhalf,
     &bdiurn,juliepun2,obltrig,spectracold,spectrahot,flagsolar,
     &nir_irwinonly,write_ff,tausum,absfact,rates_write,seas_ch4,
     &ftot_wn,fluxtot_write,hit_geisa,hitran,nir_hitirwin,nir_hitgeisa,
     &kcheck,skipwrite,firsttime
       EXTERNAL disort
!--------------------------------------------------------------------
       Real*8 mixh6fact,mixh8fact,mixh2fact,mix10fact,
     & col(0:nlay),mix11fact,mixch4fact,mix13fact,mixnh3fact,pxmr1(nr),
     & pxmr2(nr),pxmr3(nr),pxmr4(nr),pxmr5(nr),pxmr6(nr),pxmr7(nr),
     & pxmr8(nr),wvnmlo(nwn),mix_lo,c2h2_jt(81),mixch4_jt(81),
     & wvnmhi(nwn),p_jt(81),d_jt(81),c2h6_jt(81),mixc2h2_jt(81),
     & mixc2h6_jt(81),p_jtr(81),mixc2h2_jtr(81),mixc2h6_jtr(81),
     & mix_hi,slope,mixch4_jtr(81),ch4_jt(81),bjt,tauabove(nlay,nwn)
       Integer qtrtime,k,qtr,i,j,ij,jj,h,kl,ll,kj,nunwn,
     & m,callme,stepsize,saturnday,l,latint,nyear,tilt,steps,bin,
     & ed,inc,latch,stepsizenew,today,tomorrow,count,numrun,
     & flag,startstep,endstep,hh,kk,q,s,satday,n,y,qq,nnr

       Character*1 sqt1,stgk1,latst,tilts
       Character*2 sqt,stgk2,latstri,tilts2
       Character*3 planet,stgkj,stgk3,latstr,tiltstr
       Character*4 userdir
       Character*5 stgday
       Character*40 char1,label
       Character*28 char2
       CHARACTER*34 workdir
       Character*132 ptfilech4,sfile


*------------------------------------------------------------------------
*                Logical statements may be altered here
*------------------------------------------------------------------------

* Flags for Hot-version, OR using different temp profile
      res=3000000.
      resnir=1.52e6

      skipwrite=.false.  !fast: writes out mix ratios as func. of season
                         !no radiative calcs.
      ftot_wn=.false.   !write F(nu)
      fluxtot_write=.false. !write out total flux(layer) for day (heat+cool), collapsed over wno.
      nir_hitgeisa = .false.  !GEISA/HITRAN NIR.
      hit_geisa = .false. !nir is hit+geisa
      hitran = .false.  !nir hitran only OR hitran+Irwin
      nir_hitirwin = .false.
      nir_irwinonly=.true. 
      kcheck=.true.
      flagsolar=.true.

*-------------------------------------------------------------------------
*                    Mole fraction/mixing
*-------------------------------------------------------------------------

      write_ff=.false.   !write out mix ratios
      tausum=.false.
      rates_write = .false.
      irwinhalf = .false. !false is NO k/2. FALSE if absfact = true!
      absfact = .false. !true = use Bruno's modifications
      juliepun2=.false.    !TMV
      seas_ch4=.true.     !tied to juliepun2, but with seas. dep. CH4
      juliehot = .false.  !julie's test temps must be true for nom cases
      tempread = .false.  !use previously derived temp prof
              !.false. uses Julie's to and po from lbyl to start 
              !.true. requires the file 'pt_start.out' in workdir
              !it is assumed this is a previously generated 82_*.out file
      bruno = .false. !use bruno constant mixing ratio;juliepun2 is false 

*--------------------------------------------------------------------------
*                                extra
*--------------------------------------------------------------------------

      cloud = .true. !Tommy's cloudgrey opacity
      aerosol = .false.  !aerosol absorption
      scattering = .false. !when false Tomsvisabsorp.f; true disort
      obltrig=.true. !true tilt=0., lsub=0.
      bdiurn=.false.     !true --> use Bruno's diurnal approx. qtr=1

!      bin = 2 !bin=2 no binning, bin=1 binning

      print*,'Using updated TP profile?',tempread

*--------------------------------------------------------------------------
*		Begin the number of runs loop
*--------------------------------------------------------------------------
*       print*,'We are in planetocentric coords'
*       print*,'Workdirectory (4 characters):'
*       read*,userdir 

	print*,'Please input the script file'
        read*,sfile
	open(27,file=sfile,status='old')
        read(27,*) numrun
        read(27,*) obltrig
        read(27,*) nyear
        read(27,*) userdir

      do nnr=1,numrun
	read(27,*) lat
       
*      lat=80.   ! latitude in degrees
*      nyear=4 !number of Saturn years to run (currently 4 is max can run)

      obl=26.7
      if (obltrig) then
        obl = 0.   ! obliquity in degrees
      endif
       
      print*,'obliquity = ', obl,'latitude = ',lat
      pi = acos(-1.)
      pi2=2.*pi
      pi180 = pi/180.
      req = 60268.    ! equatorial radius in km for Saturn
      rpole = 54364.  ! polar radius in km for Sature
c     Note that dec = obl*sin(xlsdeg*pi180) where xlsdeg is L_s in
c     degrees.  dec = 0 at spring, fall equinox, and dec = obl at summer
c     solstice, etc.  
      latint = int(lat)
      tilt = int(obl)

      do l=1,nwn
        do ll=1,nlay
          tau_aeronew(ll,l)=0.
        enddo
      enddo

*-------------------------------------------------------------------------

      do l=1,maxcly
        do i=1,nwn
          dfwn(l,i)=0. 
        enddo
        dfwnc(l)=0.
      enddo 
*      do l=1,nlay
*        do i=1,nwn
*          fluxtott(l,i)=0.
*        enddo
*      enddo       

!---------------------------DISORT Paramters-----------------------------

      nstr=4
      NLYR=maxcly
      nmom=nstr
      do i=1,maxcly
        ssalb(i) = 0.
        pmom(0,i)=1.
        do l=1,nmom
          pmom(l,i)=0.
        enddo
      enddo
      ntau=1            !number of user optical depths 
      utau(1) =0.
      usrtau = .false.
      usrang = .true.
      numu = 1
      umu(1) =1.
      nphi =  1
      phi(1) = 0.0
      ibcnd = 0
      phi0 = 0.0
      fisot = 0.0
      lamber = .true.
      alb = 0.
      temis = 0.0
***test by tom, I don't think plank needs to be true since in the
***near IR and vis Saturn is sooooo cold no thermal emission
      plank = .false.
      onlyfl = .true.
      accur = 0.0001
      do ij=1,5
        prnt(ij) = .false.
      enddo
      header = ''

*-----------------------------------------------------------------------------
*                   create exponential tables
*-----------------------------------------------------------------------------

	call exptable(extab3,extab4)

*------------------------------------------------------------------------------
*     Searched output from pscale to determine what layer = 10mbar level
*------------------------------------------------------------------------------

      xm=0.002226 !from sat_info.dat ![Kg/mol]
      qh2=0.877    !0.u63                   ![mix ratio]  
      qhe=0.118
      ryd=8.3144           ![universal gas constant --> J mol-1 K-1]
      CpR=3.0              ![Wallace (1980); 3.0 "normal" @140K]
      Cp=((CpR)*qh2+2.5*(1.-qh2))*ryd
      print*,'Cp',Cp
      open(31,file='other/var_script_shay',status='old')
      read(31,*) mixh2fact			!ratio's to multiply mixing ratio's by for testing
      read(31,*) mixch4fact
      read(31,*) mixh6fact
      read(31,*) mixh8fact
      read(31,*) amfact
      read(31,*) mix10fact
      read(31,*) mix11fact
      read(31,*) mix13fact
      read(31,*) mixnh3fact
      close(31)

*----------------------------------------------------------------------------
*         read in T and P values corresponding to opacity file inputs
*---------------------------------------------------------------------------

	open(31,file='/data/seasonaltom/code/opacity/data/tp.dat'
     &,status='old')
	do ij=1,ntemp
	read(31,*) ttau(ij)  !in increasing order (ie. 80,90,100...)
	end do
	do ij=1,npres
	read(31,*) ptaulay(ij)  !in increasing press order (ie. 1e-8,2e-8..)
************ptaulay is the pressures the opacities are calculated for!
	end do
	do ij=1,npres
	read(31,*) ptau(ij)  !in increasing press order (ie. 1e-8,2e-8..)
	end do
        close(31)
!-----------------------------------------------------------------------------
!                            Lat file strings
!----------------------------------------------------------------------------- 

      if (latint.lt.0.and.latint.ge.-9) then
        latch=abs(latint)
        write(latst,'(I1)'),latch
        latstr='m0'//latst
      elseif (latint.le.-10) then
        latch=abs(latint)
        write(latstri,'(I2)'),latch
        latstr='m'//latstri
      elseif (latint.le.9) then
        write(latst,'(I1)'),latint
        latstr='p0'//latst
      else
        write(latstri,'(I2)'),latint
        latstr='p'//latstri
      endif

       print*,'lat dir/ is:',latstr//'degrees'

!----------------------------------------------------------------------------
!                          obl file strings
!----------------------------------------------------------------------------  

      if (tilt.le.9) then
        write(tilts,'(I1)'),tilt
        tiltstr='t0'//tilts
      else
        write(tilts2,'(I2)'),tilt
        tiltstr='t'//tilts2
      endif

!----------------------------------------------------------------------------
!                         work directory
!----------------------------------------------------------------------------

      workdir='SAT/nsubs/LAT/'//tiltstr//'/'
     &   //latstr//'degrees/'//userdir//'/'
      print*,'Working dir is:'
      print*,workdir
      open(111,file=workdir//'params.inp',status='unknown')
      print*,workdir//'params.inp'
      write(111,*) nir_hitgeisa 
      write(111,*) hit_geisa 
      write(111,*) hitran 
      write(111,*) nir_hitirwin 
      write(111,*) cloud,'cloudgrey'
      write(111,*) aerosol,'aerosol'
      write(111,*) bruno,'bruno adbunances'
      write(111,*) scattering,'scattering'
      write(111,*) tempread,'new TP profile'
      write(111,*) workdir,'workdir'
      write(111,*) nyear,'nyear'
      write(111,*) obl,'obl'
      write(111,*) lat,'lat'
      write(111,*) obltrig,'obltrig'
      write(111,*) bdiurn,'bdiurn'
      write(111,*) juliepun2,'juliepun2'
      write(111,*) irwinhalf,'irwinhalf: f is no k/2'
      write(111,*) juliepun2,'adapt. seasonal mixing?'
      write(111,*) juliehot,'juliehot? (should usu. be yes)'
      write(111,*) absfact,'Using Bruno modified k-coeff? (absfact)'


!----------------------------------------------------------------------------
!      Open up temperature output file
!----------------------------------------------------------------------------

        open(82,file=workdir//'temp.out',
     &         form='unformatted',access='direct',recl=4*(nlay+1))

!--------------------------------------------------------------------------
!                   setup gravity at nominal and new latitudes
!--------------------------------------------------------------------------

      planet='SAT'
      call gravsub(planet,lat,g)     !new latitude gravity
      g=g/100.
      latr=lat*pi180
      glatr=atan((tan(latr))*((req/rpole)**2)) !new lat !see Icarus 147,405 (2000) def'n of planetographic - planetocentric

*------------------------------------------------------------------------
*                          Pressure Scale
*------------------------------------------------------------------------

*      p(1)    =       1.02E-7          ! pressure units [bars]
*      scale   =       0.2273 
*      call pscale(p(1),scale,p,nlay)   
*      do i=1,nlay
*        pjt(i)=p(i)
*      end do
*	print*,p(nlay)

*-------------------------------------------------------------------------------
*a slight change I'm reading in play from the opacity creation program
*it should be just the same as what is done above
*-------------------------------------------------------------------------------

	do i=1,nlay
	  p(i)=ptau(i)
	  play(i)=ptaulay(i)
	  pjt(i)=ptau(i)
	end do

        write(82,rec=1) (real(0.),(real(p(ll)),ll=1,nlay))
        close(82)
        count=0

*--------------------------------------------------------------------------
*                                Wavelength Range
*-------------------------------------------------------------------------

!FIR
      wfirmin=0.
      wfirmax=599.
!MidIR
      wmidmin=599.5
      wmidmax=1600.
!NIR
      wn1n=2000.                       !5 um
      wnfn=9500.                       !1.06 um
      wn2n=4790.
!VIS
      wnfvlmx=3.3289E4           !300.4nm, 0.3 um
      wn1vlmx=1.9231E4           !520.0
      wnfvhm=1.9227E4            !520.1
      wn1vhm=1.0050E4            !995.0
****there is an issue here
      wnfvlm=1.0048E4            !995.2
      wn1vlm=9523.81                   !1050.0, 1.05 um
!UV
      wn1uv=40000.0
      wnfuv=100000.0

*-----------------------------------------------------------------------------
*   Here I'm taking care of Irwin data:
*      1) irwin_rewrite.f: rewrite the Irwin data into a usable format
*           this is only done once for ENTIRE run
*      2) irwin.f (in nir8.f subroutine): locates appropriate values (4)
*           of t_irwin and p_irwin and calculates transmission
*      3) polin2.f (in nir8.f): interpolation in 2-D for 4 values of t_irwin 
*           and p_irwin to find our desired transmission at OUR tlay and play
*           value.
*-----------------------------------------------------------------------------

      if (nir_hitirwin.or.nir_irwinonly) then  !modified to only use 4790 - 9500 cm-1
        call irwin_rewrite(pnum,tnum,info_arr,delg,pirw,tirw,wnirw,
     &wnirwl,wnirwh)
      endif

      
*-----------------------------------------------------------------------------
*                       Initializing initial T and P data
*-----------------------------------------------------------------------------

      open(13,file='other/mixsat.dat',status='old')
      do i=1,70
        read(13,*) p_real(71-i),t_o_real(71-i)
	p_real(71-i)=p_real(71-i)/1000.
      end do
      close(13)
      call tominterp(p_real,t_o_real,nlay,p,t_o,nlay)     !get T_o at level

c        lbylb writes out the LEVEL P,T', not the layers (tlay,play)
c        opacities are solved for based on tlay and play

      if (tempread) then !on exisiting 82_*.out file
        open(13,file=workdir//'hg80_48254_360.out',form='unformatted',
     &        access='direct',recl=4*nlay)

!	 write(111,*) 'pt_0903_48334.out','pt file'

        read(13,rec=3) (tinit(i),i=1,nlay)
        close(13)
        close(111)
      else
        do i=1,nlay
          tinit(i)=t_o(i)
        end do
      endif 
      do l=1,nlay
        t(l)=t_o(l)            !level
	tnews(l)=tinit(l)
      enddo
      play(1)=0.5*p(1)        !layer
      play_o(1)=play(1)
      do i=2,nlay
        play(i)=(p(i-1)+(p(i)))/2.
        play_o(i)=play(i)
      enddo
      open(12,file=workdir//'delp.out',form='unformatted',
     &  access='direct',recl=4*maxcly)
      do l=1,maxcly
        delp(l)=(p(l+1)-p(l))
      enddo
      write(12,rec=1) (real(delp(l)),l=1,maxcly)
      write(12,rec=2) (real(play(l)),l=1,maxcly)
      close(12)
      call tominterp(p,t_o,nlay,play_o,tlay_o,nlay)     !get T_o at layer

*---------------------------------------------------------------------------
*                            Mixing Ratios
*--------------------------------------------------------------------------

      open(15,file='other/mixsat.dat',status='old')
 181      format(10e10.3)
 188      format(10e12.5)   !format for new mix files
      do i=1,ntp        ! read T-P profile
        read(15,181) pr(ntp+1-i),tr(ntp+1-i),xmr1(ntp+1-i),
     &       xmr2(ntp+1-i),xmr3(ntp+1-i),xmr4(ntp+1-i),
     &       xmr5(ntp+1-i),xmr6(ntp+1-i),xmr7(ntp+1-i),
     &       xmr8(ntp+1-i)

        pr(ntp+1-i)=pr(ntp+1-i)/1000.
        tr(ntp+1-i)=tr(ntp+1-i)-50.
        pxmr1(ntp+1-i)=pr(ntp+1-i)
        pxmr2(ntp+1-i)=pr(ntp+1-i)
        pxmr3(ntp+1-i)=pr(ntp+1-i)
        pxmr4(ntp+1-i)=pr(ntp+1-i)
        pxmr5(ntp+1-i)=pr(ntp+1-i)
        pxmr6(ntp+1-i)=pr(ntp+1-i)
        pxmr7(ntp+1-i)=pr(ntp+1-i)
        pxmr8(ntp+1-i)=pr(ntp+1-i)
        xmr1(ntp+1-i)=xmr1(ntp+1-i)
        xmr2(ntp+1-i)=xmr2(ntp+1-i)
        xmr3(ntp+1-i)=xmr3(ntp+1-i)
        xmr4(ntp+1-i)=xmr4(ntp+1-i)
        xmr5(ntp+1-i)=xmr5(ntp+1-i)
        xmr6(ntp+1-i)=xmr6(ntp+1-i)
        xmr7(ntp+1-i)=xmr7(ntp+1-i)
        xmr8(ntp+1-i)=xmr8(ntp+1-i)
      enddo
      close(15)

      call mixr(nlay,1,1,nr,ntp,pxmr1,xmr1,play_o,ff_ch4)
      call mixr(nlay,1,1,nr,ntp,pxmr2,xmr2,play_o,ff_c2h2)
      call mixr(nlay,1,1,nr,ntp,pxmr3,xmr3,play_o,ff_c2h6)
      call mixr(nlay,1,1,nr,ntp,pxmr6,xmr6,play_o,ff_c2h4)
      call mixr(nlay,1,1,nr,ntp,pxmr4,xmr4,play_o,ff_c3h8)

*----------------------------------------------------------------------------
*                  Julie's Lsubs and Ephem data
*----------------------------------------------------------------------------

      if (obltrig) then
        do l=1,24148
          lsubs(l)=0.   !always same L_s and dist. from Sun
          sdis(l)=9.5
          plandec(l)=0.
        enddo
      else
        open(87,file='nsubs/lsubs.inp',status='old')
        do l=1,24148
          read(87,109) plandec(l),sdis(l),lsubs(l),day(l)
        enddo
        close(87) 
      endif

!plandec in pcentric,sdis in AU, lsubs in deg, and day in saturn days.
!L_s=90 is summer. L_s=0 is spring, L_s=180 is fall, L_s=270 is winter in the
!Northern hemisphere.

 109  format(f11.6,1x,f10.4,1x,f11.4,1x,f9.3)

*--------------------------------------------------------------------------
*                     Use a modified Julie TP-profile
*--------------------------------------------------------------------------

      if (juliehot) then 
        open(14,file='SAT/nsubs/goodavekin24sls270short2.pun',
     &  status='old')
        read(14,*)  char1
        read(14,'(8(1pe10.3))') (p_jt(i), i = 1,81)
        read(14,*) char1
        read(14,'(8(1pe10.3))') (d_jt(i), i = 1,81)
        read(14,*) char1
        read(14,'(8(1pe10.3))') (ch4_jt(i), i = 1,81)
        read(14,*) char1
        read(14,'(8(1pe10.3))') (c2h2_jt(i), i = 1,81)
        read(14,*) char1
        read(14,'(8(1pe10.3))') (c2h6_jt(i), i = 1,81)
        close(14)
        do l=1,81
          mixc2h2_jt(l)=c2h2_jt(l)/d_jt(l)  !concentration(cm-3)/density(cm-3)
          mixc2h6_jt(l)=c2h6_jt(l)/d_jt(l)
          mixch4_jt(l)=ch4_jt(l)/d_jt(l)
        enddo   
        do l=1,81
           p_jtr(82-l)=p_jt(l)/1000.
           mixc2h2_jtr(82-l)=mixc2h2_jt(l)
           mixc2h6_jtr(82-l)=mixc2h6_jt(l)
           mixch4_jtr(82-l)=mixch4_jt(l)
        enddo
        do l=1,nlay
           call locate(p_jtr,81,pjt(l),j)
           mix_lo=mixc2h2_jtr(j)
           mix_hi=mixc2h2_jtr(j+1)
           slope=(log(mix_lo)-log(mix_hi))/
     &     (log(p_jtr(j))-log(p_jtr(j+1)))
           bjt=(log(mix_lo)-slope*log(p_jtr(j)))
           ff_c2h2(1,l)=exp(slope*log(pjt(l))+bjt)
           mix_lo=mixc2h6_jtr(j)
           mix_hi=mixc2h6_jtr(j+1)
           slope=(log(mix_lo)-log(mix_hi))/
     &     (log(p_jtr(j))-log(p_jtr(j+1)))
           bjt=(log(mix_lo)-slope*log(p_jtr(j)))
           ff_c2h6(1,l)=exp(slope*log(pjt(l))+bjt)
           mix_lo=mixch4_jtr(j)
           mix_hi=mixch4_jtr(j+1)
           slope=(log(mix_lo)-log(mix_hi))/
     &     (log(p_jtr(j))-log(p_jtr(j+1)))
           bjt=(log(mix_lo)-slope*log(p_jtr(j)))
           ff_ch4(1,l)=exp(slope*log(pjt(l))+bjt)
        enddo
        If (bruno) then
          do l=1,nlay
            ff_c2h2(1,l)=(0.00000021)*qh2/0.88
            ff_c2h6(1,l)=0.000003*qh2/0.88
!wrong wasn't brunos model using something like 3.5x10-3 for Ch4
            ff_ch4(1,l)=(4.5E-3)*qh2/0.88  !Bruno fig 1 5*Solar Lambert value
          enddo   
        endif
	do l=1,nlay
	  ff_ch4(1,l)=ff_ch4(1,l)*mixch4fact
          ff_c2h2(1,l)=ff_c2h2(1,l)*mixh2fact
          ff_c2h6(1,l)=ff_c2h6(1,l)*mixh6fact
          ff_c2h4(1,l)=ff_c2h4(1,l)*mix11fact
	enddo 

!write out the mixing ratios used if we are not using temporal mix variations
	

        open(14,file='nsubs/mix_julie_new.dat',form='unformatted',
     &  access='direct',recl=4*nlay)
        write(14,rec=1) (real(ff_ch4(1,l)),l=1,nlay)
        write(14,rec=2) (real(ff_c2h2(1,l)),l=1,nlay)
        write(14,rec=3) (real(ff_c2h6(1,l)),l=1,nlay)
        write(14,rec=4) (real(ff_c3h8(1,l)),l=1,nlay)
        write(14,rec=5) (real(ff_c2h6(1,l)),l=1,nlay)
        close(14)

        if (juliepun2) then
          j=0
          open(14,file='SAT/nsubs/kin2mbar1e10_shay.pun',
     &    status='old')
          last = 81
          nmol = 7
          nn=0.
          do l=1,lsbs
             ljsbs(l)=0.+nn
             nn=nn+10.
          enddo
          nn=0.
          do l=1,nlt
            latxj(l)=-72.+nn
            nn=nn+16.
          enddo
1155    format(a)
          do l=1,lsbs
            do k=1,nlt
              read(14,1155)  label   ! label is character*40
              read(14,1155)  char2   ! char2 is character*28
              read(14,'(8(1pe10.3))') (press(i), i = 1,last)  ! last = 81 actually this is altitude, but we don't use altitude so I read it into pressure then read in the actual pressure next.
              read(14,1155)  char2
!!wrong check are all the pressures at each latitude equal?
              read(14,'(8(1pe10.3))') (press(i), i = 1,last)
              read(14,1155)  char2
              read(14,'(8(1pe10.3))') (temp(i,k), i = 1,last)
              read(14,1155)  char2
              read(14,'(8(1pe10.3))') (den(i,k), i = 1,last)
              do j=1,7
                read(14,1155)  char2
                read(14,'(8(1pe10.3))') (mix(i,j,k,l),i = 1,last)
              enddo
            enddo
          enddo
          close(14)
          do l=1,lsbs
            do k=1,nlt
              do j=1,7
                do i = 1,last
                  mixj(82-i,j,k,l)=mix(i,j,k,l)/den(i,k)
                enddo
              enddo
            enddo
          enddo
          open(14,file=workdir//'jpun2_p.out',
     &     form='unformatted', access='direct',recl=4*81)
          open(15,file=workdir//'jp2_pmixc2h2.out',
     &     form='unformatted', access='direct',recl=4*10)
          open(16,file=workdir//'jp2_pmixc2h6.out',
     &     form='unformatted', access='direct',recl=4*10)
          open(17,file=workdir//'jp2_pmixch4.out',
     &     form='unformatted', access='direct',recl=4*10)
          do l=1,81
            pxx(82-l)=press(l)/1000.
          enddo
          j=0
          do l=1,lsbs
            do k=1,nlt
              do i=1,nlay
                call locate(pxx,81,pjt(i),j)
                if(seas_ch4)then
                  slope=(log(mixj(j,3,k,l))-log(mixj(j+1,3,k,l)))/    !3=ch4
     &                (log(pxx(j))-log(pxx(j+1)))
                  bjt=log(mixj(j,3,k,l))-slope*log(pxx(j))
                  ffxch4(i,k,l)=exp(slope*log(pjt(i))+bjt)
                endif
                slope=(log(mixj(j,4,k,l))-log(mixj(j+1,4,k,l)))/    !4=c2h2
     &                (log(pxx(j))-log(pxx(j+1)))
                bjt=log(mixj(j,4,k,l))-slope*log(pxx(j))
                ffxc2h2(i,k,l)=exp(slope*log(pjt(i))+bjt)
                slope=(log(mixj(j,5,k,l))-log(mixj(j+1,5,k,l)))/    !5=c2h6
     &                (log(pxx(j))-log(pxx(j+1)))
                bjt=log(mixj(j,5,k,l))-slope*log(pxx(j))
                ffxc2h6(i,k,l)=exp(slope*log(pjt(i))+bjt)
              enddo
            enddo
          enddo
          h=0
          do l=1,nlay
            write(15,rec=l+h) (real(ffxc2h2(l,k,1)),k=1,10)
            if(seas_ch4)then
              write(17,rec=l+h) (real(ffxch4(l,k,1)),k=1,10)
            endif
          enddo
          close(15)
          close(17)
          h=h+nlay
          hh=0
          do l=1,nlay
            write(16,rec=l+hh) (real(ffxc2h6(l,k,1)),k=1,10)
          enddo
          close(16)
          hh=hh+nlay    
          q=0
          if (latxj(1).gt.lat) then
            q=1
	  else
            call locate(latxj,10,lat,q)  !find bracketting Lat's
          endif
          if (latxj(q)*latxj(q+1).eq.0.or. 
     &    (abs(latxj(q+1))-abs(latxj(q))).eq.0) then
            do l=1,nlay
              do jj=1,lsbs
                if(seas_ch4)then
                  slope=(ffxch4(l,q,jj)-ffxch4(l,q+1,jj))/
     &           (latxj(q)-latxj(q+1))
                  bjt=ffxch4(l,q,jj)-slope*latxj(q)
                  fflatch4(l,jj)=slope*lat+bjt
                endif
                slope=(ffxc2h2(l,q,jj)-ffxc2h2(l,q+1,jj))/
     &           (latxj(q)-latxj(q+1))
                bjt=ffxc2h2(l,q,jj)-slope*latxj(q)
                fflatc2h2(l,jj)=slope*lat+bjt
                slope=(ffxc2h6(l,q,jj)-ffxc2h6(l,q+1,jj))/
     &        (latxj(q)-latxj(q+1))
                bjt=ffxc2h6(l,q,jj)-slope*latxj(q)
                fflatc2h6(l,jj)=slope*lat+bjt
              enddo
            enddo
          else 
            do l=1,nlay
              do jj=1,lsbs
                if(seas_ch4)then
!wrong check this absolute value idea? I don't think it is right
                  slope=(log(ffxch4(l,q,jj))-log(ffxch4(l,q+1,jj)))/
     &              (log(abs(latxj(q)))-log(abs(latxj(q+1))))
                  bjt=log(ffxch4(l,q,jj))-slope*log(abs(latxj(q)))
                  fflatch4(l,jj)=exp(slope*log(abs(lat))+bjt)*mixch4fact
                endif
                slope=(log(ffxc2h2(l,q,jj))-log(ffxc2h2(l,q+1,jj)))/
     &              (log(abs(latxj(q)))-log(abs(latxj(q+1))))
                bjt=log(ffxc2h2(l,q,jj))-(slope)*log(abs(latxj(q)))
                fflatc2h2(l,jj)=exp(slope*log(abs(lat))+bjt)*mixh2fact
                slope=(log(ffxc2h6(l,q,jj))-log(ffxc2h6(l,q+1,jj)))/
     &          (log(abs(latxj(q)))-log(abs(latxj(q+1))))
                bjt=(log(ffxc2h6(l,q,jj)))-slope*log(abs(latxj(q)))
                fflatc2h6(l,jj)=exp(slope*log(abs(lat))+bjt)*mixh6fact
              enddo
            enddo
          ENDIF
          write(14,rec=1) (real(pxx(l)),l=1,nlay)
          close(14)
        endif  !end juliepun2
 187      format(10e10.3)
 189      format(10e12.5)   !format for new mix files
      endif  !juliehot
          
!#########################-----------------------------!! AEROSOL tau's

      if (aerosol) then
        open(11,file='SAT/nsubs/testtauabs.out',status='old')
 173       format(1p4e12.4)
        do l = 1, 191       !um wavelength
          do n = 1, 25        !deg latitude
            do ll = 1, 70     !mbar (mine?) pressure (J claims mbar)
              read(11,173) lambda_aero(l,n,ll),lat_aero(l,n,ll),
     & p_aero(l,n,ll),tau_aeroabs(l,n,ll)
              wno_aero(l,n,ll)=1./(lambda_aero(l,n,ll)*0.0001) !invcm
              pbar_aero(l,n,ll)=p_aero(l,n,ll)*0.001   !bar
            enddo
          enddo
        enddo
!!check to make sure lat exists in file. might want to later include interp.
!!  scheme to be able to use inbetween latitudes.
        do ll=1,nlay
          do n=1,25
            lat_chk(n)=lat_aero(1,n,ll)
          enddo
        enddo
        call locate(lat_chk,25,lat,j)
        if (lat.eq.0.) then
          j=25
        endif
        print*,'lat:',lat,'lat_aero:',lat_chk(j),j
        if(j.eq.0 .and. lat_chk(j).ne.lat)then
          print*,'Error: Need to interpolate lat in aerosol file'
          stop
        endif
        do l = 1, 191       !um wavelength
!         do n = 1, 25        !deg latitude
          do ll = 1, nlay
            if(lat_aero(l,j,ll).eq.lat_chk(j)) then
              wno_aerol(l)=wno_aero(l,j,ll)
              tau_aeroabsl(ll,l)=tau_aeroabs(l,j,ll)
              pbar_aerol(ll)=pbar_aero(l,j,ll)
            endif
          enddo
!         enddo
        enddo
        do l=0,190
          wno_aerolr(l+1)=wno_aerol(191-l)
          do ll = 1, nlay
            tau_aeroabslr(ll,l+1)=tau_aeroabsl(ll,191-l)
          enddo
        enddo
      endif 

!##########################__________________________________________

      h=0
      satsec=3600.*10.66    !secs in 1 satday
      do l=1,maxcly
        dtdt(l)=0.
      enddo
      do l=1,nlay
        cloudgrey(l)=0.
        cloudgrey2(l)=0.
      enddo
     
      tomorrow=1
      stepsizenew=1    !nominal stepsize  (1 day)
      totyear=nyear*24148

*-----------------------------------------------------------------------------
*                          BIG LOOP START
*-----------------------------------------------------------------------------

 989  CONTINUE       !AND now it is TODAY 
      print*,'step',stepsizenew
      print*,'tomorrow',tomorrow
      today=tomorrow  
      stepsize=stepsizenew
      do saturnday=today,totyear,stepsize    !totyear,stepsize
        print*,'runtime',totyear,'timeleft',totyear-saturnday
        firsttime=.true.
*        do l=1,maxcly
*          fluxtottheat(l)=0.
*          fluxtottcool(l)=0.
*          do ll=1,nwn 
*            fluxtottheatwn(l,ll)=0.
*            fluxtottcoolwn(l,ll)=0.
*          enddo
*        enddo

	call numtostring5(saturnday,stgday)

*-----------------------------------------------------------------

        ed=int(saturnday/24148)   !number of file repetitions
        satday=saturnday-ed*24148

!wrong could speed things up by creating decr outside of the loop and making an array
        dec=plandec(satday)
        decr=dec*pi180
        lsub=lsubs(satday)

!wrong could comment out the print below
        print*,'dec',dec,'lsubs',lsubs(satday)

        if (juliepun2)then
          jj=0
        !get mixx* to be only a function of layer
        ! interp lat and lsubs 
          call locate(ljsbs,lsbs,lsub,jj)  !bracketing L_s's
          do l=1,nlay
            if(seas_ch4)then
              slope=(fflatch4(l,jj)-fflatch4(l,jj+1))/
     &        (ljsbs(jj)-ljsbs(jj+1))
              bjt=fflatch4(l,jj)-slope*ljsbs(jj)
              ff_ch4(1,l)=(slope*lsub+bjt)
            endif
            slope=(fflatc2h2(l,jj)-fflatc2h2(l,jj+1))/
     &        (ljsbs(jj)-ljsbs(jj+1))
            bjt=fflatc2h2(l,jj)-slope*ljsbs(jj)
            ff_c2h2(1,l)=(slope*lsub+bjt)
            slope=(fflatc2h6(l,jj)-fflatc2h6(l,jj+1))/
     &        (ljsbs(jj)-ljsbs(jj+1))
            bjt=fflatc2h6(l,jj)-slope*ljsbs(jj)
            ff_c2h6(1,l)=(slope*lsub+bjt)
          enddo
          if (write_ff) then
            open(14,file=workdir//'ff_lat0_'//stgday//'.out',
     &     form='unformatted', access='direct',recl=4*nlay)
            write(14,rec=1) (real(ff_c2h2(1,l)),l=1,nlay)
            write(14,rec=2) (real(ff_c2h6(1,l)),l=1,nlay)
            write(14,rec=3) (real(ff_ch4(1,l)),l=1,nlay)
            close(14)
          endif
        endif
        if (skipwrite) then
!wrong check this carefully she had tomorrow statement commented out and after the 
!go to statment
          tomorrow = today+50
          print*,'wrote mxing ratio files...'
          go to 133
        endif

*-------------------------------------------------------------------------

	fluxbasetot=0.
	do l=1,10
	  fluxbaseday(l)=0.
	end do
        do l=1,maxcly
          sumdfdp(l)=0.          !saturnday step dfdp
        enddo
         wnf=0.
!        wnf=8024.77447349719     !1.1 micron
!        wnf=6501.41910541985     !1.3 micron
!         wnf=5267.24464642784   !1.7 micron
!        wnf=2494.74351330164    !3 micron
!         wnf=3493.85903309036   !2.3 micron
!        wnf=1305.08545733879
!         wnf=2000. !600. !wmidmin !0.  !wn1uv     !wn1vlm
        call tominterp(p,tnews,nlay,play,tlay,nlay)
        col(0) = 0.
        do j=1,nlay
          call column(p(j),tnews(j),xm,g,col(j),d(j))  !using new latitude 'g'
          cd(j) = col(j) - col(j-1)
        enddo 
        call tominterp(p,d,nlay,play_o,dlay,nlay)

        if(tlay(1).gt.300.)then
          tlay(1)=299.
	elseif (tlay(1).lt.75) then
          tlay(1)=76.
        endif

        if (cloud) then
	  do l=61,nlay
	    cloudgrey(l)=cd(l)
 	  end do
          cloudgrey2(60)=0.   !all are set to zero earlier outside seasonal loop
	  do l=60,nlay
	    cloudgrey2(l)=cloudgrey2(l-1)+cloudgrey(l)
	  end do
	  clouddiv=cloudgrey2(67)        !0.334 bar, tau=1
	  do l=1,nlay
	    cloudgrey2(l)=cloudgrey2(l)/clouddiv
	  end do
        endif

        if (nir_hitirwin) then  
          startstep=0
          endstep=260
        elseif (nir_hitgeisa) then
          startstep=0
          endstep=325
        elseif (nir_irwinonly) then
           startstep=0
           endstep=177
        else
          print*,'Error in var startstep: no logical switch defined'
          stop
        endif

        do kj=0,189 !startstep,endstep   
!(a) 0 FIR, 1-172 MID, 173-321 NIR, 322-324 VIS, 325 UV   !GEISA NIR
!(b) 0 FIR, 1-172 MID, 173-255 NIR, 256 NIR IRWIN, 257-259 VIS, 260 UV
!0 FIR; 1-184 MIR; 185 NIR; 186-188 VIS; 189 UV

*          print*,'kj = ', KJ

          if (kj.le.9) then
            write(stgk1,'(I1)'),kj
            stgkj='00'//stgk1
          elseif (kj.lt.100.and.kj.ge.10) then
            write(stgk2,'(I2)'),kj
            stgkj='0'//stgk2
          else
            write(stgk3,'(I3)'),kj
            stgkj=stgk3
          endif

*------------------------------------------------------------------------

          do l=1,maxcly
            daily_fm(l)=0.
            daily_fnvutot(l)=0.
*            sdfdp(l)=0.  !no need for this since it is always set = to something
          enddo
*          print*,'wnf',wnf,'fluxbasetot',fluxbasetot
!wrong work on making this iff as fast as possible
          if (wnf.le.wfirmax) then
            nunwn=599
            call fir(nlay,nwn,tlay,wvnmhi,wvnmlo,wnout,dtauc,
     &              nunwn,qh2,qhe,ff_ch4,cd,dlay)
            umu0=1.        !1=cos(0.)
            bin=2
            qtr=1
            callme=1
             s=1
*            print*,'fir'
          elseif (wnf.le.wmidmax) then
            nunwn=nwn
            call midtexpb(nlay,tlay,play,ptau,npres,kj,wvnmhi,
     &wvnmlo,wnout,dtauc,nwn,nunwn,cd,ff_ch4,ff_c2h2,ff_c2h6,qh2,
     &qhe,dlay)
            umu0=1.       !1=cos(0.)
            bin=1
            qtr=1
            callme=1
            s=1
*            print*,'mid'
          elseif (wnf.le.wnfn) then
            qtr=10
!            if (nir_hitgeisa) then  !HITRAN + GEISA ONLY
!              s=1
!              bin=2  !2=NO BINNING! --> binning in NIR not great, 
!                    !requires very small threshold 'value' to reproduce
!                    !unbinned df/dp with low %error (value = 1e-6)
!              nunwn=nwn
!!!!wrong, this is just like subroutine mid
!              call nirh(xc,nlay,tlay,kj,wnf,wvnmhi,wvnmlo,wnout,dtauc,
!     &             tlay_o,lat_ratio,lat_mixch4,
!     &             nwn,nunwn,hit_geisa,hitran,resnir)
!            elseif (nir_hitirwin) then !HITRAN + IRWIN
!              if (wnf.lt.wn2n) then
!                s=1
!                bin=2
!                nunwn=nwn
!                call nirh(xc,nlay,tlay,kj,wnf,wvnmhi,wvnmlo,wnout,dtauc,
!     &             tlay_o,lat_ratio,lat_mixch4,
!     &             nwn,nunwn,hit_geisa,hitran,resnir)
!              elseif (wnf.ge.wn2n.and.wnf.le.wnfn) then   !Irwin
!                print*,'Irwin params', wnf,kj
!                s=10
!                bin=2
!                nunwn=nnwn
!                print*,nunwn,'nunwn'
!                call nir(wnout,nlay,nwn,tlay,wnf,wvnmhi,wvnmlo,
!     &             dtauc,workdir,play_o,pnum,tnum,
!     &             info_arr,nunwn,ff_ch4,dlay,cd,irwinhalf,
!     &             absfact,kcheck,dtaucc5,delg,tirw,pirw)
!              endif
!            elseif (nir_irwinonly) then
              bin=2
              nunwn=1501
              s=10
              do l=1,1501
                wnout(l)=wnirw(l)
	        wvnmhi(l)=wnirwh(l)
                wvnmlo(l)=wnirwl(l)
	      end do
              call nir_IR(nlay,nwn,tlay,play_o,pnum,tnum,
     &             info_arr,nunwn,ff_ch4,cd,
     &             dtaucc5,tirw,pirw)
!            endif   !endif logical
            callme=2
*            print*,'near-ir'
          elseif (wnf.le.wnfvlm) then
            qtr=10
            nunwn=138
            bin=2
            s=1
            call vis1(nlay,nwn,tlay,wnf,wvnmhi,wvnmlo,wnout,
     &             dtauc,ff_ch4,cd,nunwn,workdir)
            callme=2 
*            print*,'vis1'
          elseif (wnf.le.wnfvhm) then
            qtr=10
            nunwn=4750
            bin=2
            call vis2(nlay,nwn,tlay,wnf,wvnmhi,wvnmlo,wnout,
     &             dtauc,ff_ch4,cd,nunwn,workdir)
            callme=2
            s=1
*            print*,'vis2'
          elseif (wnf.le.wnfvlmx) then
            bin=2
            qtr=10
            nunwn=550
            call vis3(nlay,nwn,tlay,wnf,wvnmhi,wvnmlo,wnout,
     &             dtauc,ff_ch4,cd,nunwn,workdir)
            callme=2
            s=1
*            print*,'vis3'
          else
            qtr=10
            nunwn=36
            bin=2
            call uv(nlay,nwn,tlay,wnf,wvnmhi,wvnmlo,wnout,
     &          dtauc,workdir,cd,nunwn,ff_c2h4,ff_c2h6,ff_c2h2,
     &          ff_ch4)
            callme=2
            s=1
*            print*,'uv'
          endif

*------------------------------------------------------------------------
*add in cloud opacity if needed
*------------------------------------------------------------------------
          if (cloud.and.callme.eq.2) then
            if (s.eq.1) then
              do l=1,nunwn
                do ll=1,nlay
                  dtauc(ll,l)=dtauc(ll,l)+cloudgrey2(ll)
                end do
              end do
            else
              do k=1,10
                do ll=1,nlay
                  do l=1,nunwn
                    dtaucc5(k,ll,l)=dtaucc5(k,ll,l)+cloudgrey2(ll)
                  enddo
                enddo
              enddo
            end if
          end if

*--------------------------------------------------------------------------
*                    Sum up to find unit optical depth
*--------------------------------------------------------------------------

          if (tausum) then
            if (wnf.ge.wmidmin.and.wnf.le.wmidmax) then
              open(14,file=workdir//'tausum2_'//stgday//'_'//stgkj//
     &        '.out', form='unformatted', access='direct',
     &        recl=4*nunwn)
              do ll=1,nunwn 
                dtau_sum(ll)=0.
              enddo
              flag=0
              do ll=1,nunwn 
                do l=1,nlay-1
                  dtau_sum(ll)=dtauc(l,ll)+dtau_sum(ll)
                  if (dtau_sum(ll).lt.1.) then
                    dtsum(ll) = dtau_sum(ll)
                    psave(ll) = play(l)
                  elseif(dtau_sum(ll) .ge. 1. .and. flag.eq.0)then
                    dtsum(ll) = dtau_sum(ll)
                    psave(ll) = play(l)
                    flag=1
                  elseif(dtau_sum(ll) .ge. 1.0 .and. flag.ne.0)then
                    goto 111
                  endif
!                 dta(l)=dtau_sum(l+1,ll)+dta(l)
!                 enddo
!                 print*,dta(l)
 111              continue
                enddo
              enddo
              !write(14,rec=1) (real(dtsum(l)),l=1,nunwn)
              write(14,rec=1) (real(psave(l)),l=1,nunwn)
              close(14)
            else 
              print*,'here'
              do ll=1,nunwn
                do j=1,10
!                   dtau_sum(l,ll)=0.
                  dtau_sum2(j,ll)=0.
                enddo
              enddo
              flag=0
              do ll=1,nunwn
                do j=1,s
                  do l=1,maxcly
                    if (s.ne.10) then
                      dtau_sum2(j,ll)=dtauc(l,ll)+dtau_sum2(j,ll)
                    else
                      dtau_sum2(j,ll)=dtaucc5(j,l,ll)+dtau_sum2(j,ll)
                    endif
                    if (dtau_sum2(j,ll).lt.1.) then
                      dtsum2(j,ll) = dtau_sum2(j,ll)
                      psave2(j,ll) = play(l)
                    elseif (dtau_sum2(j,ll).ge.1..and.flag.eq.0) then
                      dtsum2(j,ll) = dtau_sum2(j,ll)
                      psave2(j,ll) = play(l)
                      flag=1
                    elseif (dtau_sum2(j,ll).ge.1..and.flag.ne.0) then
                      goto 112
                    endif
  112               continue
                  enddo 
                enddo
              enddo
              print*,'here2'
              open(14,file=workdir//'psave2_'//stgday//'_'//stgkj//
     &            '.out', form='unformatted', access='direct',
     &            recl=4*nunwn)      
              do j=1,10
                write(14,rec=j) (real(psave2(j,l)),l=1,nunwn)
              enddo
              close(14)
!               open(13,file=workdir//'tausum1_'//stgday//'_'//stgkj//
!     &            '.out', form='unformatted', access='direct',
!     &            recl=4*nunwn)
!               do ll=1,maxcly
!               write(13,rec=ll) (real(dtau_sum(ll,l)),l=1,nunwn)
!               enddo
!               close(13)
            endif
          endif

*---------------------------------------------------------------------------
*                              Solar Flux
*---------------------------------------------------------------------------

         if (flagsolar.and.callme.eq.2) then
           sumflux=0.
           open(66,file='nsubs/2005tsi2.txt',
     &           status='old')
 300       format(4x,f16.12,4x,f15.11,1x,f20.16)
           do l=1,1816
             read(66,300) wnsol(l),einit(l),dum
             wnsolc(1817-l)=1./(wnsol(l)*0.0001)  !um --> cm^1
           end do
           close(66)
           do l=1,1815
             einitc(1816-l)=einit(l+1)*(wnsol(l+1)-wnsol(l))
           enddo
           open(12,file=workdir//'sunflux_invcm.dat',
     &           form='unformatted',access='direct',recl=4*1816)
           do l=1,1815
             binsize(l)=wnsolc(l+1)-wnsolc(l)
           enddo
           inc=1
           do l=1,nunwn
             fluxtot(l)=0.
             do jj=inc,1815
               if (wvnmlo(l).gt.wnsolc(jj).and.wvnmlo(l).le.
     &          wnsolc(jj+1)) then
                 if (wvnmhi(l).le.wnsolc(jj+1)) then
                   fluxtot(l)=(einitc(jj)/binsize(jj))*(wvnmhi(l)-
     &             wvnmlo(l))
                 else
                 fluxtot(l)=(einitc(jj)/binsize(jj))*(wnsolc(jj+1)-
     &           wvnmlo(l))
                 if (wvnmhi(l).le.wnsolc(jj+2)) then
                   fluxtot(l)=fluxtot(l)+(einitc(jj+1)/binsize(jj+1))*
     &             (wvnmhi(l)-wnsolc(jj+1))
                 else
                   fluxtot(l)=fluxtot(l)+(einitc(jj+1)/binsize(jj+1))*
     &	           (wnsolc(jj+2)-wnsolc(jj+1))
                 if (wvnmhi(l).le.wnsolc(jj+3)) then
                   fluxtot(l)=fluxtot(l)+(einitc(jj+2)/binsize(jj+2))*
     &             (wvnmhi(l)-wnsolc(jj+2))
                 else
                   fluxtot(l)=fluxtot(l)+(einitc(jj+2)/binsize(jj+2))*
     &             (wnsolc(jj+3)-wnsolc(jj+2))
                 if (wvnmhi(l).le.wnsolc(jj+4)) then
                   fluxtot(l)=fluxtot(l)+(einitc(jj+3)/binsize(jj+3))*
     &             (wvnmhi(l)-wnsolc(jj+3))
                 else
                   fluxtot(l)=fluxtot(l)+(einitc(jj+3)/binsize(jj+3))*
     &             (wnsolc(jj+4)-wnsolc(jj+3)) 
                 if (wvnmhi(l).le.wnsolc(jj+5)) then
	           fluxtot(l)=fluxtot(l)+(einitc(jj+4)/binsize(jj+4))*
     &             (wvnmhi(l)-wnsolc(jj+4))
                 else
                   fluxtot(l)=fluxtot(l)+(einitc(jj+4)/binsize(jj+4))*
     &             (wnsolc(jj+5)-wnsolc(jj+4))
                 if (wvnmhi(l).le.wnsolc(jj+6)) then
                   fluxtot(l)=fluxtot(l)+(einitc(jj+5)/binsize(jj+5))*
     &             (wvnmhi(l)-wnsolc(jj+5))
                 else
                   fluxtot(l)=fluxtot(l)+(einitc(jj+5)/binsize(jj+5))*
     &             (wnsolc(jj+6)-wnsolc(jj+5))
                 if (wvnmhi(l).le.wnsolc(jj+7)) then
                   fluxtot(l)=fluxtot(l)+(einitc(jj+6)/binsize(jj+6))*
     &             (wvnmhi(l)-wnsolc(jj+6))
                 else
                 print*,'Inform Tommy and he can add to this routine'
                 end if
                 end if
                 end if
                 end if
                 end if
                 end if
                 end if
                 goto 210
               end if
             end do
 210         inc=jj
             if (wnf.ge.wn1n) then
             if (wnf.le.wnfn) then
               if (nir_hitgeisa) then
                 fluxtot_nirh(l)=fluxtot(l)
                 sumflux=fluxtot_nirh(l)+sumflux
               elseif (nir_hitirwin) then
               if (wnf.lt.wn2n) then
                 fluxtot_nirhi(l)=fluxtot(l)
                 sumflux=fluxtot_nirhi(l)+sumflux
               else
                 fluxtot_niri(l)=fluxtot(l)
                 sumflux=fluxtot_niri(l)+sumflux
               endif
               elseif (nir_irwinonly) then
                 fluxtot_niri(l)=fluxtot(l)
                 sumflux=fluxtot_niri(l)+sumflux
               endif
             elseif (wnf.le.wnfvlm) then
               fluxtot_vis1(l)=fluxtot(l)
               sumflux=fluxtot_vis1(l)+sumflux
             elseif (wnf.le.wnfvhm) then
               fluxtot_vis2(l)=fluxtot(l)
               sumflux=fluxtot_vis2(l)+sumflux
             elseif (wnf.le.wnfvlmx) then
               fluxtot_vis3(l)=fluxtot(l)
               sumflux=fluxtot_vis3(l)+sumflux
             else
               fluxtot_uv(l)=fluxtot(l)
               sumflux=fluxtot_uv(l)+sumflux
             endif 
             endif
           enddo
           print*,'sumflux(nunwn)',sumflux,'sumflux at TOA',sumflux*
     &        (1./sdis(satday))**2
           close(12)
         endif	      !flagsolar


*--------------------------------------------------------------------------
	  angle=1.        !stays this way if callme eq 1
          qtr_use(1)=0.01 !stays this way if callme eq 1
          if (callme.eq.2.and.firsttime) then
            angler=((-sin(glatr)*sin(decr))/(cos(glatr)*cos(decr)))
            if (abs(angler).ge.1.) then
              if (dec*lat.gt.0.) then  !you're in perpetual daytime
                angle=180.
              else                     !you're in perpetual night
                angle=-1.
              end if
            else
              angle=acos(angler)/pi180
            end if
            if (angle.ge.0..and.angle.le.180.) then      
             !found angle that corresponds to dawn. 
             !Otherwise perpetual darkness
             ! airmass is infinite at sunrise, cos(angle at dawn)=0.            
              qtr_use(1)=angle-0.1      !dawn varies
              qtr_use(qtr)=0.01         !noon always at 0 degrees
              do qtrtime=2,qtr-1
                qtr_use(qtrtime)=qtr_use(qtrtime-1)-
     &                   (qtr_use(1)-qtr_use(qtr))/(qtr-1)
              enddo
              do qtrtime=1,qtr
                  cosz(qtrtime)=(sin(glatr)*sin(decr)+cos(glatr)*
     &             cos(decr)*cos(qtr_use(qtrtime)*pi180))
              enddo
	    end if
          endif   !callme.eq.2
*          print*,'angle',angle,qtr

          if (angle.ge.0..and.angle.le.180.) then      
            if (bdiurn) then
              qtr=1
            endif
            if (obltrig) then               
              lsub=0.   
              ringfrac=1.
            else
              call shadow(lat,lsub,ringfrac)
            endif
            ringfrac=ringfrac/sdis(satday)**2
            do qtrtime=1,qtr
              if (callme.eq.1) then   !callme1  BIN MIDir
                do l=1,maxcly
                  sumdf(l)=0.
                enddo
                if (bin.eq.1) then                   !bin!
                  starttau=dtauc(51,1)
                  startpix=1
                  value= 0.001
                  do l=1,nunwn-1
                    if (abs(starttau-dtauc(51,l+1)).gt.value) then 
                      do j=1,nlay
                        dtauc_tom(j)=dtauc(j,(l+startpix)/2)
                      enddo
*                      call tomsfastdfdtau(nlay,pi,dtauc_tom,
*     &wnout((l+startpix)/2),tnews,p,wvnmlo(startpix), 
*     &wvnmhi(l),tfluxdiff,fluxbase,fluxtottom,fflux)  
                      call tomsfastdfdtau(nlay,pi2,dtauc_tom,
     &extab3,extab4,wnout((l+startpix)/2),tnews,p,wvnmlo(startpix), 
     &wvnmhi(l),tfluxdiff,fluxbase,fluxtottom,fflux)  
                      do j=1,maxcly
                        sumdf(j)=tfluxdiff(j)+sumdf(j)
*                        sumdf2(j,l)=fluxtottom(j)
                      enddo
                     fluxbasetot=fluxbase+fluxbasetot
*                      do j=1,nlay
*                        fluxtott(j,l)=fluxtottom(j)
*                        sfflux(j,l)=fflux(j)
*                      enddo
                      starttau=dtauc(51,l+1)
                      startpix=l+1
                    endif 
                  enddo
                  if (startpix.lt.nunwn) then
                    do j=1,nlay
                      dtauc_tom(j)=dtauc(j,((nunwn-1)+startpix)/2) 
                    enddo
*                    call tomsfastdfdtau(nlay,pi,dtauc_tom,
*     &wnout(((nunwn-1)+startpix)/2),tnews,p,
*     &wvnmlo(startpix),wvnmhi(nunwn-1),tfluxdiff,fluxbase,
*     &fluxtottom,fflux)
                    call tomsfastdfdtau(nlay,pi2,dtauc_tom,
     &extab3,extab4,wnout(((nunwn-1)+startpix)/2),tnews,p,
     &wvnmlo(startpix),wvnmhi(nunwn-1),tfluxdiff,fluxbase,
     &fluxtottom,fflux)
                    do j=1,maxcly
                      sumdf(j)=tfluxdiff(j)+sumdf(j)
*                      sumdf2(j,l)=fluxtottom(j)
                    enddo
                    fluxbasetot=fluxbase+fluxbasetot
*                    do j=1,nlay
*                      fluxtott(j,((nunwn-1)+startpix)/2)=fluxtottom(j)
*                      sfflux(j,((nunwn-1)+startpix)/2)=fflux(j)  
*                    enddo
                  endif
!                 do j=1,maxcly
!                 do l=2,nunwn-1
!                   fluxtottcool(j)=fluxtott(j,l)
!                   fluxtottcoolwn(j,l)=sumdf2(j,l)
!                 enddo
!                 enddo
                else                     !no binning(bin=2) MIDir
                  do l=1,nunwn
                    do j=1,nlay
                      dtauc_tom(j)=dtauc(j,l)
*	            print*,l,'l',dtauc_tom(j),dtauc(j,l),nunwn
                    enddo
*	print*,dtauc_tom(65),wnout(l)
*        print*,tnews(65),p(65),wvnmlo(l),wvnmhi(l)
*        print*,tfluxdiff(65),fluxbase,fluxtottom(65)
*        print*,fluxtottom(65),fflux(65)
*                    call tomsfastdfdtau(nlay,pi,dtauc_tom,
*     &wnout(l),tnews,p,wvnmlo(l),wvnmhi(l)
*     &,tfluxdiff,fluxbase,fluxtottom,fflux)
                    call tomsfastdfdtau(nlay,pi2,dtauc_tom,
     &extab3,extab4,wnout(l),tnews,p,wvnmlo(l),wvnmhi(l)
     &,tfluxdiff,fluxbase,fluxtottom,fflux)
                    do j=1,maxcly
                      sumdf(j)=tfluxdiff(j)+sumdf(j)
*                      sumdf2(j,l)=fluxtottom(j)
                    enddo
                    fluxbasetot=fluxbase+fluxbasetot
*                    do j=1,nlay
*                      fluxtott(j,l)=fluxtottom(j)
*                      sfflux(j,l)=fflux(j)
*                    enddo
                  enddo
*                  do j=1,maxcly
*                    do l=1,nunwn
*                      fluxtottcool(j)=fluxtott(j,l)
*                      fluxtottcoolwn(j,l)=sumdf2(j,l)
*                    enddo
*                  enddo
                endif   !end binning
              elseif (callme.eq.2) then      !callme=2 (heating) 
                do l=1,maxcly
                  sumdf(l)=0.
                enddo
                do l=1,nlay
                  dfb(l)=0.
                enddo
                do i=1,nunwn
                  fluxtot(i)=0.
                enddo
                if (wnf.ge.wn1n) then
                if (wnf.le.wnfn) then
                    if (nir_hitgeisa) then
                      do i=1,nunwn
                        fluxtot(i)=fluxtot_nirh(i)
                      enddo
                    elseif (nir_hitirwin) then
                      if (wnf.lt.wn2n) then
	                do i=1,nunwn
                          fluxtot(i)=fluxtot_nirhi(i)
                        enddo
                      elseif (wnf.ge.wn2n) then
	                do i=1,nunwn
                          fluxtot(i)=fluxtot_niri(i)
                        enddo
                      endif
                    elseif (nir_irwinonly) then
                      do i=1,nunwn
                      fluxtot(i)=fluxtot_niri(i)
                      end do
                    endif
	        elseif (wnf.le.wnfvlm) then
                  do i=1,nunwn
	              fluxtot(i)=fluxtot_vis1(i)
                  end do
	        elseif (wnf.le.wnfvhm) then
                  do i=1,nunwn
	            fluxtot(i)=fluxtot_vis2(i)
	          end do
	        elseif (wnf.le.wnfvlmx) then
                  do i=1,nunwn
	            fluxtot(i)=fluxtot_vis3(i)
                  end do
	        else
                  do i=1,nunwn
	            fluxtot(i)=fluxtot_uv(i)
                  end do
                endif
                endif
*                fluxtotsum=0.
*                do i=1,nunwn
*                  fluxtotsum=fluxtot(i)+fluxtotsum 
*                enddo 
*                print*,fluxtotsum,'fluxtotsum'

*****************

                if (qtrtime.le.9)then
                  write(sqt1,'(I1)'),qtrtime
                  sqt='0'//sqt1
                else
                  write(sqt,'(I2)'),qtrtime
                endif

*****************

**Initialize Bruno Bdiurn:

                if (bdiurn) then
                  call brunodiurn(lat,req,rpole,dec,sdis(satday),
     &               cos_psi,lsub,avecos_day,avecos)
                endif

!!! AEROSOL tau's########################################################
!call tominterp(p,t_o,nlay,play_o,tlay_o,nlay) 

                if (aerosol) then
                  do ll=1,nunwn
                    do l=1,nlay
                      call locate(wno_aerolr,191,wnout(ll),j)
                      slope=(tau_aeroabslr(l,j+1)-tau_aeroabslr(l,j))/
     &                   (wno_aerolr(j+1)-wno_aerolr(j))
                      bb=tau_aeroabslr(l,j)-slope*wno_aerolr(j)
                      tau_aeronew(l,ll)=slope*wnout(ll)+bb
                    enddo
                  enddo
                endif

!------------------------------------------------------------------

!To test Irwin bands:
!i=1,275      3um
!i=300,500    2.3um
!i=575,900    1.7um
!i=900,1175   1.3um
!i=1180,1500  1.1um

!2,nunwn-1
!1,342 (1.7um)
!342,642 (1.3 um)
!642,943 (1.1 um)      
!2,nunwn-1      !###DISORT LOOP###

                do i=1,nunwn

                  fbeam=real(ringfrac*fluxtot(i))   !ringfrac is now ringfrac/sdis(satday)**2
                      do l=1,maxcly
                        if (s.gt.1) then
                        do j=1,s
!!wrong question the l+1 below
                          dtauc_dis(j,l)=real(dtaucc5(j,l+1,i))
!     &+real(tau_aeronew(l+1,i))
                        end do
	                else
                          dtauc_dis(1,l)=real(dtauc(l+1,i))
!     &+real(tau_aeronew(l+1,i))
	                endif
                      enddo
                  if (bdiurn) then
                    if (kcheck) then
                      if (s.eq.10) then
                        do j=1,s
                          delg(j)=delg(j)
                        enddo
                      elseif (s.eq.1) then
                        do j=1,s
                          delg(j)=1.
                        enddo
                      else
                        print*,'error in delg(j)'
                      endif
                    else
                      do j=1,s
                        delg(j)=1.
                      enddo
                    endif    !kcheck
                    call callbruno(avecos_day,fbeam,maxcly,dtauc_dis,
     &               avecos,nlay,dfb,Fbinc,etauFb,fluxbase,
     &               fluxbaseday,qtrtime,s,delg)
                    do l=2,nlay
                      sumdf(l-1)=dfb(l)+sumdf(l-1)   !sum wno to wno
                    enddo
                    do l=2,nlay
                      dfwn(l-1,i)=dfb(l)/delp(l-1)
                      dfwnc(l-1)=dfb(l)/delp(l-1)
                    enddo
                  else  !bdiurn
                    umu0=real(cosz(qtrtime))
*                    if (scattering) then
*                      dislim=5.e-6
*                    else
                      dislim=5.e-12
*                    endif
                    do l=1,maxcly
                      do j=1,s
                        if (dtauc_dis(j,l)/umu0.lt.dislim) then
                          dtauc_dis(j,l)=0.
                        endif
*!!wrong should I change the max value?
                        if (dtauc_dis(j,l).gt.1500.) then
                          dtauc_dis(j,l)=1500.
                        endif
                      enddo
                    enddo
!                  do l=1,nlay
!                    fluxtottheat(l)=Fbinc(l)+fluxtottheat(l)
!                   fluxtottheatwn(l,i)=Fbinc(l)+fluxtottheatwn(l,i)
!                  enddo
                  endif !bdiurn
  
!                  if (scattering) then
!                    do k=0,maxcly
!                      temper(k)=real(tnews(k+1))
!                    enddo
!                    btemp = temper(maxcly)
!                    ttemp = temper(0)
!                    wvnmh=wvnmhi(i)
!                    wvnml=wvnmlo(i)
!                    if (nstr.eq.4) then  !this stuff is to make sure 
!                      cmu1=0.211325 !DIS beam angle .ne. comp. angle
!                      cmu2=0.788675
!                      if (abs(umu0-cmu1).lt.0.01) UMU0=cmu1+0.01
!                      if (abs(umu0-cmu2).lt.0.01) UMU0=cmu2+0.01
!                    elseif (nstr.eq.6) then
!                      cmu1=0.112702
!                      cmu2=0.500000
!                      cmu3=0.887298
!                      if (abs(umu0-cmu1).lt.0.01) UMU0=cmu1+0.01
!                      if (abs(umu0-cmu2).lt.0.01) UMU0=cmu2+0.01
!                      if (abs(umu0-cmu3).lt.0.01) UMU0=cmu3+0.01
!                    else 
!                      print*,'nstr not account. for:potential problem'
!                    endif
*                    call DISORT(NLYR,dtauc_dis,ssalb,nmom,pmom,temper,
*     &                  wvnml,wvnmh,usrtau,ntau,utau,nstr,usrang,numu,
*     &                  umu,nphi,phi,ibcnd,fbeam,umu0,PHI0,
*     &                  FISOT,lamber,alb,btemp,ttemp,temis,
*     &                  plank,onlyfl,accur,prnt,header,maxcly,
*     &                  maxulv,maxumu,MAXPHI,maxmom, RFLDIR, RFLDN,
*     &                  FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )
*                    do l=2,nlay
*                      sumdf(l-1)=DFDT(l)*dtauc_dis(l-1)+sumdf(l-1) 
*                      !sum wno to wno
*                    enddo
!                         fluxbaseday(qtrtime)=flup(nlay)+
!     &                                        fluxbaseday(qtrtime)
!                  elseif (bdiurn) then
!                    do l=1,maxcly
!                      sumdf(l)=sumdf(l)
!                    enddo
!                  else !scattering .false., no bruno diurn ave.
                    call tomsvisabsorp(nlay,dtauc_dis,umu0,
     &                      fbeam,dfdt,fluxbase,fluxtottv,s,delg)
                    fluxbaseday(qtrtime)=fluxbase+
     &                                      fluxbaseday(qtrtime)
!wrong Are these layers correct?
                    do l=2,nlay
                      sumdf(l-1)=dfdt(l)+sumdf(l-1)   !sum wno to wno
                    enddo
                    do l=2,nlay
                      dfwn(l-1,i)=dfdt(l)/delp(l-1)
                      dfwnc(l-1)=dfdt(l)/delp(l-1)
                    enddo
!                  endif !scattering loop
                enddo             !i (nunwn) ###DISORT LOOP###
!             ENDIF !binning heat
              else 
                print*,'Problem: did not make it into *callme* loop'
              endif        !callme

              do l=1,maxcly
                df_big(l,qtrtime)=sumdf(l)
              enddo
            enddo            !qtrtime
                  wnf=wnout(nunwn)

            if (callme.eq.1) then
              do l=1,maxcly
                daily_fm(l)=sumdf(l)*satsec+daily_fm(l)
                ! Instantaneous flux? df/dfp * sec/day *1/sec?
                !summing DF over all wnos W/m^2
              enddo
            else 
              if (bdiurn) then
                do l=1,maxcly
                  daily_fnvutot(l)=sumdf(l)*satsec+daily_fnvutot(l)
                enddo
                fluxbasetot=-1.*fluxbase+fluxbasetot
              else 
                call diurnaltest(qtr,stgkj,stgday,
     &               daily_fnvu,qtr_use,maxcly,nlay,latstr,workdir,
     &               DF_big,fluxbaseday,fluxbase)
                do l=1,maxcly
                  daily_fnvutot(l)=daily_fnvu(l)+daily_fnvutot(l)
                enddo      !come out of qtr and accumulate flux
              ! This is daily DF flux accumulated.
                fluxbasetot=-1.*fluxbase/satsec+fluxbasetot
              endif
!*****Both daily_fnvutot & daily_fm should be in units of
!************  df/dp/day
            endif  !callme
          else   !angle else
            wnf=wnout(nunwn)  
            do l=1,maxcly
              daily_fnvutot(l)=0.            !check all sums and make sure
            end do
            print*,'Night, saturnday=',saturnday
          endif  !end angle if 

!-----------------------------------------------------------------------

          if (callme.eq.2) then
            do l=1,maxcly   !level=2,nlay NOT layer
              sdfdp(l)=(daily_fnvutot(l))/(delp(l))
*              print*,daily_fnvutot(l),delp(l),sdfdp(l)
            end do
          else
            do l=1,maxcly   !level=2,nlay NOT layer
              sdfdp(l)=daily_fm(l) 
            end do
          endif
!         sdfdp in Units of df/dp/day 


          if (rates_write) then
            open(81,file=workdir//'81_'//stgday//'_'//stgkj//'.out',
     &             form='unformatted',access='direct',recl=4*maxcly)
            write(81,rec=1) (real(sdfdp(ll)),ll=1,maxcly)
            write(81,rec=2) (real(play(ll)),ll=1,maxcly)
            close(81)
          endif

*          if (ftot_wn) then
*            open(85,file=workdir//'ftwnh_'//stgday//'_'//stgkj//'.out',
*     &             form='unformatted',access='direct',recl=4*nunwn)
*            open(86,file=workdir//'ftwnc_'//stgday//'_'//stgkj//'.out',
*     &             form='unformatted',access='direct',recl=4*nunwn)
*            do l=1,maxcly
*              write(85,rec=l) (real(fluxtottheatwn(l,ll)),ll=1,nunwn)
*              write(86,rec=l) (real(fluxtottcoolwn(l,ll)),ll=1,nunwn)
*            enddo
*            open(87,file=workdir//'res'//stgkj//'.out',
*     &             form='unformatted',access='direct',recl=4*nunwn)
*            write(87,rec=1) (real(wnout(ll)),ll=1,nunwn)
*            close(87)
*            close(85)
*            close(86)
*          endif 
          do l=1,maxcly   !now on same grid (layer grid)
            sumdfdp(l)=sdfdp(l)+sumdfdp(l)  !flux/satday/bar
          enddo

!----------------------[Where to go next]------------------------------
****This seems way to complicated

*          print*,'wnf=',wnf ,'this is last wno of this step'
          if (wnf.eq.wfirmax) then
            wnf=wmidmin
*            print*,'going to...',wmidmin,'(midIR)'
          elseif (wnf.ge.wmidmax.and.wnf.lt.wn1n) then
            wnf=wn1n
*            print*,'going to...',wn1n,'(NIR min)'
          elseif (nir_hitirwin) then
            if (wnf.gt.4791.179.and.wnf.lt.wn2n) then
              wnf=wn2n
            endif
          elseif (wnf.ge.wnfn.and.wnf.lt.wn1vlm)then
            wnf=wn1vlm        
*            print*,'going to...',wn1vlm,'(wn1vlm)'
          elseif (wnf.ge.wnfvlm.and.wnf.lt.wn1vhm) then
            wnf=wn1vhm
*            print*,'going to...',wn1vhm,'(wn1vhm)'
          elseif (wnf.ge.wnfvhm.and.wnf.lt.wn1vlmx) then
            wnf=wn1vlmx
*            print*,'going to...',wn1vlmx,'(wn1vlmx)'
          elseif (wnf.ge.33288.0.and.
     &        wnf.lt.39999.0)then
            wnf=wn1uv
*            print*,'going to...',wn1uv,'(wn1uv)'
          elseif (wnf.ge.98039.) then
*            print*,'End UV; end of spectrum'
          endif



        enddo                   !kj file loop
           flagsolar=.false.
	   firsttime=.false.

          do l=1,maxcly
            dtdt(l)=(sumdfdp(l)*(xm*g))*(stepsize)/((Cp*1.e5))
              !changed to divide by satsec to get Flux/sec     
              !convert bars --> Pa for units to work?
              !unit of time is already converted to per day
              !rather then per second
          enddo

        tnews(nlay)=sqrt(sqrt((4.7+(4.7-fluxbasetot))/5.6704e-8))
*        print*,'tnews flux',tnews(nlay)
        do j=1,maxcly
          pm1(j)=p(j+1)
          playm1(j)=play(j+1)
        enddo
        call tominterplin(playm1,dtdt,maxcly,p,dtdt_lay,nlay)
        stepsizenew=stepsize           
*        print*,'stepsizenew before loop',stepsizenew
cc need to remember what initial 'stepsize' is for calcs. later on. 
cc Don't change 'stepsize'
 9998   continue
        do l=1,maxcly
          if (abs(dtdt_lay(l)).ge.2.) then
            if (stepsize.ge.2) then
	      do kl=1,maxcly
                dtdt_lay(kl) = dtdt_lay(kl)/2.
	      end do
              stepsizenew=stepsizenew/2
	      goto 9998
            else
              print*,'problem with starting step, step.lt.1; critical'
            endif
	  end if
	end do
        if (stepsizenew.eq.stepsize) then
 	  do l=1,maxcly
            if (abs(dtdt_lay(l)).ge.0.8) then
   	      goto 9995
	    endif
	  end do
	  do kl=1,maxcly
            dtdt_lay(kl) = dtdt_lay(kl)*2 
	  end do
          stepsizenew=(stepsizenew*2)
	end if
*        print*,'stepsizenew',stepsizenew
! subtract off very 1st stepsize and add NEW stepsize:
        today=(saturnday-stepsize)+stepsizenew     !current day 
        tomorrow=today+stepsizenew              !day to start next
*        print*,'saturnday at start',saturnday
*        print*,'today adjusted (file write)',today
*        print*,'tomorrow adjust',tomorrow
        do l=1,maxcly
          tnews(l)=tnews(l)+dtdt_lay(l)
*          print*,'tnew',tnews(l),'change in t',dtdt_lay(l)
        enddo
*        print*,'tnewbase',tnews(nlay)
        do l=1,nlay
          if (tnews(l).gt.300.) then
            tnews(l)=299.
          elseif (tnews(l).lt.75.) then
            tnews(l)=76.
          endif
        enddo

	call numtostring5(today,stgday)

        open(82,file=workdir//'temp.out',
     &         form='unformatted',access='direct',recl=4*(nlay+1))
        open(83,file=workdir//'tempinfo.out',
     &         form='unformatted',access='direct',recl=4)
        count=count+1
        write(82,rec=1+count) (real(today),
     &(real(tnews(ll)),ll=1,nlay))
        write(83,rec=1) real(count)
        close(82)
        close(83)

*        if (fluxtot_write) then
*          open(85,file=workdir//'ftot_'//stgday//'.out',
*     &             form='unformatted',access='direct',recl=4*nlay)
*          write(85,rec=1) (real(fluxtottheat(ll)),ll=1,nlay)
*          write(85,rec=2) (real(fluxtottcool(ll)),ll=1,nlay)
*          close(85)
*        endif 

        goto 989                  
ccc either leave main loop above, or if middle ELSEIF condition met,
ccc continue from here and don't leave main saturnday loop.
 9995   CONTINUE
*        print*,'satday',saturnday
        do l=1,maxcly
          tnews(l)=tnews(l)+dtdt_lay(l)
*          print*,'tnew',tnews(l),'change in t',dtdt_lay(l)
        enddo
*        print*,'tnewbase',tnews(nlay)
        do l=1,nlay
          if (tnews(l).gt.300.) then
            tnews(l)=299.
          elseif (tnews(l).lt.75.) then
            tnews(l)=76.
          endif
        enddo

	call numtostring5(saturnday,stgday)
       
        open(82,file=workdir//'temp.out',
     &         form='unformatted',access='direct',recl=4*(nlay+1))
        open(83,file=workdir//'tempinfo.out',
     &         form='unformatted',access='direct',recl=4)
        count=count+1
        write(82,rec=1+count) (real(saturnday),
     &(real(tnews(ll)),ll=1,nlay))
        write(83,rec=1) real(count)
        close(82)
        close(83)
 133    continue
      enddo                   !saturnday loop
      enddo 		      !nnr
      close(27)

      end 
