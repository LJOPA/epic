           SUBROUTINE uv(nlay,nwn,t,wnf,wvnmhi_uv,wvnmlo_uv,
     &           wnot,dtaucc,workdir,cd,nunwn,ff_c2h4,ff_c2h6,
     &           ff_c2h2,ff_ch4)
           implicit none
           integer nunwn,nlay,nwn
           real*8 wnf,n_stp,wn(nunwn),wvnmhi_uv(nwn),
     &          wvnmlo_uv(nwn),wnot(nwn),tau_c2h2(nlay,nunwn),
     &          t(nlay),tau_ch4(nlay,nunwn),cx_c2h4(36),     
     &          cmm,angm,cd(nlay),dtaucc(nlay,nwn),cx_c2h4c(36),
     &          tau_c2h6(nlay,nunwn),uvwn_c2h4c(36),   
     &          tau_c2h4(nlay,nunwn),
     &          cx_ch4(36),cx_c2h2(36),cx_c2h6(36),uvwn_c2h6(36),
     &          cx_ch4c(36),cx_c2h2c(36),cx_c2h6c(36),uvwn_c2h4(36),
     &          uvwn_ch4(36),uvwn_c2h2(36),uvwn_ch4c(36),uvwn_c2h2c(36),
     &          uvwn_c2h6c(36),ff_ch4(1,nlay),ff_c2h2(1,nlay),
     &          ff_c2h6(1,nlay),ff_c2h4(1,nlay)
           character*34 workdir
           integer ll,k,l,i,ij

!---------------------------------------------------------------------
        print*,'UV'
        wn(1)=wnf
        n_stp=2.68E19  !Loschmidt's number in molecule/cm^3/amagat

        open(60,
     &    file='nsubs/cross_ch4.inp',
     &      status='old')
        open(20,
     &    file='nsubs/cross_c2h2.inp',
     &       status='old')
        open(29,
     &    file='nsubs/cross_c2h4.inp',
     &      status='old')
        open(38,
     &    file='nsubs/cross_c2h6.inp',
     &       status='old') 
 343    format(f6.1,1x,f28.26)
        angm=0.0000000001           !(m/ang)
        cmm=0.01                    !(m/cm)
        do l=1,nunwn
           read(60,343) uvwn_ch4(l),cx_ch4(l)
           read(20,343) uvwn_c2h2(l),cx_c2h2(l)
           read(29,343) uvwn_c2h4(l),cx_c2h4(l)
           read(38,343) uvwn_c2h6(l),cx_c2h6(l)
        enddo
        close(60)
        close(20)
        close(29)
        close(38)
        do l=0,nunwn-1
            uvwn_ch4c(l+1)=1./(uvwn_ch4(36-l)*angm/cmm)
            uvwn_c2h2c(l+1)=1./(uvwn_c2h2(36-l)*angm/cmm)
            uvwn_c2h4c(l+1)=1./(uvwn_c2h4(36-l)*angm/cmm)
            uvwn_c2h6c(l+1)=1./(uvwn_c2h6(36-l)*angm/cmm)
            cx_ch4c(l+1)=cx_ch4(36-l)
            cx_c2h2c(l+1)=cx_c2h2(36-l)
            cx_c2h4c(l+1)=cx_c2h4(36-l)
            cx_c2h6c(l+1)=cx_c2h6(36-l)
        enddo
        do i=1,nunwn
          wn(i)=uvwn_ch4c(i)
          wnot(i)=wn(i)
        enddo
        do i=2,nunwn-1
          wvnmhi_uv(i)=(wn(i)+wn(i+1))/2.
          wvnmlo_uv(i)=(wn(i)+wn(i-1))/2.
        enddo
! tau = sigma(cm^2/molec)*nSTP(molec/cm^3/amagat)*cd(km*amagat)*1E5(cm/km)
! tau = unitless!
*        open(14,file=workdir//'dtauc_uv.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
*        open(61,file=workdir//'tauuv_ch4.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
*        open(62,file=workdir//'tauuv_c2h2.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
*        open(63,file=workdir//'tauuv_c2h6.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)
*        open(64,file=workdir//'tauuv_c2h4.out',form='unformatted',
*     &     access='direct',recl=4*nunwn)

        do k=1,nlay
         do l=1,nunwn
         tau_ch4(k,l) = cd(k)*n_stp*(1.E5)*cx_ch4c(l)*ff_ch4(1,k)
         tau_c2h2(k,l) = cd(k)*n_stp*(1.E5)*cx_c2h2c(l)*ff_c2h2(1,k)
         tau_c2h4(k,l) = cd(k)*n_stp*(1.E5)*cx_c2h4c(l)*ff_c2h4(1,k)
         tau_c2h6(k,l) = cd(k)*n_stp*(1.E5)*cx_c2h6c(l) *ff_c2h6(1,k)
         enddo
*        write(61,rec=k+h) (real(tau_ch4(k,ll)),ll=1,nunwn)
*        write(62,rec=k+h) (real(tau_c2h2(k,ll)),ll=1,nunwn)
*        write(63,rec=k+h) (real(tau_c2h6(k,ll)),ll=1,nunwn)
*        write(64,rec=k+h) (real(tau_c2h4(k,ll)),ll=1,nunwn)
        enddo
      
*        close(61)
*        close(62)
*        close(63)
*        close(64)
        do l=1,nlay
          do ij=1,nunwn
             dtaucc(l,ij)=tau_ch4(l,ij)+tau_c2h2(l,ij)+
     &             tau_c2h6(l,ij)+tau_c2h4(l,ij)
*	print*,dtaucc(l,ij),tau_ch4(l,ij),tau_c2h2(l,ij),
*     &tau_c2h6(l,ij),tau_c2h4(l,ij),cd(l),ff_ch4(1,l)
          enddo
*          write(14,rec=l) (real(dtaucc(l,i)),i=1,nunwn)
        enddo
*        close(14)
        return
        END
 
