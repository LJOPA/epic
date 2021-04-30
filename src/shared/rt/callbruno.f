        Subroutine callbruno(avecos_day,FBEAM,maxcly,dtauc_tom,avecos,
     &   nlay,dFb2,Fbinc2,etauFb,fluxbase,fluxbaseday,qtrtime,s,delg)
        Implicit None
        integer s,j,maxcly,l,nlay,qtrtime
        real*8 Fb(nlay),Fbinc(10,nlay),dFb(10,nlay),Fbinc2(nlay),
     &         etautot,wnf(1),delg(10),dFb2(nlay),etauFb2,
     &         avecos_day,avecos,etauFb(10),fluxbase,fluxbaseday(10)
        real   dtauc_tom(10,nlay-1),fbeam
        logical  kcheck
         do j=1,s
        do l=1,maxcly
           dFb2(l+1)=0.
           dFb(j,l+1)=0.
           Fbinc(j,l+1)=0.
        enddo
        enddo
         IF(s.eq.1)then
            Fbinc2(1)=FBEAM*avecos_day
            do l=1,maxcly
             dFb2(l+1)=0.
             etauFb2=exp(-dtauc_tom(1,l)/(avecos))
             dFb2(l+1)=Fbinc2(l)*(1.-etauFb2)
             Fbinc2(l+1)=Fbinc2(l)*etauFb2
            enddo
            fluxbase=Fbinc2(nlay)
            fluxbaseday(qtrtime)=fluxbase+fluxbaseday(qtrtime)
        ELSEif(s.eq.10)then
         do j=1,s
          dFb2(l+1)=0.
          Fbinc2(l+1)=0.
          Fbinc(j,1)=FBEAM*avecos_day*delg(j)
          do l=1,maxcly
           dFb(j,l+1)=0.
           etauFb(j)=exp(-dtauc_tom(j,l)/(avecos))
           dFb(j,l+1)=Fbinc(j,l)*(1.-etauFb(j))  !*delg(j)
!           if(l.eq.1)then
!             Fbinc(j,l+1)=Fbinc(j,l)*delg(j)*etauFb(j)
!           else
             Fbinc(j,l+1)=Fbinc(j,l)*etauFb(j)
         !     print*,Fbinc(j,l+1),j,l+1,'f(j)'
!           endif
           dFb2(l+1)=dFb(j,l+1)+dFb2(l+1)
           Fbinc2(l+1)=Fbinc2(l)+Fbinc(j,l)
          enddo
         enddo
         fluxbase=Fbinc2(nlay)+fluxbase
         fluxbaseday(qtrtime)=fluxbase+fluxbaseday(qtrtime)
        ELSE
         print*,'problem with callbruno: s is undefined'
         stop
        ENDIF
!        do l=1,maxcly
!           dFb2(l+1)=0.
!           Fbinc2(l+1)=0.
!        enddo
!        do j=1,s
!         do l=1,maxcly
!          dFb2(l+1)=dFb(j,l+1)*delg(j)+dFb2(l+1)
!          Fbinc2(l+1)=Fbinc(j,l+1)*delg(j)+Fbinc2(l+1)
!         enddo
!        enddo
        return
        end
