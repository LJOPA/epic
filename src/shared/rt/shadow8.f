*******Subject: ring shadowing with no sind's and no NaN's allowed

	subroutine shadow(latin,lsubss,ringfrac)

      Implicit None

      REAL*8 aintc, aextc, dec, req, rpole, rcint, rcext,
     $   rb1int, rb1ext, rb2int, rb2ext, rcdint, rcdext, raint,
     $   raext, tauc, taub1, taub2, taucd, taua, pi, lat,
     $   pi180, obl, aintb1, aintb2, aextb1, aextb2, aintcd,
     $   aextcd, ainta, aexta, phiinta, phiintb1, phiintb2,
     $   phiintcd, phiintc, phiexta, phiextb1, phiextb2,
     $   phiextcd, phiextc, alpha, phie, denom, ringfrac,
     $   lsubss, latin
c
c	Saturn modification: ring shadowing
c
      pi=acos(-1.)
      pi180=pi/180.
      obl = 26.73*pi180   ! pi180 is pi/180 (converts deg to radians)
      lat=latin*pi180    !converts deg to radians
      req = 60268.
      rpole = 54364.
c
      rcint = 74655.
      rcext = 91975.
      rb1int = 91975.
      rb1ext = 98500.
      rb2int = 98500.
      rb2ext = 117510.
      rcdint = 117510.
      rcdext = 122340.
      raint = 122340.
      raext = 136780.
      tauc = 0.099
      taub1 = 0.92
      taub2 = 2.07
      taucd = 0.14
      taua = 0.51
c
c     Note: I had all angles in degrees, but my f90 does not support 
c     sind, cosd, tand, so I switched to radians.  This code below
c     has not been fully checked or tested now that I had to switch 
c     to radians.  I changed all cases of sind, etc.  Also, f90 does
c     not allow NaNs or Infty, so I needed to add a bunch of if
c     statements to avoid divide by zeros.  This portion of the 
c     code has survived limited testing/checking.  Looks okay, but
c     may still be some divide by zeros for some of the parameters 
c     below under untested conditions.
c
      dec = obl*sin(lsubss*pi180)    ! in radians

         alpha = (-1.*(req/rpole)**2.)*tan(lat)*tan(dec)
         phie = acos(alpha)
c
        if (dec .lt. 1.e-10 .and. dec .gt. -1.e-10) then
	  aintc = 1.0e+38
        else if (lat .lt. 1.e-10 .and. lat .gt. -1.0e-10) then
	  aintc = 1.0e+38
        else
          aintc = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rcint/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
        endif
        if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	  aextc = -1.0e+38
        else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	  aextc = -1.0e+38
        else
          aextc = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rcext/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
        endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aintb1 = 1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aintb1 = 1.0e+38
      else
        aintb1 = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rb1int/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aextb1 = -1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aextb1 = -1.0e+38
      else
        aextb1 = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rb1ext/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aintb2 = 1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aintb2 = 1.0e+38
      else
        aintb2 = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rb2int/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aextb2 = -1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aextb2 = -1.0e+38
      else
        aextb2 = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rb2ext/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aintcd = 1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aintcd = 1.0e+38
      else
        aintcd = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rcdint/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aextcd = -1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aextcd = -1.0e+38
      else
        aextcd = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(rcdext/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	ainta = 1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	ainta = 1.0e+38
      else
        ainta = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(raint/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (dec .lt. 1.0e-10 .and. dec .gt. -1.0e-10) then
	aexta = -1.0e+38
      else if (lat .lt. 1.0e-10 .and. lat .gt. -1.0e-10) then
	aexta = -1.0e+38
      else
        aexta = ( cos(lat)*cos(lat) + 
     $       sin(lat)*sin(lat)/(tan(dec)*tan(dec)) - 
     $       ( cos(lat)*cos(lat) + sin(lat)*sin(lat)*
     $       (req/rpole)**2 )*(raext/req)**2 )/( 2.*sin(lat)*
     $       cos(lat)/tan(dec) )
      endif
      if (aintc .gt. 1.) aintc = 1.
      if (aintc .gt. -1.) phiintc = acos(aintc)
      if (aintc .le. -1.) phiintc = pi
      if (aextc .ge. 1.) aextc = 1.
      if (aextc .ge. -1.) phiextc = acos(aextc)
      if (aextc .lt. -1.) phiextc = 1.0e+38
      if (aintb1 .gt. 1.) aintb1 = 1.
      if (aintb1 .gt. -1.) phiintb1 = acos(aintb1)
      if (aintb1 .le. -1.) phiintb1 = pi
      if (aextb1 .ge. 1.) aextb1 = 1.
      if (aextb1 .ge. -1.) phiextb1 = acos(aextb1)
      if (aextb1 .lt. -1.) phiextb1 = 1.0e+38
      if (aintb2 .gt. 1.) aintb2 = 1.
      if (aintb2 .gt. -1.) phiintb2 = acos(aintb2)
      if (aintb2 .le. -1.) phiintb2 = pi
      if (aextb2 .ge. 1.) aextb2 = 1.
      if (aextb2 .ge. -1.) phiextb2 = acos(aextb2)
      if (aextb2 .lt. -1.) phiextb2 = 1.0e+38
      if (aintcd .gt. 1.) aintcd = 1.
      if (aintcd .gt. -1.) phiintcd = acos(aintcd)
      if (aintcd .le. -1.) phiintcd = pi
      if (aextcd .ge. 1.) aextcd = 1.
      if (aextcd .ge. -1.) phiextcd = acos(aextcd)
      if (aextcd .lt. -1.) phiextcd = 1.0e+38
      if (ainta .gt. 1.) ainta = 1.
      if (ainta .gt. -1.) phiinta = acos(ainta)
      if (ainta .le. -1.) phiinta = pi
      if (aexta .ge. 1.) aexta = 1.
      if (aexta .ge. -1.) phiexta = acos(aexta)
      if (aexta .lt. -1.) phiexta = 1.0e+38
      ringfrac = 1.0
       if (phie .ge. phiexta .and. phie .le. phiintc) then
	 denom = cos(lat)*cos(dec)*sin(phie) +
     $     sin(lat)*sin(dec)*phie*(req/rpole)**2.
	 if (phie .ge. phiexta .and. phie .le. phiinta) then  ! A only
            ringfrac = cos(lat)*cos(dec)*( sin(phiexta) +
     $        (exp(-taua/cos(dec)))*(sin(phie) - sin(phiexta)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiexta + 
     $        (exp(-taua/cos(dec)))*(phie - phiexta) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextcd .and. ainta .le. 1.0 .and. phie .le. 
     $    phiintcd) then   ! A + CD only
            ringfrac = cos(lat)*cos(dec)*( sin(phiexta) +
     $        (exp(-taua/cos(dec)))*(sin(phiinta) - sin(phiexta)) + 
     $        (exp(-taucd/cos(dec)))*(sin(phie) - sin(phiextcd)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiexta + 
     $        (exp(-taua/cos(dec)))*(phiinta - phiexta)  + 
     $        (exp(-taucd/cos(dec)))*(phie - phiextcd) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextcd .and. ainta .gt. 1.0 .and. phie .le. 
     $    phiintcd) then   ! CD only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextcd) +
     $        (exp(-taucd/cos(dec)))*(sin(phie) - sin(phiextcd)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextcd + 
     $        (exp(-taucd/cos(dec)))*(phie - phiextcd) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb2 .and. ainta .le. 1.0 .and. aintcd .le. 
     $    1.0 .and. phie .le. phiintb2) then   ! A + CD + B2 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiexta) +
     $        (exp(-taua/cos(dec)))*(sin(phiinta) - sin(phiexta)) + 
     $        (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $       + (exp(-taub2/cos(dec)))*(sin(phie) - sin(phiextb2)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiexta + 
     $        (exp(-taua/cos(dec)))*(phiinta - phiexta)  + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phie - phiextb2) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb2 .and. ainta .gt. 1.0 .and. aintcd .le. 
     $    1.0 .and. phie .le. phiintb2) then   ! CD + B2 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextcd) 
     $       + (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $       + (exp(-taub2/cos(dec)))*(sin(phie) - sin(phiextb2)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextcd + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phie - phiextb2) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb2 .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. phie .le. phiintb2) then   ! B2 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextb2) 
     $       + (exp(-taub2/cos(dec)))*(sin(phie) - sin(phiextb2)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextb2 + 
     $        (exp(-taub2/cos(dec)))*(phie - phiextb2) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb1 .and. ainta .le. 1.0 .and. aintcd .le. 
     $    1.0 .and. aintb2 .le. 1.0 .and. phie .le. phiintb1) 
     $    then   ! A + CD + B2 + B1 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiexta) 
     $      + (exp(-taua/cos(dec)))*(sin(phiinta) - sin(phiexta))  
     $      + (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phie) - sin(phiextb1)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiexta + 
     $        (exp(-taua/cos(dec)))*(phiinta - phiexta)  + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phie - phiextb1) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb1 .and. ainta .gt. 1.0 .and. aintcd .le. 
     $    1.0 .and. aintb2 .le. 1.0 .and. phie .le. phiintb1) 
     $    then   ! CD + B2 + B1 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextcd) 
     $      + (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phie) - sin(phiextb1)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextcd + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phie - phiextb1) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb1 .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. aintb2 .le. 1.0 .and. phie .le. phiintb1) 
     $    then   ! B2 + B1 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextb2) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phie) - sin(phiextb1)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextb2 + 
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phie - phiextb1) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextb1 .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. aintb2 .gt. 1.0 .and. phie .le. phiintb1) 
     $    then   ! B1 only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextb1) 
     $        + (exp(-taub1/cos(dec)))*(sin(phie) - sin(phiextb1)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextb1 + 
     $        (exp(-taub1/cos(dec)))*(phie - phiextb1) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextc .and. ainta .le. 1.0 .and. aintcd .le. 
     $    1.0 .and. aintb2 .le. 1.0 .and. aintb1 .le. 1.0 .and. 
     $    phie .le. phiintc) then    ! A + CD + B2 + B1 + C
            ringfrac = cos(lat)*cos(dec)*( sin(phiexta) +
     $        (exp(-taua/cos(dec)))*(sin(phiinta) - sin(phiexta)) 
     $      + (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phiintb1) - sin(phiextb1)) 
     $      + (exp(-tauc/cos(dec)))*(sin(phie) - sin(phiextc)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiexta + 
     $        (exp(-taua/cos(dec)))*(phiinta - phiexta)  + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phiintb1 - phiextb1) +
     $        (exp(-tauc/cos(dec)))*(phie - phiextc) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextc .and. ainta .gt. 1.0 .and. aintcd .le. 
     $    1.0 .and. aintb2 .le. 1.0 .and. aintb1 .le. 1.0 .and. 
     $    phie .le. phiintc) then    ! CD + B2 + B1 + C only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextcd) 
     $      + (exp(-taucd/cos(dec)))*(sin(phiintcd) - sin(phiextcd)) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phiintb1) - sin(phiextb1)) 
     $      + (exp(-tauc/cos(dec)))*(sin(phie) - sin(phiextc)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextcd + 
     $        (exp(-taucd/cos(dec)))*(phiintcd - phiextcd) +
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phiintb1 - phiextb1) +
     $        (exp(-tauc/cos(dec)))*(phie - phiextc) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextc .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. aintb2 .le. 1.0 .and. aintb1 .le. 1.0 .and. 
     $    phie .le. phiintc) then    ! B2 + B1 + C only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextb2) 
     $      + (exp(-taub2/cos(dec)))*(sin(phiintb2) - sin(phiextb2))  
     $      + (exp(-taub1/cos(dec)))*(sin(phiintb1) - sin(phiextb1)) 
     $      + (exp(-tauc/cos(dec)))*(sin(phie) - sin(phiextc)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextb2 + 
     $        (exp(-taub2/cos(dec)))*(phiintb2 - phiextb2) +
     $        (exp(-taub1/cos(dec)))*(phiintb1 - phiextb1) +
     $        (exp(-tauc/cos(dec)))*(phie - phiextc) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextc .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. aintb2 .gt. 1.0 .and. aintb1 .le. 1.0 .and. 
     $    phie .le. phiintc) then    ! B1 + C only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextb1) 
     $      + (exp(-taub1/cos(dec)))*(sin(phiintb1) - sin(phiextb1)) 
     $      + (exp(-tauc/cos(dec)))*(sin(phie) - sin(phiextc)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextb1 + 
     $        (exp(-taub1/cos(dec)))*(phiintb1 - phiextb1) +
     $        (exp(-tauc/cos(dec)))*(phie - phiextc) )
	    ringfrac = ringfrac/denom
         endif	     
	 if (phie .ge. phiextc .and. ainta .gt. 1.0 .and. aintcd .gt. 
     $    1.0 .and. aintb2 .gt. 1.0 .and. aintb1 .gt. 1.0 .and. 
     $    phie .le. phiintc) then    ! C only
            ringfrac = cos(lat)*cos(dec)*( sin(phiextc) 
     $      + (exp(-tauc/cos(dec)))*(sin(phie) - sin(phiextc)) ) 
     $        + ((req/rpole)**2.)*sin(lat)*sin(dec)*( phiextc + 
     $        (exp(-tauc/cos(dec)))*(phie - phiextc) )
	    ringfrac = ringfrac/denom
         endif	     
       endif

	return
	end


