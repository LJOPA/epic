*******************************************************************************
***********************	Blackbody Function ************************************
******************************************************************************	



	function B(wavnum,t)		! units are erg/s/cm2/sr/(cm-1)
					! wavnum=wavenumber(cm-1)
	implicit none			! c*B(v) p.29 Goody

	Real*8 nom,plancknu,t,wavnum,B,hc2,hcok

*	h = 6.6260755e-34		! J.s
*	c = 299792458.0e0		! m
*	k = 1.380658e-23		! J/K
	nom = wavnum * 1.0e2		! m^-1
	hcok=1.4388e-2  !h*c/k
	hc2=1.191e-16	!2*h*c^2
*old	plancknu = 2.0e0 * h * c * c * nom ** 3 /
*old     &		 ( exp( h * c * nom / ( k * t ) ) - 1.0e0 )
*old	B = plancknu * 1.0e-4 * 1.0e2 * 1.0e+7
	plancknu =hc2*nom**3/(exp( hcok*nom/ t)-1.0e0)
	B = plancknu * 1.0e5   !1.0e5=1e-4*1e2*1e7
	return
	end
