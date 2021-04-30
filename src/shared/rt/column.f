	Subroutine column(p,t,xm,g,col,d)
	
	Real*8 p,t,xm,g,col,R,d
	
* ---------------------------------------------------------------------
*** Input: 	p  : 	Pressure (bars), Temperature (K)
*		xm :    Mean molecular weight (kg) 
*		g  :    gravity (m s-2)
*** Output:	col: 	column abundance of gas above this level (km amagat)
*		d  :    density of atmosphere (amagat) 
* ---------------------------------------------------------------------	

	R   = 8.3145			! gas constant (J/(mol K))
	To  = 273.15			! ref temp at which den(am)=P(To/T)

	col = (p/1.01325) * (R*To)/(g*xm)		! column abund in metre amagat
	col = col * 1.E-3 		! in km amagat 
	
	d   = (p/1.01325) * (To/t)  	! amagat  (1atm=1.01325bar) 
*********************************************************
****the density (d) is correct.  Ive checked it and compared the
****result with den=(pres*1e5)/(k*temp*2.68e25)
****   pres in bars
****  1bar=1e5 pascals (pasals are kg m-1 s-2)
****  k= boltzmanns constant= 1.38 e-23 Kg m2 s-2 K-1
****  2.68e25 = loschmidt number in m-3
	
	return 
	end 
	
