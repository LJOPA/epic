***********************************************************************
***********************Gravity Subroutine*******************************
**************copied from Bruno 2003***********************

	subroutine gravsub(pla,lat,g_1bar)

	character planet(4)*3,pla*3
	real*8 a_eq(4),g_eq(4),oblat(4),j2(4),dpi,lat,g_1bar,r_1bar
	real*8 q(4)
	integer i,ipla,nplanet


	data planet/'JUP','SAT','URA','NEP'/,a_eq/71492.,60271.,25559.,
     &24766./,g_eq/2313.,908.,865.,1091./,oblat/0.0649,0.0982,0.0229,
     &0.0171/,q/0.0891,0.1548,0.0295,0.0261/,j2/0.01469,0.01633,
     &3.35d-03,3.41d-03/,dpi/1.745329252d-02/

	nplanet=4
C**		Determination of the planet and gravity at 1 bar
	do 1 i=1,nplanet
	if(pla.ne.planet(i)) goto 1
	ipla=i
	goto 2
  1	continue
	print 200,pla
  200	format(' Unknown planet: ',a3)
	stop 11
  2	sin2=dsin(dpi*lat)**2
	ar=dsqrt(1.+sin2*oblat(ipla)*(2.-oblat(ipla))/
     &(1.-oblat(ipla))**2)
	g_1bar=g_eq(ipla)*ar*ar*(1.-1.5*j2(ipla)*(3.*sin2-1.)*ar*ar-
     &q(ipla)*(1.-sin2)/(ar*ar*ar))/(1.+1.5*j2(ipla)-q(ipla))
	r_1bar=a_eq(ipla)/ar
	write(6,201) pla,lat,g_1bar,r_1bar
  201	format(' Planet: ',a3,'  Latitude:',f6.2,'   G_1bar:',f7.1,
     &'  cm sec-2','  R_1bar:',f8.1,' km')

	return
	end
	
