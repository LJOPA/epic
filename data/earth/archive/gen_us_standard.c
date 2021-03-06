/*
 * Generates U.S. Standard Atmosphere T(p) and writes
 * as t_vs_p.earth.us_standard.
 *
 * Formulas from http://scipp.ucsc.edu/outreach/balloon/atmos/1976%20Standard%20Atmosphere.htm
 *
 * To compile: cc -lm -o gen_us_standard gen_us_standard.c
 *
 * T. Dowling, 7/27/07
 */

#include <stdio.h>
#include <math.h>

main() {
  double
    p,t,h,
    p0   = 1013.250,
    rho0 = 1.225,
    t0   = 288.15;
  int
    k,
    nk = 143;
  FILE
    *out;

  out = fopen("t_vs_p.earth.us_standard","w");
  /*
   * Output header.
   */
  fprintf(out,"Temperature vs. pressure for 1976 U.S. Standard Atmosphere\n");
  fprintf(out,"Generated by the program epic/data/earth/archive/gen_us_standard.c\n");
  fprintf(out,"Formulas from http://scipp.ucsc.edu/outreach/balloon/atmos/1976%%20Standard%%20Atmosphere.htm\n");
  fprintf(out,"\n\n");
  fprintf(out,"#  p[mbar]   T[K]\n");
  fprintf(out,"%d\n",nk);

  for (k = 0; k < nk; k++) {
    /*
     * Height field [m]
     */
    h = 71000.*(double)(nk-1-k)/(nk-1);
    if (h <= 11000.) {
      /* Region 1 */
      t = t0*(1.-h/44329.);
      p = p0*pow(1.-h/44329,5.255876);
    }
    else if (h <= 20000.) {
      /* Region 2 */
      t = t0*.751865;
      p = p0*.223361*exp((10999.-h)/6341.4);
    }
    else if (h <= 32000.) {
      /* Region 3 */
      t = t0*(0.682457+h/288136.);
      p = p0*pow(0.988626+h/198903.,-34.16319);
    }
    else if (h <= 47000.) {
      /* Region 4 */
      t = t0*(0.482561+h/102906.);
      p = p0*pow(0.898309+h/55280,-12.20114);
    }
    else if (h <= 51000.) {
      /* Region 5 */
      t = t0*.939268;
      p = p0*0.00109456*exp((46998-h)/7922);
    }
    else {
      /* Region 6 */
      t = t0*(1.434843-h/102906.);
      p = p0*pow(0.838263-h/176142.,12.20114);
    }
    fprintf(out,"%10.4f  %6.2f\n",p,t);
  }

  fclose(out);
}
