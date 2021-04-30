/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC4_PATH/License.txt                                       *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
 * Boston, MA 02110-1301, USA.                                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * epic_funcs_astron.c * * * * * * * * * * * * * * * * 
 *                                                                 *
 *  This file includes functions used to calculate astronomical    *
 *  quantities such as orbital positions.                          * 
 *                                                                 *
 *       solar_longitude()                                         *
 *       eccentric_anomaly()                                       *
 *       true_anomaly()                                            *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Variables global to this file.
 */
double
  Kepler_M,
  Kepler_e;

/*
 * Function prototypes.
 */
double kepler_zero(double E);

/*====================== solar_longitude() =========================*/

/*
 * Converts calendar date/time [C data type time_t, sec] into L_s [deg], 
 * the planetocentric longitude of the Sun.
 */

double solar_longitude(planetspec *planet,
                           time_t  date)
{
  double
    tJ2000,
    M,Mrad,E,nu,
    alpha_fms,
    l_s,
    tmp;
  static time_t
    J2000;
  struct tm
    epoch2000;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="solar_longitude";

  if (!initialized) {
    /* Trigger timezone evaluation. */
    epoch2000.tm_year  =  0;
    epoch2000.tm_mon   =  0;
    epoch2000.tm_mday  =  0;
    epoch2000.tm_hour  =  0;
    epoch2000.tm_min   =  0;
    epoch2000.tm_sec   =  0;
    epoch2000.tm_isdst = -1;
    mktime(&epoch2000);

    /* Set the standard epoch J2000.0 in terms of UTC. */
    epoch2000.tm_year = 2000-1900;
    epoch2000.tm_mon  = 0;
    epoch2000.tm_mday = 1;
    epoch2000.tm_hour = 11-timezone/3600;
    epoch2000.tm_min  = 58;
    epoch2000.tm_sec  = 56;

    J2000 = mktime(&epoch2000);

    initialized = TRUE;
  }

  /*
   * Elapsed time since the J2000 epoch [day].
   */
  tJ2000 = difftime(date,J2000)/86400.;

  switch(planet->index) {
    case MARS_INDEX:
      /*
       * Allison M, 1997, Accurate analytical representations of solar time and seasons on Mars
       * with applications to the Pathfinder/Surveyor missions, Geophys. Res. Lett. 24, 1967-1970.
       */
      Mrad      = (19.41+0.5240212*tJ2000)*DEG;
      alpha_fms = 270.39+0.5240384*tJ2000;
      l_s       = alpha_fms+(10.691+3.7e-7*tJ2000)*sin(   Mrad)
                                            +0.623*sin(2.*Mrad)
                                            +0.050*sin(3.*Mrad)
                                            +0.005*sin(4.*Mrad);
    break;
    default:
      /*
       * Use Keplerian elements to determine L_s.
       */

      /* Calculate the mean anomaly [deg] */
      M = (planet->mean_lon-planet->lon_perihelion)+(360./planet->orbit_period)*(tJ2000/365.25);

      /* Solve Kepler's equation for the eccentric anomaly, E [deg]. */
      E = eccentric_anomaly(M,planet->e);

      /* Solve for the true anomaly, nu [deg]. */
      nu = true_anomaly(E,planet->e);

      /* Calculate L_s */
      l_s = nu-planet->vernal_equinox_anomaly;
    break;
  }

  /* Map l_s to [0.,360.]. */
  l_s = 360.*modf(1.+modf(l_s/360.,&tmp),&tmp);

  return l_s;
}


/*====================== end of solar_longitude() ==================*/

/*====================== kepler_zero() =============================*/

inline double kepler_zero(double E)
{
  return (Kepler_M-E)*DEG+Kepler_e*sin(E*DEG);
}

/*====================== end of kepler_zero() ======================*/

/*====================== eccentric_anomaly() =======================*/

/*
 * Solve Kepler's equation, M = E - e sin E, for eccentric anomaly, E,
 * given the mean anomaly, M, and the orbital eccentricity, e.
 */

double eccentric_anomaly(double M,
                         double e)
{
  double
    E,
    tol = 1.e-6,
    tmp;
  int
    error_flag;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="eccentric_anomaly";

  /* Map M to [0.,360.]. */
  M = 360.*modf(1.+modf(M/360.,&tmp),&tmp);

  /*
   * Handle special cases.
   */
  if (fcmp(M,  0.) == 0) return   0.;
  if (fcmp(M,180.) == 0) return 180.;
  if (fcmp(M,360.) == 0) return 360.;

  /* Set global variable to communicate with function kepler_zero(). */
  Kepler_M = M;
  Kepler_e = e;

  error_flag = find_root(M-180.,M+180.,tol,&E,kepler_zero);
  if (error_flag) {
    if (error_flag) {
      sprintf(Message,"Error solving Kepler's equation, find_root(): error_flag=%d",error_flag);
      epic_error(dbmsname,Message);
    }
  }

  /* Map E to [0.,360.]. */
  E = 360.*modf(1.+modf(E/360.,&tmp),&tmp);

  return E;
}

/*====================== end of eccentric_anomaly() ================*/

/*====================== true_anomaly() ============================*/

/*
 * Returns true anomaly, nu [deg], given eccentric anomaly, E [deg],
 * and orbital eccentricity, e.
 */

double true_anomaly(double E,
                    double e)
{
  double
    Erad_2,
    nu,
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="true_anomaly";

  /* Map E to [0.,360.]. */
  E = 360.*modf(1.+modf(E/360.,&tmp),&tmp);

  /*
   * Handle special cases.
   */
  if (fcmp(E,  0.) == 0) return   0.;
  if (fcmp(E,180.) == 0) return 180.;
  if (fcmp(E,360.) == 0) return 360.;

  if (fcmp(e,1.) < 0) {
    if (E < 180.) {
      nu = 2.*atan(sqrt((1.+e)/(1.-e))*tan(E*.5*DEG))/DEG;
    }
    else {
      nu = 2.*atan(sqrt((1.+e)/(1.-e))*tan((E-360.)*.5*DEG))/DEG+360.;
    }
  }
  else {
    sprintf(Message,"e=%g >= 1.",e);
    epic_error(dbmsname,Message);
  }

  return nu;
}

/*====================== end of true_anomaly() =====================*/

/* * * * * * * * * * * * end of epic_funcs_astron.c * * * * * * * * */
