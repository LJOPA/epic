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

/* * * * * * * * epic_sensible_heating.c * * * * * * * * * * * * * * 
 *                                                                 *
 *  This file includes functions used to calculate the             *
 *  sensible heating (temperature changing heating) contributions  *
 *  to the interface heating rate, HEAT3(K,J,I).                   *
 *  The heating units are [J/kg/s].                                *
 *                                                                 *
 *       newtonian_cooling()                                       *
 *       perturbation_heating()                                    *
 *       t_rad()                                                   *
 *       temp_eq()                                                 *
 *       solar_insolation()					   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*====================== newtonian_cooling() ================================*/

/*
 * Force temperature profile to the given radiative-equilibrium profile
 * at an appropriate rate.  Keep the layer-average heating zero if requested.
 *
 *
 * NOTE: Lateral boundary conditions, bc_lateral(), are not applied to 
 *       HEAT3 array here.
 */

void newtonian_cooling(planetspec *planet)
{
  register int
    K,J,I;
  static int
    initialized = FALSE;
  register EPIC_FLOAT
    lat,
    pressure,
    temperature,
    fpara,
    cp_over_time,
    da,heat_avg;
  static EPIC_FLOAT
   *buffji;
  EPIC_FLOAT
    t_eq,
    area,heat_area,
    tmp;

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="newtonian_cooling";

  if (!initialized) {
    /* 
     * Allocate memory.
     */
    buffji = fvector(0,Nelem2d-1,dbmsname);

    initialized = TRUE;
  }

  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      lat = grid.lat[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        pressure     = P3(K,J,I);
        temperature  = T3(K,J,I);
        if (var.fpara.on) {
          fpara = FPARA(K,J,I);
        }
        else {
          fpara = 0.25;
        }
        cp_over_time = return_cp(planet,fpara,pressure,temperature)/t_rad(planet,K,J,I);
        /*
         * Determine T_eq:
         */
        switch(planet->index) {
          case HELD_SUAREZ_INDEX:
            t_eq = temp_eq(planet,lat,pressure,0.);
          break;
          case VENUS_LLR05_INDEX:
            t_eq = temp_eq(planet,lat,pressure,0.);
          break;
          case TITAN_INDEX:
            t_eq = temp_eq(planet,lat,pressure,0.);
          break;
          default:
            /*
             * Use t_vs_p data for radiative equilibrium profile.
             */
            get_sounding(planet,pressure,"temperature",&t_eq);
          break;
        }
        BUFFJI(J,I) = -cp_over_time*(temperature-t_eq);
      }
    }

    /*
     * Make layer average in sponge zero, and do this for the rest
     * of the model if requested.
     */
    if (K <= grid.k_sponge || grid.newt_cool_adjust) {
      /*
       * Ensure that layer average is zero.
       */
      area      = 0.;
      heat_area = 0.;
      for (J = JLO; J <= JHI; J++) {
        da = 1./grid.mn[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          area      += da;
          heat_area += BUFFJI(J,I)*da;
        }
      }

#if defined(EPIC_MPI)
      tmp = area;
      MPI_Allreduce(&tmp,&area,1,float_type,MPI_SUM,para.comm);
      tmp = heat_area;
      MPI_Allreduce(&tmp,&heat_area,1,float_type,MPI_SUM,para.comm);
#endif

      heat_avg = heat_area/area;

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          BUFFJI(J,I) -= heat_avg;
        }
      }
    }

    /*
     * Apply result to HEAT3 array.
     */
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        HEAT3(K,J,I) += BUFFJI(J,I); 
      }
    }
  }

  return;
}

/*====================== end of newtonian_cooling() =========================*/

/*====================== perturbation_heating() =============================*/

/* 
 * Apply an initial heating perturbation. The user should modify
 * the heating pattern as required.
 */

void perturbation_heating(planetspec *planet)

{
  int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="perturbation_heating";

  return;
}


/*====================== end of perturbation_heating() ======================*/

/*======================= t_rad() ===========================================*/

/*
 * Return radiative cooling time [s], given position K,J,I.
 * The value is for the bottom interface of layer K.
 */

EPIC_FLOAT t_rad(planetspec *planet,
                 int         K,
                 int         J,
                 int         I)
{
  register int
    ki;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    ka,ks;
  EPIC_FLOAT
    pressure,
    t_cool,
    p_t_cool_d,   
    neglogp;

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="t_rad";

  /*
   * Check that Newtonian cooling is turned on.
   */
  if (grid.newt_cool_on != TRUE) {
    sprintf(Message,"newt_cool_on == %d",grid.newt_cool_on);
    epic_error(dbmsname,Message);
  }

  if (!initialized) {
    if (strcmp(planet->name,"held_suarez") == 0) {
      /*
       * Set parameters ka, ks:
       */
      ka = 1./(40.*60.*60.*24.);
      ks = 1./( 4.*60.*60.*24.);
    }
    else {
      spline_pchip(var.n_t_cool,var.t_cool_table);
    }
    initialized = TRUE;
  }
  /* End of initialization. */

  /*
   * Use interface pressure.
   */
  pressure = P3(K,J,I);

  if (strcmp(planet->name,"held_suarez") == 0) {
    EPIC_FLOAT
      kt,amp,
      pbot,cos_lat;

    kt   = ka;
    pbot = P3(grid.nk,J,I);
    amp  = (pressure/pbot-0.7)/(1.-0.7);
    if (amp > 0.) {
      cos_lat = cos(grid.lat[2*J+1]*DEG);
      kt     += (ks-ka)*pow(cos_lat,4.)*amp;
    }
    t_cool = 1./kt;
  }
  else {
    /*
     *  Interpolate to get Newtonian cooling time:
     */
    neglogp = -log(pressure);
    if (neglogp >= var.t_cool_table[var.n_t_cool-1].x) {
      /* past range of t_cool data */
      t_cool = var.t_cool_table[var.n_t_cool-1].y;
    }
    else if (neglogp <= var.t_cool_table[0].x) {
      /* past range of t_cool data */
      t_cool = var.t_cool_table[0].y;
    }
    else {
      ki     = find_place_in_table(var.n_t_cool,var.t_cool_table,neglogp,&p_t_cool_d);
      t_cool = splint_pchip(neglogp,var.t_cool_table+ki,p_t_cool_d);
    }
  }

  return t_cool;
}

/*======================= end of t_rad() ====================================*/

/*======================= temp_eq() =========================================*/

/*
 * Inputs: pressure [Pa], latitude [deg], time [sec].
 *
 * time = 0.0 corresponds to northern spring equinox.
 */
EPIC_FLOAT temp_eq(planetspec *planet,
                   EPIC_FLOAT  latitude,
                   EPIC_FLOAT  pressure,
                   EPIC_FLOAT  time)
{
  int
    ki;
  static int
    initialized=FALSE;
  register EPIC_FLOAT
    t_tp,dt_tp;
  EPIC_FLOAT
    ans,
    neglogp,p_d,
    sin_y_j,sin_y;
  static EPIC_FLOAT
    dt_y,dth_z,p0;
  static float_triplet
    *t_table,
    *dt_table;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="temp_eq";

  /*
   * Check that Newtonian cooling is turned on.
   */
  if (grid.newt_cool_on != TRUE) {
    sprintf(Message,"newt_cool_on=%d",grid.newt_cool_on);
    epic_error(dbmsname,Message);
  }

  if (!initialized) {
    if (strcmp(planet->name,"held_suarez") == 0) {
      /*
       * Set parameters for the held_suarez test case.
       */
      dt_y  = 60.;
      dth_z = 10.;
      p0    = 1000.*100.;
    }
    else {
      if (var.ntp <= 0) {
        sprintf(Message,"var.ntp=%d",var.ntp);
        epic_error(dbmsname,Message);
      }

      /* Allocate memory */
      t_table  = ftriplet(0,var.ntp-1,dbmsname);
      dt_table = ftriplet(0,var.ntp-1,dbmsname);

      /* Assign table values */
      for (ki = 0; ki < var.ntp; ki++) { 
        t_table[ ki].x = dt_table[ ki].x = -log(var.pdat[ki]);
        t_table[ ki].y = var.tdat[ ki];
        dt_table[ki].y = var.dtdat[ki];
      }
      spline_pchip(var.ntp, t_table);
      spline_pchip(var.ntp,dt_table); 
    }

    initialized = TRUE;
  }
  /* end of initialization */

  if (strcmp(planet->name,"held_suarez") == 0) {
    EPIC_FLOAT
      cos2_lat,
      sin2_lat,
      p_p0;

    cos2_lat  = cos(latitude*DEG);
    cos2_lat *= cos2_lat;
    sin2_lat  = 1.-cos2_lat;
    p_p0      = pressure/p0;
    ans       = (315.-dt_y*sin2_lat-dth_z*log(p_p0)*cos2_lat)*pow(p_p0,planet->kappa);
    ans       = MAX(200.,ans);
  }
  else if (strcmp(planet->name,"venus_llr05") == 0) {
    /*
     *  Interpolate on data table:
     */
    neglogp = -log(pressure);
    if (neglogp >= t_table[var.ntp-1].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[ var.ntp-1].y;
      dt_tp = dt_table[var.ntp-1].y;
    }
    else if (neglogp <= t_table[0].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[0].y;
      dt_tp = dt_table[0].y;
    }
    else {
      ki   = find_place_in_table(var.ntp,t_table,neglogp,&p_d);
      t_tp  = splint_pchip(neglogp, t_table+ki,p_d);
      dt_tp = splint_pchip(neglogp,dt_table+ki,p_d);
    }
    /*
     * Venus forcing of Lee, Lewis, and Read (2005, Adv. Space Res. 36, 2142-2145).
     */
    ans = t_tp+dt_tp*(cos(latitude*DEG)-0.642);
  }
  else if (strcmp(planet->name,"titan") == 0) {
    /*
     *  Interpolate on data table:
     */
    neglogp = -log(pressure);
    if (neglogp >= t_table[var.ntp-1].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[ var.ntp-1].y;
      dt_tp = dt_table[var.ntp-1].y;
    }
    else if (neglogp <= t_table[0].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[0].y;
      dt_tp = dt_table[0].y;
    }
    else {
      ki   = find_place_in_table(var.ntp,t_table,neglogp,&p_d);
      t_tp  = splint_pchip(neglogp, t_table+ki,p_d);
      dt_tp = splint_pchip(neglogp,dt_table+ki,p_d);
    }
    /*
     * Titan forcing as in Flasar et al. (1981, Nature 292, 693-698).
     */
    ans = t_tp+dt_tp*cos(latitude*DEG);
  }
  else {
    sprintf(Message,"not yet implemented for %s",planet->name);
    epic_error(dbmsname,Message);
  }
   
  return ans;
}

/*======================= end of temp_eq() ==================================*/

/*======================= solar_insolation() ================================*/

/*
 * Wrapper to call planet-specific solar insolation function.
 */
void solar_insolation(planetspec *planet)
{
  if (strcmp(planet->name,"uranus") == 0) {
    uranus_solar_insolation(planet);
    return;
  }
  else {
    return;
  }
}

/*======================= end solar_insolation() ============================*/

/*======================= uranus_solar_insolation() =========================*/

/*
 * Apply solar heating to HEAT3.
 * Latitude in degrees.
 *
 * Michael Sussman.
 */
 
void uranus_solar_insolation(planetspec *planet)
{
  register int
    K,J,I;
  static int
    initialized=FALSE,
    nk, nj, ni;
  register EPIC_FLOAT
    lat,
    lon,
    add_heat;
  EPIC_FLOAT
    subsol_lat,
    subsol_lon,
    cos_theta,
    deltap,
    scale_height,
    this_scale_height,
    unit_conversion,
    absorbed,
    absorbed_heat,
    total_heat;


#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="uranus_solar_insolation";

  /*
   * Set arbitrary pressure level of max absorption = g/alpha_nu
   */
  scale_height=2000.*100.; /* 1000 hPa */
 
  /*
   * Find the subsolar point in latitude and longitude as f(t), in degrees.
   * Subsolar latitude is currently hardcoded.
   * Floor function is there to keep subsol_lon between 0 and 360.
   */

  subsol_lat = 90.;
  subsol_lon = 0.; /*NOTE: need subsol_lon */

  /*
   * Loop over our grid with k as the inner loop since we
   * propagate solar heating at a given i,j through each layer
   */

  for (I = ILO; I <= IHI; I++) {
    for (J = JLO; J <= JHI; J++) {
      lat = grid.lat[2*J+1];
      lon = grid.lon[2*I+1];
    /*
     * Find cos_theta, the cosine of the angle between this point, and
     * the subsolar point.
     */
   
      cos_theta = (cos((subsol_lon-lon)*DEG)*cos(lat*DEG)*cos(subsol_lat*DEG)) + (sin(lat*DEG)*sin(subsol_lat*DEG));

    /* 
     *  Set up total heat propagated through each i,j column. At Uranus,
     *  solar constant is 1.515 W/m2.  That must multiplied by the
     *  cosine of theta.  Greatest heating is at sub-solar point, 
     *  zero at theta >= 90 degrees (i.e. cos(theta) < 0).
     *  Additionally, the scale height (pressure of max absorption)
     *  rises as cos_theta, as well.
     */
      total_heat        = 1.515*cos_theta;
      this_scale_height = scale_height*cos_theta;

      if (cos_theta >= 0.0) { 
        for (K = KLO; K <= KHI; K++) {
          deltap = P3(K,J,I)-P3(K-1,J,I);
     
    /* Insolation equations:
     * absorbed = % light absorbed by this layer as a funtion of dP
     *          = 1 - e^(-dP/H)
     * Account for J/kg s units, then remove absorbed light from total_heat 
     */
          absorbed         = 1.0-exp(-deltap/this_scale_height);
          unit_conversion  = 1.0/(H2(K,J,I)*grid.dsgth[2*K]);
          absorbed_heat    = total_heat*absorbed;
          total_heat      -= absorbed_heat;
          HEAT3(K,J,I)    += absorbed_heat*unit_conversion;
        }
      }
    }
  }

  return;
}
/*======================= end of uranus_solar_insolation() ==================*/

/************************ end of epic_sensible_heating.c *********************/
