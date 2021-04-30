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

/* * * * * * * * * * epic_timestep.c * * * * * * * * * * * * * * * * 
 *                                                                 *
 * Integrate the prognostic variables ahead one timestep.          *
 *                                                                 *
 * This file includes the following:                               *
 *                                                                 *
 *           timestep()                                            *
 *           adams_bashforth_step()                                *
 *           leapfrog_step()                                       *
 *           uv_core()                                             *
 *           uv_pgrad()                                            *
 *           uv_drag()                                             *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>
#include <epic_pv_schemes.h>

/*
 * Local function prototypes.
 */
void uv_pgrad_traditional(  planetspec *planet);
void uv_pgrad_finite_volume(planetspec *planet);


/*======================= timestep() ========================================*/

/*
 * March the prognostic variables ahead one timestep using the 
 * specified schemes.
 *
 * If action == STEP_PROGS_AND_UPDATE_DIAGS, we assume the diagnostic variables
 * are all up-to-date upon entry.
 *
 * If action == UPDATE_DIAGS_ONLY, we bring the diagnostic variables up-to-date
 * with the prognostic variables, but do not advance the progostic variables.
 * This is used to prime the pump at startup.
 */

void timestep(planetspec  *planet,
              int          action,
              EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    i,itmp,
    is,ip,iq;
  register EPIC_FLOAT
    tmp;
  static unsigned long
    nbytes_2d,
    nbytes_3d;
  static EPIC_FLOAT
    *gz;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="timestep";

  if (!initialized) {
    nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);
    nbytes_3d = Nelem3d*sizeof(EPIC_FLOAT);

    gz = fvector(0,Nelem3d-1,dbmsname);

    if (grid.nmt_physics_on == 1) {
      /*
       * Allocate memory for nmt_physics.
       */
      nmt_allocate();
    }

    initialized = TRUE;
  }

  if (action == STEP_PROGS_AND_SYNC_DIAGS) {
    /*
     *  Clear u,v tendencies for current time.
     */
    memset(var.u.tendency+grid.it_uv_tend*Nelem3d,0,nbytes_3d);
    memset(var.v.tendency+grid.it_uv_tend*Nelem3d,0,nbytes_3d);

    if (grid.nmt_physics_on == 1) {
     /*
      * Pass diagnostics to NMT cumulus parameterization.
      *
      * Question: nmt_cumulus(): Do we need to separate out sync'ing of diagnostic variables and advancing of prognostic ones?
      */
      nmt_cumulus(planet);
    }

    /*
     * Apply sources and sinks to h-grid variables.
     */
    source_sink(planet,Buff2D);

    /*
     * The HEAT3 array is now complete with the current latent heating.
     * Update the hybrid vertical velocity, W3, on layer interfaces.
     */
    calc_w(planet);

    /*
     * Apply sources and sinks to turbulence variables.
     */
    if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
      source_sink_turb(planet,Buff2D);
    }
    else if (strcmp(grid.turbulence_scheme,"off") == 0) {
      ;
    }
    else {
      sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
      epic_error(dbmsname,Message);
    }

    /*
     * Calculate Coriolis and advection tendencies for (u,v).
     */
    uv_core(planet,Buff2D);

    /*
     * Add Rayleigh drag terms to (u,v) tendencies.
     * For the leapfrog timestep, lagging is used for stability.
     */
    uv_drag(planet,Buff2D);

    /*
     * Apply subgrid-scale model to scalar variables.
     */
    scalar_horizontal_subgrid(planet,Buff2D);
    scalar_vertical_subgrid(planet,Buff2D);

    /*
     * Advect scalar prognostic variables forward one timestep.
     */
    advection(planet,Buff2D);

    /*
     * Update GZ3.
     */
    store_pgrad_vars(planet);

    /*
     * Form GZ, a high-latitude filtered version of GZ3.
     * Without this step, noise leaks around the filter when GZ is used
     * to subtract a function of z from the variable.
     */
    memcpy(gz,var.gz3.value,Nelem3d*sizeof(EPIC_FLOAT));
    for (K = KHI-1; K >= KLO-1; K--) {
      /*
       * Use TWODIM, which does not subtract a function of GZ from the variable.
       */
      zonal_filter(GZ3_INDEX,gz+(K-Kshift)*Nelem2d,NULL,TWODIM);
      /*
       * Ensure that GZ is monotonic.
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          tmp = GZ(K+1,J,I);
          if (fcmp(tmp,0.) == 0) {
            GZ(K,J,I) = 1.e-3;
          }
          else {
            GZ(K,J,I) = MAX(GZ(K,J,I),GZ(K+1,J,I)+fabs(GZ(K+1,J,I))*1.e-5);
          }
        }  
      }
    }
    
    /*
     * Apply high-latitude filter to DUDT and DVDT before uv_pgrad(), to
     * avoid disturbing the irrotational property of the pressure gradient force.
     * Without this filter, noise tends to accumulate near the northern and southern
     * lateral boundaries.
     */
    zonal_filter(U_INDEX,var.u.tendency+grid.it_uv_tend*Nelem3d,gz,THREEDIM);
    zonal_filter(V_INDEX,var.v.tendency+grid.it_uv_tend*Nelem3d,gz,THREEDIM);

    /*
     * Apply zonal_filter() to mass variables and potential temperature.
     */
    if (var.hdry.on) {
      zonal_filter(HDRY_INDEX,var.hdry.value,gz,THREEDIM);
      restore_mass(planet,HDRY_INDEX,NO_PHASE);
    }

    for (iq = 0; iq < grid.nq; iq++) {
      zonal_filter(grid.is[iq],var.species[grid.is[iq]].phase[grid.ip[iq]].q,gz,THREEDIM);
      restore_mass(planet,grid.is[iq],grid.ip[iq]);
    }

    if (var.theta.on) {
      zonal_filter(THETA_INDEX,var.theta.value,gz,THREEDIM);
    }

    if (var.fpara.on) {
      zonal_filter(FPARA_INDEX,var.fpara.value,gz,THREEDIM);
    }

    if (var.nu_turb.on) {
      zonal_filter(NU_TURB_INDEX,var.nu_turb.value,gz,THREEDIM);
      restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
    }

    /*
     * Update the diagnostic variables needed for uv_pgrad().
     */
    set_p2etc(planet,UPDATE_THETA);
    store_pgrad_vars(planet);

    uv_pgrad(planet);

    /*
     * Add horizontal subgrid-scale model to wind tendencies.
     * The vertical effects are done implicitly on (u,v) below.
     */
    uv_horizontal_subgrid(planet,Buff2D);

    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      /*
       * March u,v forward one timestep.
       */
      adams_bashforth_step(planet,U_INDEX);
      adams_bashforth_step(planet,V_INDEX);

      /*
       *  Cycle time index backwards:
       */
      itmp      = IT_MINUS2;
      IT_MINUS2 = IT_MINUS1;
      IT_MINUS1 = IT_ZERO;
      IT_ZERO   = itmp;
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      /*
       * March u,v forward one timestep.
       */
      leapfrog_step(planet,U_INDEX);
      leapfrog_step(planet,V_INDEX);

      /* 
       * Cycle time index backwards.
       */
      itmp      = IT_MINUS1;
      IT_MINUS1 = IT_ZERO;
      IT_ZERO   = itmp;
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme = %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }

    timeplane_bookkeeping();

    /*
     * Apply vertical subgrid-scale model to U and V.
     * This is done implicitly to handle thin layers.
     */
    uv_vertical_subgrid(planet,Buff2D);

    /*
     * Apply horizontal hyperviscosity to U and V.
     *
     * NOTE: We tried applying hyperviscosity to the tendences, DUDT and DVDT,
     *       but grid-scale computational modes (stripes and checkerboards)
     *       emerged in U and V anyway.
     */
    meridional_hyperviscosity(planet,U_INDEX,var.u.value+grid.it_uv*Nelem3d,THREEDIM,Buff2D);
    meridional_hyperviscosity(planet,V_INDEX,var.v.value+grid.it_uv*Nelem3d,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,U_INDEX,var.u.value+grid.it_uv*Nelem3d,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,V_INDEX,var.v.value+grid.it_uv*Nelem3d,THREEDIM,Buff2D);

    /*
     *  Advance time:
     */
    var.model_time += (time_t)grid.dt;

    /*
     * Update solar longitude, L_s.
     */
    L_s = solar_longitude(planet,var.model_time);

    /*
     * Store commonly used diagnostic variables not calculated elsewhere.
     */
    store_diag(planet);

    /*
     * Start the calculation of HEAT3, in W/kg, on the layer interfaces (all except
     * the latent heating).
     */
    calc_heating(planet);

    /*
     * NOTE: The hybrid vertical velocity, W3, is not updated here and hence still refers to the previous
     *       timestep, whereas the other diagnostic variables now refer to the new timestep.  The problem
     *       is that the latent heating is hard to isolate from the microphysical phase changes.  The 
     *       consequence is that the W3 extracted to extract.nc files lags behind the other variables by
     *       one timestep; we believe this will not cause problems in most cases.
     */
  } 
  else if (action == SYNC_DIAGS_ONLY) {
    /*
     * NOTE: These are updated in the STEP_PROGS_AND_SYNC_DIAGS case as part of uv_pgrad(),
     *       to facilitate the "poor-man's implicit timestep."
     */
    set_p2etc(planet,UPDATE_THETA);
    store_pgrad_vars(planet);

    /*
     * Store commonly used diagnostic variables not calculated elsewhere.
     */
    store_diag(planet);

    /*
     * Calculate HEAT3, in W/kg, on the layer interfaces.
     *
     * NOTE: This currently does not include latent heating from cloud microphysics.
     */
    calc_heating(planet);

    /*
     * Calculate the hybrid vertical velocity, W3, on layer interfaces.
     */
    calc_w(planet);
  }
  else {
    sprintf(Message,"unrecognized action=%d",action);
    epic_error(dbmsname,Message);
  }

  /*
   * Put below functions common to all values of the "action" argument that 
   * may be safely executed at the end of the timestep.
   */

  if (grid.nmt_physics_on == 1) {
   /*
    * Calculate entropies, since the core diagnostics have just been updated.
    */
    nmt_entropy();
  }

  /*
   * Map Q <= Q_MIN -> 0.0, so that output files hold zero
   * instead of Q_MIN for voids.
   */
  for (iq = 0; iq < grid.nq; iq++) {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          if (fcmp(Q(grid.is[iq],grid.ip[iq],K,J,I),Q_MIN) <= 0) {
            Q(grid.is[iq],grid.ip[iq],K,J,I) = 0.;
          }
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  return;
}

/*======================= end of timestep() =================================*/

/*======================= adams_bashforth_step() ============================*/

/*
 * Use the 3rd-order Adams-Bashforth timestep.
 * This timestep is discussed by D. Durran (1991, MWR 119, 702-720).
 * It is appropriate for dissipative terms as well as for
 * conservative terms, does not suffer from the time-splitting
 * numerical instability of the leap-frog timestep, and is more
 * accurate than the leapfrog timestep.  Its main drawback is that
 * it requires more memory because it uses two previous time derivatives.
 */

void adams_bashforth_step(planetspec *planet,
                          int         index)
{
  register int
    K,J,I,
    jlo;
  wind_variable
    *wind;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="adams_bashforth_step";

  /*
   * The index jlo, as opposed to JLO, arises because U and V are on different grids
   * in the Arakawa C-grid scheme.
   */
  if (index == U_INDEX) {
    wind = &var.u;
    jlo  = JLO;
  }
  else if (index == V_INDEX) {
    wind = &var.v;
    jlo  = JFIRST;
  }
  else {
    sprintf(Message,"not set up to handle index=%d",index);
    epic_error(dbmsname,Message);
  }

  /* 
   * Specify Adams-Bashforth coefficients.
   * The FLOAT_MAX flags are set in epic_initial.c. 
   *
   * NOTE: Use jlo rather than JLO for this test.
   */
  if (DWINDDT(wind,IT_MINUS2,KLO,jlo,ILO) == FLOAT_MAX) {
    if (DWINDDT(wind,IT_MINUS1,KLO,jlo,ILO) == FLOAT_MAX) {
      /* 
       * Use 1st-order Adams-Bashforth, aka forward (Euler) difference,
       * for the initial step.
       */
      grid.ab[0] = 1.;
      grid.ab[1] = 0.;
      grid.ab[2] = 0.;
    }
    else {
      /*
       * Use 2nd-order Adams-Bashforth for the second step.
       */
      grid.ab[0] =  3./2.;
      grid.ab[1] = -1./2.;
      grid.ab[2] =     0.;
    }
  }
  else {
   /*
    * Use 3rd-order Adams-Bashforth for the third step and beyond.
    */
    grid.ab[0] =  23./12.;
    grid.ab[1] = -16./12.;
    grid.ab[2] =   5./12.;
  }

  /*
   * Advance variable.
   */
  for (K = KLO; K < KHI; K++) {
    for (J = jlo; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        WIND(wind,grid.it_uv,K,J,I) += DT*( grid.ab[0]*DWINDDT(wind,IT_ZERO,  K,J,I)
                                           +grid.ab[1]*DWINDDT(wind,IT_MINUS1,K,J,I)
                                           +grid.ab[2]*DWINDDT(wind,IT_MINUS2,K,J,I) );
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  if (index == U_INDEX) {
    bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  }
  else if (index == V_INDEX) {
    bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);
  }
  else {
    sprintf(Message,"unrecognized index=%d",index);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of adams_bashforth_step() =====================*/

/*======================= leapfrog_step() ===================================*/

/*
 * Apply the leapfrog step. 
 * Use Asselin filtering to remove the odd-even instability.
 *
 * Durran (1991) cites Williamson (1983) for the value gamma = 0.06
 * for the Asselin filter parameter.
 */

#define GAMMA_ASSELIN 0.06

void leapfrog_step(planetspec *planet,
                   int         index)
{
  register int
    K,J,I,
    jlo;
  register EPIC_FLOAT
    twodt,
    old_filtered,
    present,
    new_leap;
  wind_variable
    *wind;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="leapfrog_step";

  /*
   * The index jlo, as opposed to JLO, arises because U and V are on different grids
   * in the Arakawa C-grid scheme.
   */
  if (index == U_INDEX) {
    wind = &var.u;
    jlo  = JLO;
  }
  else if (index == V_INDEX) {
    wind = &var.v;
    jlo  = JFIRST;
  }
  else {
    sprintf(Message,"not set up to handle index=%d",index);
    epic_error(dbmsname,Message);
  }

  /* 
   * Use a forward (Euler) difference for the initial step.
   */

  /*
   * NOTE: Need jlo, not JLO for this test, otherwise the flagged V=FLOAT_MAX will be V=0.
   */
  if (WIND(wind,IT_MINUS1,KLO,jlo,ILO) == FLOAT_MAX) {
    for (K = KLO; K < KHI; K++) {
      for (J = jlo; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          present                    = WIND(wind,IT_ZERO,K,J,I);
          WIND(wind,IT_MINUS1,K,J,I) = present+DT*DWINDDT(wind,grid.it_uv_tend,K,J,I);
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    if (index == U_INDEX) {
      bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
    }
    else if (index == V_INDEX) {
      bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized index=%d",index);
      epic_error(dbmsname,Message);
    }

    return;
  }

  twodt = 2.*DT;

  /*
   * Take leapfrog step, store in old timeframe. 
   * Asselin-filter the current step.
   */
  for (K = KLO; K < KHI; K++) {
    for (J = jlo; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        old_filtered               = WIND(wind,IT_MINUS1,K,J,I);
        present                    = WIND(wind,IT_ZERO,  K,J,I);
        new_leap                   = old_filtered+twodt*DWINDDT(wind,grid.it_uv_tend,K,J,I);
        WIND(wind,IT_MINUS1,K,J,I) = new_leap;
        WIND(wind,IT_ZERO,  K,J,I) = present*(1.-2.*GAMMA_ASSELIN)+GAMMA_ASSELIN*(old_filtered+new_leap);
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  if (index == U_INDEX) {
    bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
    bc_lateral(var.u.value+IT_ZERO*Nelem3d,  THREEDIM);
  }
  else if (index == V_INDEX) {
    bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
    bc_lateral(var.v.value+IT_ZERO*Nelem3d,  THREEDIM);
  }
  else {
    sprintf(Message,"unrecognized index=%d",index);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of leapfrog_step() ============================*/

/*======================= uv_core() =========================================*/

/*
 * Calculate core tendencies for u,v for the current state.
 * These include advection and Coriolis acceleration.
 * The the pressure-gradient terms are calculated separately as part
 * of the economical explicit scheme. 
 */

void uv_core(planetspec  *planet,
             EPIC_FLOAT **Buff2D)
{
  register int    
    K,J,I,
    kk;
  unsigned long
    nbytes_2d;
  static int 
    initialized = FALSE;
  register EPIC_FLOAT
    al, be,       /* Used in AL_U, BE_U, etc. macros.         */
    ga, de,       /* See Arakawa and Lamb (1981) eqn. (3.34)  */
    ep1,ep2,      /*       "                     "            */
    ph1,ph2,      /*       "                     "            */
    m_2jp1,
    m_2j_inv,
    n_2j,n_2jp1_inv,
    havg,
    d1,d2,d1d2,davg;
  EPIC_FLOAT
    *uh,*vh,*kin;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_core";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

 /*
  * Calculate the horizontal advection and Coriolis terms for U and V.
  *
  * NOTE: In order to get the correct Coriolis terms, use the same type of 
  *       hybrid density to form UH and VH here as is used in the 
  *       potential vorticity.
  */ 
  for (K = KLO; K < KHI; K++) {
    memset(Buff2D[0],0,nbytes_2d);
    memset(Buff2D[1],0,nbytes_2d);
    memset(Buff2D[2],0,nbytes_2d);
    uh  = Buff2D[0];      
    vh  = Buff2D[1];      
    kin = Buff2D[2];  

    for (J = JLO; J <= JHI; J++) {
      n_2jp1_inv = 1./(grid.n)[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        havg    = .5*(H3(K,J,I)+H3(K,J,I-1));
        UH(J,I) = U(grid.it_uv,K,J,I)*havg*n_2jp1_inv;
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(uh,TWODIM);

    for (J = JFIRST; J <= JHI; J++) {
      m_2j_inv  = 1./(grid.m)[2*J];
      for (I = ILO; I <= IHI; I++) {
        havg    = .5*(H3(K,J,I)+H3(K,J-1,I));
        VH(J,I) = V(grid.it_uv,K,J,I)*havg*m_2j_inv;
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(vh,TWODIM);

    /*
     *  Calculate kinetic energy per unit mass, kin:
     */
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        KIN(J,I) = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                  var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(kin,TWODIM);

    /*
     * Add pv*uh, pv*vh and kin terms to dvdt:
     */
    for (J = JFIRST; J <= JHI; J++) {
      n_2j = (grid.n)[2*J];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) += ((-GA_V*UH(J  ,I+1)-DE_V*UH(J  ,I  )
		                         -AL_V*UH(J-1,I  )-BE_V*UH(J-1,I+1)
                                         +PH2*VH( J-1,I  )-PH1*VH( J+1,I  ))*PV_COEF
                                       +KIN(J-1,I)-KIN(J,I))*n_2j;
      }
    }

    /*
     * Apply pole filter to kin under gradient operator before
     * zonal component, to control numerical instability of dK/dx.
     */
    zonal_filter(HDRY_INDEX,kin,NULL,TWODIM);

    /*
     * Add pv*uh, pv*vh and kin terms to dudt.
     * The coefficients AL_U, BE_U, etc. are linear combinations
     * of PV3(K,J,I) and are defined in $EPIC4_PATH/include/epic_pv_schemes.h.
     */
    for (J = JLO; J <= JHI; J++) {
      m_2jp1 = (grid.m)[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) += ((AL_U*VH(J+1,I  )+BE_U*VH(J+1,I-1)
                                        +GA_U*VH(J  ,I-1)+DE_U*VH(J,  I  )
		                        +EP2*UH( J,  I-1)-EP1*UH( J,  I+1))*PV_COEF
                                      +KIN(J,I-1)-KIN(J,I))*m_2jp1;
      }
    }
  }

  /*
   * Vertical advection of U, V.
   */
  for (K = KLO; K < KHI; K++) {
    kk = 2*K+1;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= .5*(W3(K,J,I)+W3(K,J,I-1))*(U(grid.it_uv,K-1,J,I)-U(grid.it_uv,K+1,J,I))/
                                                                 (grid.sigmatheta[kk-2]-grid.sigmatheta[kk+2]);
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= .5*(W3(K,J,I)+W3(K,J-1,I))*(V(grid.it_uv,K-1,J,I)-V(grid.it_uv,K+1,J,I))/
                                                                  (grid.sigmatheta[kk-2]-grid.sigmatheta[kk+2]);
      }
    }
  }

  return;
}

/*======================= end of uv_core() ======================================*/

/*======================= uv_pgrad() ============================================*/

/*
 * A choice of algorithms is available for the horizontal pressure-gradient
 * force, and is set here by specifying PGRAD_TYPE.
 */
#define PGRAD_TRADITIONAL   0
#define PGRAD_FINITE_VOLUME 1
#define PGRAD_TYPE          PGRAD_FINITE_VOLUME

/*
 * Calculate pressure-gradient tendencies for (u,v).
 * A choice of algorithms is available and is set by PGRAD_TYPE above.
 */

void uv_pgrad(planetspec *planet)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_pgrad";

#if PGRAD_TYPE == PGRAD_TRADITIONAL
  uv_pgrad_traditional(planet);
#elif PGRAD_TYPE == PGRAD_FINITE_VOLUME
  uv_pgrad_finite_volume(planet);
#else
  sprintf(Message,"Unknown PGRAD_TYPE=%d",PGRAD_TYPE);
  epic_error(dbmsname,Message);
#endif
  
}

/*====================== end of uv_pgrad() ========================================*/

/*======================= uv_pgrad_traditional() ==================================*/

void uv_pgrad_traditional(planetspec *planet)
{
  register int    
    K,J,I;
  register EPIC_FLOAT
    m_2jp1,n_2j,
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_pgrad_traditional";

  for (K = KLO; K <= grid.k_sigma-2; K++) {
    for (J = JLO; J <= JHI; J++) {
      m_2jp1  = grid.m[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= m_2jp1*(MONT3(K,J,I)-MONT3(K,J,I-1)-
                                              .5*(EXNER3(K,J,I)+EXNER3(K,J,I-1))*(THETA(K,J,I)-THETA(K,J,I-1)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB3(K,J,I)+FGIBB3(K,J,I-1))
                                           *(FPARA( K,J,I)-FPARA( K,J,I-1))*m_2jp1;
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      n_2j = grid.n[2*J];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= n_2j*(MONT3(K,J,I)-MONT3(K,J-1,I)-
                                            .5*(EXNER3(K,J,I)+EXNER3(K,J-1,I))*(THETA(K,J,I)-THETA(K,J-1,I)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DVDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB3(K,J,I)+FGIBB3(K,J-1,I))
                                           *(FPARA( K,J,I)-FPARA( K,J-1,I))*n_2j;
        }
      }
    }
  }

  for (K = grid.k_sigma-1; K < KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      m_2jp1  = grid.m[2*J+1];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= m_2jp1*( (GZ3(K,J,I)*H3(K,J,I)-GZ3(K,J,I-1)*H3(K,J,I-1))*2./(H3(K,J,I)+H3(K,J,I-1))-
                                               ((GZ3(K+1,J,I)+GZ3(K+1,J,I-1))*(P3(K+1,J,I)-P3(K+1,J,I-1))
                                               -(GZ3(K-1,J,I)+GZ3(K-1,J,I-1))*(P3(K-1,J,I)-P3(K-1,J,I-1)))/
                                               (P3(K+1,J,I)+P3(K+1,J,I-1)-P3(K-1,J,I)-P3(K-1,J,I-1)) );
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB3(K,J,I)+FGIBB3(K,J,I-1))
                                           *(FPARA( K,J,I)-FPARA( K,J,I-1))*m_2jp1;
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      n_2j = grid.n[2*J];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= n_2j*( (GZ3(K,J,I)*H3(K,J,I)-GZ3(K,J-1,I)*H3(K,J-1,I))*2./(H3(K,J,I)+H3(K,J-1,I))-
                                             ((GZ3(K+1,J,I)+GZ3(K+1,J-1,I))*(P3(K+1,J,I)-P3(K+1,J-1,I))
                                             -(GZ3(K-1,J,I)+GZ3(K-1,J-1,I))*(P3(K-1,J,I)-P3(K-1,J-1,I)))/
                                             (P3(K+1,J,I)+P3(K+1,J-1,I)-P3(K-1,J,I)-P3(K-1,J-1,I)) );
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DVDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB3(K,J,I)+FGIBB3(K,J-1,I))
                                           *(FPARA( K,J,I)-FPARA( K,J-1,I))*n_2j;
        }
      }
    }
  }

  return;
}

/*======================= end of uv_pgrad_traditional() ============================*/

/*======================= uv_pgrad_finite_volume() =================================*/

/*
 * Calculate pressure-gradient tendencies for (u,v).
 *
 * Use a finite-volume approach with a 3D control volume (CV) around each u point
 * and each v point.  The 8 corners of each CV can have different GZ values, and the 
 * surface forces and mass are computed accordingly. This method expands the 2D CV
 * approach of Chu and Fan (2003, JGR 108) to 3D CVs.
 *
 * The Mathematica notebooks used to calculate the formulas are included in
 * epic/tools/mathematica.
 */

/*
 * DLOGZ_TOL is the tolerance for dz/z under which we treat the top or
 * bottom face as horizontal, thereby avoiding dividing by small dz.
 */
#undef  DLOGZ_TOL
#define DLOGZ_TOL (1.e-2)

void uv_pgrad_finite_volume(planetspec *planet)
{
  register int    
    K,J,I;
  register EPIC_FLOAT
    rho,
    lnw,lne,lts,ltn,
    rlt,rlns,rlnn,
    pgflon,pgflat,
    zbsw,zbse,zbnw,zbne,zbw,zbe,
    ztsw,ztse,ztnw,ztne,ztw,zte,
    zbs,zbn,zts,ztn,zb,zt,
    pbw,pbe,ptw,pte,
    pbs,pbn,pts,ptn,
    top,bottom,east,west,north,south,
    terms1,terms2,terms3,terms4,
    dzn,dzs,
    tmp,a1,a2,b1,b2,b3,b4;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_pgrad_finite_volume";

  /*
   * Set up boundary conditions for GZ3.  These are necessary because of the
   * extent of the finite volumes around the U and V points.
   */
  if (JLO == grid.jlo) {
    /*
     * Processor contains the southern edge of the model.
     */
    if (IS_SPOLE) {
      /*
       * For south pole, use average of southern interior edge of GZ.
       */
      for (K = 0; K <= KHI; K++) {
        tmp = 0.;
        for (I = ILO; I <= IHI; I++) {
          tmp += GZ3(K,JLO,I);
        }
        tmp /= grid.ni;
        for (I = ILOPAD; I <= IHIPAD; I++) {
          GZ3(K,JLO-1,I) = tmp;
        }
      }
    }
    else {
      /*
       * For channel edge, use zero-slope b.c.
       */
      for (K = 0; K <= KHI; K++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          GZ3(K,JLO-1,I) = GZ3(K,JLO,I);
        }
      }
    }
  }
  if (JHI == grid.nj) {
    /*
     * Processor contains the northern edge of the model.
     */
    if (IS_NPOLE) {
      /*
       * For north pole, use average of northern interior edge of GZ.
       */
      for (K = 0; K <= KHI; K++) {
        tmp = 0.;
        for (I = ILO; I <= IHI; I++) {
          tmp += GZ3(K,JHI,I);
        }
        tmp /= grid.ni;
        for (I = ILOPAD; I <= IHIPAD; I++) {
          GZ3(K,JHI+1,I) = tmp;
        }
      }
    }
    else {
      /*
       * For channel edge, use zero-slope b.c.
       */
      for (K = 0; K <= KHI; K++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          GZ3(K,JHI+1,I) = GZ3(K,JHI,I);
        }
      }
    }
  }

  /*
   * Longitudinal component
   */
  for (K = KLO; K < KHI; K++) {
    for (J = JLO; J <= JHI; J ++) {
      rlt  = grid.rlt[2*J+1];
      rlns = grid.rln[2*J  ];
      rlnn = grid.rln[2*J+2];

      I = ILO-1;

      lne  = grid.lon[2*I+1]*DEG;

      zbse = .25*(GZ3(K,J,I)+GZ3(K,J-1,I)+GZ3(K+1,J,I)+GZ3(K+1,J-1,I));
      zbne = .25*(GZ3(K,J,I)+GZ3(K,J+1,I)+GZ3(K+1,J,I)+GZ3(K+1,J+1,I));
      ztse = .25*(GZ3(K,J,I)+GZ3(K,J-1,I)+GZ3(K-1,J,I)+GZ3(K-1,J-1,I));
      ztne = .25*(GZ3(K,J,I)+GZ3(K,J+1,I)+GZ3(K-1,J,I)+GZ3(K-1,J+1,I));

      zbe  = .5*(GZ3(K,J,I)+GZ3(K+1,J,I));
      zte  = .5*(GZ3(K,J,I)+GZ3(K-1,J,I));
      pbe  = P2(K+1,J,I);
      pte  = P2(K,  J,I);

      a1     = ztne+ztse;
      b1     = ztne*ztse;
      a2     = zbne+zbse;
      b2     = zbne*zbse;
      terms1 = a1*a1-b1-(a2*a2-b2);
      terms2 = 3.*(a1-a2);
      east   = (-terms1*(pbe-pte)+terms2*(zte*pbe-zbe*pte))/(zte-zbe);

      for (I = ILO; I <= IHI; I++) {
        rho  = .5*(RHO3(K,J,I)+RHO3(K,J,I-1));

        lnw  = lne;

        zbsw = zbse;
        zbnw = zbne;
        zbw  = zbe;
        ztsw = ztse;
        ztnw = ztne;
        ztw  = zte;
        pbw  = pbe;
        ptw  = pte;

        /* West */
        west = east;

        lne  = grid.lon[2*I+1]*DEG;

        zbse = .25*(GZ3(K,J,I)+GZ3(K,J-1,I)+GZ3(K+1,J,I)+GZ3(K+1,J-1,I));
        zbne = .25*(GZ3(K,J,I)+GZ3(K,J+1,I)+GZ3(K+1,J,I)+GZ3(K+1,J+1,I));
        ztse = .25*(GZ3(K,J,I)+GZ3(K,J-1,I)+GZ3(K-1,J,I)+GZ3(K-1,J-1,I));
        ztne = .25*(GZ3(K,J,I)+GZ3(K,J+1,I)+GZ3(K-1,J,I)+GZ3(K-1,J+1,I));

        zbe  = .5*(GZ3(K,J,I)+GZ3(K+1,J,I));
        zte  = .5*(GZ3(K,J,I)+GZ3(K-1,J,I));
        pbe  = P2(K+1,J,I);
        pte  = P2(K,  J,I);

        /* East */
        a1     = ztne+ztse;
        b1     = ztne*ztse;
        a2     = zbne+zbse;
        b2     = zbne*zbse;
        terms1 = a1*a1-b1-(a2*a2-b2);
        terms2 = 3.*(a1-a2);
        east   = (-terms1*(pbe-pte)+terms2*(zte*pbe-zbe*pte))/(zte-zbe);

        /* Top */
        a1     = ztne+ztse;
        a2     = ztnw+ztsw;
        terms2 = 3.*(a1-a2);
        if (fabs(2.*(zte-ztw)/(zte+ztw-zbe-zbw)) > DLOGZ_TOL) {
          b1     = ztne*ztse;
          b2     = ztnw*ztsw;
          terms1 = a1*a1-b1-(a2*a2-b2);
          top    = (-terms1*(ptw-pte)+terms2*(zte*ptw-ztw*pte))/(zte-ztw);
        }
        else {
          /*
           * Need to assume ptw = pte to be consistent with the p = p(z) model.
           */
          top = .5*(ptw+pte)*terms2;
        }

        /* Bottom */
        a1     = zbne+zbse;
        a2     = zbnw+zbsw;
        terms2 = 3.*(a1-a2);
        if (fabs(2.*(zbe-zbw)/(zte+ztw-zbe-zbw)) > DLOGZ_TOL) {
          b1     = zbne*zbse;
          b2     = zbnw*zbsw;
          terms1 = a1*a1-b1-(a2*a2-b2);
          bottom = (-terms1*(pbw-pbe)+terms2*(zbe*pbw-zbw*pbe))/(zbe-zbw);
        }
        else {
          bottom = .5*(pbw+pbe)*terms2;
        }

        /* North */
        north   = 0.;

        /* South */
        south   = 0.;

        pgflon  = top-bottom+south-north+west-east;

        dzn     = ztne-zbne+ztnw-zbnw;
        dzs     = ztse-zbse+ztsw-zbsw;
        pgflon /= (lne-lnw)*rho*(rlns*(dzs+.5*dzn)+rlnn*(dzn+.5*dzs));
        DUDT(grid.it_uv_tend,K,J,I) += pgflon;
      }
    }
  }


  /*
   * Latitudinal component
   */
  for (K = KLO; K < KHI; K++) {
    for (I = ILO; I <= IHI; I++) {
      J = JFIRST-1;

      rlnn  = grid.rln[2*J+1];
      ltn   = grid.lat[2*J+1]*DEG;

      zbnw  = .25*(GZ3(K,J,I)+GZ3(K,J,I-1)+GZ3(K+1,J,I)+GZ3(K+1,J,I-1));
      zbne  = .25*(GZ3(K,J,I)+GZ3(K,J,I+1)+GZ3(K+1,J,I)+GZ3(K+1,J,I+1));
      ztnw  = .25*(GZ3(K,J,I)+GZ3(K,J,I-1)+GZ3(K-1,J,I)+GZ3(K-1,J,I-1));
      ztne  = .25*(GZ3(K,J,I)+GZ3(K,J,I+1)+GZ3(K-1,J,I)+GZ3(K-1,J,I+1));

      zbn   = .5*(GZ3(K,J,I)+GZ3(K+1,J,I));
      ztn   = .5*(GZ3(K,J,I)+GZ3(K-1,J,I));
      pbn   = P2(K+1,J,I);
      ptn   = P2(K,  J,I);

      a1      = ztne+ztnw;
      b1      = ztne*ztnw;
      a2      = zbne+zbnw;
      b2      = zbne*zbnw;
      terms1  = a1*a1-b1-(a2*a2-b2);
      terms2  = 3.*(a1-a2);
      north   = rlnn*(-terms1*(pbn-ptn)+terms2*(ztn*pbn-zbn*ptn))/(ztn-zbn);

      for (J = JFIRST; J <= JHI; J++) {
        rho   = .5*(RHO3(K,J,I)+RHO3(K,J-1,I));
        rlt   = grid.rlt[2*J];

        rlns  = rlnn;
        lts   = ltn;

        zbsw  = zbnw;
        zbse  = zbne;
        ztsw  = ztnw;
        ztse  = ztne;

        south = north;

        zbs   = zbn;
        zts   = ztn;
        pbs   = pbn;
        pts   = ptn;

        rlnn  = grid.rln[2*J+1];
        ltn   = grid.lat[2*J+1]*DEG;

        zbnw  = .25*(GZ3(K,J,I)+GZ3(K,J,I-1)+GZ3(K+1,J,I)+GZ3(K+1,J,I-1));
        zbne  = .25*(GZ3(K,J,I)+GZ3(K,J,I+1)+GZ3(K+1,J,I)+GZ3(K+1,J,I+1));
        ztnw  = .25*(GZ3(K,J,I)+GZ3(K,J,I-1)+GZ3(K-1,J,I)+GZ3(K-1,J,I-1));
        ztne  = .25*(GZ3(K,J,I)+GZ3(K,J,I+1)+GZ3(K-1,J,I)+GZ3(K-1,J,I+1));

        zbn   = .5*(GZ3(K,J,I)+GZ3(K+1,J,I));
        ztn   = .5*(GZ3(K,J,I)+GZ3(K-1,J,I));
        pbn   = P2(K+1,J,I);
        ptn   = P2(K,  J,I);

        /* North */
        a1     = ztne+ztnw;
        b1     = ztne*ztnw;
        a2     = zbne+zbnw;
        b2     = zbne*zbnw;
        terms1 = a1*a1-b1-(a2*a2-b2);
        terms2 = 3.*(a1-a2);
        north  = rlnn*(-terms1*(pbn-ptn)+terms2*(ztn*pbn-zbn*ptn))/(ztn-zbn);

        /* South */
        a1     = ztse+ztsw;
        b1     = ztse*ztsw;
        a2     = zbse+zbsw;
        b2     = zbse*zbsw;
        terms1 = a1*a1-b1-(a2*a2-b2);
        terms2 = 3.*(a1-a2);
        south  = rlns*(-terms1*(pbs-pts)+terms2*(zts*pbs-zbs*pts))/(zts-zbs);

        zbe    = .125*(GZ3(K,  J,I)+GZ3(K,  J,I+1)+GZ3(K,  J-1,I)+GZ3(K,  J-1,I+1)
                      +GZ3(K+1,J,I)+GZ3(K+1,J,I+1)+GZ3(K+1,J-1,I)+GZ3(K+1,J-1,I+1));
        zte    = .125*(GZ3(K-1,J,I)+GZ3(K-1,J,I+1)+GZ3(K-1,J-1,I)+GZ3(K-1,J-1,I+1)
                      +GZ3(K,  J,I)+GZ3(K,  J,I+1)+GZ3(K,  J-1,I)+GZ3(K,  J-1,I+1));
        pbe    =  .25*(P2( K+1,J,I)+P2( K+1,J,I+1)+P2( K+1,J-1,I)+P2( K+1,J-1,I+1));
        pte    =  .25*(P2( K,  J,I)+P2( K,  J,I+1)+P2( K,  J-1,I)+P2( K,  J-1,I+1));   

        /* East */
        a1     = ztne+ztse;
        b1     = ztne*ztse;
        a2     = zbne+zbse;
        b2     = zbne*zbse;
        terms1 = a1*a1-b1-(a2*a2-b2);
        terms2 = 3.*(a1-a2);
        east   = .5*(rlns-rlnn)*(-terms1*(pbe-pte)+terms2*(zte*pbe-zbe*pte))/(zte-zbe);

        zbw    = .125*(GZ3(K,  J,I)+GZ3(K,  J,I-1)+GZ3(K,  J-1,I)+GZ3(K,  J-1,I-1)
                      +GZ3(K+1,J,I)+GZ3(K+1,J,I-1)+GZ3(K+1,J-1,I)+GZ3(K+1,J-1,I-1));
        ztw    = .125*(GZ3(K-1,J,I)+GZ3(K-1,J,I-1)+GZ3(K-1,J-1,I)+GZ3(K-1,J-1,I-1)
                      +GZ3(K,  J,I)+GZ3(K,  J,I-1)+GZ3(K,  J-1,I)+GZ3(K,  J-1,I-1));
        pbw    =  .25*(P2( K+1,J,I)+P2( K+1,J,I-1)+P2( K+1,J-1,I)+P2( K+1,J-1,I-1));
        ptw    =  .25*(P2( K,  J,I)+P2( K,  J,I-1)+P2( K,  J-1,I)+P2( K,  J-1,I-1));

        /* West */
        a1     = ztnw+ztsw;
        b1     = ztnw*ztsw;
        a2     = zbnw+zbsw;
        b2     = zbnw*zbsw;
        terms1 = a1*a1-b1-(a2*a2-b2);
        terms2 = 3.*(a1-a2);
        west   = .5*(rlns-rlnn)*(-terms1*(pbw-ptw)+terms2*(ztw*pbw-zbw*ptw))/(ztw-zbw);

        /* Top */
        a1     = ztne+ztnw;
        a2     = ztse+ztsw;
        terms2 = 3.*(a1-a2);
        if (fabs(2.*(ztn-zts)/(ztn+zts-zbn-zbs)) > DLOGZ_TOL) {
          b1     = ztne*ztnw;
          b2     = ztse*ztsw;
          b3     = ztne*ztse+ztnw*ztsw;
          terms1 = a1*a1-b1-(a2*a2-b2);
          top    = .5*((ptn-pts)*((rlnn+rlns)*terms1+(rlnn-rlns)*(b1+b2-b3))+(rlnn+rlns)*terms2*(ztn*pts-zts*ptn))/(ztn-zts);
        }
        else {
          top = .25*(rlnn+rlns)*(pts+ptn)*terms2;
        }

        /* Bottom */
        a1     = zbne+zbnw;
        a2     = zbse+zbsw;
        terms2 = 3.*(a1-a2);
        if (fabs(2.*(zbn-zbs)/(ztn+zts-zbn-zbs)) > DLOGZ_TOL) {
          b1     = zbne*zbnw;
          b2     = zbse*zbsw;
          b3     = zbne*zbse+zbnw*zbsw;
          terms1 = a1*a1-b1-(a2*a2-b2);
          bottom = .5*((pbn-pbs)*((rlnn+rlns)*terms1+(rlnn-rlns)*(b1+b2-b3))+(rlnn+rlns)*terms2*(zbn*pbs-zbs*pbn))/(zbn-zbs);
        }
        else {
          bottom = .25*(rlnn+rlns)*(pbs+pbn)*terms2;
        }

        /* 
         *  Note: West and east have the same sign here; in the northern hemisphere they could
         *        legitimately be called "northwest" and "northeast," and likewise in the southern hemisphere.
         */
        pgflat  = top-bottom+south-(north+west+east);

        dzn     = ztne-zbne+ztnw-zbnw;
        dzs     = ztse-zbse+ztsw-zbsw;
        pgflat /= (ltn-lts)*rho*rlt*(rlns*(dzs+.5*dzn)+rlnn*(dzn+.5*dzs));

        DVDT(grid.it_uv_tend,K,J,I) += pgflat;
      }
    }
  }

  return;
}

/*======================= end of uv_pgrad_finite_volume() =======================*/

/*======================= uv_drag() =============================================*/

/*
 * Calculate Rayleigh drag tendencies.
 *
 * For the Adams-Bashforth timestep, use the current state.
 * For the leapfrog timestep, use the previous state (IT_MINUS1), which is
 * called lagging (using the current state is numerically unstable).
 *
 *  Aaron Herrnstein. 03-29-06
 *  Added option for applying drag to either 
 *     1.) each individual U(K,J,I) or 
 *     2.) the zonal average of U(K,J,:).
 */

#define HELD_SUAREZ_PBL_NU0 (1./(24.*60.*60.))

/*
 * Set LATERAL_SPONGE_DEPTH to -1 for no effect
 * (setting it to 0 can produce an incorrect warning about dividing by zero).
 * Otherwise, it is the number of J levels on each of the northern and
 * southern boundaries associated with lateral sponges. The value 5
 * provides a reasonable effect.
 */
#define LATERAL_SPONGE_DEPTH -1

void uv_drag(planetspec  *planet,
             EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    uavg,vavg,
    uray,vray,
    p3,pbot,amp,
    nu0,
    tmp;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *j_t_sponge_inv;
  wind_variable
    *u,*v;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_drag";

  if (!initialized) {
    if (LATERAL_SPONGE_DEPTH > 0) {
      /* 
       * Minimum relaxation time (strongest effect) in units of timesteps.
       * NOTE: Use >= 5. to avoid numerical instability.
       */
      EPIC_FLOAT
        min_relax_time = 30.,
        t0_inv,
        tmp;

      /*
       * Allocate memory.
       */
      j_t_sponge_inv = fvector(grid.jlo,grid.nj,dbmsname);

      /*
       * Assign values to lateral sponge dissipation rate.
       */
      t0_inv = 1./(min_relax_time*(EPIC_FLOAT)grid.dt);
      for (J = 0; J < LATERAL_SPONGE_DEPTH; J++) {
        tmp                        = (EPIC_FLOAT)(LATERAL_SPONGE_DEPTH-J)/LATERAL_SPONGE_DEPTH;
        j_t_sponge_inv[grid.jlo+J] = t0_inv*.5*(1.-cos(M_PI*tmp));
        j_t_sponge_inv[grid.nj -J] = t0_inv*.5*(1.-cos(M_PI*tmp));
      }
    }

    initialized = TRUE;
  }

  /*
   *  Add Rayleigh friction term for U and V to damp gravity-wave 
   *  reflections at the model's top.  
   *  This is the "sponge."
   */
  for (K = KLO-1; K <= grid.k_sponge; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= (grid.t_sponge_inv[K])*
                                       (U(grid.it_uv_dis,K,J,I)-U_SPINUP(K,J,I));
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= (grid.t_sponge_inv[K])*
                                       (V(grid.it_uv_dis,K,J,I)-0.);
      }
    }
  }
  /*
   * Add lateral sponges if requested. 
   */
  if (LATERAL_SPONGE_DEPTH > 0) {
    for (K = grid.k_sponge+1; K < KHI; K++) {
      for (I = ILO; I <= IHI; I++) {
        for (J = JLO; J <= JHI; J++) {
          if (J < grid.jlo+LATERAL_SPONGE_DEPTH ||
              J > grid.nj -LATERAL_SPONGE_DEPTH   ) {
            DUDT(grid.it_uv_tend,K,J,I) -= (j_t_sponge_inv[J])*
                                           (U(grid.it_uv_dis,K,J,I)-U_SPINUP(K,J,I));
          }
        }
        for (J = JFIRST; J <= JHI; J++) {
          if (J < grid.jlo +LATERAL_SPONGE_DEPTH ||
              J > grid.nj+1-LATERAL_SPONGE_DEPTH   ) {
            DVDT(grid.it_uv_tend,K,J,I) -= (j_t_sponge_inv[J])*
                                           (V(grid.it_uv_dis,K,J,I)-0.);
          }
        }
      }
    }
  }

  if ((grid.newt_cool_on && grid.prandtl > 0.      ) ||
      (grid.tau_drag > 0. && grid.tau_drag < 1.e+20)   ) {
    /* 
     * Add Rayleigh drag.
     */
    for (K = grid.k_sponge+1; K < KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        /* compute zonal average of u */
        uavg = 0.;
        for (I = ILO; I <= IHI; I++) {
          /* NOTE: Only valid if I direction isn't cut */
          uavg += WIND(&var.u,grid.it_uv_dis,K,J,I);
        }
        uavg /= (grid.ni);
        if (grid.newt_cool_on == TRUE) {
          for (I = ILO; I <= IHI; I++) {
            /*
             * Average t_rad() onto U grid.
             */
            uray                         = (  grid.drag_zonal_avg)*uavg 
                                         + (1-grid.drag_zonal_avg)*WIND(&var.u,grid.it_uv_dis,K,J,I);
            tmp                          = .5*(t_rad(planet,K,J,I)+t_rad(planet,K,J,I-1));
            nu0                          = grid.prandtl/tmp;
            DUDT(grid.it_uv_tend,K,J,I) -= nu0*(uray-U_SPINUP(K,J,I));
          }
        }
        else {
          for (I = ILO; I <= IHI; I++) {
            uray                         = (  grid.drag_zonal_avg)*uavg 
                                         + (1-grid.drag_zonal_avg)*WIND(&var.u,grid.it_uv_dis,K,J,I);
            nu0                          = 1./grid.tau_drag;
            DUDT(grid.it_uv_tend,K,J,I) -= nu0*(uray-U_SPINUP(K,J,I));
          }
        }
      }
      if (grid.drag_v) {
        for (J = JFIRST; J <= JHI; J++) {
          /* compute zonal average of v */
          vavg = 0.;
          for (I = ILO; I <= IHI; I++) {
            vavg += WIND(&var.v,grid.it_uv_dis,K,J,I);
          }
          vavg /= grid.ni;
          if (grid.newt_cool_on == TRUE) {
            for (I = ILO; I <= IHI; I++) {
              /*
               * Average t_rad() onto V grid.
               */
              vray                         = (  grid.drag_zonal_avg)*vavg 
                                           + (1-grid.drag_zonal_avg)*WIND(&var.v,grid.it_uv_dis,K,J,I);
              tmp                          = .5*(t_rad(planet,K,J,I)+t_rad(planet,K,J-1,I));
              nu0                          = grid.prandtl/tmp;
              DVDT(grid.it_uv_tend,K,J,I) -= nu0*vray;
            }
          }
          else {
            for (I = ILO; I <= IHI; I++) {
              vray                         = (  grid.drag_zonal_avg)*vavg 
                                           + (1-grid.drag_zonal_avg)*WIND(&var.v,grid.it_uv_dis,K,J,I);
              nu0                          = 1./grid.tau_drag;
              DVDT(grid.it_uv_tend,K,J,I) -= nu0*vray;
            }
          }
        }
      }
    } /* end of loop over K */
  }

  /*
   * For the Held-Suarez test case, add planetary boundary layer (PBL) Rayleigh drag.
   * Otherwise, the PBL is handled by epic_subgrid.c subroutines.
   */
  if (strcmp(planet->name,"held_suarez") == 0) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pbot  = 0.5*(P3(grid.nk,J,I)+P3(grid.nk,J,I-1));
        for (K = grid.nk-1; K >= 1; K--) {
          p3  = .5*(P3(K,J,I)+P3(K,J,I-1));
           amp = (p3/pbot-0.7)/(1.-0.7);
          if (amp < 0.) {
            break;
          }
          else {
            nu0                          = HELD_SUAREZ_PBL_NU0*amp;
            DUDT(grid.it_uv_tend,K,J,I) -= nu0*WIND(&var.u,grid.it_uv_dis,K,J,I);
          }
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pbot  = 0.5*(P3(grid.nk,J,I)+P3(grid.nk,J-1,I));
        for (K = grid.nk-1; K >= 1; K--) {
          p3  = .5*(P3(K,J,I)+P3(K,J-1,I));
          amp = (p3/pbot-0.7)/(1.-0.7);
          if (amp < 0.) {
            break;
          }
          else {
            nu0                          = HELD_SUAREZ_PBL_NU0*amp;
            DVDT(grid.it_uv_tend,K,J,I) -= nu0*WIND(&var.v,grid.it_uv_dis,K,J,I);
          }
        }
      }
    }
  }

  return;
}

/*======================= end of uv_drag() ======================================*/

/* * * * * * * * * * * * end of epic_timestep.c * * * * * * * * * * * * * * * */
