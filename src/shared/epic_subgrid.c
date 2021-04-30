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

/* * * * * * * * * * epic_subgrid.c  * * * * * * * * * * * * * * * *
 *                                                                 *
 * Vimal Kumar Parimi, Raymond P. LeBeau, T. Dowling               *
 *                                                                 *
 * Subgrid-scale closure subroutines for the EPIC model.           *
 *                                                                 *
 * This file contains the following subroutines:                   *
 *                                                                 *
 *  ____Called from outside___    _______Internal____________      *
 *                                                                 *
 *  set_max_nu()                                                   *
 *  init_viscosity()                                               *
 *                                                                 *
 *  scalar_horizontal_subgrid()   scalar_horizontal_diffusion()    *
 *                                laplacian_h()                    *
 *                                                                 *
 *                                zonal_hyperviscosity()           *
 *                                meridional_hyperviscosity()      *
 *                                d2dy2()                          *
 *                                                                 *
 *  scalar_vertical_subgrid()     scalar_vertical_diffusion()      *
 *                                                                 *
 *  uv_horizontal_subgrid()       uv_horizontal_diffusion()        *
 *                                divergence_damping()             *
 *                                                                 *
 *  uv_vertical_subgrid()         uv_vertical_diffusion()          *
 *                                                                 *
 *  make_arrays_subgrid()                                          *
 *  free_arrays_subgrid()                                          *
 *  init_subgrid()                                                 *
 *  set_diffusion_coef()                                           *
 *  stability_factor()                                             *
 *                                                                 *
 *  source_sink_turb()            source_sink_SA()                 *
 *                                                                 *
 *                                dwall_SA()                       *
 *                                delta_SA()                       *
 *                                law_of_the_wall()                *
 *                                tau_surface()                    *
 *                                func_utau()                      *
 *                                invert_fv1()                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Variables with scope of epic_subgrid.c.
 */
EPIC_FLOAT
  *d_wall;

/*=============== set_max_nu() ====================================*/

/*
 * NOTE: max_nu[] can be large, so use double precision and
 * order the arithmetic wisely.
 */

void set_max_nu(double *max_nu)
{
  int
    K;
  register double
    dt,dy0,dy2,
    t2,dz,min_dz;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_max_nu";

  dt  = (double)(grid.dt);
  if (grid.nj == 0 && grid.ni == 1) {
    /*
     * One-dimesional, vertical case. Use dz instead of dy for the grid spacing.
     * The vertical diffusion is done implicitly, so the issue
     * is not numerical stability but accuracy.
     */
    min_dz = FLOAT_MAX;
    for (K = 1; K <= grid.nk; K++) {
      /*
       * Estimate layer temperature.
       */
      t2 = return_temp(planet,0.25,grid.p_ref[2*K],grid.theta_ref[2*K]);
      /*
       * Estimate dz and find minimum value.
       */
      dz     = log(grid.p_ref[2*K+1]/grid.p_ref[2*K-1])*planet->rgas*t2/grid.g[2*grid.jtp+1];
      min_dz = MIN(min_dz,dz);   
    }
    dy0 = min_dz;
  }
  else {
    dy0 = grid.dy0;         
  }
  dy2 = dy0*dy0;

  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    /* 
     * Ray LeBeau's estimates of safe max nu's (Dowling et al, 1998) for
     * the 3rd-Order Adams-Bashforth timestep.
     */
    max_nu[0] = (1./3.    )/dt;
    max_nu[2] = (dy2/30.  )/dt;
    max_nu[4] = (dy2/240. )*(dy2/dt);
    max_nu[6] = (dy2/800. )*(dy2/dt)*dy2;
    max_nu[8] = (dy2/2400.)*(dy2/dt)*dy2*dy2;
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    /*
     * In the leapfrog-timestep case, we lag the viscosity terms with
     * a long forward step.
     *
     * The factor 1/2 for max_nu[0] is mentioned in Durran (1991, eq(21)).
     * The factor 1/8 for max_nu[2] is derived in Tannehill, Anderson, and
     * Pletcher (1997, p. 137).
     *
     * We assume (1./8.)^(n/2) in general.
     */
    max_nu[0] = .5/dt;
    max_nu[2] = dy2/(8.*dt);
    max_nu[4] = max_nu[2]*dy2/8.;
    max_nu[6] = max_nu[4]*dy2/8.;
    max_nu[8] = max_nu[6]*dy2/8.;
  }
  else {
    sprintf(Message,"unrecognized uv_timestep_scheme: %s \n",
                    grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============== end of set_max_nu() ==============================*/

/*============== init_viscosity() =================================*/

/* 
 * Initialize hyperviscosity coefficients. 
 * See Dowling et al (1998) for details.
 */

void init_viscosity(planetspec *planet)
{
  register int
    ii;
  double  
    max_nu[MAX_NU_ORDER+1];  /* NOTE: declared as double, not EPIC_FLOAT */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_viscosity";

  set_max_nu(max_nu);

  grid.nu[2] = MIN(planet->kinvisc,max_nu[2]);

  sprintf(Message,"Divergence damping coeff, fraction of max\n");
  grid.nudiv_nondim = input_float(Message,grid.nudiv_nondim);

  for (ii = 4; ii <= MAX_NU_ORDER; ii+=2) {
    if (max_nu[ii] > 0.) {
      sprintf(Message,"nu[%d], fraction of max\n",ii);
      grid.nu_nondim[ii] = input_float(Message,grid.nu_nondim[ii]);
    }
    grid.nu[ii] = grid.nu_nondim[ii]*max_nu[ii];
  }

  return;
}

/*============== end of init_viscosity() ==========================*/

/*============== scalar_horizontal_subgrid() ======================*/
/*
 * Wrapper function.
 */
void scalar_horizontal_subgrid(planetspec  *planet,
                               EPIC_FLOAT **Buff2D)
{
  register int
    iq;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_horizontal_subgrid";

  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    if (grid.diffusion_direction == HORIZONTAL_AND_VERTICAL ||
        grid.diffusion_direction == JUST_HORIZONTAL) {
      scalar_horizontal_diffusion(planet,Buff2D);
    }
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * Apply hyperviscosity in meridional and zonal directions.
   */
  if (var.hdry.on) {
    meridional_hyperviscosity(planet,HDRY_INDEX,var.hdry.value,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,HDRY_INDEX,var.hdry.value,THREEDIM,Buff2D);
    restore_mass(planet,HDRY_INDEX,NO_PHASE);
  }

  if (var.theta.on) {
    meridional_hyperviscosity(planet,THETA_INDEX,var.theta.value,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,THETA_INDEX,var.theta.value,THREEDIM,Buff2D);
  }

  if (var.nu_turb.on) {
    meridional_hyperviscosity(planet,NU_TURB_INDEX,var.nu_turb.value,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,NU_TURB_INDEX,var.nu_turb.value,THREEDIM,Buff2D);
    restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
  }

  if (var.fpara.on) {
    meridional_hyperviscosity(planet,FPARA_INDEX,var.fpara.value,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,FPARA_INDEX,var.fpara.value,THREEDIM,Buff2D); 
  }

  /*
   * The iq loop is set up to cover only the species/phase fields that have been turned on.
   */ 
  for (iq = 0; iq < grid.nq; iq++) {
    meridional_hyperviscosity(planet,grid.is[iq],var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM,Buff2D);
    zonal_hyperviscosity(     planet,grid.is[iq],var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM,Buff2D);
    restore_mass(planet,grid.is[iq],grid.ip[iq]);
  }

  return;
}

/*============== end of scalar_horizontal_subgrid() ===============*/

/*============== scalar_horizontal_diffusion ======================*/

void scalar_horizontal_diffusion(planetspec  *planet,
                                 EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    iq,is,ip;
  register EPIC_FLOAT
    dt;
  const EPIC_FLOAT
    sigma_inv = 3./2.;
  EPIC_FLOAT
    *hh,
    *diff_coef,
    *laph;

  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_horizontal_diffusion";

  hh        = Buff2D[0];
  laph      = Buff2D[1];
  diff_coef = Buff2D[2];


  dt = (EPIC_FLOAT)grid.dt;

  /*
   * Apply diffusion to THETA.
   *
   * NOTE: We have not included molecular diffusion for THETA.
   *       If it is added, may have to deal with
   *       density weighting and units of the transport coefficient.
   */
  for (K = grid.k_sigma-1; K <= KHI; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HH(J,I)        = THETA(K,J,I);
        DIFF_COEF(J,I) = DIFFUSION_COEF_THETA(K,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
    
    laplacian_h(planet,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        THETA(K,J,I) += dt*LAPH(J,I);
      }
    }
  }

  /*
   * NOTE: Currently not applying diffusion to HDRY.
   *       If it is added, recall that DIFFUSION_COEF_MASS is carried on the layer interfaces.
   */

  /*
   * Apply diffusion to mixing ratios, Q, which are carried on the layer interfaces.
   */

  for (iq = 0; iq < grid.nq; iq++) {
    is = grid.is[iq];
    ip = grid.ip[iq];
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        if (grid.ip[iq] == VAPOR) {
          /*
           * Include molecular diffusion for vapor phase.
           */
          for (I = ILOPAD; I <= IHIPAD; I++) {
            HH(J,I)        = Q(is,ip,K,J,I);
            DIFF_COEF(J,I) = mass_diffusivity(planet,is,T3(K,J,I),P3(K,J,I))
                            +DIFFUSION_COEF_MASS(K,J,I);
          }
        }
        else {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            HH(J,I)        = Q(is,ip,K,J,I);
            DIFF_COEF(J,I) = DIFFUSION_COEF_MASS(K,J,I);
          }
        }
      }
      /* No need to apply bc_lateral() here. */
      
      laplacian_h(planet,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          Q(is,ip,K,J,I) += dt*LAPH(J,I);
        }
      }
      /* No need to apply bc_lateral() here. */
    }
    /*
     * Clean up any negative mass introduced by diffusion truncation error.
     */
    restore_mass(planet,is,ip);
  }

  /*
   * Update pressures, etc.
   */
  set_p2etc(planet,UPDATE_THETA);

  /*
   * Apply diffusion to NU_TURB.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HH(J,I)        = NU_TURB(K,J,I);
        DIFF_COEF(J,I) = sigma_inv*(planet->kinvisc+NU_TURB(K,J,I));
      }
    }
    /* No need to apply bc_lateral() here. */

    laplacian_h(planet,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        NU_TURB(K,J,I) += dt*LAPH(J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  restore_mass(planet,NU_TURB_INDEX,NO_PHASE);

  return;
}

/*============== end of scalar_horizontal_diffusion ===============*/

/*============== laplacian_h() ====================================*/

/*
 *  Calculates the 2D Laplacian of an h-grid scalar field,
 *  multiplied by the input diffusion coefficient, DIFF_COEF(J,I).
 *
 *  The pointer hh should point to the appropriate JI plane.
 *
 *  Pointers to memory for two working JI-plane buffers
 *  are passed in as buff1 and buff2.
 */

void laplacian_h(planetspec *planet,
                 EPIC_FLOAT *hh,
                 EPIC_FLOAT *diff_coef,
                 EPIC_FLOAT *lph,
                 EPIC_FLOAT *buff1,
                 EPIC_FLOAT *buff2)
{
  int 
    J,I;
  EPIC_FLOAT 
    h_edge,
    m_2j,n_2j,m_2jp1,n_2jp1,
    m_2j_inv,m_2jp2_inv;
  EPIC_FLOAT
    *gh1,*gh2;
  static EPIC_FLOAT
    *buffji;
  static int
    initialized = FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="laplacian_h";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    buffji = fvector(0,Nelem2d-1,dbmsname);

    initialized = TRUE;
  }

  /* 
   * Check that the pointer hh is not NULL.
   */
  if (!hh) {
    sprintf(Message,"input pointer hh=NULL");
    epic_error(dbmsname,Message);
  }

  /* Zero working buffers: */
  memset(buff1,0,Nelem2d*sizeof(EPIC_FLOAT));
  memset(buff2,0,Nelem2d*sizeof(EPIC_FLOAT));

  gh1 = buff1;
  gh2 = buff2;

  /* 
   * Compute grad h = (gh1,gh2).
   *
   * Apply zonal_filter() to a copy of HH before d/dx, to control numerical instability.
   */
  memcpy(buffji,hh,Nelem2d*sizeof(EPIC_FLOAT));
  zonal_filter(H2_INDEX,buffji,NULL,TWODIM);
  for (J = JLO; J <= JHI; J++) {
    m_2jp1 = (grid.m)[2*J+1];
    for (I = ILO; I <= IHI; I++) {
      GH1(J,I) = m_2jp1*(BUFFJI(J,I)-BUFFJI(J,I-1));    
    }
  }
  /* update gh1 edges below */
  for (J = JFIRST; J <= JHI; J++) {
    n_2j = (grid.n)[2*J];
    for (I = ILO; I <= IHI; I++) {
      GH2(J,I) = n_2j*(HH(J,I)-HH(J-1,I));
    }
  }

  /* 
   * Fill in gh2 for top and bottom channel boundaries. 
   * NOTE: Code only valid if i-direction is not cut. 
   */
  if (JLO == grid.jlo && !IS_SPOLE) {
    /* southern edge */
    h_edge = 0.;
    for (I = ILO; I <= IHI; I++) {
      h_edge += HH(JLO,I);
    }
    h_edge /= (IHI-ILO+1);
    n_2j = (grid.n)[2*JLO];
    for (I = ILOPAD; I <= IHIPAD; I++) {
      GH2(JLO,I) = n_2j*(HH(JLO,I)-h_edge);
    }
  }
  if (JHI == grid.nj && !IS_NPOLE) {
    /* northern edge */
    h_edge = 0.;
    for (I = ILO; I <= IHI; I++) {
      h_edge += HH(JHI,I);
    }
    h_edge /= (IHI-ILO+1);
    n_2j = (grid.n)[2*(JHI+1)];
    for (I = ILOPAD; I <= IHIPAD; I++) {
      GH2(JHI+1,I) = n_2j*(h_edge-HH(JHI,I));
    }
  }
  /* update gh2 edges below */

  /* Update edges for gh1, gh2: */
  bc_lateral(gh1,TWODIM);
  bc_lateral(gh2,TWODIM);

  /*
   * Multiply by DIFF_COEF.
   */
  /* GH1 is on the U-grid. */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      GH1(J,I) *= .5*(DIFF_COEF(J,I)+DIFF_COEF(J,I-1));
    }
  }
  bc_lateral(gh1,TWODIM);

  /* GH2 is on the V-grid. */
  for (J = JFIRST; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      GH2(J,I) *= .5*(DIFF_COEF(J,I)+DIFF_COEF(J-1,I));
    }
  }
  bc_lateral(gh2,TWODIM);

  /* 
   * Compute laplacian.
   *
   * Apply zonal_filter to GH1 to control numerical instability of d/dx.
   */
  zonal_filter(U_INDEX,gh1,NULL,TWODIM);
  for (J = JLO; J <= JHI; J++) {
    m_2jp1 = (grid.m)[2*J+1];
    n_2jp1 = (grid.n)[2*J+1];
    if (J == grid.jlo && IS_SPOLE) {
      m_2j_inv = 0.;
    }
    else {
      m_2j_inv = 1./(grid.m)[2*J];
    }
    if (J == grid.nj && IS_NPOLE) {
      m_2jp2_inv = 0.;
    }
    else {
      m_2jp2_inv = 1./(grid.m)[2*J+2];
    }
    for (I = ILO; I <= IHI; I++) {
      LPH(J,I) = m_2jp1*( (GH1(J,I+1)-GH1(J,I))
                  +n_2jp1*(GH2(J+1,I)*m_2jp2_inv-GH2(J,I)*m_2j_inv) );
    }
  }
  bc_lateral(lph,TWODIM);

  return;
}

/*============== end of laplacian_h() =============================*/

/*============== zonal_hyperviscosity() ===========================*/

/*
 * Apply hyperviscosity in zonal direction.
 * Use a 3-point stencil for each second derivative and
 * taper the hyperviscosity coefficient for |lat| > LAT0 to maintain
 * numerical stability.
 */

void zonal_hyperviscosity(planetspec  *planet,
                          int          is,
                          EPIC_FLOAT  *variable,
                          int          dim,
                          EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kstart,kend,
    jbot,
    itmp,itmp2,
    jhalfstep,
    sign;
  static int
    initialized = FALSE;
  register double
    tmp,rln,
    diff_coef;
  static double
    m0,
    coef[  MAX_NU_ORDER+1];
  double
    max_nu[MAX_NU_ORDER+1];
  EPIC_FLOAT
    *a,
    *hh,
    *d2dx2var,
    *ptmp,
     false_mean;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="zonal_hyperviscosity";

  if (!initialized) {
    set_max_nu(max_nu);

    /*
     * grid.nu[2] should be zero, since regular viscosity
     * is handled elsewhere.
     */
    grid.nu[2] = 0.;

    for (itmp2 = 2; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
      /*
       * Factor the viscosity coefficients to
       * help prevent floating point overflow/underflow.
       */
      if (grid.nu[itmp2] > 0.) {
        coef[itmp2] = grid.nu[itmp2];
      }
      else {
        coef[itmp2] = 0.01*max_nu[itmp2];
      }
      for (itmp = 2; itmp < itmp2; itmp+=2) {
        coef[itmp2] /= coef[itmp];
      }
    }

    /*
     * Set m0=1/dx0 corresponding to LAT0, the latitude above which measures are
     * taken to deal with the vanishing size of dx on the lon-lat grid.
     */
    if (strcmp(grid.geometry,"globe") == 0) {
      rln  = planet->re/sqrt(1.+ pow(planet->rp/planet->re*tan(LAT0*DEG),2.) );
      m0   = 1./(rln*grid.dln*DEG);
    }
    else if (strcmp(grid.geometry,"f-plane")  == 0) { 
      if (strcmp(grid.f_plane_map,"polar") == 0) {
        rln = .5*grid.f_plane_half_width;
        m0  = 1./(rln*grid.dln*DEG);
      }
      else {
        m0 = 1./(grid.f_plane_half_width/grid.ni);
      }
    }
    else {
      sprintf(Message,"unrecognized grid.geometry=%s");
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  /*
   * No need to apply zonal hyperviscosity if ni < 2;
   */
  if (grid.ni < 2) return;

  tmp = 0.;
  for (itmp2 = 4; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
    tmp += grid.nu[itmp2];
  }
  if (tmp == 0.) {
    /* 
     * Return if all hyperviscosity coefficients are zero: 
     */
    return;
  }

  if (dim == TWODIM) {
    kstart = Kshift;
    kend   = Kshift;  
  }
  else if (dim == THREEDIM) {
    kstart = KLO;
    kend   = KHI;
  }
  else {
    sprintf(Message,"dim=%d not recognized");
    epic_error(dbmsname,Message);
  }

  switch(is) {
    case U_INDEX:
      jbot      = JLO;
      jhalfstep = 1;
      if (dim == THREEDIM) {
        kend = KHI-1;
      }
    break;
    case V_INDEX:
      jbot      = JFIRST;
      jhalfstep = 0;
      if (dim == THREEDIM) {
        kend = KHI-1;
      }
    break;
    case PV3_INDEX:
      jbot      = JFIRST;
      jhalfstep = 0;
    break;
    case THETA_INDEX:
      jbot      = JLO;
      jhalfstep = 1;
      if (dim == THREEDIM) {
        kstart = grid.k_sigma-1;
      }
    break;
    default:
      jbot      = JLO;
      jhalfstep = 1;
    break;
  }

  for (K = kstart; K <= kend; K++) {
    a = variable+(K-Kshift)*Nelem2d;

    hh       = Buff2D[0];
    d2dx2var = Buff2D[1];

    for (J = JLOPAD; J <= JHIPADPV; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HH(J,I) = A(J,I);
      }
    }

    sign = 1;
    itmp = 2;
    tmp  = 1.;
    while (tmp > 0.) {
      /*
       * Calculate second derivative.
       */
      for (J = jbot; J <= JHI; J++) {
        if (fabs(grid.lat[2*J+jhalfstep]) > LAT0) {
          /*
           * Taper hyperviscosity coefficient at high latitudes to maintain numerical stability.
           */
          diff_coef = coef[itmp]*m0*m0;
        }
        else {
          diff_coef = coef[itmp]*grid.m[2*J+jhalfstep]*grid.m[2*J+jhalfstep];
        }

        for (I = ILO; I <= IHI; I++) {
          D2DX2VAR(J,I) = diff_coef*(HH(J,I-1)-2.*HH(J,I)+HH(J,I+1));
        }
        /* Apply periodicity. */
        D2DX2VAR(J,grid.ni+1) = D2DX2VAR(J,1      );
        D2DX2VAR(J,0        ) = D2DX2VAR(J,grid.ni);

        /*
         * Ensure that D2DX2VAR has zero zonal-mean value, consistent with the 
         * periodic zonal boundary conditions.
         */
        false_mean = 0.;
        for (I = ILO; I <= IHI; I++) {
          false_mean += D2DX2VAR(J,I);
        }
        false_mean /= grid.ni;
        for (I = ILOPAD; I <= IHIPAD; I++) {
          D2DX2VAR(J,I) -= false_mean;
        }

        if (grid.nu[itmp] > 0.) {
          /*
           * Apply hyperviscosity.
           */
          for (I = ILO; I <= IHI; I++) {
            A(J,I) += DT*(EPIC_FLOAT)sign*D2DX2VAR(J,I);
          }
        }
      }
      if (grid.nu[itmp] > 0.) {
        /* Need to apply bc_lateral() here. */
        bc_lateral(a,TWODIM);
      }

      sign     *= -1;
      itmp     +=  2;
      ptmp      = hh;
      hh        = d2dx2var;
      d2dx2var  = ptmp;
      /* 
       * Set up test for whether all remaining hyperviscosity coefficients are zero.
       */
      tmp = 0.;
      for (itmp2 = itmp; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
        tmp += grid.nu[itmp2];
      }
    }
  } /* K loop */

  return;
}

/*============== end of zonal_hyperviscosity() ====================*/

/*============== meridional_hyperviscosity() ======================*/

/*
 * Apply explicit hyperviscosity in meridional direction.
 * The derivatives are done by calls to d2dy2().
 */

void meridional_hyperviscosity(planetspec  *planet,
                               int          is,
                               EPIC_FLOAT  *variable,
                               int          dim,
                               EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kstart,kend,
    jbot,
    itmp,itmp2,
    sign;
  static int
    initialized = FALSE;
  register double
    tmp,
    diff_coef;
  static double
    coef[  MAX_NU_ORDER+1];
  double
    max_nu[MAX_NU_ORDER+1];
  EPIC_FLOAT
    *a,
    *hh,
    *d2dy2var,
    *ptmp; 
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="meridional_hyperviscosity";

  if (!initialized) {
    set_max_nu(max_nu);

    /*
     * grid.nu[2] should be zero, since regular viscosity
     * is handled elsewhere.
     */
    grid.nu[2] = 0.;

    for (itmp2 = 2; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
      /*
       * Factor the viscosity coefficients to
       * help prevent floating point overflow/underflow.
       */
      if (grid.nu[itmp2] > 0.) {
        coef[itmp2] = grid.nu[itmp2];
      }
      else {
        coef[itmp2] = 0.01*max_nu[itmp2];
      }
      for (itmp = 2; itmp < itmp2; itmp+=2) {
        coef[itmp2] /= coef[itmp];
      }
    }

    initialized = TRUE;
  }

  tmp = 0.;
  for (itmp2 = 4; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
    tmp += grid.nu[itmp2];
  }
  if (tmp == 0.) {
    /* 
     * Return if all hyperviscosity coefficients are zero: 
     */
    return;
  }

  if (dim == TWODIM) {
    kstart = Kshift;
    kend   = Kshift;  
  }
  else if (dim == THREEDIM) {
    kstart = KLO;
    kend   = KHI;
  }
  else {
    sprintf(Message,"dim=%d not recognized");
    epic_error(dbmsname,Message);
  }

  switch(is) {
    case U_INDEX:
      jbot = JLO;
      if (dim == THREEDIM) {
        kend = KHI-1;
      }
    break;
    case V_INDEX:
      jbot = JFIRST;
      if (dim == THREEDIM) {
        kend = KHI-1;
      }
    break;
    case PV3_INDEX:
      jbot = JFIRST;
    break;
    case THETA_INDEX:
      jbot = JLO;
      if (dim == THREEDIM) {
        kstart = grid.k_sigma-1;
      }
    break;
    default:
      jbot = JLO;
    break;
  }

  for (K = kstart; K <= kend; K++) {
    a = variable+(K-Kshift)*Nelem2d;

    hh       = Buff2D[0];
    d2dy2var = Buff2D[1];

    for (J = JLOPAD; J <= JHIPADPV; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        HH(J,I) = A(J,I);
      }
    }

    sign = 1;
    itmp = 2;
    tmp  = 1.;
    while (tmp > 0.) {
      d2dy2(planet,is,hh,coef[itmp],d2dy2var);

      if (grid.nu[itmp] > 0.) {
        /* 
         * Apply hyperviscosity. Use a forward (Euler) step.
         * The hyperviscosity coefficient is handled by the coef[]
         * factoring to help prevent floating point overflow/underflow.
         *
         * Apply zonal_filter() to D2DY2VAR before applying hyperviscosity
         * to avoid adding zonal noise to the polar regions.
         */
        zonal_filter(is,d2dy2var,NULL,TWODIM);

        if (is == V_INDEX || is == PV3_INDEX) {
          for (J = JFIRST; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              A(J,I) += DT*(EPIC_FLOAT)sign*D2DY2VAR(J,I);
            }
          }
        }
        else {
          for (J = JLO; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              A(J,I) += DT*(EPIC_FLOAT)sign*D2DY2VAR(J,I);
            }
          }
        }

        /* Need to apply bc_lateral() here. */
        bc_lateral(a,TWODIM);
      }

      sign     *= -1;
      itmp     +=  2;
      ptmp      = hh;
      hh        = d2dy2var;
      d2dy2var  = ptmp;
      /* 
       * Set up test for whether all remaining viscosity coefficients are zero.
       */
      tmp = 0.;
      for (itmp2 = itmp; itmp2 <= MAX_NU_ORDER; itmp2+=2) {
        tmp += grid.nu[itmp2];
      }
    }
  } /* K loop */

  return;
}

/*============== end of meridional_hyperviscosity =================*/

/*============== d2dy2() ==========================================*/

/*
 * Calculates the meridional component of the Laplacian, 
 * multiplied by the constant input diffusion coefficient diff_coef.
 *
 * We create the differencing scheme by fitting the first
 * 3 Legendre polynomials to the 3-pt stencil and taking
 * the Laplacian on the result to generate the coefficients.
 *
 *  The pointer hh should point to the appropriate JI plane.
 */

void d2dy2(planetspec *planet,
           int         index,
           EPIC_FLOAT *hh,
           EPIC_FLOAT  diff_coef,
           EPIC_FLOAT *d2dy2var)
{
  int 
    J,I;
  static int
    initialized = FALSE;
  static float_triplet
    *h_coef,
    *v_coef;
  register double
    r2_inv,
    sum,
    lat[3],
    a[3];
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="d2dy2";

  if (!initialized) {
    h_coef = ftriplet(0,grid.nj,  dbmsname);
    v_coef = ftriplet(0,grid.nj+1,dbmsname);
    /*
     * Calculate coefficients for the H-grid/U-grid latitude positions,
     * based on applying the first 3 Legendre polynomials to the 3-pt stencil.
     * The southern and northern boundaries are handled with
     * an off-center scheme, and otherwise we use a centered scheme.
     * The math can be seen in the Mathematica notebook epic/tools/mathematica/merid_lap.nb.
     */
    J = grid.jlo;
    r2_inv      = (double)grid.rlt[2*J+1];
    r2_inv      = 1./(r2_inv*r2_inv);
    lat[0]      = (double)grid.lat[2*J+1]*DEG;
    lat[1]      = (double)grid.lat[2*J+3]*DEG;
    lat[2]      = (double)grid.lat[2*J+5]*DEG;
    h_coef[J].x =  r2_inv*(1.+sin(lat[0])*(-3.*sin(lat[0])+sin(lat[1])+sin(lat[2])))/
                   ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
    h_coef[J].y = -r2_inv*2.*(cos(2.*lat[0])+sin(lat[0])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    h_coef[J].z =  r2_inv*(cos(2.*lat[0])+sin(lat[0])*sin(lat[1]))/
                   (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));

    for (J = grid.jlo+1; J < grid.nj; J++) {
      r2_inv      = (double)grid.rlt[2*J+1];
      r2_inv      = 1./(r2_inv*r2_inv);
      lat[0]      = (double)grid.lat[2*J-1]*DEG;
      lat[1]      = (double)grid.lat[2*J+1]*DEG;
      lat[2]      = (double)grid.lat[2*J+3]*DEG;
      h_coef[J].x =  r2_inv*(cos(2.*lat[1])+sin(lat[1])*sin(lat[2]))/
                     ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
      h_coef[J].y = -r2_inv*2.*(1.+sin(lat[1])*(sin(lat[0])-3.*sin(lat[1])+sin(lat[2])))/
                     ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
      h_coef[J].z =  r2_inv*(cos(2.*lat[1])+sin(lat[0])*sin(lat[1]))/
                     (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    }
    J = grid.nj;
    r2_inv      = (double)grid.rlt[2*J+1];
    r2_inv      = 1./(r2_inv*r2_inv);
    lat[0]      = (double)grid.lat[2*J-3]*DEG;
    lat[1]      = (double)grid.lat[2*J-1]*DEG;
    lat[2]      = (double)grid.lat[2*J+1]*DEG;
    h_coef[J].x =  r2_inv*(cos(2.*lat[2])+sin(lat[1])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
    h_coef[J].y = -r2_inv*2.*(cos(2.*lat[2])+sin(lat[0])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    h_coef[J].z =  r2_inv*(1.+sin(lat[2])*(sin(lat[0])+sin(lat[1])-3.*sin(lat[2])))/
                   (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));

    /*
     * Calculate coefficients for the V-grid/PV-grid latitude positions.
     */
    J = grid.jlo;
    r2_inv      = (double)grid.rlt[2*J];
    r2_inv      = 1./(r2_inv*r2_inv);
    lat[0]      = (double)grid.lat[2*J  ]*DEG;
    lat[1]      = (double)grid.lat[2*J+2]*DEG;
    lat[2]      = (double)grid.lat[2*J+4]*DEG;
    v_coef[J].x =  r2_inv*(1.+sin(lat[0])*(-3.*sin(lat[0])+sin(lat[1])+sin(lat[2])))/
                   ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
    v_coef[J].y = -r2_inv*2.*(cos(2.*lat[0])+sin(lat[0])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    v_coef[J].z =  r2_inv*(cos(2.*lat[0])+sin(lat[0])*sin(lat[1]))/
                   (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));

    for (J = grid.jlo+1; J <= grid.nj; J++) {
      r2_inv      = (double)grid.rlt[2*J];
      r2_inv      = 1./(r2_inv*r2_inv);
      lat[0]      = (double)grid.lat[2*J-2]*DEG;
      lat[1]      = (double)grid.lat[2*J  ]*DEG;
      lat[2]      = (double)grid.lat[2*J+2]*DEG;
      v_coef[J].x =  r2_inv*(cos(2.*lat[1])+sin(lat[1])*sin(lat[2]))/
                     ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
      v_coef[J].y = -r2_inv*2.*(1.+sin(lat[1])*(sin(lat[0])-3.*sin(lat[1])+sin(lat[2])))/
                     ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
      v_coef[J].z =  r2_inv*(cos(2.*lat[1])+sin(lat[0])*sin(lat[1]))/
                     (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    }
    J = grid.nj+1;
    r2_inv      = (double)grid.rlt[2*J];
    r2_inv      = 1./(r2_inv*r2_inv);
    lat[0]      = (double)grid.lat[2*J-4]*DEG;
    lat[1]      = (double)grid.lat[2*J-2]*DEG;
    lat[2]      = (double)grid.lat[2*J  ]*DEG;
    v_coef[J].x =  r2_inv*(cos(2.*lat[2])+sin(lat[1])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*.5*(sin(lat[2])-sin(lat[0])));
    v_coef[J].y = -r2_inv*2.*(cos(2.*lat[2])+sin(lat[0])*sin(lat[2]))/
                   ((sin(lat[1])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));
    v_coef[J].z =  r2_inv*(1.+sin(lat[2])*(sin(lat[0])+sin(lat[1])-3.*sin(lat[2])))/
                   (.5*(sin(lat[2])-sin(lat[0]))*(sin(lat[2])-sin(lat[1])));

    initialized = TRUE;
  }
  /* end of initialization */

  /* 
   * Check that the pointer hh is not NULL.
   */
  if (!hh) {
    sprintf(Message,"input pointer hh=NULL");
    epic_error(dbmsname,Message);
  }

  if (index == V_INDEX || index == PV3_INDEX) {
    if (JLO == 0) {
      /*
       * South pole uses an uncentered scheme.
       */
      J = 0;
      for (I = ILO; I <= IHI; I++) {
        /*
         * Multiply the coefficients first to avoid underflows.
         */
        a[0] = diff_coef*v_coef[J].x;
        a[1] = diff_coef*v_coef[J].y;
        a[2] = diff_coef*v_coef[J].z;
        D2DY2VAR(J,I) = a[0]*HH(J,I)+a[1]*HH(J+1,I)+a[2]*HH(J+2,I);
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      /*
       * Interior points use a centered scheme.
       */
      for (I = ILO; I <= IHI; I++) {
        a[0] = diff_coef*v_coef[J].x;
        a[1] = diff_coef*v_coef[J].y;
        a[2] = diff_coef*v_coef[J].z;
        D2DY2VAR(J,I) = a[0]*HH(J-1,I)+a[1]*HH(J,I)+a[2]*HH(J+1,I);
      }
    }
    if (JHI == grid.nj) {
      /*
       * North pole uses an uncentered scheme.
       */
      J = grid.nj+1;
      for (I = ILO; I <= IHI; I++) {
        a[0] = diff_coef*v_coef[J].x;
        a[1] = diff_coef*v_coef[J].y;
        a[2] = diff_coef*v_coef[J].z;
        D2DY2VAR(J,I) = a[0]*HH(J-2,I)+a[1]*HH(J-1,I)+a[2]*HH(J,I);
      }
    }
  }
  else {
    for (J = JLO; J <= JHI; J++) {
      if (J == grid.jlo) {
        /*
         * Southern edge uses an uncentered scheme.
         */
        for (I = ILO; I <= IHI; I++) {
          a[0] = diff_coef*h_coef[J].x;
          a[1] = diff_coef*h_coef[J].y;
          a[2] = diff_coef*h_coef[J].z;
          D2DY2VAR(J,I) = a[0]*HH(J,I)+a[1]*HH(J+1,I)+a[2]*HH(J+2,I);
        }
      }
      else if (J == grid.nj) {
        /*
         * Northern edge uses an uncentered scheme.
         */
        for (I = ILO; I <= IHI; I++) {
          a[0] = diff_coef*h_coef[J].x;
          a[1] = diff_coef*h_coef[J].y;
          a[2] = diff_coef*h_coef[J].z;
          D2DY2VAR(J,I) = a[0]*HH(J-2,I)+a[1]*HH(J-1,I)+a[2]*HH(J,I);
        }
      }
      else {
        /*
         * Interior points use a centered scheme.
         */
        for (I = ILO; I <= IHI; I++) {
          a[0] = diff_coef*h_coef[J].x;
          a[1] = diff_coef*h_coef[J].y;
          a[2] = diff_coef*h_coef[J].z;
          D2DY2VAR(J,I) = a[0]*HH(J-1,I)+a[1]*HH(J,I)+a[2]*HH(J+1,I);
        }
      }
    }
  }
  bc_lateral(d2dy2var,TWODIM);

  return;
}

/*============== end of d2dy2() ===================================*/

/*============== scalar_vertical_subgrid() ========================*/

/*
 * Wrapper function.
 */
void scalar_vertical_subgrid(planetspec  *planet,
                             EPIC_FLOAT **Buff2D)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_vertical_subgrid";

  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    if (grid.diffusion_direction == HORIZONTAL_AND_VERTICAL ||
        grid.diffusion_direction == JUST_VERTICAL) {
      scalar_vertical_diffusion(planet,Buff2D);
    }
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============== end of scalar_vertical_subgrid() =================*/

/*============== scalar_vertical_diffusion() ======================*/

void scalar_vertical_diffusion(planetspec  *planet,
                               EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kay,
    iq;
  static int
    nnk,
    initialized = FALSE;
  EPIC_FLOAT
    diffusion_coeff,
    *tau_wall;
  const EPIC_FLOAT
    sigma_inv = 3./2.;
  static EPIC_FLOAT
    *stab_factor,
    *zee,
    *aaa,
    *dee,
    *ans;
  unsigned long
    nbytes_2d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_vertical_diffusion";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  memset(Buff2D[0],0,nbytes_2d);
  tau_wall = Buff2D[0];

  if (!initialized) {
    nnk = KHI-KLO+1;

    /* Allocate memory: */
    stab_factor = fvector(0,2*KHI+1,dbmsname);
    zee         = fvector(0,KHI+2,  dbmsname);
    aaa         = fvector(0,KHI+2,  dbmsname);
    dee         = fvector(0,KHI+2,  dbmsname);
    ans         = fvector(0,KHI+2,  dbmsname);

    initialized = TRUE;
  }

  /*
   * Apply vertical diffusion to THETA, 
   * which is carried on the layer interfaces.
   *
   * NOTE: Have not yet included molecular diffusion for THETA.
   *
   * Do not apply diffusion to THETA outside of pure-sigma region of model,
   * since there it is a diagnostic variable. 
   */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      stability_factor(THETA_INDEX,J,I,stab_factor);
      /*
       * Use no-flux boundary condition at bottom.
       */
      aaa[0] = 1.e+20;
      zee[0] = GZ3(KHI,J,I)/grid.g[2*J+1];
      for (K = KHI; K >= grid.k_sigma-2; K--) {
        /*
         * Load in positive-z direction (which unfortunately
         * fights against the top-down K numbering of layers).
         */
        kay      = KHI-K+1;
        aaa[kay] = THETA(K,J,I);
        zee[kay] = GZ3(K,J,I)/grid.g[2*J+1];
      }
      for (K = KHI; K > grid.k_sigma-2; K--) {
        kay      = KHI-K+1;
        dee[kay] = stab_factor[2*K]*.5*(DIFFUSION_COEF_THETA(K,J,I)+DIFFUSION_COEF_THETA(K-1,J,I));
      }

      crank_nicolson(KHI-grid.k_sigma+2,DT,zee,aaa,dee,NULL,ans);

      for (K = KHI; K >= grid.k_sigma-1; K--) {
        kay          = KHI-K+1;
        THETA(K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.theta.value,THREEDIM);

  /*
   * NOTE: Not applying diffusion to HDRY.
   */

  /*
   * Apply vertical diffusion to mixing ratios, Q,
   * which are carried on the layer interfaces.
   *
   * Loop over all activated species/phase variables.
   */
  for (iq = 0; iq < grid.nq; iq++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        stability_factor(grid.is[iq],J,I,stab_factor);

        /*
         * Use no-flux boundary conditions at top and bottom.
         */
        aaa[0    ] = 1.e+20;
        aaa[KHI+2] = 1.e+20;
        zee[0]     = (GZ3(KHI,J,I)-(GZ3(KHI-1,J,I)-GZ3(KHI,J,I)))/grid.g[2*J+1];
        zee[KHI+2] = (GZ3(0,  J,I)+(GZ3(0,    J,I)-GZ3(1,  J,I)))/grid.g[2*J+1];
        for (K = KHI; K >= KLO-1; K--) {
          kay      = KHI-K+1;
          /*
           * Apply diffusion to mixing ratio.
           */
          aaa[kay] = Q(grid.is[iq],grid.ip[iq],K,J,I);
          zee[kay] = GZ3(K,J,I)/grid.g[2*J+1];
        }

        /*
         * NOTE: The dee array is on a staggered grid.
         */
        for (K = KHI; K > KLO-1; K--) {
          kay      = KHI-K+1;
          dee[kay] = stab_factor[2*K]*.5*(DIFFUSION_COEF_MASS(K,J,I)+DIFFUSION_COEF_MASS(K-1,J,I));
        }
        if (grid.ip[iq] == VAPOR) {
          /*
           * Include molecular mass diffusivity for vapor phase.
           */
          for (K = KHI; K > KLO-1; K--) {
            kay       = KHI-K+1;
            dee[kay] += mass_diffusivity(planet,grid.is[iq],T2(K,J,I),P2(K,J,I));
          }
        }

        crank_nicolson(KHI+1,DT,zee,aaa,dee,NULL,ans);

        for (K = KHI; K >= KLO-1; K--) {
          kay                              = KHI-K+1;
          Q(grid.is[iq],grid.ip[iq],K,J,I) = ans[kay];
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);
    /*
     * Clean up any negative mass introduced by diffusion truncation errors.
     */
    restore_mass(planet,grid.is[iq],grid.ip[iq]);
  }

  /*
   * Update P2, etc.
   */
  set_p2etc(planet,UPDATE_THETA);

  /*
   * Apply vertical diffusion to NU_TURB, which is carried on the layer interfaces.
   *
   * NOTE: The function tau_surface() currently sets TAU_WALL to zero for gas-giant case.
   */
  tau_surface(planet,NU_TURB_INDEX,tau_wall,Buff2D[1]);

  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      /*
       * Use no-flux boundary condition at top.
       */
      aaa[KHI+2] = 1.e+20;
      zee[KHI+2] = (GZ3(0,J,I)+(GZ3(0,J,I)-GZ3(1,J,I)))/grid.g[2*J+1];
      for (K = KHI; K >= KLO-1; K--) {
        /*
         * Load in positive-z direction (which unfortunately
         * fights against the top-down K numbering of layers).
         */
        kay      = KHI-K+1;
        aaa[kay] = NU_TURB(K,J,I);
        zee[kay] = GZ3(K,J,I)/grid.g[2*J+1];
      }
      /*
       * NOTE: The diffusion coefficient is on a staggered grid.
       */
      for (K = KHI; K > KLO-1; K--) {
        kay      = KHI-K+1;
        dee[kay] = sigma_inv*(planet->kinvisc+.5*(NU_TURB(K,J,I)+NU_TURB(K-1,J,I)));
      }
      /*
       * Set bottom boundary condition consistent with TAU_WALL.
       */
      zee[0] = (GZ3(KHI,J,I)-(GZ3(KHI-1,J,I)-GZ3(KHI,J,I)))/grid.g[2*J+1];
      dee[0] = dee[1];
      aaa[0] = aaa[1]-TAU_WALL(J,I)*(zee[1]-zee[0])/dee[0];

      crank_nicolson(KHI+1,DT,zee,aaa,dee,NULL,ans);

      for (K = KHI; K >= KLO-1; K--) {
        kay            = KHI-K+1;
        NU_TURB(K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.nu_turb.value,THREEDIM);

  restore_mass(planet,NU_TURB_INDEX,NO_PHASE);

  return;
}

/*============== end of scalar_vertical_diffusion() ===============*/

/*============== uv_horizontal_subgrid() ==========================*/

/*
 * Wrapper function.
 */
void uv_horizontal_subgrid(planetspec  *planet,
                           EPIC_FLOAT **Buff2D)
{
  EPIC_FLOAT
    *pt_dudt,
    *pt_dvdt;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_horizontal_subgrid";

  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    if (grid.diffusion_direction == HORIZONTAL_AND_VERTICAL ||
        grid.diffusion_direction == JUST_HORIZONTAL) {
      uv_horizontal_diffusion(planet,Buff2D);
    }
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognzied grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  /*
   * Add divergence damping.
   */
  divergence_damping(planet);

  /*
   * NOTE: We tried applying hyperviscosity to the tendencies DUDT and DVDT, but
   *       found that grid-scale computational modes (stripes and checkerboards)
   *       emerged in U and V anyway. Hence, we apply hyperviscosity to U and V
   *       directly.
   */

  return;
}

/*============== end of uv_horizontal_subgrid() ===================*/

/*=================== uv_horizontal_diffusion() ===================*/

/*
 * We follow Tannehill, Anderson, and Pletcher (1997), p. 269
 * for the components of the viscosity in general coordinates.
 * Here,
 *   x1 = phi,         east longitude
 *   x2 = lambda,      planetographic latitude
 *   x3 = zeta,        hybrid vertical coordinate (a.k.a. sgth)
 *
 *   h1 = r(lambda),
 *   h2 = R(lambda),
 *   h3 = h/rho,       where rho is density and h = -(1/g)dp/dzeta is 
 *                     hybrid density (assumes plane-parallel geometry)
 *   u1 = u,
 *   u2 = v,
 *   u3 = 0,           shallow-atmosphere approximation
 *
 * We set all d/dzeta = 0 in this subroutine.  The vertical-gradient
 * terms are handled with an implicit timestep in uv_vertical_diffusion()
 * for numerical stability.
 *
 * NOTE: Currently not lagging the viscosity and density factors themselves in the 
 *       leapfrog case, which might possibly cause numerical stability problems, 
 *       but hasn't presented an obvious problem. We are lagging the (u,v) values 
 *       themselves via grid.it_uv_dis.
 */

#undef  COEFFD
#define COEFFD(j,i) coeffd[i+(j)*Iadim-Shift2d]
#undef  H_RHO
#define H_RHO(j,i) h_rho[i+(j)*Iadim-Shift2d]

void uv_horizontal_diffusion(planetspec  *planet,
			     EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    jj;
  register EPIC_FLOAT
    coeff,rho,nu,
    rln_inv,rlt_inv,
    h_inv,
    e11,e12,e22,e33;
  EPIC_FLOAT
    *tau11,*tau12,*tau22,*tau33,
    *uu,*vv,
    *coeffd,*h_rho;
  register EPIC_FLOAT
    m_2j,m_2jp1,
    n_2j,n_2jp1,n_2jp2;
  unsigned long
    nbytes_2d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_horizontal_diffusion"; 

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  tau11  = Buff2D[0];
  tau22  = Buff2D[1];
  tau33  = Buff2D[2];
  tau12  = Buff2D[3];
  coeffd = Buff2D[4];
  h_rho  = Buff2D[5];
  uu     = Buff2D[6];
  vv     = Buff2D[7];

  for (K = KLO; K < KHI; K++) { 
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        COEFFD(J,I) = planet->dynvisc+RHO3(K,J,I)*DIFFUSION_COEF_UV(K,J,I);
        H_RHO( J,I) = H3(K,J,I)/RHO3(K,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */

    /*
     * Copy U and V into UU and VV.
     */
    memcpy(uu,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    memcpy(vv,var.v.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));

    /*
     * Apply zonal_filter() to inputs, to control numerical instabilities at high latitudes.
     */
    zonal_filter(H2_INDEX,coeffd,NULL,TWODIM);
    zonal_filter(H2_INDEX,h_rho, NULL,TWODIM);
    zonal_filter(U_INDEX, uu,    NULL,TWODIM);
    zonal_filter(V_INDEX, vv,    NULL,TWODIM);

    /*
     * TAU11 and TAU22 adapt naturally to the h-grid,
     * whereas TAU33 is awkward (in the C-grid system).
     */
    for (J = JLO; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./grid.rln[jj+1];
      m_2jp1  = grid.m[jj+1];
      m_2j    = grid.m[jj  ];
      n_2jp2  = grid.n[jj+2];
      n_2jp1  = grid.n[jj+1];
      n_2j    = grid.n[jj  ];

      for (I = ILO; I <= IHI; I++) {
        e11    = (UU(J,I+1)-UU(J,I))*m_2jp1
	         +.5*(VV(J,I)+VV(J+1,I))*rln_inv*(grid.rln[jj+2]-grid.rln[jj])*n_2jp1;
        e22    = (VV(J+1,I)-VV(J,I))*n_2jp1;

        coeff  = .5/H_RHO(J,I);
        e33    = coeff*(UU(J,I+1)*(H_RHO(J,I+1)-H_RHO(J,I  ))*m_2jp1
                       +UU(J,I  )*(H_RHO(J,I  )-H_RHO(J,I-1))*m_2jp1);
        if (VV(J+1,I) != 0.) {
          e33 += coeff*VV(J+1,I)*(H_RHO(J+1,I)-H_RHO(J,  I))*n_2jp2;
        }
        if (VV(J,I) != 0.) {
          e33 += coeff*VV(J,  I)*(H_RHO(J,  I)-H_RHO(J-1,I))*n_2j;
        }

	/*
	 * TAU11 and TAU22 are calculated using e's and weighted by H_RHO.
         * The same for TAU33, except it is not weighted by H_RHO.
	 */
        coeff      = H_RHO(J,I)*(2./3.)*COEFFD(J,I);
	TAU11(J,I) = coeff*(2.*e11-e22-e33);
	TAU22(J,I) = coeff*(2.*e22-e11-e33);
	TAU33(J,I) = (coeff/H_RHO(J,I))*(2.*e33-e11-e22);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(tau11,TWODIM);
    bc_lateral(tau22,TWODIM);
    bc_lateral(tau33,TWODIM);

    /*
     * TAU12 falls naturally on the pv-grid.
     * Call vorticity() to handle the northern and southern poles or edges.
     *
     * NOTE: the call to vorticity() must be by every node, even though we
     *       are only using the results at the model's northen and southern extremes.
     *       If the special edge cases were to be broken out into a separate subroutine,
     *       this would not be necessary (this is not a pressing concern).
     */
    vorticity(planet,ON_SIGMATHETA,RELATIVE,uu,vv,NULL,tau12);

    if (JHI == grid.nj) {
      /*
       * Northern pole or edge.
       */
      coeff = 0.;
      J     = grid.nj;
      for (I = ILO; I <= IHI; I++) {
        coeff += COEFFD(J,I)*H_RHO(J,I);
      }
      coeff /= grid.ni;
      J = grid.nj+1;
      for (I = ILO; I <= IHI; I++) {
        /*
         * The minus sign comes from the fact that we need du/dy, not -du/dy.
         */
        TAU12(grid.nj+1,I) *= -coeff;
      }
    }
    if (JLO == grid.jlo) {
      /*
       * Southern pole or edge.
       */
      coeff = 0.;
      J     = 0;
      for (I = ILO; I <= IHI; I++) {
        coeff += COEFFD(J,I)*H_RHO(J,I);
      }
      coeff /= grid.ni;
      for (I = ILO; I <= IHI; I++) {
        TAU12(0,I) *= -coeff;
      }
    }
    /*
     * Fill in interior of TAU12.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj   = 2*J;
      m_2j = grid.m[jj];
      n_2j = grid.n[jj];

      for (I = ILO; I <= IHI; I++) {
        e12 = (VV(J,I)-VV(J,I-1))*m_2j 
	     +grid.rln[jj]*(UU(J,I)/grid.rln[jj+1]-UU(J-1,I)/grid.rln[jj-1])*n_2j;
	/*
	 * TAU12 is calculated using e12 and weighted by H/RHO2.
         * Average onto pv-grid.
	 */
        coeff    = .25*(COEFFD(J,  I  )*H_RHO(J,  I  )
                       +COEFFD(J,  I-1)*H_RHO(J,  I-1)
                       +COEFFD(J-1,I  )*H_RHO(J-1,I  )
                       +COEFFD(J-1,I-1)*H_RHO(J-1,I-1));
	TAU12(J,I) = coeff*e12;
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(tau12,TWODIM);
    
    /*
     * Horizontal-gradient diffusion of U.
     */
    for (J = JLO; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./(grid.rln)[jj+1];
      m_2jp1  = grid.m[jj+1];
      n_2jp1  = grid.n[jj+1];
      for (I = ILO; I <= IHI; I++) {
        h_inv                        = 2./(H2(K,J,I)+H2(K,J,I-1));
    	DUDT(grid.it_uv_tend,K,J,I) += h_inv*(m_2jp1*(TAU11(J,I)-TAU11(J,I-1))
                                             +n_2jp1*rln_inv*(grid.rln[jj+2]*TAU12(J+1,I)-grid.rln[jj]*TAU12(J,I)
                                                             +.5*(TAU12(J+1,I)+TAU12(J,I))*(grid.rln[jj+2]-grid.rln[jj]))
                                              -.5*(TAU33(J,I)+TAU33(J,I-1))*m_2jp1*(H_RHO(J,I)-H_RHO(J,I-1)));
      }
    }

    /*
     * Horizontal-gradient diffusion of V.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./(grid.rln)[jj];
      m_2j    = grid.m[jj];
      n_2j    = grid.n[jj];
      for (I = ILO; I <= IHI; I++) {
        h_inv                        = 2./(H2(K,J,I)+H2(K,J-1,I));
    	DVDT(grid.it_uv_tend,K,J,I) += h_inv*(m_2j*(TAU12(J,I)-TAU12(J,I-1))
                                             -n_2j*rln_inv*(grid.rln[jj+1]*TAU22(J,I)-grid.rln[jj-1]*TAU22(J-1,I)
                                                             -.5*(TAU11(J,I)+TAU11(J-1,I))*(grid.rln[jj+1]-grid.rln[jj-1]))
                                              -.5*(TAU33(J,I)+TAU33(J-1,I))*n_2j*(H_RHO(J,I)-H_RHO(J-1,I)));
      }
    }
  }
 
  return;
}

/*================== end of uv_horizontal_diffusion() =============*/

/*================== divergence_damping() =========================*/

/*
 * Add artificial damping of horizontal divergence to control
 * numerical instabilities associated with gravity waves.
 * See Skamarock and Klemp (1992, Mon. Wea. Rev. 120, 2109-2127).
 */

void divergence_damping(planetspec *planet)
{
  register int
    K,J,I;
  static int
    initialized = FALSE;
  static double
    max_nu[MAX_NU_ORDER+1];
  static EPIC_FLOAT
    nudiv,m0,
   *div,
   *buff2d;
  register double
    coef,rln;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="divergence_damping";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    buff2d = fvector(0,Nelem2d-1,dbmsname);
    div    = fvector(0,Nelem2d-1,dbmsname);

    set_max_nu(max_nu);
    nudiv = grid.nudiv_nondim*max_nu[2];
    
    /*
     * Set m0=1/dx0 corresponding to LAT0, the latitude above which measures are
     * taken to deal with the vanishing size of dx on the lon-lat grid.
     */
    if (strcmp(grid.geometry,"globe") == 0) {
      rln  = planet->re/sqrt(1.+ pow(planet->rp/planet->re*tan(LAT0*DEG),2.) );
      m0   = 1./(rln*grid.dln*DEG);
    }
    else if (strcmp(grid.geometry,"f-plane")  == 0) { 
      if (strcmp(grid.f_plane_map,"polar") == 0) {
        rln = .5*grid.f_plane_half_width;
        m0  = 1./(rln*grid.dln*DEG);
      }
      else {
        m0 = 1./(grid.f_plane_half_width/grid.ni);
      }
    }
    else {
      sprintf(Message,"unrecognized grid.geometry=%s");
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  for (K = KLO; K < KHI; K++) {
    /*
     * We calculate the horizontal divergence locally, instead of using
     * DIV_UV3, because we need to lag it for numerical stability.
     */
    memcpy(buff2d,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    zonal_filter(U_INDEX,buff2d,NULL,TWODIM);
    divergence(buff2d,
               var.v.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,
               div);
    for (J = JLO; J <= JHI; J++) {
      if (grid.lat[2*J+1] > LAT0) {
        /*
         * Taper damping coefficient at high latitudes to maintain numerical stability.
         */
        coef = nudiv*m0;
      }
      else {
        coef = nudiv*grid.m[2*J+1];
      }
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) += coef*(DIV(J,I)-DIV(J,I-1));
      }
    }
    /* No need to call bc_lateral() here. */

    for (J = JFIRST; J <= JHI; J++) {
      coef = nudiv*grid.n[2*J];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) += coef*(DIV(J,I)-DIV(J-1,I));
      }
    }
    /* No need to call bc_lateral() here. */
  }

  return;
}

/*================== end of divergence_damping() ==================*/

/*============== uv_vertical_subgrid() ============================*/

/*
 * Wrapper function.
 */
void uv_vertical_subgrid(planetspec  *planet,
                         EPIC_FLOAT **Buff2D)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_vertical_subgrid";

  if (strcmp(grid.turbulence_scheme,"Spalart-Allmaras DES") == 0) {
    if (grid.diffusion_direction == HORIZONTAL_AND_VERTICAL ||
        grid.diffusion_direction == JUST_VERTICAL) {
      uv_vertical_diffusion(planet,Buff2D);
    }
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }
  
  return;
}

/*============== end of uv_vertical_subgrid() =====================*/

/*============== uv_vertical_diffusion() ==========================*/

/*
 * The vertical-gradient diffusion terms for (u,v) are handled here using 
 * an implicit timestep to avoid numerical instabilities arising from
 * relatively small dz values.  The horizontal-gradient terms are
 * handled in uv_horizontal_diffusion().
 *
 * The shallow-atmosphere approximation has been applied, which eliminates
 * terms involving the vertical velocity.
 *
 * NOTE: The timeplane for U and V is grid.it_uv rather than 
 *       grid.it_uv_dis, because the implicit timestep does not need
 *       to be lagged.
 */

void uv_vertical_diffusion(planetspec  *planet,
                           EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kay;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    *tau_wall;
  static EPIC_FLOAT
    *zee,
    *aaa,
    *mu,
    *rho,
    *ans,
    *stab_factor0,
    *stab_factor2;
  unsigned long
    nbytes_2d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_vertical_diffusion";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);
  memset(Buff2D[0],0,nbytes_2d);
  tau_wall = Buff2D[0];

  if (!initialized) {
    /* Allocate memory: */
    zee          = fvector(0,KHI,    dbmsname);
    aaa          = fvector(0,KHI,    dbmsname);
    mu           = fvector(0,KHI,    dbmsname);
    rho          = fvector(0,KHI,    dbmsname);
    ans          = fvector(0,KHI,    dbmsname);
    stab_factor0 = fvector(0,2*KHI+1,dbmsname);
    stab_factor2 = fvector(0,2*KHI+1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Apply vertical diffusion to U.
   */
  tau_surface(planet,U_INDEX,tau_wall,Buff2D[1]);
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      stability_factor(U_INDEX,J,I,  stab_factor2);
      stability_factor(U_INDEX,J,I-1,stab_factor0);
      /*
       * Use no-flux boundary condition at top.
       */ 
      aaa[KHI] = 1.e+20;
      zee[KHI] = .5*(GZ3(0,J,I)+GZ3(0,J,I-1))/grid.g[2*J+1];
      for (K = KHI-1; K >= KLO; K--) {
        kay      = KHI-K;
        aaa[kay] = U(grid.it_uv,K,J,I);
        zee[kay] = .5*(GZ3(K,J,I)+GZ3(K,J,I-1))/grid.g[2*J+1];
        rho[kay] = .5*(RHO3(K,J,I)+RHO3(K,J,I-1));
      }

      for (K = KHI-1; K >= KLO; K--) {
        kay     = KHI-K;
        mu[kay] = planet->dynvisc
                 +.5*(stab_factor2[2*K]*RHO2(K,J,I  )*.5*(DIFFUSION_COEF_UV(K-1,J,I  )+DIFFUSION_COEF_UV(K,J,I  ))
                     +stab_factor0[2*K]*RHO2(K,J,I-1)*.5*(DIFFUSION_COEF_UV(K-1,J,I-1)+DIFFUSION_COEF_UV(K,J,I-1)));
      }
      /*
       * Set bottom boundary condition consistent with TAU_WALL.
       */
      zee[0] = .5*(GZ3(KHI,J,I)+GZ3(KHI,J,I-1))/grid.g[2*J+1];
      mu[ 0] = mu[1];
      aaa[0] = aaa[1]-TAU_WALL(J,I)*(zee[1]-zee[0])/mu[0];

      crank_nicolson(KHI-1,DT,zee,aaa,mu,rho,ans);

      for (K = KHI-1; K >= KLO; K--) {
        kay                 = KHI-K;
        U(grid.it_uv,K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

  /*
   * Apply vertical diffusion to V.
   */
  tau_surface(planet,V_INDEX,tau_wall,Buff2D[1]);
  for (J = JFIRST; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      stability_factor(V_INDEX,J,  I,stab_factor2);
      stability_factor(V_INDEX,J-1,I,stab_factor0);
      /*
       * Use no-flux boundary condition at top.
       */
      aaa[KHI] = 1.e+20;
      zee[KHI] = .5*(GZ3(0,J,I)+GZ3(0,J-1,I))/grid.g[2*J];
      for (K = KHI-1; K >= KLO; K--) {
        kay      = KHI-K;
        aaa[kay] = V(grid.it_uv,K,J,I);
        zee[kay] = .5*(GZ3(K,J,I)+GZ3(K,J-1,I))/grid.g[2*J];
        rho[kay] = .5*(RHO3(K,J,I)+RHO3(K,J-1,I));
      }
      for (K = KHI-1; K >= KLO; K--) {
        kay     = KHI-K;
        mu[kay] = planet->dynvisc+
                  .5*(stab_factor2[2*K]*RHO2(K,J,  I)*.5*(DIFFUSION_COEF_UV(K-1,J,  I)+DIFFUSION_COEF_UV(K,J,  I))
                     +stab_factor0[2*K]*RHO2(K,J-1,I)*.5*(DIFFUSION_COEF_UV(K-1,J-1,I)+DIFFUSION_COEF_UV(K,J-1,I)));
      }
      /*
       * Set bottom boundary condition consistent with TAU_WALL.
       *
       * NOTE: TAU_WALL is zero for gas giants.
       */
      zee[0] = .5*(GZ3(KHI,J,I)+GZ3(KHI,J-1,I))/grid.g[2*J];
      mu[ 0] = mu[1];
      aaa[0] = aaa[1]-TAU_WALL(J,I)*(zee[1]-zee[0])/mu[0];

      crank_nicolson(KHI-1,DT,zee,aaa,mu,rho,ans);

      for (K = KHI-1; K >= KLO; K--) {
        kay                 = KHI-K;
        V(grid.it_uv,K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

  return;
}

/*============== end of uv_vertical_diffusion() ====================*/

/*======================= make_arrays_subgrid() ====================*/

void make_arrays_subgrid(void)
{
  register int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="make_arrays_subgrid";

  /*
   * Allocate memory.
   */
  d_wall = fvector(0,Nelem3d-1,dbmsname);

  return;
}

/*======================= end of make_arrays_subgrid() =============*/

/*======================= free_arrays_subgrid() ====================*/

void free_arrays_subgrid(void)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_arrays_subgrid";

  /*
   * Free allocated memory.
   */
  free_fvector(d_wall,0,Nelem3d-1,dbmsname);

  return;
}

/*======================= end of free_arrays_subgrid() ============*/

/*======================= init_subgrid() ==========================*/

void init_subgrid(planetspec *planet)
{
  register int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_subgrid";

  if (!var.nu_turb.on) {
    sprintf(Message,"var.nu_turb.on is off",var.nu_turb.on);
    epic_error(dbmsname,Message);
  }

  /*
   * Initialize arrays.
   */  
  for (K = KLOPAD; K <= KHIPAD; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        NU_TURB(K,J,I) = planet->kinvisc*20.;
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  return;
}

/*======================= end of init_subgrid() ====================*/

/*======================= set_diffusion_coef() =====================*/

/*
  * Calculate turbulent diffusion coefficients.
  * We carry DIFFUSION_COEF_UV, DIFFUSION_COEF_THETA, and DIFFUSION_COEF_MASS
  * on the p3-grid.
  *
  * Molecular diffusion should be accounted for elsewhere.
  */

void set_diffusion_coef(planetspec *planet)
{
  register int
    K,J,I,
    kend;
  register EPIC_FLOAT
    chi3,fv1,nu_turb,turb,
    u_tan,kin,
    tmp,max_nu2;
  EPIC_FLOAT
    *u2d,*v2d;
  const EPIC_FLOAT
    cv1          = 7.1,
    cv1_3        = cv1*cv1*cv1;
  static double
    max_nu[MAX_NU_ORDER+1];
  static int
    initialized=FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_diffusion_coef";

  if (!initialized) {
    set_max_nu(max_nu);

    initialized = TRUE;
  }

  if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * For terrestrial planets, the bottom layer, K = KHI,
     * is treated separately below.
     */
    kend = KHI-1;
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    kend = KHI;
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  for (K = KLO; K <= kend; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      /*
       * Taper maximum viscosity to avoid numerical instability.
       */
      max_nu2 = max_nu[2]*MIN(1.,pow(grid.dy0*grid.m[2*J+1],-2.));

      for (I = ILOPAD; I <= IHIPAD; I++) {
        chi3  = NU_TURB(K,J,I)/planet->kinvisc;
        chi3 *= chi3*chi3;
        fv1   = chi3/(chi3+cv1_3);
        turb  = fv1*NU_TURB(K,J,I);
         
        tmp                         = turb;
        tmp                         = LIMIT_RANGE(0.,tmp,max_nu2);

        /*
         * Take the turbulent Prandtl number to be unity,
         * following Collins et al. (2004, NCAR/TN-464+STR).
         */
        tmp                         = turb;
        tmp                         = LIMIT_RANGE(0.,tmp,max_nu2);
        DIFFUSION_COEF_THETA(K,J,I) = tmp;

        /*
         * Mass diffusivity.
         * Currently using the same value for turbulent mass diffusivity
         * as for temperature.
         */
        tmp                        = turb;
        tmp                        = LIMIT_RANGE(0.,tmp,max_nu2);
        DIFFUSION_COEF_MASS(K,J,I) = tmp;
      }
    }
    /* No need to apply g() here. */
  }

  
  if (strcmp(planet->type,"terrestrial") == 0) {
    dwall_SA(planet,d_wall);

    K   = KHI-1;
    u2d = var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    v2d = var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    for (J = JLO; J <= JHI; J++) {
      /*
       * Taper maximum viscosity to avoid numerical instability.
       */
      max_nu2 = max_nu[2]*MIN(1.,pow(grid.dy0*grid.m[2*J+1],-2.));

      for (I = ILO; I <= IHI; I++) {
        nu_turb = NU_TURB(K,J,I);
        kin     = get_kin(planet,u2d,v2d,J,I);
        u_tan   = sqrt(2.*kin);
        turb    = law_of_the_wall(planet,K,J,I,NU_TURB_INDEX,nu_turb,u_tan);
         
        tmp                           = turb+planet->kinvisc;
        tmp                           = LIMIT_RANGE(0.,tmp,max_nu2);
        DIFFUSION_COEF_UV(KHI,J,I)    = tmp;

        tmp                           = turb;
        tmp                           = LIMIT_RANGE(0.,tmp,max_nu2);
        DIFFUSION_COEF_THETA(KHI,J,I) = tmp;

        tmp                           = turb;
        tmp                           = LIMIT_RANGE(0.,tmp,max_nu2);
        DIFFUSION_COEF_MASS(KHI,J,I)  = tmp;
      }
    }
  }
   
  /* Need to apply bc_lateral() here, because get_kin() above doesn't work on pads. */
  bc_lateral(var.diffusion_coef_uv.value,THREEDIM);
  bc_lateral(var.diffusion_coef_theta.value,THREEDIM);
  bc_lateral(var.diffusion_coef_mass.value,THREEDIM);

  return;
}

/*======================= end of set_diffusion_coef() ==============*/

/*======================= stability_factor() =======================*/
/*
 * Stability factor for vertical turbulent diffusion.
 * This function fills in the entire kk column of stab_factor[],
 * which should include the range [1,2*nk+1] as in the declaration 
 * of Ri[] below.
 */
#undef  MAX_RI
#define MAX_RI 100.

void stability_factor(int         index,
                      int         J,
                      int         I,
                      EPIC_FLOAT *stab_factor)
{
  register int
    K,kk;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *Ri;
  EPIC_FLOAT
    slope;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="stability_factor";

  if (!initialized) {
    /* Allocate memory. */
    Ri = fvector(0,2*KHI+1,dbmsname);

    initialized = TRUE;
  }

  if (strcmp(grid.stability_factor,"off") == 0) {
    for (kk = 1; kk <= 2*KHI+1; kk++) {
      stab_factor[kk] = 1.;
    }
    return;
  }
  else if (strcmp(grid.stability_factor,"Collins et al. (2004)") == 0) {
    /*
     * Set Richardson-number vector, limited by MAX_RI.
     */
    K = 0;
    Ri[1] = LIMIT_RANGE(-MAX_RI,RI2(K+1,J,I),MAX_RI);
    for (K = KLO; K < KHI; K++) {
      kk = 2*K;
      Ri[kk  ] = LIMIT_RANGE(-MAX_RI,RI2(K,J,I),MAX_RI);
      Ri[kk+1] = LIMIT_RANGE(-MAX_RI,.5*(RI2(K,J,I)+RI2(K+1,J,I)),MAX_RI);
    }
    K  = KHI;
    kk = 2*K;
    Ri[kk  ] = LIMIT_RANGE(-MAX_RI,RI2(K,J,I),MAX_RI);
    Ri[kk+1] = Ri[kk];

    /*
     * Limit the derivative dRi/dkk to a reasonable range,
     * to remove computational oscillations
     * in low N^2, low (du/dz)^2 regions.
     */
    for (kk = 1; kk <= 2*KHI; kk++) {
      slope = Ri[kk+1]-Ri[kk];
      if (-slope > .9*MAX_RI) {
        Ri[kk+1] = 0.;
      }
/********This code causes problems with Venus_LLR05 runs.
      if (-slope > .5*MAX_RI) {
        Ri[kk+1] = MAX(-MAX_RI,.5*(-MAX_RI+Ri[kk]));
      }
      else if (slope > .5*MAX_RI) {
        Ri[kk+1] = MIN( MAX_RI,.5*( MAX_RI+Ri[kk]));
      }
********/
    }

    switch (index) {
      case U_INDEX:
      case V_INDEX:
      case THETA_INDEX:
      case HDRY_INDEX:
      case HDRY3_INDEX:
      default:
        for (kk = 1; kk <= 2*KHI+1; kk++) {
          /*
           * From the NCAR Community Atmosphere Model (CAM 3.0),
           * as described by Collins et al. (2004, NCAR/TN-464+STR).
           * Evaluated in the layer.  Assumes RI2(K,J,I) has been calculated.
           */ 
          if (Ri[kk] > 0.) {
            stab_factor[kk] = 1./(1.+10.*Ri[kk]*(1.+8.*Ri[kk]));
          }
          else {
            stab_factor[kk] = sqrt(1.-18.*Ri[kk]);
          }
        }
        return;
      break;
    }
  }
  else {
    sprintf(Message,"unrecognized grid.stability_factor=%s",grid.stability_factor);
    epic_error(dbmsname,Message);
  }
}

/*======================= end of stability_factor() ================*/

/*======================= source_sink_turb() =======================*/

/*
 * Wrapper for turbulence source-sink subroutine.
 * This allows flexibility in what turbulence model is used without
 * having to change the hook to the rest of the EPIC model.
 *
 * NOTE: There is a name conflict with "source_sink_subgrid" and 
 *       LAM MPI, in the file lam_config_file.h.
 */
inline void source_sink_turb(planetspec  *planet,
                             EPIC_FLOAT **Buff2D)
{

  source_sink_SA(planet,Buff2D);

  return;
}

/*======================= end of source_sink_turb() ================*/

/*======================= source_sink_SA() =========================*/

void source_sink_SA(planetspec  *planet,
                    EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kk,jj;
  const EPIC_FLOAT    
    cw2             = 0.3,
    cw3_6           = pow(2.,6.),
    cb1             = 0.1355,
    cb2             = 0.622,
    cv1             = 7.1,
    cv1_3           = cv1*cv1*cv1,
    sigma           = 2./3.,
    kappa           = 0.41,
    kappa_2         = kappa*kappa,
    cw1             = (cb1/kappa_2)+((1.+cb2)/sigma),
    C_DES           = 0.65,
    nu_turb_limiter = 3000.;
  EPIC_FLOAT
    ptop,pbot,
    var1,var3,S,
    u_var1,u_var3,
    v_var1,v_var3,
    w_var1,w_var3,
    s11,s22,s33,s12,s13,s23,
    r_sa,g_sa,chi,fv1,fv2,fw,d_tilda,S_tilda,
    sig_z_conv,delta,
    tmp1,tmp2,dz_inv,
    nu_turb,u_tan,kin,chi3,
    brunt2,brunt,
    dnudt,
    h_rho;
  EPIC_FLOAT
    *u2d,*v2d;
  register EPIC_FLOAT
    m_2jp1,n_2j,n_2jp1,n_2jp2,
    m_2j_inv,n_2jp1_inv,
    mn_2jm1_inv,mn_2j_inv,
    mn_2jp1_inv,mn_2jp2_inv,
    mn_u,mn_v;
  static unsigned long
    nbytes_2d;
  static int
    solid_surface,
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="source_sink_SA";

  if (!initialized) {
    nbytes_2d = (unsigned long)(Nelem2d*sizeof(EPIC_FLOAT));

    if (strcmp(planet->type,"terrestrial") == 0) {
      solid_surface = TRUE;
    }
    else if (strcmp(planet->type,"gas-giant") == 0) {
      solid_surface = FALSE;
    }
    else {
      sprintf(Message,"unrecognized planet->type=%s",planet->type);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;
  }

  /*
   * D_WALL(K,J,I) is the distance to the wall.
   * For gas giants (no wall), dwall_SA sets this to FLOAT_MAX.
   */
  dwall_SA(planet,d_wall);

  for (K = KLO; K < KHI; K++) {
    kk = 2*K+1;
    /*
     * Assign zero'd memory to u,v:
     */
    memset(Buff2D[0],0,nbytes_2d);
    memset(Buff2D[1],0,nbytes_2d);
    u2d = Buff2D[0];
    v2d = Buff2D[1];
    /*
     * Averaging U values onto V grid and storing in U2D array.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;

      mn_2jp1_inv = 1./(grid.mn)[jj+1];
      mn_2jm1_inv = 1./(grid.mn)[jj-1];
      mn_v        = .5/(mn_2jm1_inv+mn_2jp1_inv);
      for (I = ILO; I <= IHI; I++) {
        U2D(J,I) = ((U(grid.it_uv_dis,K,J,  I)+U(grid.it_uv_dis,K,  J,I+1))*mn_2jp1_inv
                   +(U(grid.it_uv_dis,K,J-1,I)+U(grid.it_uv_dis,K,J-1,I+1))*mn_2jm1_inv)*mn_v;
      }
    }
    if (JLO == grid.jlo) {
      J      = grid.jlo;
      u_var1 = 0.;
      if (!IS_SPOLE) {
        /* 
         * Need U on southern-edge of channel.
         */
        for (I = ILO; I <= IHI; I++) {
          u_var1 += U(grid.it_uv_dis,K,J,I);
        }
        u_var1 /= grid.ni;
      }
      for (I = ILO; I <= IHI; I++) {
        U2D(J,I) = u_var1;
      }
    }
    if (JHI == grid.nj) {
      J      = grid.nj;
      u_var3 = 0.;
      if (!IS_NPOLE) {
        /* 
         * Need U on northern-edge of channel.
         */
        for (I = ILO; I <= IHI; I++) {
          u_var3 += U(grid.it_uv_dis,K,J,I);
        }
        u_var3 /= grid.ni;
      }
      for (I = ILO; I <= IHI; I++) {
        U2D(J+1,I) = u_var3;
      }
    }
    bc_lateral(u2d,TWODIM);

    /*
     * Averaging V values onto U grid and 
     * storing in V2D array.
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J;

      mn_2j_inv   = 1./(grid.mn)[jj];
      mn_2jp2_inv = 1./(grid.mn)[jj+2];
      mn_u        = .5/(mn_2j_inv+mn_2jp2_inv);
      for (I = ILO; I <= IHI; I++) {
	    V2D(J,I) = ((V(grid.it_uv_dis,K,J,  I)+V(grid.it_uv_dis,K,J,  I-1))*mn_2j_inv
		       +(V(grid.it_uv_dis,K,J+1,I)+V(grid.it_uv_dis,K,J+1,I-1))*mn_2jp2_inv)*mn_u;
      }
    }
    bc_lateral(v2d,TWODIM);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J;

      m_2jp1 = grid.m[jj+1];
      n_2jp1 = grid.n[jj+1];
      n_2jp2 = grid.n[jj+2];
      n_2j   = grid.n[jj  ];
      for (I = ILO; I <= IHI; I++) {
        delta   = delta_SA(planet,K,J,I);
        d_tilda = MIN(D_WALL(K,J,I),C_DES*delta);
        dz_inv  = grid.g[jj+1]/(GZ3(K-1,J,I)-GZ3(K+1,J,I));
        chi     = NU_TURB(K,J,I)/planet->kinvisc;
        fv1     = pow(chi,3.)/(pow(chi,3.)+cv1_3);
        fv2     = 1.-(chi/(1.+chi*fv1));
        /*
         *  Turbulence eddy viscosity production
         */
        s11 = (U(grid.it_uv_dis,K,J,I+1)-U(grid.it_uv_dis,K,J,I))*m_2jp1;
        s22 = (V(grid.it_uv_dis,K,J+1,I)-V(grid.it_uv_dis,K,J,I))*n_2jp1;
        s33 = (W3(K-1,J,I)-W3(K+1,J,I))/(grid.sigmatheta[kk-2]-grid.sigmatheta[kk+2]); 

        s12 = .5*((U2D(J+1,I  )-U2D(J,I))*n_2jp1
                 +(V2D(J,  I+1)-V2D(J,I))*m_2jp1);

        u_var3 = 0.5*(U(grid.it_uv_dis,K+1,J,I)+U(grid.it_uv_dis,K+1,J,I+1));
        u_var1 = 0.5*(U(grid.it_uv_dis,K-1,J,I)+U(grid.it_uv_dis,K-1,J,I+1));

	w_var3 = 0.5*(W3(K,J,I)+W3(K,J,I+1));
	w_var1 = 0.5*(W3(K,J,I)+W3(K,J,I-1));

        /* h*dzeta/dt = rho*dz/dt */
        h_rho  = H3(K,J,I)/RHO3(K,J,I);

	s13    = .5*((u_var1-u_var3)*dz_inv+h_rho*(w_var3-w_var1)*m_2jp1);
	
        v_var3 = 0.5*(V(grid.it_uv_dis,K+1,J,I)+V(grid.it_uv_dis,K+1,J+1,I));
        v_var1 = 0.5*(V(grid.it_uv_dis,K-1,J,I)+V(grid.it_uv_dis,K-1,J+1,I));

	w_var3 = 0.5*(W3(K,J,I)+W3(K,J+1,I));
	w_var1 = 0.5*(W3(K,J,I)+W3(K,J-1,I));

	s23    = .5*((v_var1-v_var3)*dz_inv+h_rho*(w_var3-w_var1)*n_2jp1);

        /*
         * 2.*s12 = s12+s21, etc. 
         */
	S = sqrt(2.*(s11*s11+s22*s22+s33*s33+2.*(s12*s12+s13*s13+s23*s23)));

        tmp1    = kappa_2*d_tilda*d_tilda;

        S_tilda = S+(NU_TURB(K,J,I)*fv2)/tmp1;
        tmp2    = S_tilda*tmp1;

        dnudt  = cb1*S_tilda*NU_TURB(K,J,I);

        dnudt += (cb2/sigma)*pow(((NU_TURB(K,J,I+1)-NU_TURB(K,J,I-1))*m_2jp1*.5),2.);
        dnudt += (cb2/sigma)*pow(.5*((NU_TURB(K,J+1,I)-NU_TURB(K,  J,I))*n_2jp2
                                    +(NU_TURB(K,J,  I)-NU_TURB(K,J-1,I))*n_2j  ),2.);
        dnudt += (cb2/sigma)*pow(((NU_TURB(K-1,J,I)-NU_TURB(K+1,J,I))*dz_inv),2.);

        /*
         *  SA-model turbulence destruction term.
         */
        if (tmp2 != 0.) {
          r_sa = NU_TURB(K,J,I)/tmp2;
          if (r_sa < 0.) {
            sprintf(Message,"r_sa=%g < 0., NU_TURB(%2d,%2d,%2d)=%g",r_sa,K,J,I,NU_TURB(K,J,I));
            epic_error(dbmsname,Message);
          }
          else if (r_sa < 10.) {
            g_sa = r_sa+cw2*(pow(r_sa,6.)-r_sa);
            fw   = g_sa*pow(((1.+cw3_6)/(pow(g_sa,6.)+cw3_6)),(1./6.));
          }
          else {
            fw = 1.;
          }
        }
        else {
          fw = 1.;
        }
 
        dnudt -= cw1*fw*pow((NU_TURB(K,J,I)/d_tilda),2.);

        /*
         * Buoyancy term.  
         *
         * NOTE: We need better physics and better rationale here.
         *
         * The parameter nu_turb_limiter prevents NU_TURB
         * from growing exponentially in an unbounded manner
         * in N^2 < 0 regions.
         */
        brunt2 = get_brunt2(planet,kk,J,I);
        brunt  = sqrt(fabs(brunt2));
        if (brunt2 < 0.) {
          dnudt -= NU_TURB(K,J,I)*exp(-NU_TURB(K,J,I)/nu_turb_limiter)*brunt2/brunt;
        }

	NU_TURB(K,J,I) += DT*dnudt;
        NU_TURB(K,J,I)  = MAX(1.e-4*planet->kinvisc,NU_TURB(K,J,I));
      }
    }
  }
  
  if (solid_surface) {
    K   = KHI-1;
    u2d = var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    v2d = var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        nu_turb          = NU_TURB(K,J,I);
        kin              = get_kin(planet,u2d,v2d,J,I);
        u_tan            = sqrt(2.*kin);
        nu_turb          = law_of_the_wall(planet,K,J,I,NU_TURB_INDEX,nu_turb,u_tan);
        nu_turb          = invert_fv1(planet,nu_turb);
        NU_TURB(KHI,J,I) = nu_turb;
        NU_TURB(KHI,J,I) = MAX(1.e-4*planet->kinvisc,NU_TURB(K,J,I));
      }
    }
  }
  
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.nu_turb.value,THREEDIM);

  return;
}

/*======================= end of source_sink_SA() ===========================*/

/*======================= dwall_SA() ========================================*/

void dwall_SA(planetspec *planet,
              EPIC_FLOAT *d_wall)
{
  register int
    K,J,I;
  
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dwall_SA";

  if (strcmp(planet->type,"gas-giant") == 0) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD;I++) {
          /*
           * Set distance to "wall" to FLOAT_MAX.
           */
          D_WALL(K,J,I) = FLOAT_MAX;
        }
      }
    }
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    for (K = KLO-1; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD;I++) {
          /*
           * Calculating distance from the wall.
           */
          D_WALL(K,J,I) = (GZ3(K,J,I)-grid.topo_scale*GZ_SURFACE(J,I))/grid.g[2*J+1];
        }
      }
    }
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  return;
}
      

/*======================= end of dwall_SA() ====================================*/

/*======================= delta_SA() ===========================================*/

/*
 * Compute delta in (5.9) of Dowling et al (2006) for the Spalart-Allmaras
 * turbulence model.
 *
 * NOTE: We have changed the definition from the Dowling et al (2006) paper.
 */

EPIC_FLOAT delta_SA(planetspec *planet,
                    int         K,
                    int         J,
                    int         I)
{
  register int
    jj = 2*J+1;
  EPIC_FLOAT
    dx,dy,dz,
    delta;
  
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="delta_SA";

  dx    = 1./grid.m[jj];
  dy    = 1./grid.n[jj];
  if (K < grid.nk) {
    dz = (GZ3(K-1,J,I)-GZ3(K+1,J,I))/grid.g[jj];
  }
  else {
    dz = (GZ3(K-1,J,I)-GZ3(grid.nk,J,I))/grid.g[jj];
  }

  delta = MIN(dx,dy);
  delta = MIN(delta,dz);

  return delta;
}

/*======================= end of delta_SA() =================================*/

/*======================= tau_surface() =====================================*/

/*
 * Calculate the shear stress at the surface.
 *
 * NOTE: Because of the call to get_kin, do not call this function
 *       from a pad position (like J=JLOPAD or I=IHIPAD).
 */

EPIC_FLOAT tau_surface(planetspec  *planet,
                       int          index,
                       EPIC_FLOAT  *tau_wall,
                       EPIC_FLOAT  *buffji)
{
  register int
    K,J,I;
  EPIC_FLOAT
    chi3,fv1,
    rho,diffusion_coef,nu_turb,dz,
    kin,u_tan,theta0,
    *u2d,
    *v2d;
  EPIC_FLOAT
    *kie;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1;
  static int
    nbytes_2d,
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="tau_surface";

  if (!initialized) {
    nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

    initialized = TRUE;
  }

  K   = grid.nk-1;
  u2d = var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d;
  v2d = var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d;

  memset(buffji,0,nbytes_2d);
  kie = buffji;

  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      KIE(J,I) = get_kin(planet,u2d,v2d,J,I);
    }
  }
  bc_lateral(kie,TWODIM);

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * NOTE: Placeholder, need better tau values for bottom of gas giant atmosphere.
     */
    if (index == U_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          TAU_WALL(J,I) = 0.;
        }
      }
    }
    else if (index == V_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          TAU_WALL(J,I) = 0.;
        }
      }
      /* No need to apply bc_lateral() here */
    }
    else if (index == THETA_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          TAU_WALL(J,I) = 0.;
        }
      }
    }
    else if (index == NU_TURB_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          TAU_WALL(J,I) = 0.;
        }
      }
    }
    else {
      sprintf(Message,"unrecognized index=%d",index);
      epic_warning(dbmsname,Message);
    }
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    dwall_SA(planet,d_wall);

    if (index == U_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          rho            = .5*(RHO3(K,J,I)+RHO3(K,J,I-1));
          nu_turb        = .5*(NU_TURB(K,J,I)+NU_TURB(K,J,I-1));
          dz             = .5*(GZ3(K,J,I)+GZ3(K,J,I-1)-GZ3(K+1,J,I)-GZ3(K+1,J,I-1))/grid.g[2*J+1];
          kin            = .5*(KIE(J,I)+KIE(J,I-1));
          u_tan          = sqrt(2.*kin);
          diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
          TAU_WALL(J,I)  = (planet->dynvisc+rho*diffusion_coef)*(U(grid.it_uv,K,J,I)-0.)/dz;	
        }
      }	
    }
    else if (index == V_INDEX) {
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          rho            = .5*(RHO3(K,J,I)+RHO3(K,J-1,I));
          nu_turb        = .5*(NU_TURB(K,J,I)+NU_TURB(K,J-1,I));
          dz             = .5*(GZ3(K,J,I)+GZ3(K,J-1,I)-GZ3(K+1,J,I)-GZ3(K+1,J-1,I))/grid.g[2*J];
          kin            = .5*(KIE(J,I)+KIE(J-1,I));
          u_tan          = sqrt(2.*kin);
          diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
          TAU_WALL(J,I)  = (planet->dynvisc+rho*diffusion_coef)*(V(grid.it_uv,K,J,I)-0.)/dz;
        }
      }
    }
    else if (index == NU_TURB_INDEX) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          rho            = RHO3(K,J,I);
          nu_turb        = NU_TURB(K,J,I);
          dz             = (GZ3(K,J,I)-GZ3(K+1,J,I))/grid.g[2*J+1];
          kin            = get_kin(planet,u2d,v2d,J,I);
          u_tan          = sqrt(2.*kin);
          diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
          nu_turb        = invert_fv1(planet,diffusion_coef);
          TAU_WALL(J,I)  = (planet->dynvisc+rho*nu_turb)*(u_tan-0.)/dz;
        }
      }
    }
    else {
      sprintf(Message,"unrecognized index=%d",index);
      epic_error(dbmsname,Message);
    }
  }
  else {
    sprintf(Message,"unrecognzied planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

#if EPIC_CHECK == 1
  /*
   * Screen for nan.
   */
  for (J = JLO; J <= JHI; J++) {
     for (I = ILO; I <= IHI; I++) {
      if (!finite(TAU_WALL(J,I))) {
        sprintf(Message,"index=%d, TAU_WALL(%d,%d)=%g",index,J,I,TAU_WALL(J,I));
        epic_error(dbmsname,Message);
      }
    }
  }
#endif

  /* Need to apply bc_lateral() here. */
  bc_lateral(tau_wall,TWODIM);

  return;
}

/*======================= end of tau_surface() ====================*/

/*======================= law_of_the_wall() =======================*/

/*
 * Raymond P. LeBeau
 *
 * Please put some comments here. -T.D.
 * Yes, please do! -A.H.
 *
 * Modified by Aaron Herrnstein on 03-13-06.  The extra argument "index"
 * is used to specify the h-grid, u-grid, or v-grid.
 *
 * Modified by Tim Dowling on 07-14-08. Added K argument.
 */

EPIC_FLOAT law_of_the_wall(planetspec *planet,
                           int         K,
                           int         J,
                           int         I,
                           int         index,
			   EPIC_FLOAT  t_vis,
			   EPIC_FLOAT  u_tan)
{
  EPIC_FLOAT
    t_vis_new,dwall,dwall_inv,
    u_tau,y_plus,u_plus,
    x1,x2,xl,dx,
    fl,f,
    rts,swap,
    chi3,fv1;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1,
    tol        = 1.e-6;
  int
    iter,
    itmax = 30;
  static int
    warned = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="law_of_the_wall";

  if (strcmp(planet->type,"gas-giant") == 0) {
    return 1.0e-4*planet->kinvisc;
  }
  else {
    chi3   = t_vis/planet->kinvisc;
    chi3  *= chi3*chi3;
    fv1    = chi3/(chi3+cv1_3);
    t_vis *= fv1;

    if (index == U_INDEX) {
      dwall = 0.5*(D_WALL(K,J,I)+D_WALL(K,J,I-1));
    }
    else if (index == V_INDEX  ||  index == PV3_INDEX) {
      dwall = 0.5*grid.m[2*J]*( D_WALL(K,J  ,I)/grid.m[2*J+1] 
                               +D_WALL(K,J-1,I)/grid.m[2*J-1] );
    }
    else {
      /*
       * Default is the h-grid.
       */
      dwall = D_WALL(K,J,I);
    }

    if (dwall*u_tan < 0.00023) {
      return 0.01*planet->kinvisc;
    }

    if (dwall > 0.) {
      dwall_inv = 1./dwall;
    }
    else {
      sprintf(Message,"JI=%d %d, dwall = %g",J,I,dwall);
      epic_error(dbmsname,Message);
    }

    x1 = sqrt(planet->kinvisc*(1.+ 0.01)*u_tan*dwall_inv);
    x2 = sqrt(planet->kinvisc*(1.+20.  )*u_tan*dwall_inv);
    fl = func_utau(x1,u_tan,dwall);
    f  = func_utau(x2,u_tan,dwall);

    if(fabs(fl) < fabs(f)) {
      rts  = x1;
      xl   = x2;
      swap = fl;
      fl   = f;
      f    = swap;
    }
    else {
      xl  = x1;
      rts = x2;
    }

    for (iter = 0; iter < itmax; iter++) {
      dx   = (xl-rts)*f/(f-fl);
      xl   = rts;
      fl   = f;
      rts += dx;
      f    = func_utau(rts,u_tan,dwall);
      if (f > .99e+20) {
        t_vis_new = 0.01*planet->kinvisc;
        return t_vis_new;
      }
      if (fabs(dx) < tol || f == 0.0) {
        u_tau     = rts;
        u_plus    = u_tan/u_tau;
        y_plus    = u_tau*dwall/planet->kinvisc;
        t_vis_new = MAX((u_tau*u_tau*dwall/u_tan)-planet->kinvisc,0.01*planet->kinvisc);
        break;
      }
    }
    if (iter < itmax) {
      t_vis_new = MAX((u_tau*u_tau*dwall/u_tan)-planet->kinvisc,0.01*planet->kinvisc);
    }
    else {
      t_vis_new = 20.*planet->kinvisc;
    }
  }

  return t_vis_new;
} 

/*=================== end of law_of_the_wall() ====================*/

/*=================== func_utau() =================================*/

EPIC_FLOAT func_utau(EPIC_FLOAT u_tau, 
                     EPIC_FLOAT u_tan, 
                     EPIC_FLOAT dwall)
{
  register EPIC_FLOAT
    res,u_plus,y_plus,
    yy;
  const EPIC_FLOAT
    C              = 5.,
    const1         = 0.127,
    const2         = 1./.41,
    const3         = 0.41*1.43e-3,
    one_third      = 1./3.,
    sqrt_one_third = sqrt(1./3.),
    pi_six         = M_PI/6.;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="func_utau";
  
  if (u_tau != 0.) {
    u_plus = u_tan/u_tau;
  }
  else {
    sprintf(Message,"u_tau=%g",u_tau);
    epic_error(dbmsname,Message);
  }
  y_plus = fabs(u_tau*dwall/planet->kinvisc);

  if (u_plus < 0.1) {
    return 1.e+20;
  }
  
  if (y_plus < 5.) {
    res = y_plus-u_plus;
  }
  else if (y_plus < 70.) {
    yy  = const1*y_plus;
    res = (      one_third*log((1.+yy)/(sqrt(yy*(yy-1.)+1.)))
           +sqrt_one_third*(atan(sqrt_one_third*(2.*yy-1.))+pi_six) )/const1
         +.25*const2*log(1.+const3*pow(y_plus,4.))-u_plus;
  }
  else {
    res = const2*log(y_plus)+C-u_plus;
  }
  
  return res;
}

/*=================== end of func_utau() ==========================*/

/*======================= invert_fv1() ============================*/

/* 
 *Calculate nu_turb (nu_tilde) from DIFF_COEF (nu_t).
 */
 
EPIC_FLOAT invert_fv1(planetspec *planet,
                      EPIC_FLOAT  t_vis)

{
  EPIC_FLOAT
    nu_turb,
    chi3,fv1,chi;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1,
    pi_factor  = M_PI*0.5/30.0;
  int
    iter,
    it_max=100;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="invert_fv1";

  /* The following fits are good to about +/-5%, will try to improve later -rpl */

  /* note chi here is t_vis/planet->kinvisc, not the high Re version */

  if (t_vis < 0.) {
    sprintf(Message,"t_vis=%g<0.",t_vis);
    epic_error(dbmsname,Message);
  }

  chi = t_vis/planet->kinvisc;
  if (chi > 30) {
    /* fvl ~ 1.0 */
    nu_turb = t_vis;
  }
  else if (chi < 0.4) {
    /* chi3 in denominator is small vs. 7.1^3  */
    nu_turb = planet->kinvisc*pow(chi*cv1_3,.25);
  }
  else {
    /* mid-range, fit */
    fv1 = cos(pi_factor*pow(chi,0.35));
    fv1 = 1.0-fv1*fv1;
    nu_turb = planet->kinvisc*chi/fv1;
  }

  return nu_turb;
} 

/*=================== end of invert_fv1() =========================*/

/* * * * * * * * * * end of epic_subgrid.c * * * * * * * * * * * * */




