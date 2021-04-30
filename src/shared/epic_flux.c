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

/*
 * Functions contained in this file:
 *    mass_flux_top()
 *    mass_flux_bot()
 *    calc_w()
 *    advection()
 *    hsu_advection()
 *    akima_advection()
 */

#include <epic.h>

/*======================= mass_flux_top() ==================================*/

/*
 * Hybrid mass flux for variable at the top of the model.
 *
 * For index == HDRY3_INDEX, this is WH.
 * For index == H_2O_INDEX, this is WH*(mixing ratio of H_2O), etc.
 *
 * Returns as an argument the hybrid vertical velocity at the
 * top of the model, wtop.
 */

EPIC_FLOAT mass_flux_top(planetspec *planet,
                         int         index,
                         int         J,
                         int         I,
                         EPIC_FLOAT *wtop)
{
  EPIC_FLOAT
    flux;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="mass_flux_top";

  flux  = 0.;
  if (wtop) {
   *wtop = 0.;
  }

  return flux;
}

/*======================= end of mass_flux_top() ===========================*/

/*======================= mass_flux_bot() ==================================*/

/*
 * Hybrid mass flux for variable at the bottom of the model.
 *
 * For index == HDRY3_INDEX, this is WH.
 * For index == H_2O_INDEX, this is WH*(mixing ratio of phase), etc.
 *
 * Returns as an argument the hybrid vertical velocity at the
 * bottom of the model, wbot.
 */

EPIC_FLOAT mass_flux_bot(planetspec *planet,
                         int         index,
                         int         J,
                         int         I,
                         EPIC_FLOAT *wbot)
{
  EPIC_FLOAT
    flux;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="mass_flux_bot";

  flux  = 0.;
  if (wbot) {
    *wbot = 0.;
  }

  return flux;
}

/*======================= end of mass_flux_bot() ===========================*/

/*======================= calc_w() =========================================*/

/*
 * Calculate the hybrid vertical velocity, W3 = d(sigmatheta)/dt,
 * at all the layer interfaces.
 */

#undef  SUM_DIV
#define SUM_DIV(k,j,i)               sum_div[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#undef  HORIZONTAL_DIV
#define HORIZONTAL_DIV(k,j,i) horizontal_div[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void calc_w(planetspec *planet)
{
  register int
    K,J,I,
    kk,m;
  register EPIC_FLOAT
    pbot,ptop,psigma,
    sigma,
    dsgth1,dsgth3_inv,
    dFdp,dFdpbot,dFdptop,dFdth,
    theta3,
    p2,p3,p4,dp,
    pbot2,pbot3,pbot4,dpbot,
    epsilon,
    dlnp2,dlnp4,
    dz10,dz21,dz32,dz20,dz30,dz31,
    df0,df1,df2,df3,
    dsumdiv,
    tmp;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    sgth,sgth_d;           
  static EPIC_FLOAT
    *sum_div;

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
    dbmsname[]="calc_w";

  if (!initialized) {
    /* Allocate memory */
    sum_div     = fvector(0,Nelem3d-1,dbmsname);
    initialized = TRUE;
  }

  /* 
   * Assign the hybrid vertical velocity, W3, for the
   * top and bottom of the model.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      mass_flux_top(planet,HDRY3_INDEX,J,I,&W3(  0,J,I));
      mass_flux_bot(planet,HDRY3_INDEX,J,I,&W3(KHI,J,I));
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Calculate the heating contribution to W3 for the
   * model's interior interfaces.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      ptop = P3(  0,J,I);
      pbot = P3(KHI,J,I);
      for (K = KLO; K < grid.k_sigma-1; K++) {
        sigma     = get_sigma(pbot,P3(K,J,I),ptop);
        dFdth     = g_sigma(sigma);
        W3(K,J,I) = dFdth*HEAT3(K,J,I)/EXNER3(K,J,I);
      }
      for (K = grid.k_sigma-1; K < KHI; K++) {
        W3(K,J,I) = 0.;
      }
    }
  }

  /*
   * NOTE: Gravity, either g or planet->g, is not to be
   *       confused with the function g(sigma) that is used to 
   *       define the hybrid coordinate.
   */

  /*
   * Calculate and store the summed horizontal-divergence terms.
   *
   * This accurate integration scheme is based on the one described and illustrated in Fig. 8(d) of 
   *   Leslie LM,  Purser RJ, 1992,  A comparative study of the performance of various
   *   vertical discretization schemes, Meteorol. Atmos. Phys. 50, 61-73.
   */

  /* Zero-out SUM_DIV array. */
  memset(sum_div,0,Nelem3d*sizeof(EPIC_FLOAT));

  /*
   * Start with top layer.
   */
  K = KLO;
  SUM_DIV(K,J,I) = .5*(P3(K-1,J,I)*DIV_UV3(K-1,J,I)
                      +P3(K,  J,I)*DIV_UV3(K,  J,I))*log(P3(K,J,I)/P3(K-1,J,I));

  /*
   * Continue down with K loop.
   */
  for (K = KLO+1; K < KHI; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        dz10 = log(P3(K-1,J,I)/P3(K-2,J,I));
        dz21 = log(P3(K,  J,I)/P3(K-1,J,I));
        dz32 = log(P3(K+1,J,I)/P3(K,  J,I));
        dz20 = dz21+dz10;
        dz30 = dz32+dz20;
        dz31 = dz32+dz21;

        df0  = P3(K-2,J,I)*DIV_UV3(K-2,J,I);
        df1  = P3(K-1,J,I)*DIV_UV3(K-1,J,I);
        df2  = P3(K,  J,I)*DIV_UV3(K,  J,I);
        df3  = P3(K+1,J,I)*DIV_UV3(K+1,J,I);

        dsumdiv = -(df0*(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)
                   +df1*(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)
                   +df2*(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)
                   +df3*(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)               )/12.;

        SUM_DIV(K,J,I) = SUM_DIV(K-1,J,I)+dsumdiv;
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  /*
   * Finish with bottom layer.
   */
  K = KHI;
  SUM_DIV(K,J,I) = SUM_DIV(K-1,J,I)+.5*(P3(K-1,J,I)*DIV_UV3(K-1,J,I)
                                       +P3(K,  J,I)*DIV_UV3(K,  J,I))*log(P3(K,J,I)/P3(K-1,J,I));

  /*
   * Add remaining terms to W3.
   */
  epsilon = pow(machine_epsilon(),1./3.);

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      ptop   = P3(0,             J,I);
      psigma = P3(grid.k_sigma-1,J,I);
      pbot   = P3(KHI,           J,I);

      for (K = KLO; K < KHI; K++) {
        theta3 = THETA(K,J,I);

        p3     = P3(K,J,I);
        pbot3  = pbot;

        dp     = p3*epsilon;
        dpbot  = pbot3*epsilon;

        p2     = MAX(p3-dp*.5,ptop);
        p4     = MIN(p3+dp*.5,pbot);

        /*
         * NOTE: Care must be taken when calculating vertical derivatives of zeta = F(theta,p,pbot)
         *       near psigma, because of the cusp in g_sigma(). Avoid the region above sigma_sigma.
         *       Sigma decreases with increasing p.
         */
        if (p2 < psigma && p4 >= psigma) {
          p2 = p3;
        }

        pbot2  = MAX(pbot3-dpbot*.5,p3);
        pbot4  = MAX(pbot3+dpbot*.5,p3);

        /*
         * NOTE: Care must be taken when calculating vertical derivatives of zeta = F(theta,p,pbot)
         *       near psigma, because of the cusp in g_sigma(). Avoid the region above sigma_sigma.
         *       Sigma increases with increasing pbot.
         */
        if (get_sigma(pbot4,p3,ptop) > grid.sigma_sigma && get_sigma(pbot2,p3,ptop) <= grid.sigma_sigma) {
          pbot4 = pbot3;
        }

        if (p2 < p4) {
          dFdp = (return_sigmatheta(theta3,p2,pbot3,ptop)
                 -return_sigmatheta(theta3,p4,pbot3,ptop))/(p2-p4);
        }
        else {
          sprintf(Message,"p2=%g p4=%g",p2,p4);
          epic_error(dbmsname,Message);
        }

        if (pbot2 < pbot4) {
          dFdpbot = (return_sigmatheta(theta3,p3,pbot2,ptop)
                    -return_sigmatheta(theta3,p3,pbot4,ptop))/(pbot2-pbot4);
        }
        else {
          sprintf(Message,"pbot2=%g pbot4=%g",pbot2,pbot4);
          epic_error(dbmsname,Message);
        }

        /*
         * Add the dFdp and dFdpbot terms to W3.
         */
        W3(K,J,I) -= SUM_DIV(K,  J,I)*dFdp
                    +SUM_DIV(KHI,J,I)*dFdpbot;
      }
    }
  }

  return;
}

/*======================= end of calc_w() ==================================*/

/*======================= advection() ======================================*/

/*
 * Advect scalar prognostic variables forward one timestep using 
 * the specified scheme(s).  Ensure non-negative status as necessary.
 */

void advection(planetspec  *planet,
               EPIC_FLOAT **Buff2D)

{
  register int
    is,ip,iq,
    K,J,I;
  static EPIC_FLOAT
    *unity,
    *old;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="advection";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    unity = fvector(0,Nelem3d-1,dbmsname);
    old   = fvector(0,Nelem3d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Advect the UNITY array to calculate the divergence appropriate for
   * the P3 grid, and use it in the advection calculations for
   * the non-mass variables carried on that grid (THETA, FPARA).
   */
  for (K = KLOPAD; K <= KHIPAD; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        UNITY(K,J,I) = 1.;
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.theta.advection_scheme) == 0) {
    hsu_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
  }
  else if (strcmp("Monotonized, 4th-order centered (Akima)",var.theta.advection_scheme) == 0) {
    akima_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
  }
  else {
    sprintf(Message,"unrecognized var.theta.advection_scheme=%s",var.theta.advection_scheme);
    epic_error(dbmsname,Message);
    exit(1);
  }

  if (var.theta.on) {
    /*
     * Advect THETA, which is carried on the interfaces.
     * Recall that THETA is carried on the interfaces and is only predicted in
     * the sigma region of the model, for K >= grid.k_sigma-1;
     */ 
    for (K = grid.k_sigma-1; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          OLD(K,J,I) = THETA(K,J,I);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
    if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.theta.advection_scheme) == 0) {
      hsu_advection(planet,THETA_INDEX,NO_PHASE,var.theta.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp("Monotonized, 4th-order centered (Akima)",var.theta.advection_scheme) == 0) {
      akima_advection(planet,THETA_INDEX,NO_PHASE,var.theta.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.theta.advection_scheme=%s",var.theta.advection_scheme);
      epic_error(dbmsname,Message);
      exit(1);
    }
    for (K = grid.k_sigma-1; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          THETA(K,J,I) -= OLD(K,J,I)*(UNITY(K,J,I)-1.);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  if (var.fpara.on) {
    /*
     * Advect FPARA, which is carried on the interfaces.
     * We assume here the same advection scheme is used for THETA and FPARA.
     */
    if (strcmp(var.fpara.advection_scheme,var.theta.advection_scheme) != 0) {
      sprintf(Message,"var.fpara.advection_scheme=%s != var.theta_advection_scheme=%s",
                       var.fpara.advection_scheme,var.theta.advection_scheme);
      epic_error(dbmsname,Message);
    }

    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          OLD(K,J,I) = FPARA(K,J,I);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
    if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.fpara.advection_scheme) == 0) {
      hsu_advection(planet,FPARA_INDEX,NO_PHASE,var.fpara.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp("Monotonized, 4th-order centered (Akima)",var.fpara.advection_scheme) == 0) {
      akima_advection(planet,FPARA_INDEX,NO_PHASE,var.fpara.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.fpara.advection_scheme=%s",var.fpara.advection_scheme);
      epic_error(dbmsname,Message);
      exit(1);
    }
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          FPARA(K,J,I) -= OLD(K,J,I)*(UNITY(K,J,I)-1.);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * The condensables are carried on the P3 grid, but may have a different
   * advection scheme than above, so update the UNITY array if necessary.
   */
  if (grid.nq > 0 &&
      strcmp(var.hdry.advection_scheme,var.theta.advection_scheme) != 0) {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          UNITY(K,J,I) = 1.;
        }
      }
    }
    /* No need to apply bc_lateral() here. */

    if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.hdry.advection_scheme) == 0) {
      hsu_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp("Monotonized, 4th-order centered (Akima)",var.hdry.advection_scheme) == 0) {
      akima_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.hdry.advection_scheme=%s",var.hdry.advection_scheme);
      epic_error(dbmsname,Message);
      exit(1);
    }
  }

  for (iq = 0; iq < grid.nq; iq++) {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          /*
           * Set Q = MAX(Q,Q_MIN), store value in OLD.
           */
          OLD(K,J,I) = Q(grid.is[iq],grid.ip[iq],K,J,I) = MAX(Q(grid.is[iq],grid.ip[iq],K,J,I),Q_MIN);
        }
      }
    }
    if (strcmp(var.species[grid.is[iq]].advection_scheme,"Predictor-corrector (Hsu, Konor, and Arakawa)") == 0) {
      hsu_advection(planet,grid.is[iq],grid.ip[iq],
                    var.species[grid.is[iq]].phase[grid.ip[iq]].q,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp(var.species[grid.is[iq]].advection_scheme,"Monotonized, 4th-order centered (Akima)") == 0) {
      akima_advection(planet,grid.is[iq],grid.ip[iq],
                      var.species[grid.is[iq]].phase[grid.ip[iq]].q,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.species[%d].advection_scheme=%s",
                      grid.is[iq],var.species[grid.is[iq]].advection_scheme);
      epic_error(dbmsname,Message);
    }

    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          Q(grid.is[iq],grid.ip[iq],K,J,I) -= OLD(K,J,I)*(UNITY(K,J,I)-1.);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  if (var.nu_turb.on) {
    /*
     * NU_TURB is carried on the P3 grid, but may have a different
     * advection scheme than above, so update the UNITY array if necessary.
     */
    if (strcmp(var.nu_turb.advection_scheme,var.hdry.advection_scheme) != 0) {
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            UNITY(K,J,I) = 1.;
          }
        }
      }
      /* No need to apply bc_lateral() here. */

      if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.nu_turb.advection_scheme) == 0) {
        hsu_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else if (strcmp("Monotonized, 4th-order centered (Akima)",var.nu_turb.advection_scheme) == 0) {
        akima_advection(planet,P3_INDEX,NO_PHASE,unity,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else {
        sprintf(Message,"unrecognized var.nu_turb.advection_scheme=%s",var.nu_turb.advection_scheme);
        epic_error(dbmsname,Message);
        exit(1);
      }
    }

    /*
     * Advect NU_TURB, which is carried on the layer interfaces.
     */
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          OLD(K,J,I) = NU_TURB(K,J,I);
        }
      }
    }
    /* No need to apply bc_lateral() here. */
    if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.nu_turb.advection_scheme) == 0) {
      hsu_advection(planet,NU_TURB_INDEX,NO_PHASE,var.nu_turb.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp("Monotonized, 4th-order centered (Akima)",var.nu_turb.advection_scheme) == 0) {
      akima_advection(planet,NU_TURB_INDEX,NO_PHASE,var.nu_turb.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.nu_turb.advection_scheme=%s",var.nu_turb.advection_scheme);
      epic_error(dbmsname,Message);
      exit(1);
    }
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          NU_TURB(K,J,I) -= OLD(K,J,I)*(UNITY(K,J,I)-1.);
          /*
           * Make sure NU_TURB is positive definite.
           */
          NU_TURB(K,J,I) = MAX(1.e-4*planet->kinvisc,NU_TURB(K,J,I));
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * All variables carried on the interfaces should now be advected. *
   * The following variables are carried in the layers.              *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Advect HDRY, which is carried in the layers.
   */
  if (var.hdry.on) {
    if (strcmp(var.hdry.advection_scheme,"Predictor-corrector (Hsu, Konor, and Arakawa)") == 0) {
      hsu_advection(planet,HDRY_INDEX,NO_PHASE,var.hdry.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else if (strcmp(var.hdry.advection_scheme,"Monotonized, 4th-order centered (Akima)") == 0) {
      akima_advection(planet,HDRY_INDEX,NO_PHASE,var.hdry.value,HORIZONTAL_AND_VERTICAL,Buff2D);
    }
    else {
      sprintf(Message,"unrecognized var.hdry.advection_scheme=%s",var.hdry.advection_scheme);
      epic_error(dbmsname,Message);
    }
    restore_mass(planet,HDRY_INDEX,NO_PHASE);
    /*
     * Need to update pressures, etc.
     */
    set_p2etc(planet,DONT_UPDATE_THETA);
  }

  /*
   * NOTE: A call to set_p2etc(planet,UPDATE_THETA) should be made to update
   *       theta = theta_diag, which is a function of pressure, 
   *       as soon as possible after returning from this subroutine.   
   */

  return;
}

/*======================= end of advection() ===============================*/

/*======================= hsu_advection() ==================================*/

/*
 *  Hsu and Arakawa's predictor-corrector advection, as modified 
 *  by eqs. (B.5) and (B.6) of Konor and Arakawa (1997).
 *  We apply it to each direction sequentially.
 *
 *  This is a positive-definite, divergence-form advection scheme
 *  that handles steep gradients.
 */

void hsu_advection(planetspec  *planet,
                   int          is,
                   int          ip,
                   EPIC_FLOAT  *buff3d,
                   int          direction,
                   EPIC_FLOAT **Buff2D)
{
  int   
    K,J,I,
    kk,k_first;
  static int
    initialized = FALSE;
  unsigned long
    nbytes_2d;
  register EPIC_FLOAT
    uhatp,uhatm,
    vhatp,vhatm,
    whatp,whatm,
    mu,vertvel,
    al,gap,gam,bep,bem,behatp,behatm,g_hsu,
    mn_2jp1,mn_2j,m_2j_inv,n_2jp1_inv,tmp,
    ep,a_min,
    da;
  EPIC_FLOAT
    *a,*a_diff1,*a_diff2,*a_old,
    *um,*up,*u2d,*vm,*vp,*v2d,
    *ff,*a_pred,*aaa;
  static EPIC_FLOAT
    *wm,*wp,
    *w_a,*w_a_pred,
    *w_a_min,*w_a_diff1,*w_a_diff2,
    *w_ff,*w_aaa;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="hsu_advection";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    wm        = fvector(0,KHI,  dbmsname);
    wp        = fvector(0,KHI,  dbmsname);
    w_a       = fvector(0,KHI,  dbmsname);
    w_a_pred  = fvector(0,KHI,  dbmsname);
    w_a_min   = fvector(0,KHI,  dbmsname);
    w_ff      = fvector(0,KHI+1,dbmsname);
    w_a_diff1 = fvector(0,KHI,  dbmsname);
    w_a_diff2 = fvector(0,KHI,  dbmsname);
    w_aaa     = fvector(0,KHI,  dbmsname);

    initialized = TRUE;
  }

  /*
   * Small constant ep (epsilon) introduced in Konor and Arakawa's 
   * eqs. (B.5,B.6) as a modification to Hsu and Arakawa's eqs. (6.15) and (6.16):
   */
  ep = 1.e-10;

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_VERTICAL) {
    /***********************
     * Vertical advection. *
     ***********************/

    /*
     * NOTE: For variables on the interface, in the bottom layer we are moving material
     *       down to the bottom interface and up from it, but in the top layer we are 
     *       not currently doing the same for the top interface of the model.
     */

    switch(is) {
      case HDRY_INDEX:
        if (!var.hdry.on) return;
        k_first = KLO;
        for (K = k_first; K <= KHI; K++) {
          w_a_min[K] = grid.h_min[K];
        }
      break;
      case P3_INDEX:
        k_first = KLO;
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = 0.;
        }
      break;
      case THETA_INDEX:
        if (!var.theta.on) return;
        k_first = grid.k_sigma-1;
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = 0.;
        }
      break;
      case NU_TURB_INDEX:
        if (!var.nu_turb.on) return;
        k_first = KLO;
        for (K = k_first; K <= KHI; K++) {
          w_a_min[K] = NU_TURB_MIN;
        }
      break;
      case FPARA_INDEX:
        if (!var.fpara.on) return;
        k_first = KLO;
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = Q_MIN;
        }
      break;
      default:
        if (is < FIRST_SPECIES || is > LAST_SPECIES) {
          sprintf(Message,"case is=%d not recognized",is);
          epic_error(dbmsname,Message);
        }
        if (!var.species[is].phase[ip].on) return;
        k_first = KLO;
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = Q_MIN;
        }
      break;
    }

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        switch(is) {
          case HDRY_INDEX:
            for (K = k_first; K <= KHI; K++) {
              w_a[K] = BUFF3D(K,J,I);
            }
          break;
          case THETA_INDEX:
          case FPARA_INDEX:
          case P3_INDEX:
          case NU_TURB_INDEX:
          default:
            for (K = k_first-1; K <= KHI; K++) {
              w_a[K] = BUFF3D(K,J,I);
            }
          break;
        }

        if (is == HDRY_INDEX) {
          /*
           * Variable is carried in layer.
           */
          for (K = k_first-1; K <= KHI; K++) {
            vertvel = W3(K,J,I);
            wm[K]   = MIN(0.,vertvel);
            wp[K]   = MAX(0.,vertvel);
          }
        }
        else {
          /*
           * Variable is carried on layer interfaces.
           */
          if (is >= FIRST_SPECIES && is <= LAST_SPECIES &&
              ip >= FIRST_PRECIP  && ip <= LAST_PRECIP    ) {
            for (K = k_first; K <= KHI; K++) {
              if (K == KLO) {
                /*
                 * Assume zero vertical velocity in top half of 
                 * topmost layer.
                 */
                vertvel = 0.;
              }
              else if (K == KHI) {
                /*
                 * Assume zero vertical velocity in bottom half of
                 * bottommost layer.
                 *
                 * NOTE: Doing otherwise leads to spurious advection
                 *       off the bottom of the model.
                 */
                vertvel = 0.;
                /*
                 * Include terminal velocity for precipitation.
                 */
                vertvel -= fabs(terminal_velocity(is,ip,2*K,J,I))*RHO2(K,J,I)/H2(K,J,I);
              }
              else {
                vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));
                /*
                 * Include terminal velocity for precipitation.
                 */
                vertvel -= fabs(terminal_velocity(is,ip,2*K,J,I))*RHO2(K,J,I)/H2(K,J,I);
              }

              wm[K] = MIN(0.,vertvel);
              wp[K] = MAX(0.,vertvel);
            }
          }
          else {
            for (K = k_first; K <= KHI; K++) {
              if (K == KLO) {
                /*
                 * Assume zero vertical velocity in top half of 
                 * topmost layer.
                 */
                vertvel = 0.;
              }
              else if (K == KHI) {
                /*
                 * Assume zero vertical velocity in bottom half of
                 * bottommost layer.
                 *
                 * NOTE: Doing otherwise leads to spurious advection
                 *       off the bottom of the model.
                 */
                vertvel = 0.;
              }
              else {
                vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));
              }
              wm[K]   = MIN(0.,vertvel);
              wp[K]   = MAX(0.,vertvel);
            }
          }
        }

        /* 
         * Predictor: 
         */
        if (is == HDRY_INDEX || is == NU_TURB_INDEX) {
          /*
           *  Here, w_a[K] is located above wp[K],wm[K].
           */
          if (k_first == 1) {
            w_ff[k_first-1] = mass_flux_top(planet,is,J,I,NULL);
          }
          else {
            w_ff[k_first-1] = 0.;
          }
          w_ff[KHI] = mass_flux_bot(planet,is,J,I,NULL);
          for (K = k_first; K < KHI; K++) {
            w_ff[K] = wp[K]*w_a[K+1]+wm[K]*w_a[K];
          }
          for (K = k_first; K <= KHI; K++) {
            da          = (w_ff[K]-w_ff[K-1])*grid.dsgth_inv[2*K]*DT;
            w_a_pred[K] = w_a[K]+da;
          }
        }
        else {
          /*
           *  Here, w_a[K] is located below wp[K],wm[K].
           */
          for (K = k_first; K <= KHI; K++) {
            w_ff[K] = wp[K]*w_a[K]+wm[K]*w_a[K-1];
          }
          w_ff[KHI+1] = 0.;
          w_a_pred[k_first-1] = w_a[k_first-1];
          for (K = k_first; K <= KHI; K++) {
            da          = (w_ff[K+1]-w_ff[K])*grid.dsgth_inv[2*K+1]*DT;
            w_a_pred[K] = w_a[K]+da;
          }
        }

        /* 
         * Corrector: 
         */
        if (is == HDRY_INDEX || is == NU_TURB_INDEX) {
          for (K = k_first; K < KHI; K++) {
            w_ff[K] = .5*(wp[K]*(w_a_pred[K  ]+w_a[K+1])
                         +wm[K]*(w_a_pred[K+1]+w_a[K  ]));
          }
          for (K = k_first; K <= KHI; K++) {
            if (K == KLO) {
              w_a_diff1[K] = w_a_pred[K]-w_a[K+1];
              w_a_diff2[K] = 0.;
              tmp          = fabs(w_a[K+1]-w_a[K])+ep;
              w_aaa[K]     = tmp*tmp;
            }
            else if (K == KHI) {
              w_a_diff1[K] = 0.;
              w_a_diff2[K] = w_a_pred[K]-w_a[K-1];
              tmp          = fabs(-w_a[K]+w_a[K-1])+ep;
              w_aaa[K]     = tmp*tmp;
            }
            else {
              w_a_diff1[K] = w_a_pred[K]-w_a[K+1];
              w_a_diff2[K] = w_a_pred[K]-w_a[K-1];
              tmp          = fabs(w_a[K+1]-2.*w_a[K]+w_a[K-1])+ep;
              w_aaa[K]     = tmp*tmp;
            }
          }

          for (K = k_first; K < KHI; K++) {
            whatp =  sqrt(wp[K]*wp[K+1]);
            whatm = -sqrt(wm[K]*wm[K-1]);
            /* Courant number, mu */
            mu     = (wp[K]-wm[K])*grid.dsgth_inv[2*K+1]*DT;
            al     = (1.+mu)/6.;
            gap    = w_aaa[K+1]/(w_aaa[K+1]+(w_a[K+1]-w_a_min[K+1])
                                           *(w_a[K  ]-w_a_min[K  ]));
            gap    = gap*gap;
            gam    = w_aaa[K  ]/(w_aaa[K  ]+(w_a[K+1]-w_a_min[K+1])
                                           *(w_a[K  ]-w_a_min[K  ]));
            gam    = gam*gam;
            bep    = 1.+(1./(2.*al)-1.)*gap;
            bem    = 1.+(1./(2.*al)-1.)*gam;
            behatp = 1.-gap;
            behatm = 1.-gam;
            g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                         -whatp*behatp*w_a_diff1[K+1]
                            +wm[K]*bem*w_a_diff2[K+1]
                         -whatm*behatm*w_a_diff2[K  ]);
            w_ff[K] += g_hsu;
          }
        }
        else {
          for (K = k_first; K <= KHI; K++) {
            w_ff[K] = .5*(wp[K]*(w_a_pred[K-1]+w_a[K  ])
                         +wm[K]*(w_a_pred[K  ]+w_a[K-1]));
          }

          for (K = k_first-1; K <= KHI; K++) {
            if (K == 0) {
              tmp      = fabs(w_a[K+1]-w_a[K])+ep;
              w_aaa[K] = tmp*tmp;
            }
            else if (K == KHI) {
              tmp      = fabs(-w_a[K]+w_a[K-1])+ep;
              w_aaa[K] = tmp*tmp;
            }
            else {
              tmp      = fabs(w_a[K+1]-2.*w_a[K]+w_a[K-1])+ep;
              w_aaa[K] = tmp*tmp;
            }
          }

          w_a_diff2[k_first-1] = 0.;
          for (K = k_first; K <= KHI; K++) {
            w_a_diff1[K] = w_a_pred[K-1]-w_a[K  ];
            w_a_diff2[K] = w_a_pred[K  ]-w_a[K-1];
          }

          for (K = k_first; K <= KHI; K++) {
            if (K == KLO) {
              whatm = 0.;
            }
            else {
              whatm = -sqrt(wm[K]*wm[K-1]);
            }
            if (K == KHI) {
              whatp = 0.;
            }
            else {
              whatp = sqrt(wp[K]*wp[K+1]);
            }
            /* Courant number, mu */
            mu     = (wp[K]-wm[K])*grid.dsgth_inv[2*K]*DT;
            al     = (1.+mu)/6.;
            gap    = w_aaa[K  ]/(w_aaa[K  ]+(w_a[K  ]-w_a_min[K  ])
                                           *(w_a[K-1]-w_a_min[K-1]));
            gap    = gap*gap;
            gam    = w_aaa[K-1]/(w_aaa[K-1]+(w_a[K  ]-w_a_min[K  ])
                                           *(w_a[K-1]-w_a_min[K-1]));
            gam    = gam*gam;
            bep    = 1.+(1./(2.*al)-1.)*gap;
            bem    = 1.+(1./(2.*al)-1.)*gam;
            behatp = 1.-gap;
            behatm = 1.-gam;
            if (K < KHI) {
              g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                           -whatp*behatp*w_a_diff1[K+1]
                              +wm[K]*bem*w_a_diff2[K  ]
                           -whatm*behatm*w_a_diff2[K-1]);
            }
            else {
              g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                              +wm[K]*bem*w_a_diff2[K  ]
                           -whatm*behatm*w_a_diff2[K-1]);
            }
            w_ff[K] += g_hsu;
          }
          w_ff[KHI+1] = 0.;
        }

        switch(is) {
          case HDRY_INDEX:
            for (K = k_first; K <= KHI; K++) {
              da             = (w_ff[K]-w_ff[K-1])*grid.dsgth_inv[2*K]*DT;
              BUFF3D(K,J,I) += da;
            }
          break;
          case THETA_INDEX:
          case FPARA_INDEX:
          case P3_INDEX:
          case NU_TURB_INDEX:
          default:
            for (K = k_first; K <= KHI; K++) {
              da             = (w_ff[K+1]-w_ff[K])*grid.dsgth_inv[2*K+1]*DT;
              BUFF3D(K,J,I) += da;
            }
          break;
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_HORIZONTAL) {

    /*************************
     * Meridional advection. *
     *************************/

    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;

      switch(is) {
        case HDRY_INDEX:
          if (!var.hdry.on) return;
          a_min = grid.h_min[K];
        break;
        case P3_INDEX:
          a_min = 0.;
        break;
        case THETA_INDEX:
          if (!var.theta.on) return;
          if (K <= grid.k_sigma-2) continue;
          a_min = 0.;
        break;
        case NU_TURB_INDEX:
          if (!var.nu_turb.on) return;
          a_min = NU_TURB_MIN;
        break;
        case FPARA_INDEX:
          if (!var.fpara.on) return;
          a_min = Q_MIN;
        break;
        default:
          if (!var.species[is].phase[ip].on) return;
          a_min = Q_MIN;
        break;
      }

      a = buff3d+(K-Kshift)*Nelem2d;

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);

      vm  = Buff2D[0];
      vp  = Buff2D[1];
      v2d = Buff2D[2];

      if (is == HDRY_INDEX) {
        for (J = JFIRST; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            V2D(J,I) = get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(v2d,TWODIM);
      }
      else {
        for (J = JLOPAD; J <= JHIPADPV; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            V2D(J,I) = V(grid.it_uv,K,J,I);
          }
        }
        /* No need to apply bc_lateral() here. */
      }

      /*
       * Take range of J for VM and VP to be JLOPAD, JHIPADPV, thereby
       * eliminating the need to apply bc_lateral(). It is not necessary
       * to treat the I index in this manner for VM and VP.
       */
      for (J = JLOPAD; J <= JHIPADPV; J++) {
        if (fabs(grid.lat[2*J]) == 90.) {
          m_2j_inv = 0.;
        }
        else {
          m_2j_inv = 1./(grid.m)[2*J];
        }
        for (I = ILO; I <= IHI; I++) {
          VM(J,I) = MIN(0.,V2D(J,I)*m_2j_inv);
          VP(J,I) = MAX(0.,V2D(J,I)*m_2j_inv);
        }
      }
      /* 
       * Do not need to apply bc_lateral() to VM and VP. 
       * Done with v2d memory.
       */

      /* 
       * Predictor: 
       */
      memset(Buff2D[4],0,nbytes_2d);
      ff = Buff2D[4];
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = VP(J,I)*A(J-1,I)+VM(J,I)*A(J,I);
        }
      }
      bc_lateral(ff,TWODIM);

      memset(Buff2D[5],0,nbytes_2d);
      a_pred = Buff2D[5];
      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da          = (FF(J,I)-FF(J+1,I))*mn_2jp1*DT;
          A_PRED(J,I) = A(J,I)+da;
        }
      }
      bc_lateral(a_pred,TWODIM);

      /*
       * Corrector:
       */
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = .5*(VP(J,I)*(A_PRED(J,I  )+A(J-1,I))
                       +VM(J,I)*(A_PRED(J-1,I)+A(J,I  )));
        }
      }

      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);
      a_diff1 = Buff2D[2];
      a_diff2 = Buff2D[3];
      if (strcmp(grid.geometry,"f-plane") == 0 
          && strcmp(grid.f_plane_map,"cartesian") == 0) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            A_DIFF1(J,I) = A_PRED(J,I)-A(J-1,I);
            A_DIFF2(J,I) = A_PRED(J,I)-A(J+1,I);
          }
        }
      }
      else if (strcmp(grid.geometry,"globe") == 0 ||
              (strcmp(grid.geometry,"f-plane") == 0 
               && strcmp(grid.f_plane_map,"polar") == 0))  {
        /* NOTE: do interior points, ie J = JFIRST; J < JHI */
        for (J = JFIRST; J < JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            A_DIFF1(J,I) = A_PRED(J,I)-A(J-1,I);
            A_DIFF2(J,I) = A_PRED(J,I)-A(J+1,I);
          }
        }
        /* Special a_diff's at ends. */
        J = JLO;
        if (JLO == grid.jlo) {
          for (I = ILO; I <= IHI; I++) {
            A_DIFF1(J,I) = 0.;
            A_DIFF2(J,I) = A_PRED(J,I)-A(J+1,I);
          }
        }
        J = JHI;
        if (JHI == grid.nj) {
          /* edge */
          for (I = ILO; I <= IHI; I++) {
            A_DIFF1(J,I) = A_PRED(J,I)-A(J-1,I);
            A_DIFF2(J,I) = 0.;
          }
        }
        else {
          /* interior */
          for (I = ILO; I <= IHI; I++) {
            A_DIFF1(J,I) = A_PRED(J,I)-A(J-1,I);
            A_DIFF2(J,I) = A_PRED(J,I)-A(J+1,I);
          }
        }
      } 
      else {
        fprintf(stderr,"Unrecognized geometry in epic_timestep \n");
        exit(1);
      }
      bc_lateral(a_diff1,TWODIM);
      bc_lateral(a_diff2,TWODIM);

      /* Finished with a_pred memory. */
      memset(Buff2D[5],0,nbytes_2d);
      aaa = Buff2D[5];

      if (strcmp(grid.geometry,"f-plane") == 0 
          && strcmp(grid.f_plane_map,"cartesian") == 0) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            tmp      = fabs(A(J-1,I)-2.*A(J,I)+A(J+1,I))+ep;
            AAA(J,I) = tmp*tmp;
          }
        }
      }
      else if (strcmp(grid.geometry,"globe") == 0 ||
              (strcmp(grid.geometry,"f-plane") == 0 
               && strcmp(grid.f_plane_map,"polar") == 0))  {
        /* NOTE: do interior points, ie J = JFIRST; J < JHI */
        for (J = JFIRST; J < JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            tmp      = fabs(A(J-1,I)-2.*A(J,I)+A(J+1,I))+ep;
            AAA(J,I) = tmp*tmp;
          }
        }
        /* Special a_diff's at ends. */
        J = JLO;
        if (JLO == grid.jlo) {
          for (I = ILO; I <= IHI; I++) {
            tmp      = fabs(-A(J,I)+A(J+1,I))+ep;
            AAA(J,I) = tmp*tmp;
          }
        }
        J = JHI;
        if (JHI == grid.nj) {
          /* edge */
          for (I = ILO; I <= IHI; I++) {
            tmp      = fabs(A(J-1,I)-A(J,I))+ep;
            AAA(J,I) = tmp*tmp;
          }
        }
        else {
          /* interior */
          for (I = ILO; I <= IHI; I++) {
            tmp      = fabs(A(J-1,I)-2.*A(J,I)+A(J+1,I))+ep;
            AAA(J,I) = tmp*tmp;
          }
        }
      } 
      bc_lateral(aaa,TWODIM);

      for (J = JFIRST; J <= JHI; J++) {
        mn_2j = (grid.mn)[2*J];
        for (I = ILO; I <= IHI; I++) {
          vhatp  =  sqrt(VP(J,I)*VP(J-1,I));
          vhatm  = -sqrt(VM(J,I)*VM(J+1,I));
          /* Courant number, mu */
          mu     = (VP(J,I)-VM(J,I))*mn_2j*DT;
          al     = (1.+mu)/6.;
          gap    = AAA(J-1,I)/(AAA(J-1,I)+(A(J-1,I)-a_min)*(A(J,I)-a_min));
          gap    = gap*gap;
          gam    = AAA(J,I  )/(AAA(J,I  )+(A(J-1,I)-a_min)*(A(J,I)-a_min));
          gam    = gam*gam;
          bep    = 1.+(1./(2.*al)-1.)*gap;
          bem    = 1.+(1./(2.*al)-1.)*gam;
          behatp = 1.-gap;
          behatm = 1.-gam;
          g_hsu  = -al*( VP(J,I)*bep*A_DIFF1(J,  I)
                       -vhatp*behatp*A_DIFF1(J-1,I)
                        +VM(J,I)*bem*A_DIFF2(J-1,I)
                       -vhatm*behatm*A_DIFF2(J,  I));
          FF(J,I) += g_hsu;
        }  
      }
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da      = (FF(J,I)-FF(J+1,I))*mn_2jp1*DT;
          A(J,I) += da;
        }
      }
      bc_lateral(a,TWODIM);
    }

    /********************
     * Zonal advection. *
     ********************/

    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;

      switch(is) {
        case HDRY_INDEX:
          if (!var.hdry.on) return;
          a_min = grid.h_min[K];
        break;
        case P3_INDEX:
          a_min = 0.;
        break;
        case THETA_INDEX:
          if (!var.theta.on) return;
          if (K <= grid.k_sigma-2) continue;
          a_min = 0.;
        break;
        case NU_TURB_INDEX:
          if (!var.nu_turb.on) return;
          a_min = NU_TURB_MIN;
        break;
        case FPARA_INDEX:
          if (!var.fpara.on) return;
          a_min = Q_MIN;
        break;
        default:
          if (!var.species[is].phase[ip].on) return;
          a_min = Q_MIN;
        break;
      }

      /*
       * Apply zonal_filter() to a copy of variable.
       */
      a = Buff2D[6];
      memcpy(a,buff3d+(K-Kshift)*Nelem2d,Nelem2d*sizeof(EPIC_FLOAT));
      zonal_filter(is,a,NULL,TWODIM);

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);
      um  = Buff2D[0];
      up  = Buff2D[1];
      u2d = Buff2D[2];

      if (is == HDRY_INDEX) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            U2D(J,I) = get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(u2d,TWODIM);
      }
      else {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            U2D(J,I) = U(grid.it_uv,K,J,I);
          }
        }
        /* Do not need to call bc_lateral() here. */
      }

      /*
       * Take range of I for UM and UP to be ILOPAD, IHIPAD, thereby
       * eliminating the need to apply bc_lateral(). It is not necessary
       * to treat the J index in this manner for UM and UP.
       */
      for (J = JLO; J <= JHI; J++) {
        n_2jp1_inv = 1./(grid.n)[2*J+1];
        for (I = ILOPAD; I <= IHIPAD; I++) {
          UM(J,I) = MIN(0.,U2D(J,I)*n_2jp1_inv);
          UP(J,I) = MAX(0.,U2D(J,I)*n_2jp1_inv);
        }
      }
      /* No need to apply bc_lateral(). */

      /*
       * Done with u2d memory. 
       */

      /* 
       * Predictor: 
       */
      memset(Buff2D[4],0,nbytes_2d);
      ff = Buff2D[4];
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = UP(J,I)*A(J,I-1)+UM(J,I)*A(J,I);
        }
      }
      bc_lateral(ff,TWODIM);

      memset(Buff2D[5],0,nbytes_2d);
      a_pred = Buff2D[5];
      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da          = (FF(J,I)-FF(J,I+1))*mn_2jp1*DT;
          A_PRED(J,I) = A(J,I)+da;
        }
      }
      bc_lateral(a_pred,TWODIM);

      /* 
       * Corrector: 
       */
      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);
      a_diff1 = Buff2D[2];
      a_diff2 = Buff2D[3];
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = .5*(UP(J,I)*(A_PRED(J,I  )+A(J,I-1))
                       +UM(J,I)*(A_PRED(J,I-1)+A(J,I  )));
          A_DIFF1(J,I) = A_PRED(J,I)-A(J,I-1);
          A_DIFF2(J,I) = A_PRED(J,I)-A(J,I+1);
        }
      }
      bc_lateral(a_diff1,TWODIM);
      bc_lateral(a_diff2,TWODIM);

      /* Finished with a_pred memory. */

      memset(Buff2D[5],0,nbytes_2d);
      aaa = Buff2D[5];
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          tmp      = fabs(A(J,I-1)-2.*A(J,I)+A(J,I+1))+ep;
          AAA(J,I) = tmp*tmp;
        }
      }
      bc_lateral(aaa,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          uhatp  =  sqrt(UP(J,I)*UP(J,I-1));
          uhatm  = -sqrt(UM(J,I)*UM(J,I+1));
          /* Courant number, mu */
          mu     = (UP(J,I)-UM(J,I))*mn_2jp1*DT;
          al     = (1.+mu)/6.;
          gap    = AAA(J,I-1)/(AAA(J,I-1)+(A(J,I-1)-a_min)*(A(J,I)-a_min));
          gap    = gap*gap;
          gam    = AAA(J,I  )/(AAA(J,I  )+(A(J,I-1)-a_min)*(A(J,I)-a_min));
          gam    = gam*gam;
          bep    = 1.+(1./(2.*al)-1.)*gap;
          bem    = 1.+(1./(2.*al)-1.)*gam;
          behatp = 1.-gap;
          behatm = 1.-gam;
          g_hsu  = -al*( UP(J,I)*bep*A_DIFF1(J,I  )
                       -uhatp*behatp*A_DIFF1(J,I-1)
                        +UM(J,I)*bem*A_DIFF2(J,I-1)
                       -uhatm*behatm*A_DIFF2(J,I  ));
          FF(J,I) += g_hsu;
        }  
      }
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da             = (FF(J,I)-FF(J,I+1))*mn_2jp1*DT;
          BUFF3D(K,J,I) += da;
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(buff3d,THREEDIM);
  }

  return;
}

/*======================= end of hsu_advection() ===========================*/

/*======================= akima_advection() ================================*/

/*
 *  Monotonized, fourth-order, centered, divergence-form advection scheme.
 *
 *  From (4.9)-(4.12) of Shchpetkin & McWilliams (2003, JGR 108, doi:10.1029/2001JC001047).
 *  Alex Shchepetkin corresponded with Tim Dowling by email that this simple advection scheme
 *  is based on a midpoint interpolation due to Akima, and that it behaves well in practice.
 */

void akima_advection(planetspec  *planet,
                     int          is,
                     int          ip,
                     EPIC_FLOAT  *buff3d,
                     int          direction,
                     EPIC_FLOAT **Buff2D)
{
  int   
    K,J,I,
    kk,k_first;
  static int
    initialized = FALSE;
  unsigned long
    nbytes_2d;
  register EPIC_FLOAT
    aval,vertvel,
    mn_2jp1,mn_2j,m_2j_inv,n_2jp1_inv,tmp;
  EPIC_FLOAT
    *a,*deltaa,*dela,
    *u2d,*v2d,
    *ff;
  static EPIC_FLOAT
    *w_a,*w_deltaa,*w_dela,
    *w,
    *w_ff;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="akima_advection";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    w         = fvector(0,KHI,  dbmsname);
    w_a       = fvector(0,KHI,  dbmsname);
    w_deltaa  = fvector(0,KHI,  dbmsname);
    w_dela    = fvector(0,KHI,  dbmsname);
    w_ff      = fvector(0,KHI+1,dbmsname);

    initialized = TRUE;
  }

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_VERTICAL) {
    /***********************
     * Vertical advection. *
     ***********************/

    /*
     * NOTE: For variables on the interface, in the bottom layer we are moving material
     *       down to the bottom interface and up from it, but in the top layer we are 
     *       not currently doing the same for the top interface of the model.
     */

    switch(is) {
      case HDRY_INDEX:
        if (!var.hdry.on) return;
        k_first = KLO;
      break;
      case P3_INDEX:
        k_first = KLO;
      break;
      case NU_TURB_INDEX:
        if (!var.nu_turb.on) return;
        k_first = KLO;
      break;
      case THETA_INDEX:
        if (!var.theta.on) return;
        k_first = grid.k_sigma-1;
      break;
      case FPARA_INDEX:
        if (!var.fpara.on) return;
        k_first = KLO;
      break;
      default:
        if (!var.species[is].phase[ip].on) return;
        k_first = KLO;
      break;
    }

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        switch(is) {
          case HDRY_INDEX:
            for (K = k_first; K <= KHI; K++) {
              w_a[K] = BUFF3D(K,J,I);
            }
          break;
          case THETA_INDEX:
          case FPARA_INDEX:
          case P3_INDEX:
          case NU_TURB_INDEX:
          default:
            for (K = k_first-1; K <= KHI; K++) {
              w_a[K] = BUFF3D(K,J,I);
            }
          break;
        }

        if (is == HDRY_INDEX) {
          /*
           * Variable is carried in layer.
           */
          for (K = k_first-1; K <= KHI; K++) {
            vertvel = W3(K,J,I);
            w[K]    = vertvel;
          }
        }
        else {
          /*
           * Variable is carried on layer interfaces.
           */
          if (is >= FIRST_SPECIES && is <= LAST_SPECIES &&
              ip >= FIRST_PRECIP  && ip <= LAST_PRECIP    ) {
            for (K = k_first; K <= KHI; K++) {
              if (K == KLO) {
                /*
                 * Assume zero vertical velocity in top half of 
                 * topmost layer.
                 */
                vertvel = 0.;
              }
              else if (K == KHI) {
                /*
                 * Assume zero vertical velocity in bottom half of
                 * bottommost layer.
                 *
                 * NOTE: Doing otherwise leads to spurious advection
                 *       off the bottom of the model.
                 */
                vertvel = 0.;
                /*
                 * Include terminal velocity for precipitation.
                 */
                vertvel -= fabs(terminal_velocity(is,ip,2*K,J,I))*RHO2(K,J,I)/H2(K,J,I);
              }
              else {
                vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));
                /*
                 * Include terminal velocity for precipitation.
                 */
                vertvel -= fabs(terminal_velocity(is,ip,2*K,J,I))*RHO2(K,J,I)/H2(K,J,I);
              }
              w[K] = vertvel;
            }
          }
          else {
            for (K = k_first; K <= KHI; K++) {
              if (K == KLO) {
                /*
                 * Assume zero vertical velocity in top half of 
                 * topmost layer.
                 */
                vertvel = 0.;
              }
              else if (K == KHI) {
                /*
                 * Assume zero vertical velocity in bottom half of
                 * bottommost layer.
                 *
                 * NOTE: Doing otherwise leads to spurious advection
                 *       off the bottom of the model.
                 */
                vertvel = 0.;
              }
              else {
                vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));
              }
              w[K] = vertvel;
            }
          }
        }

        if (is == HDRY_INDEX || is == NU_TURB_INDEX) {
          /*
           *  Variable carried in the layer.
           *  Here, w_a[K] is located a half K step before w[K].
           */
          w_deltaa[k_first-1]= 0.;
          for (K = k_first; K < KHI; K++) {
            w_deltaa[K] = w_a[K]-w_a[K+1];
          }
          w_deltaa[KHI] = 0.;

          for (K = k_first; K <= KHI; K++) {
            tmp       = w_deltaa[K]*w_deltaa[K-1];
            w_dela[K] = (tmp > 0.) ? 2.*tmp/(w_deltaa[K]+w_deltaa[K-1]) : 0.;
          }

          if (k_first == 1) {
            w_ff[k_first-1] = mass_flux_top(planet,is,J,I,NULL);
          }
          else {
            w_ff[k_first-1] = 0.;
          }
          for (K = k_first; K < KHI; K++) {
            aval    = .5*(w_a[K]+w_a[K+1])-(w_dela[K]-w_dela[K+1])/6.;
            w_ff[K] = w[K]*aval;
          }
          w_ff[KHI] = mass_flux_bot(planet,is,J,I,NULL);
        }
        else {
          /*
           *  Variable carried on the bottom interface.
           *  Here, w_a[K] is located a half K step after w[K].
           */
          for (K = k_first; K <= KHI; K++) {
            w_deltaa[K] = w_a[K-1]-w_a[K];
          }

          w_dela[k_first-1] = 0.;
          for (K = k_first; K < KHI; K++) {
            tmp       = w_deltaa[K]*w_deltaa[K+1];
            w_dela[K] = (tmp > 0.) ? 2.*tmp/(w_deltaa[K]+w_deltaa[K+1]) : 0.;
          }
          w_dela[KHI] = 0.;

          for (K = k_first; K <= KHI; K++) {
            aval    = .5*(w_a[K-1]+w_a[K])-(w_dela[K-1]-w_dela[K])/6.;
            w_ff[K] = w[K]*aval;
          }
          w_ff[KHI+1] = 0.;
        }

        switch(is) {
          case HDRY_INDEX:
            for (K = k_first; K <= KHI; K++) {
              BUFF3D(K,J,I) += (w_ff[K]-w_ff[K-1])*grid.dsgth_inv[2*K]*DT;
            }
          break;
          case THETA_INDEX:
          case FPARA_INDEX:
          case P3_INDEX:
          case NU_TURB_INDEX:
          default:
            for (K = k_first; K <= KHI; K++) {
              BUFF3D(K,J,I) += (w_ff[K+1]-w_ff[K])*grid.dsgth_inv[2*K+1]*DT;
            }
          break;
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_HORIZONTAL) {

    /*************************
     * Meridional advection. *
     *************************/

    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;

      if (is == THETA_INDEX && K <= grid.k_sigma-2) continue;

      a = buff3d+(K-Kshift)*Nelem2d;

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);

      v2d    = Buff2D[0];
      deltaa = Buff2D[1];
      dela   = Buff2D[2];
      ff     = Buff2D[3];

      if (is == HDRY_INDEX) {
        for (J = JFIRST; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            V2D(J,I) = get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(v2d,TWODIM);
      }
      else {
        for (J = JLOPAD; J <= JHIPADPV; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            V2D(J,I) = V(grid.it_uv,K,J,I);
          }
        }
        /* No need to apply bc_lateral() here. */
      }

      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DELTAA(J,I) = A(J,I)-A(J-1,I);
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(deltaa,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          tmp       = DELTAA(J,I)*DELTAA(J+1,I);
          DELA(J,I) = (tmp > 0.) ? 2.*tmp/(DELTAA(J,I)+DELTAA(J+1,I)) : 0.;
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(dela,TWODIM);

      for (J = JFIRST; J <= JHI; J++) {
        if (fabs(grid.lat[2*J]) == 90.) {
          m_2j_inv = 0.;
        }
        else {
          m_2j_inv = 1./(grid.m)[2*J];
        }
        for (I = ILO; I <= IHI; I++) {
          aval    = .5*(A(J,I)+A(J-1,I))-(DELA(J,I)-DELA(J-1,I))/6.;
          FF(J,I) = V2D(J,I)*aval*m_2j_inv;
        }
      }
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          A(J,I) += (FF(J,I)-FF(J+1,I))*mn_2jp1*DT;
        }
      }
      bc_lateral(a,TWODIM);
    }

    /********************
     * Zonal advection. *
     ********************/

    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;

      if (is == THETA_INDEX && K <= grid.k_sigma-2) continue;

      /*
       * Apply zonal_filter() to a copy of variable.
       */
      a = Buff2D[6];
      memcpy(a,buff3d+(K-Kshift)*Nelem2d,Nelem2d*sizeof(EPIC_FLOAT));
      zonal_filter(is,a,NULL,TWODIM);

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);

      u2d    = Buff2D[0];
      deltaa = Buff2D[1];
      dela   = Buff2D[2];
      ff     = Buff2D[3];

      if (is == HDRY_INDEX) {
        /*
         * Variable is in the layer.
         */
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            U2D(J,I) = get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to call bc_lateral() here. */
        bc_lateral(u2d,TWODIM);
      }
      else {
        /*
         * Variable is on the interface.
         */
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            U2D(J,I) = U(grid.it_uv,K,J,I);
          }
        }
        /* No need to apply bc_lateral() here. */
      }

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          DELTAA(J,I) = A(J,I)-A(J,I-1);
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(deltaa,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          tmp       = DELTAA(J,I)*DELTAA(J,I+1);
          DELA(J,I) = (tmp > 0.) ? 2.*tmp/(DELTAA(J,I)+DELTAA(J,I+1)) : 0.;
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(dela,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        n_2jp1_inv = 1./(grid.n)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          aval    = .5*(A(J,I)+A(J,I-1))-(DELA(J,I)-DELA(J,I-1))/6.;
          FF(J,I) = U2D(J,I)*aval*n_2jp1_inv;
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = (grid.mn)[2*J+1];
        for (I = ILO; I <= IHI; I++) {
          BUFF3D(K,J,I) += (FF(J,I)-FF(J,I+1))*mn_2jp1*DT;
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(buff3d,THREEDIM);
  }

  return;
}

/*======================= end of akima_advection() =========================*/

/* * * * * * * * * * * * *  end of epic_flux.c  * * * * * * * * * * * * * * */
















