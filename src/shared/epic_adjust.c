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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 *  This file contains the following functions:                    *
 *      restore_mass()                                             *
 *      zonal_filter()                                             *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*=================== restore_mass() ========================================*/

/*
 * We currently just maintain a floor of grid.h_min[K] for H and 0.0
 * for mixing-ratio variables like H_2O_liquid.
 *
 * NOTE: Previously we distributed mass vertically, but we found that
 *       this tends to generate unwanted gravity waves.
 *       
 * NOTE: This function should not call set_p2etc() itself.
 */

void restore_mass(planetspec *planet,
                  int         species_index,
                  int         phase_index)
{
  register int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="restore_mass";
    
  if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
    /*
     * For mixing ratios, maintain a floor of 0.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          Q(species_index,phase_index,K,J,I) = MAX(0.,Q(species_index,phase_index,K,J,I));
        }
      }
    }
  }
  else if (species_index == HDRY_INDEX) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          HDRY(K,J,I) = MAX(grid.h_min[K],HDRY(K,J,I));
        }
      }
    }
  }
  else if (species_index == NU_TURB_INDEX) {
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          NU_TURB(K,J,I) = MAX(1.e-4*planet->kinvisc,NU_TURB(K,J,I));
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else {
    sprintf(Message,"unimplemented index=%d",species_index);
    epic_error(dbmsname,Message);
  }

  return;
}

/*=================== end of restore_mass() =================================*/

/*======================= zonal_filter() ====================================*/

void zonal_filter(int         index,
                  EPIC_FLOAT *buff3d,
                  EPIC_FLOAT *gz,
                  int         dim)
/* 
 * Filter to remove the CFL violation at the poles.
 * Use dim to indicate THREEDIM (KJI cube) or TWODIM (JI plane).
 * The TWODIM case should not subtract a function of GZ (which is assumed
 * to be carried on the P3 grid).
 *
 * NOTE: ni must be an integer power of 2, in order to use realft().
 *
 * NOTE: We used to have the zonal hyperviscosity here, but took
 *       it out because it was prone to Gibbs-effect ringing.
 */
{
  register int 
    K,J,I,
    kstart,kend,kay,
    jstart,
    Imin;
  EPIC_FLOAT
    geezee,gzmin,dgz;
  static int 
    initialized=FALSE,
    has_pole   =FALSE;
  static EPIC_FLOAT 
    *data, 
    *high_lat_h,
    *high_lat_pv,
    *varz;
  static float_triplet
    *table;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="zonal_filter";

  /* No need to filter if ni = 1. */
  if (grid.ni == 1) {
    return;
  }

  if(!initialized) {
    register int
      n,nn;
    EPIC_FLOAT 
      m0,n0,rln,rlt,lat0,
      re,rp,r,
      tmp0,tmp1,
      fac,
      al1 = 2.00,
      al2 = 2.00;
   
    data         = fvector(1,grid.ni,  dbmsname);
    high_lat_h   = fvector(0,Nelem2d-1,dbmsname);
    high_lat_pv  = fvector(0,Nelem2d-1,dbmsname);

    varz         = fvector( 0,Nelem3d-1,dbmsname);
    table        = ftriplet(0,grid.nk+1,dbmsname);

    /*
     * NOTE: The filter for the H-grid and the U-grid is the same.
     */

    if (strcmp(grid.geometry,"globe") == 0) {
      re   = planet->re;
      rp   = planet->rp;
      /* Filter is applied poleward of LAT0 */
      lat0 = LAT0*DEG;  
      rln  = re/sqrt( 1.+ pow(rp/re*tan(lat0),2.) );
      rlt  = rln/( cos(lat0)*( pow(sin(lat0),2.) +
                   pow(re/rp*cos(lat0),2.) ) );
      m0 = 1./(rln*grid.dln*DEG);
      n0 = 1./(rlt*grid.dlt*DEG);

      has_pole = TRUE;
    }
    else if (strcmp(grid.geometry,"f-plane")  == 0 &&
             strcmp(grid.f_plane_map,"polar") == 0) {
      /* nj = 2*(nj+1-1)/2 */
      m0 = grid.m[grid.nj];
      n0 = grid.n[grid.nj];

      has_pole = TRUE;
    }

    /*
     * h, u grid
     */
    for (J = JLO; J <= JHI; J++) {
      HIGH_LAT_H(J,1) = 1.;    
      if (has_pole) {   
        /* r is roughly 1.5 at a pole, 1 otherwise: */  
        r    = grid.mn[2*J+1]/(grid.m[2*J+1]*grid.n[2*J+1]);
        tmp0 = ((grid.n)[2*J+1]/n0)/((grid.m)[2*J+1]/m0);
        tmp0 = pow(tmp0,al1);
        for (I = 2; I <= grid.ni/2+1; I++) {
          tmp1  = sin((I-1)*(grid.dln)*DEG/2.);
          tmp1  = pow(tmp1,-al2);
          tmp1 *= r*tmp0;
          tmp1  = (tmp1 < 1.) ? tmp1 : 1.;
          HIGH_LAT_H(J,I) = tmp1;
        }
      }
      else {
        for (I = 2; I <= grid.ni/2+1; I++) {
          HIGH_LAT_H(J,I) = 1.;
        }
      }
    }

    /*
     * pv, v grid
     */
    for (J = JFIRST; J <= JHI; J++) {
      HIGH_LAT_PV(J,1) = 1.; 
      if (has_pole) {      
        /* r is roughly 1.5 at a pole, 1 otherwise: */  
        r    = grid.mn[2*J]/(grid.m[2*J]*grid.n[2*J]);
        tmp0 = ((grid.n)[2*J]/n0)/((grid.m)[2*J]/m0);
        tmp0 = pow(tmp0,al1);
        for (I = 2; I <= grid.ni/2+1; I++) {
          tmp1  = sin((I-1)*(grid.dln)*DEG/2.);
          tmp1  = pow(tmp1,-al2);
          tmp1 *= r*tmp0;
          tmp1  = (tmp1 < 1.) ? tmp1 : 1.;
          HIGH_LAT_PV(J,I) = tmp1;
        }
      }
      else {
        for (I = 2; I <= grid.ni/2+1; I++) {
          HIGH_LAT_PV(J,I) = 1.;
        }
      }
    }
    initialized = TRUE;
  } 
  /* End of initialization. */

  switch(index) {
    case V_INDEX:
    case PV3_INDEX:
      jstart = JFIRST;
    break;
    default:
      jstart = JLO;
    break;
  }

  if (dim == TWODIM) {
    kstart = Kshift;
    kend   = Kshift;  
  }
  else if (dim == THREEDIM) {
    /*
     * Subtract a representative profile vs GZ before applying filter.
     * The issue is that the zonal filter acts along hybrid-coordinate surfaces,
     * and hence tends to smear constant-with-z structure.
     */

    for (J = jstart; J <= JHI; J++) {
      /*
       * Find location of lowest GZ at this latitude.
       */
      K     = grid.nk;
      gzmin = FLOAT_MAX;
      for (I = ILO; I <= IHI; I++) {
        if (GZ(K,J,I) < gzmin) {
          Imin  = I;
          gzmin = GZ(K,J,I);
        }
      }
      /*
       * Set up spline for var(gz).
       */
      I = Imin;

      switch(index) {
        case U_INDEX:
          kstart = KLO;
          kend   = KHI;
          for (K = kend; K >= kstart; K--) {
            kay          = kend-K;
            geezee       = .5*(GZ(K,J,I)+GZ(K,J,I-1));
            table[kay].x = geezee;
            table[kay].y = BUFF3D(K,J,I);
          }
          spline_pchip(kend-kstart+1,table);

          for (I = ILO; I <= IHI; I++) {
            for (K = kstart; K <= kend; K++) {
              geezee         = .5*(GZ(K,J,I)+GZ(K,J,I-1));
              kay            = find_place_in_table(kend-kstart+1,table,geezee,&dgz);
              VARZ(K,J,I)    = splint_pchip(geezee,table+kay,dgz);
              BUFF3D(K,J,I) -= VARZ(K,J,I);
            }
          }
        break;
        case V_INDEX:
          kstart = KLO;
          kend   = KHI;
          for (K = kend; K >= kstart; K--) {
            kay          = kend-K;
            geezee       = .5*(GZ(K,J,I)+GZ(K,J-1,I));
            table[kay].x = geezee;
            table[kay].y = BUFF3D(K,J,I);
          }
          spline_pchip(kend-kstart+1,table);

          for (I = ILO; I <= IHI; I++) {
            for (K = kstart; K <= kend; K++) {
              geezee         = .5*(GZ(K,J,I)+GZ(K,J-1,I));
              kay            = find_place_in_table(kend-kstart+1,table,geezee,&dgz);
              VARZ(K,J,I)    = splint_pchip(geezee,table+kay,dgz);
              BUFF3D(K,J,I) -= VARZ(K,J,I);
            }
          }
        break;
        case HDRY_INDEX:
        case H2_INDEX:
          /*
           * Variable is defined in the layer.
           * Interpolate on log of variable.
           */
          kstart = KLO;
          kend   = KHI;
          for (K = kend; K >= kstart; K--) {
            kay          = kend-K;
            geezee       = .5*(GZ(K,J,I)+GZ(K-1,J,I));
            table[kay].x = geezee;
            table[kay].y = log(BUFF3D(K,J,I));
          }
          spline_pchip(kend-kstart+1,table);

          for (I = ILO; I <= IHI; I++) {
            for (K = kstart; K <= kend; K++) {
              geezee         = .5*(GZ(K,J,I)+GZ(K-1,J,I));
              kay            = find_place_in_table(kend-kstart+1,table,geezee,&dgz);
              VARZ(K,J,I)    = exp(splint_pchip(geezee,table+kay,dgz));
              BUFF3D(K,J,I) -= VARZ(K,J,I);
            }
          }
        break;
        case THETA_INDEX:
          /*
           * Variable is defined on the layer interfaces, and is only prognostic in
           * the sigma-coordinate region.
           */
          kstart = grid.k_sigma-1;
          kend   = KHI;

          for (K = kend; K >= kstart; K--) {
            kay          = kend-K;
            geezee       = GZ(K,J,I);
            table[kay].x = geezee;
            table[kay].y = BUFF3D(K,J,I);
          }
          spline_pchip(kend-kstart+1,table);

          for (I = ILO; I <= IHI; I++) {
            for (K = kstart; K <= kend; K++) {
              geezee         = GZ(K,J,I);
              kay            = find_place_in_table(kend-kstart+1,table,geezee,&dgz);
              VARZ(K,J,I)    = splint_pchip(geezee,table+kay,dgz);
              BUFF3D(K,J,I) -= VARZ(K,J,I);
            }
          }
        break;
        default:
          kstart = KLO-1;
          kend   = KHI;
          for (K = kend; K >= kstart; K--) {
            kay          = kend-K;
            geezee       = GZ(K,J,I);
            table[kay].x = geezee;
            table[kay].y = BUFF3D(K,J,I);
          }
          spline_pchip(kend-kstart+1,table);

          for (I = ILO; I <= IHI; I++) {
            for (K = kstart; K <= kend; K++) {
              geezee         = GZ(K,J,I);
              kay            = find_place_in_table(kend-kstart+1,table,geezee,&dgz);
              VARZ(K,J,I)    = splint_pchip(geezee,table+kay,dgz);
              BUFF3D(K,J,I) -= VARZ(K,J,I);
            }
          }
        break;
      }
    }
  }
  else {
    sprintf(Message,"dim=%d not recognized");
    epic_error(dbmsname,Message);
  }

  for (K = kstart; K <= kend; K++) {
    switch(index) {
      case PV3_INDEX:
      case V_INDEX:
        /*
         * Variable defined on whole-integer J grid (pv-grid or v-grid).
         */
        for (J = jstart; J <= JHI; J++) {
          if (fabs(grid.lat[2*J]) < LAT0) {
            /*
             * Only apply for lat >= |LAT0|.
             */
            continue;
          }
          for (I = 1; I <= grid.ni; I++) {
            data[I] = BUFF3D(K,J,I);
          }
          realft(data,grid.ni,1);

          data[1] *= HIGH_LAT_PV(J,1);
          for (I = 2; I <= (grid.ni)/2; I++) {
            data[2*I-1] *= HIGH_LAT_PV(J,I);
            data[2*I  ] *= HIGH_LAT_PV(J,I);
          }
          data[2] *= HIGH_LAT_PV(J,(grid.ni)/2+1);

          realft(data,grid.ni,-1);
          for (I = 1; I <= grid.ni; I++) {
            BUFF3D(K,J,I) = data[I]*(2./(grid.ni));
          }
        }
      break;
      default:
        for (J = jstart; J <= JHI; J++) {
          if (fabs(grid.lat[2*J+1]) < LAT0) {
            /*
             * Only apply for lat >= |LAT0|.
             */
            continue;
          }
          for (I = 1; I <= grid.ni; I++) {
            data[I] = BUFF3D(K,J,I);
          }
          realft(data,grid.ni,1);

          data[1] *= HIGH_LAT_H(J,1);
          for (I = 2; I <= (grid.ni)/2; I++) {
            data[2*I-1] *= HIGH_LAT_H(J,I);
            data[2*I  ] *= HIGH_LAT_H(J,I);
          }
          data[2] *= HIGH_LAT_H(J,(grid.ni)/2+1);

          realft(data,grid.ni,-1);
          for (I = 1; I <= grid.ni; I++) {
            BUFF3D(K,J,I) = data[I]*(2./(grid.ni));
          }
        }
      break;
    }
  }

  if (dim == THREEDIM) {
    /*
     * Add VARZ back in.
     */
    for (K = kstart; K <= kend; K++) {
      for (J = jstart; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          BUFF3D(K,J,I) += VARZ(K,J,I);
        }
      }
    }
  }

  /* Need to call bc_lateral() here. */
  bc_lateral(buff3d,dim);

  return;
}

/*======================= end of zonal_filter() ==========================*/

/* * * * * * * * * * * * end of epic_adjust.c * * * * * * * * * * * * * * */













