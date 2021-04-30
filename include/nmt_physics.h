/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2007-2008 David Raymond                           *
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


/* * * * * * * * * * * * * nmt_physics.h * * * * * * * * * * * * * * * * * * 
 *                                                                         *
 *                       David J. Raymond                                  *
 *                                                                         *
 *              Header file for the EPIC nmt_physics model                 *
 *                                                                         *
 *               For the latest version of the code email:                 *
 *                            mherman@nmt.edu                              *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



/* * * * * * * * * * * * * * * * * * * * * * * *
 *                                             *
 *            Defines and Includes             *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */

#include "epic_datatypes.h"

/*********************************************************************
 * set the types of the variables. this is here so that, for instance,
 * all floats can be set to doubles, etc. Also, we match the precision
 * of EPIC floats.
 *********************************************************************/
#if EPIC_PRECISION == 8
#  define TFLOAT double
#  define TDOUBLE double
#elif EPIC_PRECISION == 4
#  define TFLOAT float
#  define TDOUBLE double
#else
#  error Unrecognized value for EPIC_PRECISION environment variable.
#endif

/* these are needed for diabat3 */ 
#define TINT int
#define TLONG int


/********************************************
 * thermodynamic constants
 * 
 * NOTE: There is another constant 'EPSILON', 
 *       that is used in epic_microphysics.c
 ********************************************/
#define GEE                   9800.0
#define FREEZING              273.15
#define PREF                  100000.0
#define TREF                  300.0
#define KAPPA                 0.286
#define CP                    1005.0
#define RHOREF                1.2
#define LVAPOR                2.50e6
#define LFREEZ                3.33e5
#define ALFA                  ((LVAPOR + LFREEZ)/(CP*TREF))
#define EPSI_LON              0.622
#define RADITERMAX            5

/****************************************
 * define indices for the diabat3 outputs
 ****************************************/
#define STH                   0
#define SRT                   1
#define SU                    2
#define SV                    3
#define STHE_RAD              4
#define STHE_CU               5
#define IRFLUX                6
#define NUM_DIABAT3_OUTS      7

/***************************************************
 * define count of number of input arrays to diabat3
 ***************************************************/
#define NUM_DIABAT3_INS       15

/**************************************************
 * array shift macros
 **************************************************/
#define DRY_ENTROPY(k,j,i)        var.dry_entropy.value[       i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define MOIST_ENTROPY(k,j,i)      var.moist_entropy.value[     i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define SAT_MOIST_ENTROPY(k,j,i)  var.sat_moist_entropy.value[ i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define THE_FLUX(k,j,i)           var.the_flux.value[          i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define RT_FLUX(k,j,i)            var.rt_flux.value[           i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define U_FLUX(k,j,i)             var.u_flux.value[            i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define V_FLUX(k,j,i)             var.v_flux.value[            i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define CONVTHROT(k,j,i)          var.convthrot.value[         i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define FLUXTHROT(k,j,i)          var.fluxthrot.value[         i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define RAIN_RATE(k,j,i)          var.rain_rate.value[         i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/**************************************************
 * define a macro offset for diabat3 output arrays
 **************************************************/
#define DB_OUTS(var,k,j,i) nmt.outs[var][i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/**************************************
 * structure for holding diabat3 inputs
 **************************************/
typedef struct {
  int
    nz,
    land,
    dorad,
    frad;
  EPIC_FLOAT
    radcool,
    tpause,
    radbrk,
    cvc,
    cvs,
    cvp,
    cve,
    theslop,
    pstiff,
    pscale,
    pbltop,
    cdrag,
    wscale,
    sst,
    lfrac,
    eflux0,
    cld,
    cfract,
    dt,
   *zc,
   *xcell,
   *th,
   *rt,
   *pi,
   *u,
   *v,
   *xrho,
   *sth,
   *sthe_cu,
   *sthe_rad,
   *srt,
   *su,
   *sv,
   *irflux,
    rain,
    eflux,
    rflux,
    uflux,
    vflux,
    convthrot,
    fluxthrot,
    cvcaug,
   *outs[NUM_DIABAT3_OUTS];
} diabat3spec;

/***********************************
 * declare a diabat3spec for the run
 ***********************************/
diabat3spec
  nmt;

/*************************
 * function prototypes
 *************************/
void diabat3(TLONG nz, TLONG land, TLONG dorad, TLONG frad, TFLOAT radcool,
	     TFLOAT tpause, TFLOAT radbrk, TFLOAT cvc,
	     TFLOAT cvs, TFLOAT cvp, TFLOAT cve, TFLOAT theslop, TFLOAT pstiff,
	     TFLOAT pscale, TFLOAT pbltop, TFLOAT cdrag, TFLOAT wscale,
	     TFLOAT sst, TFLOAT lfrac, TFLOAT eflux0, TFLOAT cld,
	     TFLOAT cfract, TFLOAT dt, TFLOAT *zc, TFLOAT *xcell,
	     TFLOAT *th, TFLOAT *rt, TFLOAT *pi, TFLOAT *u, TFLOAT *v,
	     TFLOAT *xrho, TFLOAT *sth, TFLOAT *sthe_cu,
	     TFLOAT *sthe_rad, TFLOAT *srt, TFLOAT *su, TFLOAT *sv,
	     TFLOAT *irflux, TFLOAT *rain, TFLOAT *eflux, TFLOAT *rflux,
	     TFLOAT *uflux, TFLOAT *vflux, TFLOAT *convthrot,
	     TFLOAT *fluxthrot, TFLOAT *cvcaug);

TFLOAT rs(TFLOAT th, TFLOAT pi);
TFLOAT escalc(TFLOAT td);
TFLOAT thetacomp2(TFLOAT the, TFLOAT r, TFLOAT pi);
TFLOAT thetacomp(TFLOAT the, TFLOAT r, TFLOAT pi);
TFLOAT thetaecomp(TFLOAT th, TFLOAT r, TFLOAT pi);
TFLOAT thetaecomp2(TFLOAT th, TFLOAT r, TFLOAT rsat);
TFLOAT thetaescomp(TFLOAT th, TFLOAT pi);
TFLOAT gmoist(TFLOAT th, TFLOAT pi);
TFLOAT pifromp(TFLOAT p);
TFLOAT pfrompi(TFLOAT pi);
TFLOAT rhocomp(TFLOAT th, TFLOAT pi);
static double myrand(void);

void check(char *place, TFLOAT *array, int type, TLONG nelem);
void nmt_allocate();
void nmt_cumulus(planetspec *planet);
void nmt_apply_sources();
void nmt_memset();
void nmt_init_diag();
void nmt_init_water();
void nmt_randomize(int index);
void nmt_entropy();

/* not ready at this time */
int read_rt_vs_p(planetspec *planet,EPIC_FLOAT *pdat,EPIC_FLOAT *tdat,int portion);
