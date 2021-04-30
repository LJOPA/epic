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

/*******************************************************************
 * This is a concatenation of what had been three different files.
 * Any changes to those files should be placed in the proper slot
 * in this file. Later, this file will be reworked to follow the
 * style of the EPIC code. Please move all #includes to the top.
 *******************************************************************/

/* * * * * * * * * * * * * epic_cumulus.c  * * * * * * * * * * * * *
 *                         v.2.0                                   *        
 *                                                                 *
 *         For the latest version of the code email:               *
 *                      mherman@nmt.edu                            *
 *                                                                 *
 *                                                                 *
 *       Functions effecting the cumulus parameterization.         *
 *       This file includes the following functions:               *
 *                                                                 *
 *           from diabat3.c                                        *
 *           ***************                                       *
 *           floatbuff()                                           *
 *           fbuf()                                                *
 *           surfextrap()                                          *
 *           pblfld()                                              *
 *           thresh()                                              *
 *           wt()                                                  *
 *           throtcomp()                                           *
 *           sources()                                             *
 *           cuparam()                                             *
 *           elimit()                                              *
 *           zedgecomp()                                           *
 *           lcomp()                                               *
 *           kcomp()                                               *
 *           taucomp()                                             *
 *           radiation()                                           *
 *           fradiation()                                          *
 *           diabat3()                                             *
 *                                                                 *
 *           from thermo.c                                         *
 *           **************                                        *
 *           thetasat()                                            *
 *           rs()                                                  *
 *           escalc()                                              *
 *           thetacomp2()                                          *
 *           thetacomp()                                           *
 *           thetaecomp()                                          *
 *           thetaecomp2()                                         *
 *           thetaescomp()                                         *
 *           gmoist()                                              *
 *           pifromp()                                             *
 *           pfrompi()                                             *
 *           rhocomp()                                             *
 *                                                                 *
 *           from wrapper.c                                        *
 *           ***************                                       *
 *           wrapper()                                             *
 *           check()                                               *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*****************************
 * These are all the includes
 * for all three files.
 *****************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "epic.h"


/*******************************************************************
 * diabat3.c
 *******************************************************************/

/* diabat3.c -- Version 3 of Raymond convective/radiative package.
 *
 * Note: This subroutine provides a parameterization for **all**
 * diabatic effects, including latent heat release in explicit
 * vertical motion and in hard adjustment due to dry static
 * instability.  Don't include any separate explicit diabatic
 * effects or you will be double counting them!
 *
 * The model accepts input fields zc thru w in a single column at cell-
 * centered levels defined by zc.  The output profiles sth thru sv
 * have the same structure.  The values of dz for each cell
 * should be consistent with zc.  The levels must be ordered from
 * lowest to highest in height.  Tendencies for potential temperature,
 * mixing ratio, and the two horizontal wind components are returned
 * at these levels.  Note that the unit ``ks'' is kiloseconds.
 *
 * $Id: diabat3.c,v 1.24 2007/07/21 10:47:48 raymond Exp $
 */

/* floatbuff -- This routine uses malloc to allocate float buffer
 * space.  The argument nelem defines the number of floating
 * point elements.  It returns a float pointer to the space.  If malloc
 * fails, getbuff dies with an error message.
 */
static TFLOAT *floatbuff(TLONG nelem)
{
  TFLOAT *point;
  if ((point = (TFLOAT*)malloc((unsigned)(sizeof(TFLOAT)*nelem)))
       == (TFLOAT*)NULL) {
    fprintf(stderr,"floatbuff: can't get %d words\n",nelem);
    exit(1);
  }
  return(point);
}

/* fbuf -- Free a buffer obtained with floatbuff */
static void fbuf(TFLOAT *ptr)
{
  free((char*)ptr);
}

/* surfextrap -- Extrapolate a field to the surface */
static TFLOAT surfextrap(TFLOAT *fld, TFLOAT *dz)
{
  return(fld[0] + (fld[0] - fld[1])*dz[0]/(dz[1] + dz[0]));
}



/* cuparam -- The convective parameterization */

/* constants for the cumulus parameterization */
#define FTHROTSLOP 0.5
#define WTRAIN 5.
#define WTSNOW 5.
#define MAXFRAC 0.3
#define ALPHA1 0.
#define ALPHA2 -1.
#define DTMAX 1.0

/* pblfld -- PBL mean value of field */
static TFLOAT pblfld(nz,pbltop,zc,fld)
TLONG nz;                    /* number of levels */
TFLOAT pbltop;               /* height of top of PBL */
TFLOAT *zc;                  /* grid centers */
TFLOAT *fld;                 /* some field */
{
  TLONG lev1, nlevs;
  TFLOAT fldval;

  nlevs = 0;
  fldval = 0.;
  while (1) {
    fldval += fld[nlevs];
    nlevs++;
    if (nlevs >= nz - 2 || zc[nlevs] > pbltop) break;
  }
  fldval /= (TFLOAT)nlevs;
  return(fldval);
}

/* thresh -- compute threshold value of theta-e for convection --
 * this is equal to the saturated equivalent potential temperature
 * between one and two PBL depths
 */
static TFLOAT thresh(TLONG nz, TFLOAT pbltop, TFLOAT *zc, TFLOAT *th,
		     TFLOAT *rt, TFLOAT *pi)
{
  TFLOAT thes;
  TINT lev1, nlevs;

  /* find starting level for averaging threshold saturated theta-e -- 
   * must be above PBL and lifted parcel must be saturated
   */
  lev1 = 1;
  while (lev1 < nz - 2 && (zc[lev1] < pbltop ||
			   rt[0]/rs(th[0],pi[lev1]) < 1.)) {
    lev1++;
  }

  /* average saturated theta-e to twice PBL thickness */
  nlevs = 0;
  thes = 0.;
  while (1) {
    nlevs++;
    thes += thetaescomp(th[lev1],pi[lev1]);
    lev1++;
    if (lev1 >= nz - 2 || zc[lev1] > 2.*pbltop) break;
  }
  thes /= (TFLOAT)nlevs;
  return(thes);
}

/* wt -- compute particle terminal velocity */
static TFLOAT wt(TFLOAT th, TFLOAT pi)
{
  TFLOAT temp;
  temp = th*pi/CP;
  return(temp > FREEZING ? WTRAIN : WTSNOW);
}

/* throtcomp -- evaluate the throttle function */
static TFLOAT throtcomp(TFLOAT value, TFLOAT slop)
{
  TFLOAT tvalue;
  tvalue = 0.5 + value/slop;
  if (tvalue < 0.) tvalue = 0.;
  if (tvalue > 1.) tvalue = 1.;
  return(tvalue);
}

/* buoy -- compute buoyancy of a surface parcel at an arbitrary level */
static TFLOAT buoy(TLONG level, TFLOAT thesurf, TFLOAT rtsurf,
		   TFLOAT *pi, TFLOAT *th)
{
  return(thetacomp2(thesurf,rtsurf,pi[level]) - th[level]);
}

/* sources -- do the real work of the convective parameterization */

static void sources(TLONG nz, TLONG mode, TFLOAT dt, TFLOAT alpha,
		    TFLOAT cvc, TFLOAT cvs, TFLOAT cvp, TFLOAT cve,
		    TFLOAT pstiff, TFLOAT pscale,
		    TFLOAT pbltop, TFLOAT *zc, TFLOAT *dz, TFLOAT *pi,
		    TFLOAT *th, TFLOAT *the, TFLOAT *rt, TFLOAT *u,
		    TFLOAT *v, TFLOAT *rho, TFLOAT *rsat, TFLOAT eflux,
		    TFLOAT rflux, TFLOAT uflux, TFLOAT vflux,
		    TFLOAT *sthe, TFLOAT *srt, TFLOAT *su,
		    TFLOAT *sv, TFLOAT *rain, TFLOAT *cvcaug)
{
  TLONG iz,jz;               /* z indices */
  TLONG top;                 /* index just above cloud top */
  TFLOAT stratx,penex,evapx,totx;  /* for precipitation computation */
  TFLOAT precip,satvalue,rp; /* for precipitation computation */
  TFLOAT cvcactual,pwat,satpwat,sfrac,zscaled;
  TFLOAT *pshape;
  TFLOAT stab,dzri,dzrj,denom;   /* for hard adjustment */
  TFLOAT validfract;         /* fraction of cell occupied by cloud */
  TFLOAT thegrad;            /* vertical gradient of the */
  TFLOAT emean,rmean,umean,vmean;  /* stuff for mean value comps */
  TFLOAT ztfract,buoy1,buoy2; /* for fractional high top cell calc */
  TFLOAT relhum;              /* for convective precip calc */
  TFLOAT curateint;           /* mass-weighted integral of curate(z) */
  TFLOAT mass;                /* mass per unit area in a layer */
  TFLOAT *curate;             /* weighting profile for mixing */
  TFLOAT zscale;              /* scale height of convection */
  TFLOAT thepbl,rtpbl;        /* pbl values of fields */

  /* allocate space */
  curate = floatbuff(nz);
  pshape = floatbuff(nz);

  /* compute the top of the convective layer as the level above the first
   * positively buoyant layer marching from the top down -- the top layer
   * used in all the computations is actually top - 1, i. e., the last
   * positively buoyant layer
   */
  thepbl = pblfld(nz,pbltop,zc,the);
  rtpbl = pblfld(nz,pbltop,zc,rt);
  top = nz - 2;
  while (top > 0 && buoy(top,thepbl,rtpbl,pi,th) < 0.) {
    top--;
  }

  /* note that top is always > 0 */
  top++;

  /* compute the fractional position of neutral level in high top layer */
  buoy1 = buoy(top,the[0],rt[0],pi,th);
  buoy2 = buoy(top - 1,the[0],rt[0],pi,th);
  if (top == 1) {
    ztfract = 1.;
  }
  else if (buoy2 - buoy1 < 1.e-10) {
    ztfract = 0.5;
  }
  else {
    ztfract = buoy2/(buoy2 - buoy1);
  }
  if (ztfract > 1.) ztfract = 1.;
  if (ztfract < 0.) ztfract = 0.;

  if (mode == 1) {  /* deep mode */
    zscale = zc[top - 1] + ztfract*(dz[top - 1] + dz[top])/2.;
    if (zscale < pbltop) zscale = pbltop;
  }
  else {  /* pbl mode */
    zscale = pbltop;
  }

  /* get some vertical profiles */
  pwat = satpwat = 0.;
  for (iz = 0; iz < nz; iz++) {

    /* compute valid fraction of cell to include */
    if (iz >= top) {
      validfract = 0.;
    }
    else if (iz == top - 1) {
      validfract = ztfract;
    }
    else {
      validfract = 1.;
    }

    /* compute saturation fraction and augmented mixing parameter */
    mass = dz[iz]*rho[iz];
    pwat += mass*rt[iz];
    satpwat += mass*rsat[iz];

    /* compute height dependence of convective precipitation */
    zscaled = zc[iz]/pscale;
    pshape[iz] = (1. - exp(-zscaled*zscaled))*validfract;
  }
  sfrac = pwat/satpwat;
  if (sfrac < 0.) sfrac = 0.;
  cvcactual = cvc*pow(sfrac,pstiff);
  if (cvcactual*DTMAX > MAXFRAC) cvcactual = MAXFRAC/DTMAX;

  /* compute means */
  curateint = emean = rmean = umean = vmean = 0.;
  for (iz = 0; iz < nz; iz++) {

    /* compute valid fraction of cell to include */
    if (iz >= top) {
      validfract = 0.;
    }
    else if (iz == top - 1) {
      validfract = ztfract;
    }
    else {
      validfract = 1.;
    }

    /* compute height dependent convective rate */
    curate[iz] = zc[iz] < zscale ? cvcactual*
      (1. + alpha*zc[iz]/zscale)*validfract : 0.;

    /* compute the means */
    mass = dz[iz]*rho[iz];
    curateint += curate[iz]*mass;
    emean += curate[iz]*the[iz]*mass;
    rmean += curate[iz]*rt[iz]*mass;
    umean += curate[iz]*u[iz]*mass;
    vmean += curate[iz]*v[iz]*mass;
  }
  emean /= curateint;
  rmean /= curateint;
  umean /= curateint;
  vmean /= curateint;

  /* compute precip and all source terms -- march down from highest level */
  iz = nz - 1;
  precip = 0.;
  while (iz >= 0) {

    /* compute valid fraction of cell to include */
    if (iz >= top) {
      validfract = 0.;
    }
    else if (iz == top - 1) {
      validfract = ztfract;
    }
    else {
      validfract = 1.;
    }

    /* stratiform rain production -- can't occur if dthe/dz < 0 */
    /*if (iz < 2){
      fprintf(stderr,"(iz=%d) rt=%.10f, rsat=%.10f\n",iz,rt[iz],rsat[iz]);
      }*/        
    
    satvalue = rt[iz] - rsat[iz];
    if (iz >= top - 2) {
      thegrad = 10.;
    }
    else {
      thegrad = (the[iz + 1] - the[iz])/(zc[iz + 1] - zc[iz]);
    }
    stratx = satvalue > 0. ? satvalue*cvs*throtcomp(thegrad - 2.,2.) : 0.;

    /* penetrative rain production */
    relhum = rt[iz]/rsat[iz];
    if (relhum > 1.) relhum = 1.;
    if (relhum < 0.) relhum = 0.;
    penex = cvp*pow(relhum,pstiff)*pshape[iz];
    if (pstiff*penex*DTMAX > MAXFRAC*rt[iz]) {
      penex = MAXFRAC*rt[iz]/(pstiff*DTMAX);
    }

    /* evaporation of rain */
    rp = precip/(rho[iz]*wt(th[iz],pi[iz]));
    evapx = satvalue < 0. ? satvalue*rp*cve : 0.;
    /*if (iz < 2){
      fprintf(stderr,"(iz=%d) evapx=%.10f, satvalue=%.10f, rp=%.10f, cve=%.10f, precip=%.10f, rho[iz]=%.10f\n",iz,evapx,satvalue,rp,cve,precip,rho[iz]);
      }*/


    /* total precipitation effects advanced to current level */
    totx = stratx + penex + evapx;

    /* compute the precipitation, taking precaution against negative value */
    precip += totx*rho[iz]*dz[iz];
    if (precip < 0.) {
      totx -= precip/(rho[iz]*dz[iz]);
      precip = 0.;
    }
    /*if (iz < 2){
      fprintf(stderr,"(iz=%d) totx=%.10f, stratx=%.10f, penex=%.10f, evapx=%.10f, rho[iz]=%.10f, dz[iz]=%.10f\n",iz,totx,stratx,penex,evapx,rho[iz],dz[iz]);
      }*/

    
    /* source terms */
    /*if (iz < 2){
      fprintf(stderr,"rmean=%.16f, rt[%d]=%.16f, rflux=%.16f, curate[%d]=%.16f, totx=%.16f\n",rmean,iz,rt[iz],rflux,iz,curate[iz],totx);
      }*/
    sthe[iz] = (emean - the[iz] + eflux/curateint)*curate[iz];
    srt[iz] = (rmean - rt[iz] + rflux/curateint)*curate[iz] - totx;
    su[iz] = (umean - u[iz] + uflux/curateint)*curate[iz];
    sv[iz] = (vmean - v[iz] + vflux/curateint)*curate[iz];
    
    /* decrement counter */
    iz--;

  }


  /* return the rain, which is the downward integrated, mass weighted totx */
  *rain = precip;
  /* return the augmented convective mixing rate */
  *cvcaug = cvcactual;

  /* hard adjustment if dry statically unstable -- this shouldn't
   * happen, so the adjustment code is commented out */
#if 0
  for (iz = nz - 2; iz >= 0; iz--) {
    jz = iz + 1;
    stab = (th[jz] - th[iz])/(zc[jz] - zc[iz]);

    /* if unstable compute the target values of fields */
    if (stab < 0. && dt > 0.) {
      dzri = dz[iz]*rho[iz];
      dzrj = dz[jz]*rho[jz];
      denom = dzri + dzrj;
      rmean = (rt[jz]*dzrj + rt[iz]*dzri)/denom;
      emean = (the[jz]*dzrj + the[iz]*dzri)/denom;
      umean = (u[jz]*dzrj + u[iz]*dzri)/denom;
      vmean = (v[jz]*dzrj + v[iz]*dzri)/denom;

      /* increment source terms to do adjustment */
      srt[jz] += (rmean - rt[jz])/dt;
      srt[iz] += (rmean - rt[iz])/dt;
      sthe[jz] += (emean - the[jz])/dt;
      sthe[iz] += (emean - the[iz])/dt;
      su[jz] += (umean - u[jz])/dt;
      su[iz] += (umean - u[iz])/dt;
      sv[jz] += (vmean - v[jz])/dt;
      sv[iz] += (vmean - v[iz])/dt;
    }
  }
#endif

  /* free space */
  fbuf(pshape);
  fbuf(curate);
}

/* the main routine for the cumulus parameterization */
static void cuparam(TLONG nz, TLONG land, TFLOAT cvc,
		    TFLOAT cvs, TFLOAT cvp, TFLOAT cve, TFLOAT theslop,
		    TFLOAT pstiff, TFLOAT pscale, TFLOAT pbltop,
		    TFLOAT cdrag, TFLOAT wscale, TFLOAT sst, TFLOAT lfrac,
		    TFLOAT eflux0, TFLOAT dt, TFLOAT *zc, TFLOAT *dz,
		    TFLOAT *th, TFLOAT *the, TFLOAT *rt, TFLOAT *pi,
		    TFLOAT *u, TFLOAT *v, TFLOAT *rho, TFLOAT *rsat,
		    TFLOAT *sthe, TFLOAT *srt, TFLOAT *su, TFLOAT *sv,
		    TFLOAT *rain, TFLOAT *eflux, TFLOAT *rflux,
		    TFLOAT *uflux, TFLOAT *vflux, TFLOAT *convthrot,
		    TFLOAT *fluxthrot, TFLOAT *cvcaug)
{
  TLONG iz;                   /* looping index */
  TFLOAT ewind;               /* effective wind for fluxes */
  TFLOAT surfpi;              /* extrapolated Exner function at surface */
  TFLOAT surftemp;            /* extrapolated temperature at surface */
  TFLOAT surfthe;             /* extrapolated thetae at surface */
  TFLOAT surfrt;              /* extrapolated mixing ratio at surface */
  TFLOAT surfth;              /* extrapolated potential temperature at surf */
  TFLOAT surfu,surfv;         /* extrapolated wind at surface */
  TFLOAT sstheta;             /* potential temperature of sea surface */
  TFLOAT ssthe;               /* equivalent potential temperature of surface */
  TFLOAT ssrt;                /* mixing ratio of surface */
  TFLOAT fthrottle;           /* throttle on fluxes due to sea-air temp diff */
  TFLOAT *sthe1,*srt1,*su1,*sv1,rain1;  /* deep forcing terms */
  TFLOAT *sthe2,*srt2,*su2,*sv2,rain2;  /* pbl forcing terms */
  TFLOAT cvcaug1,cvcaug2;     /* augmented mixing rate */
  TFLOAT throttle,cthrottle;  /* deep convective throttle and complement */
  TFLOAT srtcorr,precipcorr;  /* correct for negative rt generation */
  TFLOAT srtint;

  /* check size of dt */
  if (dt > DTMAX) {
    fprintf(stderr,"diabat3: Sorry, must have dt < %f\n",DTMAX);
    exit(1);
  }

  /* allocate space */
  sthe1 = floatbuff(nz);
  srt1 = floatbuff(nz);
  su1 = floatbuff(nz);
  sv1 = floatbuff(nz);
  sthe2 = floatbuff(nz);
  srt2 = floatbuff(nz);
  su2 = floatbuff(nz);
  sv2 = floatbuff(nz);

  /* compute surface fluxes -- currently uses same drag coeff for
   * thermodynamics and wind
   */

  /* extrapolate things to the surface */
  /*fprintf(stderr,"dz[1]=%.15f, dz[0]=%.15f\n",dz[1],dz[0]);*/
  surfpi = surfextrap(pi,dz);
  surfth = surfextrap(th,dz);
  surfthe = surfextrap(the,dz);
  surfrt = surfextrap(rt,dz);
  surfu = surfextrap(u,dz);
  surfv = surfextrap(v,dz);
  surftemp = surfth*surfpi/CP;
  sstheta = sst*CP/surfpi;       /* sea surface theta */
  /*fprintf(stderr,"sst=%.15f, sstheta=%.15f, surfpi=%.15f\n",sst,sstheta,surfpi);*/
  ssrt = rs(sstheta,surfpi);     /* surface mixing ratio */
  ssthe = thetaecomp(sstheta,ssrt,surfpi);  /* surface thetae */
  /*fprintf(stderr,"ssthe=%.15f, sstheta=%.15f, ssrt=%.15f, surfpi=%.15f\n",ssthe,sstheta,ssrt,surfpi);*/

  /* this formula for the effective wind keeps the fluxes from going to
   * zero at zero real wind
   */
  ewind = pow(surfu*surfu + surfv*surfv + wscale*wscale,0.5);

  /* choice in fluxes of land versus ocean */
  if (land) {

    /* throttle is set to unity over land */
    fthrottle = 1.;

    /* land fluxes are specified */
    *eflux = eflux0;
    /*fprintf(stderr,"lfrac=%.15f, eflux0=%.15f, surfthe=%.15f\n",lfrac,eflux0,surfthe);*/
    *rflux = lfrac*eflux0/(ALFA*surfthe);
  }
  else {

    /* this throttle function turns off surface fluxes if the air temperature
     * at the surface is much in excess of the sea surface temperature
     */
    fthrottle = throtcomp(sst - surftemp,FTHROTSLOP);

    /* ocean fluxes are done via a bulk formula */

    /*fprintf(stderr,"fthrottle=%.15f, rho[0]=%.15f, ewind=%.15f, cdrag=%.15f, ssthe=%.15f, surfthe=%.15f\n",fthrottle,rho[0],ewind,cdrag,ssthe,surfthe);*/
    *eflux = fthrottle*rho[0]*ewind*cdrag*(ssthe - surfthe);
    /*fprintf(stderr,"in cuparam after assignment)(: eflux = %.15f\n",*eflux);*/

    *rflux = fthrottle*rho[0]*ewind*cdrag*(ssrt - surfrt);

    /* kill negative thermodynamic fluxes */
    if (*eflux < 0.) *eflux = 0.;
    if (*rflux < 0.) *rflux = 0.;
  }

  /* momentum fluxes */
  *uflux = -fthrottle*rho[0]*ewind*cdrag*surfu;
  *vflux = -fthrottle*rho[0]*ewind*cdrag*surfv;

  /* compute source terms for convection -- mode 1 --> deep, mode2 --> pbl */
  sources(nz,1,dt,ALPHA1,cvc,cvs,cvp,cve,pstiff,pscale,pbltop,zc,dz,
	  pi,th,the,rt,u,v,rho,rsat,*eflux,*rflux,*uflux,*vflux,
	  sthe1,srt1,su1,sv1,&rain1,&cvcaug1);
  sources(nz,2,dt,ALPHA2,5.*cvc,cvs,0.,cve,pstiff,pscale,pbltop,zc,dz,
	  pi,th,the,rt,u,v,rho,rsat,*eflux,*rflux,*uflux,*vflux,
	  sthe2,srt2,su2,sv2,&rain2,&cvcaug2);

  /* compute convective threshold throttle function and its complement */
  throttle = throtcomp(pblfld(nz,pbltop,zc,the) -
		       thresh(nz,pbltop,zc,th,rt,pi),theslop);
  cthrottle = 1. - throttle;
  


  /* mix deep and pbl source terms */
  srtint = 0.;
  precipcorr = 0.;
  for (iz = 0; iz < nz; iz++) {
    sthe[iz] = throttle*sthe1[iz] + cthrottle*sthe2[iz];
    srt[iz] = throttle*srt1[iz] + cthrottle*srt2[iz];
    su[iz] = throttle*su1[iz] + cthrottle*su2[iz];
    sv[iz] = throttle*sv1[iz] + cthrottle*sv2[iz];

    /* keep srt from producing negative total water */
    if (dt > 0. && rt[iz] + srt[iz]*dt < 0.) {
      srtcorr = -srt[iz] - rt[iz]/dt;
      srt[iz] += srtcorr;
      precipcorr -= srtcorr*rho[iz]*dz[iz];
    }
    /*if (iz < 2){
      fprintf(stderr,"(iz=%d): throttle=%.10f, srt1=%.10f, cthrottle=%.10f, srt2=%.10f, rt=%.10f, srt after checking:%.10f\n",iz,throttle,srt1[iz],cthrottle,srt2[iz],rt[iz],srt[iz]);
      }*/
    srtint += srt[iz]*rho[iz]*dz[iz];
  }
  *rain = throttle*rain1 + cthrottle*rain2 + precipcorr;
  *fluxthrot = fthrottle;
  *convthrot = throttle;
  *cvcaug = cvcaug1;   /* the deep mixing rate is our concern */

  /* free space */
  fbuf(sv2);
  fbuf(su2);
  fbuf(srt2);
  fbuf(sthe2);
  fbuf(sv1);
  fbuf(su1);
  fbuf(srt1);
  fbuf(sthe1);
}


/* radiation -- Parameterization of long and shortwave radiation
 * this model is tuned by making rough comparisons with NCAR's
 * CCM2 radiation model (constructed by Jeff Kiehl and associates)
 */

/* one and only one of the following radiation options must be
   turned on (on = 1, off = 0) -- 1 channel radiation is faster, but is
   just a gray body scheme -- 8 channel radiation is
   the radiation scheme used in the
   Raymond-Torres (1998) paper -- J. Atmos. Sci., 55, 1771-1790 --
   and includes a crude representation of radiation-water vapor
   interactions */
#define ONECHANNEL 0
#define MULTICHANNEL 1

/* constants for the radiation parameterization */
#if ONECHANNEL
#define NBANDS 1
#endif

#if MULTICHANNEL
#define NBANDS 8
#endif

#define CO2 0
#define CONT 1
#define RTREF 0.02
#define ZSCALE 8.
#define SBCONST 5.67e-8

/* elimit -- limit the value to +- 75 to avoid blowing up exponential fn */
static TDOUBLE elimit(val)
TFLOAT val;
{
  TFLOAT lval;
  lval = val;
  if (lval > 75.) lval = 75.;
  if (lval < -75.) lval = -75.;
  return(lval);
}

/* zedgecomp -- extrapolate cell edge array from cell centered array */
static void zedgecomp(TLONG nz, TFLOAT *dz, TFLOAT *vcell,
		      TFLOAT *vedge)
{
  TLONG iz,jz;
  for (iz = 1; iz < nz; iz++) {
    jz = iz - 1;
    vedge[iz] = (vcell[iz]*dz[jz] + vcell[jz]*dz[iz])/(dz[iz] + dz[jz]);
  }
  vedge[0] = (vcell[0]*(dz[1] + 2.*dz[0]) - vcell[1]*dz[0])/(dz[1] + dz[0]);
  vedge[nz] = (vcell[nz - 1]*(dz[nz - 2] + 2.*dz[nz - 1]) -
	       vcell[nz - 2]*dz[nz - 1])/(dz[nz - 2] + dz[nz - 1]);
}

/* lcomp -- compute the l factor from delta tau and the b values
 * this assumes that b = sigma*T^4 is linear across the grid cell
 * -- integrate analytically and watch out for small dtau
 */
static TFLOAT lcomp(TFLOAT dtau, TFLOAT b1, TFLOAT b2)
{
  TFLOAT efact,omefact,lfactor;
  efact = exp(elimit(-dtau));
  omefact = 1. - efact;
  if (dtau < 0.001) {
    lfactor = 0.5*dtau*(b1 + b2);
  }
  else {
    lfactor = b1*(omefact/dtau - efact) + b2*(1. - omefact/dtau);
  }
  return(lfactor);
}

/* kcomp -- compute the kplus and kminus profiles */
static void kcomp(TLONG nz, TFLOAT *dz, TFLOAT *tau, TFLOAT *th,
		   TFLOAT *pi, TFLOAT *kplus, TFLOAT *kminus)
{
  TLONG iz,jz,mz;
  TFLOAT temp,dtau,efact;
  TFLOAT *eth,*epi;           /* extrapolated cell edge th and pi */
  TFLOAT *b;                  /* extrapolated cell edge sigma*T^4 */

  /* make space */
  mz = nz + 1;
  eth = floatbuff(mz);
  epi = floatbuff(mz);
  b = floatbuff(mz);

  /* compute the Planck function profile on cell edges */
  zedgecomp(nz,dz,th,eth);
  zedgecomp(nz,dz,pi,epi);
  eth[0] = th[0];               /* neutral bl assumed */
  for (iz = 0; iz < mz; iz++) {
    temp = eth[iz]*epi[iz]/CP;
    b[iz] = SBCONST*temp*temp*temp*temp;
  }

  /* do the kplus profile */
  kplus[0] = 0.;
  for (iz = 1; iz < mz; iz++) {
    jz = iz - 1;
    dtau = tau[iz] - tau[jz];
    efact = exp(elimit(-dtau));
    kplus[iz] = efact*kplus[jz] + lcomp(dtau, b[jz], b[iz]);
  }

  /* do the kminus profile */
  kminus[nz] = 0.;
  for (iz = nz - 1; iz >= 0; iz--) {
    jz = iz + 1;
    dtau = tau[jz] - tau[iz];
    efact = exp(elimit(-dtau));
    kminus[iz] = efact*kminus[jz] + lcomp(dtau,b[jz],b[iz]);
  }

  /* free space */
  fbuf(b);
  fbuf(epi);
  fbuf(eth);
}

/* taucomp -- compute the path length for a band -- tau is defined on
 * cell edges and is zero at the surface
 */
static void taucomp(TLONG mz, TLONG ib, TFLOAT *dz, TFLOAT *th,
		    TFLOAT *rt, TFLOAT *pi, TFLOAT *rho, TFLOAT *rsat,
		    TFLOAT abs, TFLOAT cld, TFLOAT *tau)
{
  TLONG iz,kz;
  TFLOAT mutot,temp,pbroad,rl,rvap,mutotmax,rfract;

  /* integrate upward from surface to get tau -- iz = upper edge index,
   * kz = cell center and lower edge index
   */
  tau[0] = 0.;
  for (iz = 1; iz < mz; iz++) {
    kz = iz - 1;

    /* compute pressure broadening factor */
    temp = th[kz]*pi[kz]/CP;
    pbroad = pow(temp/TREF,0.5)*rho[kz]/RHOREF;

    /* compute vapor-condensate split */
    rl = rt[kz] - rsat[kz];
    if (rl < 0.) rl = 0.;
    rvap = rt[kz] - rl;
    rfract = rvap/RTREF;

    /* compute mutot in cell center at level kz for desired band */
    if (ib == CO2) {                 /* carbon dioxide bands */
      mutot = abs*pbroad;
    }
    else if (ib == CONT) {           /* water vapor continuum */
      mutot = abs*rfract*rfract;
    }
    else {                           /* water vapor bands */
      mutot = abs*pbroad*rfract;
    }

    /* add in cloud part */
    mutot += cld*rl*rho[kz]/(ZSCALE*RHOREF);

    /* limit mutot for numerical reasons -- larger values are simply opaque */
    mutotmax = 30./(rho[kz]*dz[kz]);
    if (mutot > mutotmax) mutot = mutotmax;

    /* compute tau at the next level */
    tau[iz] = tau[kz] + rho[kz]*dz[kz]*mutot;
  }
}

/* the main radiation routine */
static void radiation(TLONG nz, TLONG nosun, TFLOAT cld, TFLOAT sst,
		      TFLOAT *zc, TFLOAT *dz, TFLOAT *th, TFLOAT *the,
		      TFLOAT *rt, TFLOAT *pi, TFLOAT *rho, TFLOAT *rsat,
		      TFLOAT *sthe, TFLOAT *flux)
{
  TLONG ib,iz;                /* index parameters */
  TLONG mz;                   /* number of grid edge points */
  TFLOAT scos;                /* factor for solar zenith angle */
  TFLOAT bsurf;               /* upward energy flux from surface */
  TFLOAT topefact,botefact,iplus,iminus;  /* intermediate results */
  TFLOAT normfactor;          /* normalization factor */
  TFLOAT *eflux;              /* total flux on cell edges */
  TFLOAT *tau;                /* actual absorptivities */
  TFLOAT *kplus,*kminus;      /* intermediate profiles */

  /* arrays containing per-band information about radiation model */

#if ONECHANNEL
  /* this has a single band independent of water vapor */

  /* specific absorption coefficients for each band */
  static TFLOAT abs[NBANDS] = {1.};

  /* downward flux of solar radiation for each band for overhead sun */
  static TFLOAT solar[NBANDS] = {0.};

  /* fraction of Planck-weighted spectrum occupied by each band */
  static TFLOAT fract[NBANDS] = {1.0};
#endif

#if MULTICHANNEL
  /* arrays containing per-band information about radiation model
   * band 0 is CO2 while band 1 is water vapor continuum
   * bands 2-7 are water vapor lines
   * This is the radiation scheme used in Raymond and Torres (1998)
   */

  /* specific absorption coefficients for each band */
  static float abs[NBANDS] = {2.,0.5,1.,10.,100.,1000.,10000.,100000.};

  /* downward flux of solar radiation for each band for overhead sun */
  static float solar[NBANDS] = {0.,160.,60.,30.,15.,7.,0.,0.};

  /* fraction of Planck-weighted spectrum occupied by each band */
  static float fract[NBANDS] = {0.11,0.20,0.19,0.16,0.13,0.10,0.07,0.04};
#endif

  /* make space -- nz = number of cells, mz = number of cell edges */
  mz = nz + 1;
  eflux = floatbuff(mz);
  tau = floatbuff(mz);
  kplus = floatbuff(mz);
  kminus = floatbuff(mz);

  /* define solar zenith angle factor and surface Planck function
   * we include the nosun flag for testing the scheme with no
   * solar input from the top
   */
  if (nosun) {
    scos = 0.;
  }
  else {
    scos = 1./3.14159;
  }
  bsurf = SBCONST*sst*sst*sst*sst;

  /* zero edge flux array */
  for (iz = 0; iz < mz; iz++) eflux[iz] = 0.;

  /* loop on bands */
  for (ib = 0; ib < NBANDS; ib++) {

    /* compute tau, the actual absorption coefficient, for this band */
    taucomp(mz,ib,dz,th,rt,pi,rho,rsat,abs[ib],cld,tau);

    /* compute the kplus and kminus profiles */
    kcomp(nz,dz,tau,th,pi,kplus,kminus);

    /* add in the flux profile for this band */
    for (iz = 0; iz < mz; iz++) {
      botefact = exp(elimit(-tau[iz]));
      topefact = exp(elimit(tau[iz] - tau[nz]));
      iplus = fract[ib]*(bsurf*botefact + kplus[iz]);
      iminus = scos*solar[ib]*topefact + fract[ib]*kminus[iz];
      eflux[iz] += (iplus - iminus);
    }
  }

  /* compute the cell-centered radiative source and flux */
  for (iz = 0; iz < nz; iz++) {
    normfactor = the[iz]/(rho[iz]*pi[iz]*th[iz]);
    sthe[iz] = -normfactor*(eflux[iz + 1] - eflux[iz])/dz[iz];
    flux[iz] = 0.5*(eflux[iz + 1] + eflux[iz]);
  }

  /* free space */
  fbuf(kminus);
  fbuf(kplus);
  fbuf(tau);
  fbuf(eflux);
}

/* fixed radiation code */
static void fradiation(TLONG nz, TFLOAT *zc, TFLOAT radcool, TFLOAT tpause,
	   TFLOAT radbrk, TFLOAT *sthe_rad, TFLOAT *sthe_radclr,
	   TFLOAT *irflux, TFLOAT *irfluxclr)
{
  TLONG iz;
  TFLOAT z;

  for (iz = 0; iz < nz; iz++) {
    z = zc[iz];

    /* compute the fixed radiative cooling profile */
    if (z < radbrk*tpause) {
      sthe_rad[iz] = 1.;
    }
    else if (z < tpause && z > radbrk*tpause) {
      sthe_rad[iz] = 1. - (z/tpause - radbrk)/(1. - radbrk);
    }
    else {
      sthe_rad[iz] = 0.;
    }
    sthe_rad[iz] = -sthe_rad[iz]*radcool/86.4;

    /* fill in other stuff -- ignore irflux var for now */
    sthe_radclr[iz] = sthe_rad[iz];
    irflux[iz] = irfluxclr[iz] = 0.;
  }
}



/* diabat3 -- Make all the diabatic calculations
 * nz - number of vertical grid cells
 * land - land/water flag
 * dorad - do radiation flag
 * frad - flag indicating fixed radiation
 * radcool - fixed radiation cooling rate
 * tpause - height of tropopause
 * radbrk - frac of tpause for which fixed radiative cooling starts decrease
 * cvc - convective mixing rate parameter
 * cvs - stratiform rain rate parameter
 * cvp - convective rain rate parameter
 * cve - rain evaporation rate parameter
 * theslop - range of theta_e going from suppressed to full deep convection
 * pstiff - stiffness parameter in convective precip generation
 * pscale - shape parameter for convective precip production
 * pbltop - top of planetary boundary layer
 * cdrag - surface drag coefficient, used also for thermodynamics
 * wscale - gustiness parameter for surface fluxes
 * sst - sea surface temperature
 * lfrac - latent heat fraction out of total over land
 * eflux0 - thetae flux over lang
 * cld - cloud absorption parameter
 * cfract - cloud fractional area coverage
 * dt - time step
 * zc - cell center z values
 * xcell - vertical cell sizes in calling model space
 * th - potential temperature profile - cell centered, as all that follow
 * rt - total water mixing ratio profile
 * pi - exner function profile
 * u - x wind profile
 * v - y wind profile
 * xrho - mass density profile in calling model space
 * sth - theta source term profile, return variable as all below
 * sthe_cu - theta-e source term profile from convection
 * sthe_rad - theta-e source term profile from radiation
 * srt - total cloud water mixing ratio source term profile
 * su - x momentum source term profile
 * sv - y momentum source term profile
 * irflux - radiative flux profile
 * rain - rainfall rate
 * eflux - surface theta-e flux
 * rflux - surface total water flux (evaporation)
 * uflux - surface x momentum flux
 * vflux - surface y momentum flux
 * convthrot - convective throttle value
 * fluxthrot - flux throttle value
 * cvcaug - deep convective mixing rate
 */

/* this variable is used to be sure initialization is done right */
static TLONG firsttime = 1;

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
	     TFLOAT *fluxthrot, TFLOAT *cvcaug)
{
  TLONG iz;                   /* index variables */
  TFLOAT *the;                /* equivalent potential temperature */
  TFLOAT *rsat;               /* saturation mixing ratio */
  TFLOAT *sthe_radclr;        /* clear sky radiation source */
  TFLOAT *irfluxclr;          /* clear sky radiation flux */
  TFLOAT *dz;                 /* vertical cell size in geometric space */
  TFLOAT *rho;                /* geometric space density */
  TFLOAT surfpi, surfth, surftemp;   /* surface stuff for radiation */

  /* radiation routines need to be called on first time step */
  if (firsttime) {
    firsttime = 0;
    if (dorad == 0) {
      fprintf(stderr,"diabat3: dorad must be nonzero on first call\n");
      exit(1);
    }
  }
/*         fprintf(stderr,"time: %ds\n",(int)TIME); */
/*       fprintf(stderr,"nz: %d\n",nz); */
/*       fprintf(stderr,"land: %d\n",land); */
/*       fprintf(stderr,"dorad: %d\n",dorad); */
/*       fprintf(stderr,"frad: %d\n",frad); */
/*       fprintf(stderr,"radcool: %f\n",radcool); */
/*       fprintf(stderr,"tpause: %f\n",tpause); */
/*       fprintf(stderr,"radbrk: %f\n",radbrk); */
/*       fprintf(stderr,"cvc: %f\n",cvc); */
/*       fprintf(stderr,"cvs: %f\n",cvs); */
/*       fprintf(stderr,"cvp: %f\n",cvp); */
/*       fprintf(stderr,"cve: %f\n",cve); */
/*       fprintf(stderr,"theslop: %f\n",theslop); */
/*       fprintf(stderr,"pstiff: %f\n",pstiff); */
/*       fprintf(stderr,"pscale: %f\n",pscale); */
/*       fprintf(stderr,"pbltop: %f\n",pbltop); */
/*       fprintf(stderr,"cdrag: %f\n",cdrag); */
/*       fprintf(stderr,"wscale: %f\n",wscale); */
/*       fprintf(stderr,"sst: %f\n",sst); */
/*       fprintf(stderr,"lfrac: %f\n",lfrac); */
/*       fprintf(stderr,"elfux0: %f\n",eflux0); */
/*       fprintf(stderr,"cld: %f\n",cld); */
/*       fprintf(stderr,"cfract: %f\n",cfract); */
/*       fprintf(stderr,"dt: %f\n",dt); */
/*       check("post-diabat3 zc 1",zc, 0, nmt.nz); */
/*       check("post-diabat3 xcell 1",xcell, 0, nmt.nz); */
/*       check("post-diabat3 th 1",th, 0, nmt.nz); */
/*       check("post-diabat3 rt 1",rt, 0, nmt.nz); */
/*       check("post-diabat3 pi 1",pi, 0, nmt.nz); */
/*       check("inside-diabat3 u 1",u, 0, nmt.nz); */
/*       check("post-diabat3 v 1",v, 0, nmt.nz); */
/*       check("post-diabat3 xrho 1",xrho, 0, nmt.nz); */
/*       check("post-diabat3 sth 1",sth, 0, nmt.nz); */
/*       check("post-diabat3 sthe_cu 1",sthe_cu, 0, nmt.nz); */
/*       check("post-diabat3 sthe_rad 1",sthe_rad, 0, nmt.nz); */
/*       check("post-diabat3 srt 1",srt, 0, nmt.nz); */
/*       check("inside-diabat3 su 1",su, 0, nmt.nz); */
/*       check("inside-diabat3 sv 1",sv, 0, nmt.nz); */
/*       check("post-diabat3 irflux 1",irflux, 0, nmt.nz); */
/*       fprintf(stderr,"rain: %f\n",*rain); */
/*         fprintf(stderr,"in diabat3: eflux = %f\n",*eflux); */
/*       fprintf(stderr,"rflux: %f\n",*rflux); */
/*       fprintf(stderr,"uflux: %f\n",*uflux); */
/*       fprintf(stderr,"vflux: %f\n",*vflux); */
/*       fprintf(stderr,"convthrot: %f\n",*convthrot); */
/*       fprintf(stderr,"fluxthrot: %f\n",*fluxthrot); */
      
  /* make space */
  the = floatbuff(nz);
  rsat = floatbuff(nz);
  sthe_radclr = floatbuff(nz);
  irfluxclr = floatbuff(nz);
  dz = floatbuff(nz);
  rho = floatbuff(nz);

  /* compute the equivalent potential temperature and
     saturation mixing ratio */
  for (iz = 0; iz < nz; iz++) {
    rsat[iz] = rs(th[iz],pi[iz]);
    the[iz] = thetaecomp2(th[iz],rt[iz],rsat[iz]);

    /* compute cell sizes and geometric density */
    
    rho[iz] = PREF*pow(pi[iz]/CP,1./KAPPA)/(KAPPA*th[iz]*pi[iz]);
    /*fprintf(stderr,"xrho[%d]=%.15f, xcell[%d]=%.15f, rho[%d]=%.15f\n",\
      iz,xrho[iz],iz,xcell[iz],iz,rho[iz]);*/
    dz[iz] = xrho[iz]*xcell[iz]/rho[iz];
  }

  /* call the convective parameterization */

  cuparam(nz,land,cvc,cvs,cvp,cve,theslop,pstiff,pscale,pbltop,cdrag,
	  wscale,sst,lfrac,eflux0,dt,zc,dz,th,the,rt,pi,u,v,rho,rsat,
	  sthe_cu,srt,su,sv,rain,eflux,rflux,uflux,vflux,
	  convthrot,fluxthrot,cvcaug);

  /* do radiation only if specified */
  if (dorad) {

    /* use fixed radiation if frad not zero */
    if (frad) {
      fradiation(nz,zc,radcool,tpause,radbrk,sthe_rad,sthe_radclr,
		 irflux,irfluxclr);
    }

    /* otherwise call the radiative parameterization -- if over land, use the
     * extrapolated surface temperature in place of the sea surface
     * temperature */
    else {
      if (land) {
	surfpi = surfextrap(pi,dz);
	surfth = surfextrap(th,dz);
	surftemp = surfth*surfpi/1005.;
      }
      else {
	surftemp = sst;
      }

      /* cloudy atmosphere radiation */
      radiation(nz,0,cld,surftemp,zc,dz,th,the,rt,pi,rho,rsat,sthe_rad,
		irflux);

      /* clear atmosphere radiation */
      radiation(nz,0,0.,surftemp,zc,dz,th,the,rt,pi,rho,rsat,sthe_radclr,
		irfluxclr);
    }
      
    /* the actual radiation is a weighted average of cloudy and clear */
    for (iz = 0; iz < nz; iz++) {
      sthe_rad[iz] = cfract*sthe_rad[iz] + (1. - cfract)*sthe_radclr[iz];
      irflux[iz] = cfract*irflux[iz] + (1. - cfract)*irfluxclr[iz];
    }

  }

  /* compute the theta source term -- note that the formula below depends
   * on the definition of theta-e used in this program --
   * compute sth every time since sthe_cu and srt are updated every time
   */
  for (iz = 0; iz < nz; iz++) {
    if (iz < 2){
      /*fprintf(stderr,"(iz=%d) th=%.12f, sthe_cu=%.12f, sthe_rad=%.15f, the=%.15f, srt=%.15f\n",iz,th[iz],sthe_cu[iz],sthe_rad[iz],the[iz],srt[iz]);*/
    }
    sth[iz] = th[iz]*((sthe_cu[iz] + sthe_rad[iz])/the[iz] - ALFA*srt[iz]);
    fflush(stderr);
  }


  /* free space */
  fbuf(rho);
  fbuf(dz);
  fbuf(irfluxclr);
  fbuf(sthe_radclr);
  fbuf(rsat);
  fbuf(the);
}




/*******************************************************************
 * thermo.c
 *******************************************************************/

/* thermo.c -- Moist thermodynamics package -- derived from diabat2.
 * Variable naming convection:
 * th: Potential temperature (K)
 * pi: Exner function (J/kg/K)
 * td: Dewpoint temperature (K)
 * the: Equivalent potential temperature (K)
 * rt: Total cloud water (vapor plus droplets) mixing ratio (g/g)
 * rsat: Saturation mixing ratio (g/g)
 * p: Pressure (pa)
 *
 * $Id: thermo.c,v 1.3 2005/03/27 18:01:58 raymond Exp $
 */

/* thetasat -- assume saturation and calculate theta -- really
 * theta-l */
static TFLOAT thetasat(TFLOAT the,TFLOAT pi,TFLOAT thstart,TINT nloops)
{
  TFLOAT theta,rsval,expfactor,drsdth;
  TINT loop;
  theta = thstart;
  for (loop = 0; loop < nloops; loop++) {
    rsval = rs(theta,pi);
    expfactor = the*exp(-ALFA*rsval);
    drsdth = rs(theta + 1.,pi) - rsval;
    theta -= (theta - expfactor)/(1. + ALFA*expfactor*drsdth);
  }
  return(theta);
}

/* rs - compute saturation mixing ratio (g/g) as function of
 * potential temp (K) and Exner function (J/kg/K)
 * use simple formula
 */
TFLOAT rs(TFLOAT th,TFLOAT pi)
{
  TFLOAT p,es,rsat;
  p = PREF*pow(pi/CP,1./KAPPA);
  /*fprintf(stderr,"in rs():th=%.15f, pi=%.15f\n",th,pi); */
  es = escalc(th*pi/CP);
  /*fprintf(stderr,"in rs():es=%.15f\n",es);*/
  rsat = EPSI_LON*es/p;
  /*fprintf(stderr,"in rs():rsat=%.15f\n",rsat); */
  return(rsat);
}

/* escalc -- compute vapor pressure (pa) from dewpoint (K) -- Teten formula */
TFLOAT escalc(TFLOAT td)
{
  /*fprintf(stderr,"in escalc: td=%.15f\n",td);*/
  return(611.2*exp(17.67*(td - 273.15)/(td - 29.65)));
}

/* thetacomp2 -- do this computation by a lookup table */
#define ITHEMAX 401
#define IPIMAX 501
#define DTHE 1.
#define DPI 2.
#define THE0 200.
#define PI0 200.

TFLOAT thetacomp2(TFLOAT the,TFLOAT rt,TFLOAT pi)
{
  static TFLOAT v[ITHEMAX][IPIMAX];
  static TINT begin = 1;
  TINT ithe,ipi,jthe,jpi;
  TFLOAT theval,pival,athe,api,othe,opi,theta;
  
  /* fill the lookup table the first time around */
  if (begin) {
    for (ithe = 0; ithe < ITHEMAX; ithe++) {
      for (ipi = 0; ipi < IPIMAX; ipi++) {
	theval = THE0 + ithe*DTHE;
	pival = PI0 + ipi*DPI;
	v[ithe][ipi] = thetasat(theval,pival,theval*exp(-ALFA*rt),5);
      }
    }
    begin = 0;
  }
  
  /* compute test value of theta */
  theta = the*exp(-ALFA*rt);
  
  /* if unsaturated, don't do anything */
  if (rt < rs(theta,pi)) {
  }
  
  /* if saturated, we interpolate */
  else {
    ithe = (the - THE0)/DTHE;
    if (ithe < 0) ithe = 0;
    if (ithe >= ITHEMAX - 1) ithe = ITHEMAX - 2;
    jthe = ithe + 1;
    ipi = (pi - PI0)/DPI;
    if (ipi < 0) ipi = 0;
    if (ipi >= IPIMAX - 1) ipi = IPIMAX - 2;
    jpi = ipi + 1;
    athe = (the - (THE0 + DTHE*ithe))/DTHE;
    othe = 1. - athe;
    api = (pi - (PI0 + DPI*ipi))/DPI;
    opi = 1. - api;
    theta = othe*(opi*v[ithe][ipi] + api*v[ithe][jpi])
      + athe*(opi*v[jthe][ipi] + api*v[jthe][jpi]);
  }
  return(theta);
}

/* thetacomp -- compute potential temperature from equivalent
 * potential temperature, total water mixing ratio, and Exner function.
 */
TFLOAT thetacomp(TFLOAT the,TFLOAT rt,TFLOAT pi)
{
  TFLOAT theta;
  
  /* compute test value of theta */
  theta = the*exp(-ALFA*rt);
  
  /* if unsaturated, don't do anything */
  if (rt < rs(theta,pi)) {
  }
  
  /* if saturated -- we have to do a Newton's method iteration */
  else {
    theta = thetasat(the,pi,theta,5);
  }
  
  /* return with theta value */
  return(theta);
}

/* calculate theta_e */
TFLOAT thetaecomp(TFLOAT th,TFLOAT rt,TFLOAT pi)
{
  return(thetaecomp2(th,rt,rs(th,pi)));
}

/* calculate theta_e given the saturation mixing ratio */
TFLOAT thetaecomp2(TFLOAT th,TFLOAT rt,TFLOAT rsat)
{
  TFLOAT rvap;
  rvap = rt > rsat ? rsat : rt;
  return(th*exp(ALFA*rvap));
}

/* thetaescomp -- compute saturated thetae value given theta (K) and
 * pi (J/kg/K)
 */
TFLOAT thetaescomp(TFLOAT th,TFLOAT pi)
{
  return(th*exp(ALFA*rs(th,pi)));
}

/* gmoist -- compute (d theta/dz) for moist adiabatic parcel */
TFLOAT gmoist(TFLOAT th,TFLOAT pi)
{
  TFLOAT drsdpi,drsdth,rate;
  drsdpi = rs(th,pi + 1.) - rs(th,pi);
  drsdth = rs(th + 1.,pi) - rs(th,pi);
  rate = GEE*ALFA*drsdpi/(1. + th*ALFA*drsdth);
  return(rate);
}

/* pifromp -- compute exner function (J/K/kg) from pressure (Pa) */
TFLOAT pifromp(TFLOAT p)
{
  return(CP*pow(p/PREF,KAPPA));
}

/* pfrompi -- compute pressure (Pa) from exner function (J/K/kg) */
TFLOAT pfrompi(TFLOAT pi)
{
  return(PREF*pow(pi/CP,1./KAPPA));
}

/* rhocomp -- compute density (kg/m^3) from theta (K) and exner (J/K/kg) */
TFLOAT rhocomp(TFLOAT th, TFLOAT pi)
{
  return(PREF*pow(pi/CP,1./KAPPA)/(KAPPA*th*pi));
}

/* random number generator uniform on -1 to 1 */
static double myrand(void)
{
  double rvalue;
  rvalue = 2.*(((double)random())/((double)RAND_MAX)) - 1.;
  return(rvalue);
}






/*******************************************************************
 * wrapper.c
 *******************************************************************/

/*******************************************************************************
 * wrapper.c
 * 
 * This is a wrapper function to take EPIC 3D arrays, convert them into
 * 1D vertical columns, and send them to David Raymond's cumulus 
 * parametrization, diabat3().
 * 
 * The following papers are cited:
 *
 * EPICAM - "The Explicit Planetary Isentropic-Coordinate (EPIC) Model"
 *          (Dowling et al., 1998).
 *
 * EAMITF - "The EPIC atmospheric model with an isentropic/terrain-following
 *          hybrid vertical coordinate" (Dowling et al., 2006). 
 *
 * Programmer: Mike Herman (mherman@nmt.edu)
 * Date:       June 13, 2007.
 *******************************************************************************/

/*******************************************************************
 * define constant scalar model parameters
 * 
 * These could be defined in a data file, to be read by this program.
 * Later, we can pull them in, using EPIC's 'initial' program. 
 *******************************************************************/

/*******************************************************************************
 * nmt_cumulus()
 *
 * A wrapper function that culls data from EPIC and sends it to diabat3. The
 * diagnostic varaibles from EPIC are in the form of 3d arrays representing
 * the entire EPIC model domain. Here, we grab a vertical column at a time
 * and send it to diabat3, along with several scalar values determined by
 * the user. Then, diabat3 returns several scalar values and tendencies
 * in the form of 1d arrays. 
 *******************************************************************************/
void nmt_cumulus(planetspec *planet){

  /********************************
   * cumulus' own little variables 
   ********************************/
  int       kr;                     /* for vertical layers of diabat3      */ 
  int       kk;                     /* for get_var calls in EPIC           */ 
  int       I,J,K;                  /* for iterating over EPIC dimensions  */
  TFLOAT    g = 9.807;              /* define accel. of gravity            */
  TFLOAT    new_sst = 0;            /* to hold the latest value of sst     */ 

  /***************************************
   * define a static variable to count the 
   * number of iterations btw radiations
   ***************************************/
  static TLONG raditer = 5;
  
  /********************
   * set the dorad flag
   ********************/
  if (raditer++ >= RADITERMAX){
    raditer = 0;
    nmt.dorad = 1;
  }else{
    nmt.dorad = 0;
  }

  /********************************************************************
   * Get the (I,J)th column of each necessary EPIC diagnostic variable.
   * 
   * NOTE: We send in the entire EPIC column.
   *
   * nmt.nz - the number of cells in diabat3. It is given 
   *      by the number of vertical layers in EPIC.
   *      If there are nk layers in EPIC, there are 
   *      nmt.nz = nk cells in diabat3.
   *
   * kr - allows us to reference the diabat3 inputs starting from the
   *      top of the diabat3 model. We subtract 1 from kr (below), 
   *      since the diabat3 cell indices begin at 0. 
   *      kr = nmt.nz-1 = nk-2. 
   *
   * kk - kk = 2*K. This is used for  referencing any value of any 
   *      EPIC prognostic variable, and for referencing sigmatheta.
   *
   * diabat3():
   * The model accepts input fields zc thru w in a single column at cell-
   * centered levels defined by zc.  The output profiles sth thru sv
   * have the same structure.  The values of dz for each cell
   * should be consistent with zc.
   *
   ********************************************************************/
  for (I = ILO; I <= IHI; I++) {
    for (J = JLO; J <= JHI; J++) {
      kr = nmt.nz-1;
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K;

        /* get kr'th value of zc (km). 
         * 
         * GZ3 is the geopotential (m^2/s^2) defined on
         * the bottom interface of the Kth layer of the EPIC model. 
         * We divide by g to get the height in meters, 
         * then convert to km. 
         */
/*NOTE: need to convert to GZ3
        nmt.zc[kr] = GZ(K,J,I)/g/1000;
*/

        /* get the kr'th value of th (K).
         *
         * THETA is the value of theta (K) defined on the bottom 
         * interface of the Kth layer of the EPIC model.
         * Theta is actually a hybrid of diagnostic and 
         * prognostic values of used in EPIC. See EAMITF, p4.
         */
        nmt.th[kr] = THETA2(K,J,I);
        
        /* get the kr'th value of the exner function (J/kg/K).
         *
         * EXNER3 is the value of the exner function (J/kg/K)
         * on the bottom interface of the Kth layer of the EPIC model.
         */
/*NOTE: need to convert to EXNER3
        nmt.pi[kr] = EXNER(K,J,I);
*/

        /* get the xrho of the EPIC layer in model density units
         *
         * HDRY3 is the value of hybrid density (kg/m^2/K) on the
         * bottom interface of the Kth layer of the EPIC model.
         */
        nmt.xrho[kr] = HDRY(K,J,I);

        /* get kr'th value of the total cloud water mixing ratio (g/g).
         * 
         * Here, we derive the mixing ratio, rs, from theta (th)
         * and the exner funtion (pi), using a function defined
         * in thermo.c, rs().
         *
         * Alternatively, we grab the current value of total water density,
         * H_2O VAPOR.
         */
        /*rt[kr] = rs(th[kr],pi[kr]);*/
        nmt.rt[kr] = get_var(planet,H_2O_INDEX,VAPOR,grid.it_h,kk,J,I);

        /* get the kr'th xcell column height in model vertical coord. units.
         * 
         * We measure the cells as the difference in the values of 
         * the hybrid coordinate, sigmatheta, on the EPIC layers 
         * above and below the bottom interface of the Kth layer of the 
         * EPIC model. However, at the bottom cell in diabat3, we must
         * get the difference in sigmatheta btw the layer bounding the
         * cell above, and the interface immediately below that layer. 
         * This is due to the fact that there are no more layers below the 
         * interface. This also means that all the variables defined at cell
         * center in diabat3 for this bottom cell are actually defined at the
         * floor of the cell. 
         */
          nmt.xcell[kr] = grid.sigmatheta[kk-1] - grid.sigmatheta[kk+1];

        /* get the kr'th value of u (km/ks=m/s).
         * 
         * Here, we grab the longitudinal velocity, U (m/s), on the bottom 
         * interface of the Kth layer of the EPIC model. This is done using
         * the EPIC get_var() function, which returns any EPIC prognostic at
         * given value of sigmatheta.
         */
        nmt.u[kr] = U(grid.it_uv,K,J,I);
        /*fprintf(stderr,"U (%d,%d,%d) = %.16f\n",K,J,I,nmt.u[kr]);*/

        /* get the kr'th value of v (km/ks=m/s).
         * 
         * Here, we grab the meridional velocity, V (m/s), on the bottom 
         * interface of the Kth layer of the EPIC model. This is done using
         * the EPIC get_var() function, which returns any EPIC prognostic at
         * given value of sigmatheta.
         */
        nmt.v[kr] = V(grid.it_uv,K,J,I);

        /* decrement kr, since we are sending the columns in bottom-to-top */
        kr--;
      }

      /***********************************************************
       * call diabat3 after the (i,j)th column arrays are complete
       ***********************************************************/
      diabat3(nmt.nz,
              nmt.land,
              nmt.dorad,
              nmt.frad,
              nmt.radcool,
              nmt.tpause,
              nmt.radbrk,
              nmt.cvc,
              nmt.cvs,
              nmt.cvp,
              nmt.cve,
              nmt.theslop,
              nmt.pstiff,
              nmt.pscale,
              nmt.pbltop,
              nmt.cdrag,
              nmt.wscale,
              nmt.sst,
              nmt.lfrac,
              nmt.eflux0,
              nmt.cld,
              nmt.cfract,
              nmt.dt,
              nmt.zc,
              nmt.xcell,
              nmt.th,
              nmt.rt,
              nmt.pi,
              nmt.u,
              nmt.v,
              nmt.xrho,
              nmt.sth,
              nmt.sthe_cu,
              nmt.sthe_rad,
              nmt.srt,
              nmt.su,
              nmt.sv,
              nmt.irflux,
              &nmt.rain,
              &nmt.eflux,
              &nmt.rflux,
              &nmt.uflux,
              &nmt.vflux,
              &nmt.convthrot,
              &nmt.fluxthrot,
              &nmt.cvcaug);
      

      /* return tendency values to EPIC */
      kr = nmt.nz-1;
      for (K = KLO; K <= KHI; K++) {
        
        /*********************************************************
         * get source terms and convert from units/ks to units/s.
         * Since su and sv are in km/ks, we need not convert.
         *
         * NOTE: DB_OUTS(STH,2,J,I) corresponds to K = 2 in EPIC, and
         *       the top cell center in diabat3.
         *       Also, the indices range up to K = nk.
         *
         * UNITS:      [STH] = (K/s)
         *             [SRT] = (kg/kg/s)
         *              [SU] = (m/s^2)
         *              [SV] = (m/s^2)
         *        [STHE_RAD] = (K/s)
         *         [STHE_CU] = (K/s)
         *          [IRFLUX] = (W/m)
         *      
         *********************************************************/
        DB_OUTS(STH,K,J,I) = nmt.sth[kr]/1000.;
        DB_OUTS(SRT,K,J,I) = nmt.srt[kr]/1000.;
        DB_OUTS(SU,K,J,I) = nmt.su[kr];
        DB_OUTS(SV,K,J,I) = nmt.sv[kr];
        DB_OUTS(STHE_RAD,K,J,I) = nmt.sthe_rad[kr]/1000.;
        DB_OUTS(STHE_CU,K,J,I) = nmt.sthe_cu[kr]/1000.;
        DB_OUTS(IRFLUX,K,J,I) = nmt.irflux[kr];

        /* grab 2d output diagnostics for extracting */
        /*fprintf(stderr,"nmt.eflux: %.15f\n",nmt.eflux);*/
        THE_FLUX(K,J,I) = nmt.eflux;
        /*fprintf(stderr,"THE_FLUX(%d,%d,%d)=%.15f\n",K,J,I,THE_FLUX(K,J,I)); */
        RT_FLUX(K,J,I) = nmt.rflux;
        U_FLUX(K,J,I) = nmt.uflux;
        V_FLUX(K,J,I) = nmt.vflux;
        CONVTHROT(K,J,I) = nmt.convthrot;
        FLUXTHROT(K,J,I) = nmt.fluxthrot;
        RAIN_RATE(K,J,I) = nmt.rain;
        kr--;
      }
    }
  }
} /* end of nmt_cumulus() */


/************************************************************************
 * check() - a checking function to output values at certain locations
 * 
 * To check arrays, pass the array name: array_name, type = 0, size.
 * To check scalars, pass scalar's address: &scalar_name, type = 1, size. 
 ************************************************************************/
void check(char *place, TFLOAT *array, int type, TLONG nelem){
  int             i, 
                  tabcounter = 0;
  TFLOAT          value;
  /* label and format the output */
  fprintf(stderr,"\n%s:\n\t\t", place);
  /* to show array values */
  if(type == 0){
    for(i = 0; i < nelem; i++){
      /* format array elements in nice rows of 5 */
      if(tabcounter == 5){
        fprintf(stderr,"\n\t\t");
        tabcounter = 0;
      }
      value = array[i];
      fprintf(stderr,",%.15f ",value);
      tabcounter++;
    }
  /* to show scalar values */ 
  }else{
    fprintf(stderr,"%.15f ",array[0]);
  }
  fprintf(stderr,"\n\n");
}

/***********************************************************
 * nmt_allocate() - allocate memory for the diabat3 values
 *
 ***********************************************************/
void nmt_allocate(){
  if(grid.nmt_physics_on){
    int
      i;
    static char
      dbmsname[]="nmt_allocate";

    for (i = 0; i < NUM_DIABAT3_OUTS; i++) {
      nmt.outs[i] = fvector(0,Nelem3d-1,dbmsname);
    }
    nmt.zc = fvector(0,nmt.nz-1,dbmsname);
    nmt.xcell = fvector(0,nmt.nz-1,dbmsname);
    nmt.th = fvector(0,nmt.nz-1,dbmsname);
    nmt.rt = fvector(0,nmt.nz-1,dbmsname);
    nmt.pi = fvector(0,nmt.nz-1,dbmsname);
    nmt.u = fvector(0,nmt.nz-1,dbmsname);
    nmt.v = fvector(0,nmt.nz-1,dbmsname);
    nmt.xrho = fvector(0,nmt.nz-1,dbmsname);
    nmt.sth = fvector(0,nmt.nz-1,dbmsname);
    nmt.sthe_cu = fvector(0,nmt.nz-1,dbmsname);
    nmt.sthe_rad = fvector(0,nmt.nz-1,dbmsname);
    nmt.srt = fvector(0,nmt.nz-1,dbmsname);
    nmt.su = fvector(0,nmt.nz-1,dbmsname);
    nmt.sv = fvector(0,nmt.nz-1,dbmsname);
    nmt.irflux = fvector(0,nmt.nz-1,dbmsname);
    nmt_memset();
  }
}


/******************************************************************
 * nmt_apply_sources() - send diabat3 source terms to EPIC
 *
 * This is done from epic_funcs_diag.c, source_sink().
 *
 ******************************************************************/
void nmt_apply_sources(planetspec *planet){
  if(grid.nmt_physics_on == 1){
    int
      K,J,I,kk;
    EPIC_FLOAT
      sth,srt;

    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {

          /*********************************************************************
           * pass diabat3 theta source term to EPIC's HEAT3.
           * 
           * STH and sth are from the Kth EPIC layer. All other values are 
           * defined on the bottom interface of the Kth EPIC layer.
           *
           * [STH] = (K/s)
           * [EXNER3] = (J/kg/K)
           * Therefore, [STH*EXNER3] = (J/kg/s), which are the units of HEAT3.
           * 
           *********************************************************************/
          if ( K == KHI){
            /************************************************************** 
             * if at the bottom layer, extrapolate to the bottom interface, 
             * using get_var's formula 
             **************************************************************/
            sth = DB_OUTS(STH,K,J,I)+(DB_OUTS(STH,K,J,I)-DB_OUTS(STH,K-1,J,I))* \
              (grid.sigmatheta[kk+1]-grid.sigmatheta[kk])*grid.dsgth_inv[kk-1];
          }else{
            /****************************************************************** 
             * otherwise, take a simple average of the current and lower layers 
             ******************************************************************/
            sth = (DB_OUTS(STH,K,J,I) + DB_OUTS(STH,K+1,J,I))*.5;
          }
          HEAT3(K,J,I) += sth*EXNER3(K,J,I);
          /*fprintf(stderr,"nmt_apply_sources: DB_OUTS(STH,%d,%d,%d) = %.16f\n",\
            K,J,I,sth);*/

          /********************************************************************
           * pass diabat3 momentum sources to EPIC's DUDT, DVDT.
           * 
           * Here, we are adding to the EPIC momentum tendencies. This pertains
           * to eqs.4a, 4b, p223, EPICAM. Since DUDT and DVDT are advected on 
           * the EPIC layers, we can just add the diabat3 sources, which are 
           * also on the layers.
           *
           * [SU] = (m/s^2)
           * [SV] = (m/s^2)
           * Also, [DUDT and DVDT] = (m/s^2)
           *
           ********************************************************************/
          DUDT(grid.it_uv_tend,K,J,I) += DB_OUTS(SU,K,J,I); 
          DVDT(grid.it_uv_tend,K,J,I) += DB_OUTS(SV,K,J,I); 
          /*fprintf(stderr,"(%d,%d,%d) nmt_apply_sources: SU=%.16f,SV=%.16f\n", \
            K,J,I,DB_OUTS(SU,K,J,I),DB_OUTS(SV,K,J,I));*/

          /**********************************************************************
           * pass diabat3 total cloud water source term to EPIC's Q
           *
           * SRT and srt are defined on the Kth EPIC layer. All other variables
           * are defined on the bottom interface of the Kth EPIC layer. Thus,
           * we must average and extrapolate to get SRT, srt on the interfaces.
           *
           * [SRT] = [srt] = (kg/kg/s)
           * [grid.dt] = (s)
           * Therefore, [srt*grid.dt] = (kg/kg), which are the units of epic_q.
           *
           **********************************************************************/
          if ( K == KHI){
            /************************************************************** 
             * if at the bottom layer, extrapolate to the bottom interface, 
             * using get_var's formula 
             **************************************************************/
            srt = DB_OUTS(SRT,K,J,I)+(DB_OUTS(SRT,K,J,I)-DB_OUTS(SRT,K-1,J,I))* \
              (grid.sigmatheta[kk+1]-grid.sigmatheta[kk])*grid.dsgth_inv[kk-1];
          }else{
            /****************************************************************** 
             * otherwise, take a simple average of the current and lower layers 
             ******************************************************************/
            srt = (DB_OUTS(SRT,K,J,I) + DB_OUTS(SRT,K+1,J,I))*.5;
          }
          /***************************************************************
           * convert diabat3 source into an increment of total cloud water
           ***************************************************************/
          srt *= grid.dt;
          /**************************************** 
           * set new mix. ratio value to EPIC macro
           ****************************************/
          Q(H_2O_INDEX,VAPOR,K,J,I) += srt;
          /*fprintf(stderr,"(%d,%d,%d) nmt_apply_sources: srt = %.16f\n",\
            K,J,I,srt);*/
        }
      }
    }
  }
}

/**********************************************************************
 * nmt_initialize - initialize some variables for nmt physics
 **********************************************************************/
void nmt_init_diag(){
  if (grid.nmt_physics_on == 1){
    nmt_entropy();
  }
}

/**********************************************************************
 * nmt_init_water - initialize total cloud water for nmt physics
 **********************************************************************/
void nmt_init_water(){
  if (grid.nmt_physics_on == 1){
    int
      K,J,I;
    EPIC_FLOAT
      rel_hum,
      rsat, rt;
    
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          rsat = rs(THETA(K,J,I),EXNER3(K,J,I));
          rel_hum = 0.80;
          rel_hum += rel_hum*(.05*myrand());
          rt = rel_hum*rsat;
          Q(H_2O_INDEX,VAPOR,K,J,I) = rt;
          /*fprintf(stderr,"Q(%d,%d,%d) = %.16f\n",\
            K,J,I,Q(H_2O_INDEX,VAPOR,K,J,I)); */
        }
      }
    }
  }
}

/**********************************************************************
 * nmt_randomize - add random noise to a variable field
 **********************************************************************/
void nmt_randomize(int index){
  if (grid.nmt_physics_on == 1){
    int K,J,I;

    for (J = JLO; J <= JHI; J++){
      for (I = ILO; I <= IHI; I++){
        if (index == T3_INDEX){
          T3(KHI,J,I) += T3(KHI,J,I)*(.007*myrand());
        }
        for (K = KLO; K <= KHI; K++){        }
        ;
      }
    }
  }   
}  

/**********************************************************************
 * nmt_entropy - calculate entropies for each timestep
 **********************************************************************/
void nmt_entropy(){
  if (grid.nmt_physics_on == 1){
    int K,J,I,kk;
    TFLOAT    rtot;                   /* total cloud water mix. ratio        */
    TFLOAT    rsat;                   /* saturation mixing ratio for entropy */
    TFLOAT    rvap;                   /* vapor mixing ratio for entropy      */
    TFLOAT    latent = LVAPOR+LFREEZ; /* sum of latent heats for entropy     */
    
    for (K = KLO; K <= KHI; K++){
      kk = 2*K;
      for (J = JLO; J <= JHI; J++){
        for (I = ILO; I <= IHI; I++){

        /* grab some values for computing the entropies */
        rtot = get_var(planet,H_2O_INDEX,VAPOR,grid.it_h,kk,J,I);
/*NOTE: Need to convert to EXNER3
        rsat = rs(THETA2(K,J,I),EXNER(K,J,I));
*/
        rvap = MIN(rtot,rsat);
        /*fprintf(stderr,"theta2(%d,%d,%d) = %.16f\n",K,J,I,THETA2(K,J,I));*/
        
        /* grab diagnostics to save to extract file */
        DRY_ENTROPY(K,J,I) = CP*log(THETA2(K,J,I)/FREEZING);
        /*fprintf(stderr,"theta2(%d,%d,%d) = %.16f\n",K,J,I,THETA2(K,J,I));*/
        /*fprintf(stderr,"theta2(%d,%d,%d) = %.16f\n",K,J,I,THETA2(K,J,I));*/
        /*fprintf(stderr,"dry entropy(%d,%d,%d) = %.16f\n",K,J,I,DRY_ENTROPY(K,J,I));*/
        MOIST_ENTROPY(K,J,I) = DRY_ENTROPY(K,J,I)+(latent*rvap)/FREEZING;
        SAT_MOIST_ENTROPY(K,J,I) = DRY_ENTROPY(K,J,I)+(latent*rsat)/FREEZING;
        /*fprintf(stderr,"time=%d\n",(int)TIME);*/
        if ((K == KHI || K == KHI-1) && (int)TIME%86400 < 0.5){
          if (J == (int)(JHI-JLO)/2){
            if (I == (int)(IHI-ILO)/2){
              fprintf(stderr,"Sd=%.3f, S=%.3f, Ss=%.3f, rt=%.10f, rs=%.10f, rv=%.10f, theta=%.3f, CP=%.1f, FREEZ=%.2f, L=%f\n",DRY_ENTROPY(K,J,I),MOIST_ENTROPY(K,J,I),SAT_MOIST_ENTROPY(K,J,I),rtot,rsat,rvap,THETA2(K,J,I),CP,FREEZING,latent);
            }
          }
        }
        }
      }
    }
  }
}

/**********************************************************************
 * nmt_memset() - zero all the nmt output diagnostics
 **********************************************************************/
void nmt_memset(){
  if(grid.nmt_physics_on == 1){
    memset(var.dry_entropy.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.moist_entropy.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.sat_moist_entropy.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.the_flux.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.rt_flux.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.u_flux.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.v_flux.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.convthrot.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.fluxthrot.value,0,Nelem3d*sizeof(EPIC_FLOAT));
    memset(var.rain_rate.value,0,Nelem3d*sizeof(EPIC_FLOAT));
  }
}
