/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2002-2009 Timothy E. Dowling                      *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/GNU_General_Public_License.                        *
 *                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
 *                                                                 *
 * You should have received a copy of the GNU General Public       *
 * License along with this program; if not, write to the Free      *
 * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
 * Boston, MA  02111-1307, USA.                                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * gen_t_vs_p.c  * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Merge the Galileo Probe T(p) data of Seiff et al with the       *
 * orton.dat profile, to create the file t_vs_p.jupiter.orton for  *
 * use in the EPIC model.  Smooth with a boxcar average, and       *
 * record these details in the file header.                        *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* 
 * Compile with:                                             
 * cc -I$EPIC_PATH/include -DEPIC_PATH=\"$EPIC_PATH\" -lm -o gen_t_vs_p gen_t_vs_p.c
 */

#include <stdio.h>
#include <math.h>

/*
 * Boxcar smoothing in range [p/(1+ALPHA/100.), p*(1+ALPHA/100.)].
 */
#define ALPHA 50.
/*
 * Seam between data files [hPa].
 */
#define SEAM_P 0.000719

/*
 * Function prototypes:
 */
double *dvector(int nl, int nh);
void  free_dvector(double *m, int nl, int nh);

main(int   argc,
     char *argv[]) 
{
  char   
    header[128];
  int
    kk,k,
    orton_ntp,orton_ceiling,
    seiff_ntp,seiff_floor,
    ceiling,floor,
    ntp;
  double
    *orton_t,*orton_p,
    *seiff_t,*seiff_p,
    *tdat,*pdat,
    *smooth_t,
    *count,
     pup,pdn,
     a,b,
     tmp;
  FILE 
    *orton,
    *seiff,
    *t_vs_p;

  /* Open orton.dat: */
  orton = fopen(EPIC_PATH"/data/jupiter/archive/orton.dat","r");

  /* Skip over header: */
  for (kk = 0; kk < 27 ; kk++) {
    fgets(header,100,orton);  
  }
 
  /* 
   * Allocate arrays for orton.dat. There are 131 points in the original
   * file, plus we have added some T(p) points deeper than 100 bars.
   */
  orton_ntp = 131+7;
  orton_t   = dvector(0,orton_ntp-1);
  orton_p   = dvector(0,orton_ntp-1);

  /* Read data and reorder to get increasing p: */
  for (kk = orton_ntp-1; kk >= 0; kk--) {
    fscanf(orton,"%lf %lf %lf %lf %lf %lf %lf", 
           &tmp,orton_p+kk,orton_t+kk,&tmp,&tmp,&tmp,&tmp);
    /* Convert to hPa: */
    orton_p[kk] *= 1.e+3; 
  }
  fclose(orton);

  /* Open seiff.dat: */
  seiff = fopen(EPIC_PATH"/data/jupiter/archive/seiff.dat","r");

  /* Skip over header: */
  for (kk = 0; kk < 5 ; kk++) {
    fgets(header,100,seiff);  
  }

  /* Read size of file: */
  fscanf(seiff,"%d",&seiff_ntp);
  fgets(header,100,seiff);
  fgets(header,100,seiff);

  /* Allocate arrays for seiff.dat: */
  seiff_t  = dvector(0,seiff_ntp-1);
  seiff_p  = dvector(0,seiff_ntp-1);

  for (kk = 0; kk < seiff_ntp; kk++) {
    fscanf(seiff,"%lf %lf %lf %lf",
                 &tmp,seiff_p+kk,seiff_t+kk,&tmp);
    /* Seiff pressures are in hPa. */
  }
  fclose(seiff);

  /* 
   * Determine orton_ceiling, the index of the first point with p larger
   * than SEAM_P, and the corresponding seiff_floor.
   */
  orton_ceiling = 0;
  kk            = 0;
  while (orton_p[kk] < SEAM_P) {
    kk++;
  }
  orton_ceiling = kk;

  seiff_floor = seiff_ntp-1;
  kk          = seiff_ntp-1;
  while (seiff_p[kk] >= SEAM_P) {
    kk--;
  }
  seiff_floor = kk;

  /*
   * Concatenate the two data arrays into a single array.
   * Allocate memory.
   */
  ntp      = seiff_floor+1+orton_ntp-orton_ceiling;
  pdat     = dvector(0,ntp-1);
  tdat     = dvector(0,ntp-1);
  smooth_t = dvector(0,ntp-1);
  count    = dvector(0,ntp-1);

  for (kk = 0; kk <= seiff_floor; kk++) {
    pdat[kk] = seiff_p[kk];
    tdat[kk] = seiff_t[kk];
  }
  for (k = orton_ceiling; k < orton_ntp; k++) {
    pdat[kk] = orton_p[k];
    tdat[kk] = orton_t[k];
    kk++;
  }

  /* 
   * Smooth temperature data using a boxcar averaging
   * between P*(1+ALPHA/100.) and P/(1+ALPHA/100.).
   */
  ceiling = 0;
  floor   = ntp-1;
  for (kk = 0; kk < ntp; kk++) {
    smooth_t[kk] = 0.;
    count[kk]    = 0.;
    pup          = pdat[kk]/(1.+ALPHA/100.);
    pdn          = pdat[kk]*(1.+ALPHA/100.);
    if (pup < pdat[0]) {
      ceiling = (ceiling > kk) ? ceiling : kk;
    }
    else if (pdn > pdat[ntp-1]) {
      floor = (floor < kk) ? floor : kk;
    }
    else {
      for (k = 0; k < ntp; k++) {
        if (pdat[k] >= pup && pdat[k] <= pdn) {
          count[kk]    += 1.;
          smooth_t[kk] += tdat[k];
        }
      }
      smooth_t[kk] /= count[kk];
    }
  }

  /*
   * Extrapolate data at bottom end where averaging is not centered.
   */
  a = (tdat[floor-2]-tdat[floor-1])/(log(pdat[floor-2])-log(pdat[floor-1]));
  b = tdat[floor-1];
  for (kk = floor; kk < ntp; kk++) {
    smooth_t[kk] = a*(log(pdat[kk])-log(pdat[floor-1]))+b;
  }


  /* Open t_vs_p.jupiter: */
  t_vs_p = fopen(EPIC_PATH"/data/jupiter/t_vs_p.jupiter.orton","w");

  /* Write preamble: */
  fprintf(t_vs_p," Temperature versus pressure for Jupiter, generated by epic/data/jupiter/gen_t_vs_p.c.\n");
  fprintf(t_vs_p," Data for p <= %.3e hPa from Seiff et al, 1997, ``Thermal structure"
                 " of Jupiter's upper atmosphere derived from the \n",SEAM_P);
  fprintf(t_vs_p," Galileo Probe,'' Science 276, 102-104. "
                 " Data for p > %.3e hPa from orton.dat. \n",SEAM_P);
  fprintf(t_vs_p," Data have been smoothed with a boxcar average");
  fprintf(t_vs_p," in the range [p/(1+alpha),p*(1+alpha)], with alpha = %.2f.\n",
                   ALPHA/100.);
  fprintf(t_vs_p,"\n");
  fprintf(t_vs_p,"#   p[hPa]      T[K] \n");

  /* 
   * Write number of data points.
   * Cut off top end where averaging is not centered.
   */
  fprintf(t_vs_p,"%d\n",ntp-(ceiling+1));

  /* Write data: */
  for (kk = ceiling+1; kk < ntp; kk++) {
    fprintf(t_vs_p,"   %10.3e  %7.2f\n",pdat[kk],smooth_t[kk]);
  }
  fclose(t_vs_p);

  /* Free allocated memory: */
  free_dvector(smooth_t,0,ntp-1      );
  free_dvector(pdat,    0,ntp-1      );
  free_dvector(tdat,    0,ntp-1      );
  free_dvector(seiff_t, 0,seiff_ntp-1);
  free_dvector(seiff_p, 0,seiff_ntp-1);
  free_dvector(orton_p, 0,orton_ntp-1);
  free_dvector(orton_t, 0,orton_ntp-1);

  return;
}


/*======================= dvector() ============================================*/

double *dvector(int nl, int nh)
      /*
       *  Allocates memory for a 1D double array 
       *  with range [nl..nh].
       */
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  double         
    *m;

  if (nh < nl) {
    fprintf(stderr,"called dvector(%d,%d) \n",nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe - nl_safe + 1);

  m = (double *)calloc(len_safe, sizeof(double));
  if (!m) {
    fprintf(stderr, "calloc error in dvector \n");
    exit(1);
  }
  m -= nl_safe;
  return m;
}

/*======================= end of dvector() ====================================*/

/*======================= free_dvector() ======================================*/

void  free_dvector(double *m, int nl, int nh)
      /*
       *  Frees memory allocated by dvector().
       */
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);
}

/*======================= end of free_dvector() ===============================*/

/* * * * * * * * * * end of gen_t_vs_p.c * * * * * * * * * * * * * */
