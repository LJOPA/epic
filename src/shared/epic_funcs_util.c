/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * All material in this file is covered by the                     *
 * GNU General Public License.                                     *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling, except as noted.    *
 *                                                                 *
 * The fcmp() function is derived from the fcmp() function that is *
 * Copyright (c) 1998-2000 Theodore C. Belding                     *
 * University of Michigan Center for the Study of Complex Systems. *
 *                                                                 *
 * The monotonic spline, spline_pchip(), is derived from the PCHIP *
 * (Piecewise Cubic Hermite Interpolating Polynomial) code written *
 * by Fred Fritsch, which was tranlated into C by John Burkardt,   *
 * and modified for use here by Tim Dowling.                       *
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

/* * * * * * * * *  epic_funcs_util.c  * * * * * * * * * * * * * * * * * * * * 
 *                                                                           *
 *       These functions for math and utility operations should  not         *
 *       reference EPIC variables by name or index, such that they can be    *
 *       used outside the model.                                             *
 *                                                                           *
 *       NOTE: One dependency with the EPIC model is the environment         *
 *             variable EPIC_PRECISION, which specifies the floating-point   *
 *             precision.                                                    *
 *                                                                           *
 *       This file includes the following:                                   *
 *                                                                           *
 *           ivector(),free_ivector()                                        *
 *           fvector(),free_fvector()                                        *
 *           dvector(),free_dvector()                                        *
 *           ftriplet(),free_ftriplet()                                      *
 *           spline(),splint()                                               *
 *           linint()                                                        *
 *           spline_pchip(), splint_pchip()                                  *
 *           pchst()                                                         *
 *           lagrange_interp()                                               *
 *           gamma_nr()                                                      *
 *           sech2()                                                         *
 *           normed_legendre()                                               *
 *           machine_epsilon()                                               *
 *           find_root()                                                     *
 *           broyden_root()                                                  *
 *           global_step()                                                   *
 *           line_search()                                                   *
 *           dogleg_driver()                                                 *
 *           dogleg_step()                                                   *
 *           trust_region()                                                  *
 *           qr_decompose()                                                  *
 *           qr_update()                                                     *
 *           qr_rotate()                                                     *
 *           lu_decompose()                                                  *
 *           lu_backsub()                                                    *
 *           lu_improve()                                                    *
 *           find_place_in_table()                                           *
 *           tridiag()                                                       *
 *           band_decomp()                                                   *
 *           band_back_sub()                                                 *
 *           band_multiply()                                                 *
 *           band_improve()                                                  *
 *           poly_interp()                                                   *
 *           nth_trapezoidal()                                               *
 *           romberg_integral()                                              *
 *           crank_nicolson()                                                *
 *           hqr()                                                           *
 *           quicksort()                                                     *
 *           swap()                                                          *
 *           four1(),realft()                                                *
 *           c_num(),c_mult(),c_add(),c_sub()                                *
 *           c_exp(),c_abs(),c_real(),c_imag                                 *
 *           fcmp()                                                          *
 *           least_squares()                                                 *
 *           savitzky_golay()                                                *
 *           random_number()                                                 *
 *           lat_centric_to_graphic(), lat_graphic_to_centric()              *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic_funcs_util.h>

/*
 * The following macro is useful to insert into code while trying to corral a problem.
 */
#undef  DEBUG_MILESTONE
#define DEBUG_MILESTONE(comment) fprintf(stderr,"%s, %2d: "#comment"\n", \
                                         dbmsname,++idbms);fflush(stderr);


/*======================= ivector() ============================================*/
      
/*
 *  Allocates memory for a 1D int array 
 *  with range [nl..nh].
 */

#undef  DEBUG

int *ivector(int  nl, 
             int  nh,
             char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  int         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="ivector";

  if (nh < nl) {
    fprintf(stderr,"**error:%s, called by %s, range (%d,%d)\n",
                    dbmsname,calling_func,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (int *)calloc(len_safe,sizeof(int));
  if (!m) {
    fprintf(stderr,"**error:%s, called by %s\n",dbmsname,calling_func);
    exit(1);
  }
  m -= nl_safe;

#if defined(DEBUG)
  fprintf(stderr,"ivector() called by %s \n",calling_func);
#endif

  return m;
}

/*======================= end of ivector() ====================================*/

/*======================= free_ivector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by ivector().
 */

void free_ivector(int  *m, 
                  int  nl, 
                  int  nh,
                  char *calling_func)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

#if defined(DEBUG)
  fprintf(stderr,"free_ivector() called by %s \n",calling_func);
#endif

  return;
}

/*======================= end of free_ivector() ================================*/

/*======================= fvector() ============================================*/
      
/*
 *  Allocates memory for a 1D FLOAT array 
 *  with range [nl..nh].
 */

#undef  DEBUG 

FLOAT *fvector(int  nl, 
               int  nh,
               char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  FLOAT         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="fvector";

#if defined(DEBUG)
  fprintf(stderr,"fvector() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    fprintf(stderr,"**error:%s, called by %s, range (%d,%d)\n",
                   dbmsname,calling_func,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (FLOAT *)calloc(len_safe,sizeof(FLOAT));

  if (!m) {
    fprintf(stderr,"**error:%s, called by %s, nl=%d,nh=%d,len_safe=%d\n",
            dbmsname,calling_func,nl,nh,len_safe);
    exit(1);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of fvector() ====================================*/

/*======================= free_fvector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by fvector().
 */

void free_fvector(FLOAT *m, 
                  int     nl, 
                  int     nh,
                  char   *calling_func)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

#if defined(DEBUG)
  fprintf(stderr,"free_fvector() called by %s \n",calling_func);
#endif

  return;
}

/*======================= end of free_fvector() ===============================*/

/*======================= dvector() ============================================*/
      
/*
 *  Allocates memory for a 1D double array 
 *  with range [nl..nh].
 */

#undef  DEBUG 

double *dvector(int  nl, 
                int  nh,
                char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  double         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dvector";

#if defined(DEBUG)
  fprintf(stderr,"dvector() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    fprintf(stderr,"**error:%s, called by %s, range (%d,%d)\n",
                   dbmsname,calling_func,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (double *)calloc(len_safe,sizeof(double));

  if (!m) {
    fprintf(stderr,"**error:%s, called by %s, nl=%d,nh=%d,len_safe=%d\n",
            dbmsname,calling_func,nl,nh,len_safe);
    exit(1);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of dvector() ====================================*/

/*======================= free_dvector() ======================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by dvector().
 */

void free_dvector(double *m, 
                  int     nl, 
                  int     nh,
                  char   *calling_func)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

#if defined(DEBUG)
  fprintf(stderr,"free_dvector() called by %s \n",calling_func);
#endif

  return;
}

/*======================= end of free_dvector() ===============================*/

/*======================= ftriplet() ==========================================*/
      
/*
 *  Allocates memory for a 1D float_triplet array 
 *  with range [nl..nh].
 */

/*
 * Define DEBUG here to enable this state.
 */
#undef  DEBUG 

float_triplet *ftriplet(int  nl, 
                        int  nh,
                        char *calling_func)
{
  unsigned int  
    len_safe;
  int           
    nl_safe, nh_safe;
  float_triplet         
    *m;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="ftriplet";

#if defined(DEBUG)
  fprintf(stderr,"ftriplet() called by %s \n",calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    fprintf(stderr,"**error:%s, called by %s, range (%d,%d)\n",
                   dbmsname,calling_func,nl,nh);
    exit(1);
  }

  nl_safe  = (nl < 0) ? nl : 0;
  nh_safe  = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe-nl_safe+1);

  m = (float_triplet *)calloc(len_safe,sizeof(float_triplet));

  if (!m) {
    fprintf(stderr,"**error:%s, called by %s\n",dbmsname,calling_func);
    exit(1);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of ftriplet() ===================================*/

/*======================= free_ftriplet() =====================================*/

#undef  DEBUG

/*
 *  Frees memory allocated by ftriplet().
 */

void free_ftriplet(float_triplet *m, 
                  int             nl, 
                  int             nh,
                  char           *calling_func)
{
  int  
    nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

#if defined(DEBUG)
  fprintf(stderr,"free_ftriplet() called by %s \n",calling_func);
#endif

  return;
}

/*======================= end of free_ftriplet() ==============================*/

/*======================= spline() ============================================*/

/* 
 * Cubic spline routine. 
 * Adapted from Numerical Recipes in C, p. 96. 
 * The tridiagonal system is solved with pivoting
 * to avoid numerical instability.
 * Assumes zero-offset arrays.
 *
 * NOTE: We stripe the data into one array with float_triplet to get
 *       a cache-aware memory layout.
 */

void spline(int            n,
            float_triplet *table,
            FLOAT          y1_bot, 
            FLOAT          y1_top)
{
  int     
    j,jm1,jp1;
  FLOAT 
     dx_a,dx_b,dx_c,  
    *a,
    *b,
    *c,
    *r,
    *y2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="spline";

  /*
   * Check validity of n.
   */
  if (n < 3) {
    fprintf(stderr,"**error:%s, n=%d < 3\n",dbmsname,n);
    exit(1);
  }
  
  /*
   * Allocate memory.
   */
  a  = fvector(0,n-1,dbmsname);
  b  = fvector(0,n-1,dbmsname);
  c  = fvector(0,n-1,dbmsname);
  r  = fvector(0,n-1,dbmsname);
  y2 = fvector(0,n-1,dbmsname);

  for (j = 0; j < n; j++) {
    y2[j] = (table+j)->z;
  }

  for (j = 0; j < n; j++) {
    jm1 = j-1;
    jp1 = j+1;
    if (jm1 >= 0) {
      dx_a = (table+j)->x - (table+jm1)->x; 
      if (dx_a <= 0.) {
        fprintf(stderr,"**error:%s, x[%d]=%g x[%d]=%g\n",
                        dbmsname,jm1,(table+jm1)->x,j,(table+j)->x);
        exit(1);
      }
      if (jp1 < n) {
        dx_b = (table+jp1)->x - (table+jm1)->x;
        if (dx_b <= 0.) {
          fprintf(stderr,"**error:%s, x[%d]=%g x[%d]=%g\n",
                          dbmsname,jm1,(table+jm1)->x,jp1,(table+jp1)->x);
          exit(1);
        }
      }
    }
    if (jp1 < n) {
      dx_c = (table+jp1)->x - (table+j)->x;
      if (dx_c <= 0.) {
        fprintf(stderr,"**error:%s, x[%d]=%g x[%d]= %g\n",
                        dbmsname,j,(table+j)->x,jp1,(table+jp1)->x);
        exit(1);
      }
    }

    if (j == 0) {
      if (y1_bot > 0.99e+30) {
        y2[j] = 0.;
      }
      else {
        b[j] = dx_c/3.;
        c[j] = dx_c/6.;
        r[j] = ((table+jp1)->y - (table+j)->y)/dx_c-y1_bot;
      }
    }
    else if (j == n-1) {
      if (y1_top > 0.99e+30) {
        y2[j] = 0.;
      }
      else {
        a[j] = dx_a/6.;
        b[j] = dx_a/3.;
        r[j] = y1_top-((table+j)->y - (table+jm1)->y)/dx_a;
      }
    }
    else {
      a[j] = dx_a/6.;
      b[j] = dx_b/3.;
      c[j] = dx_c/6.;
      r[j] = ((table+jp1)->y - (table+j)->y)/dx_c - ((table+j)->y - (table+jm1)->y)/dx_a;
    }
  }

  if (y1_bot > 0.99e+30) {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 on both ends. */
      tridiag(n-2,a+1,b+1,c+1,r+1,y2+1,WITH_PIVOTING);
    }
    else {
      /* y2 = 0 at start. */
      tridiag(n-1,a+1,b+1,c+1,r+1,y2+1,WITH_PIVOTING);
    }
  }
  else {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 at end. */
      tridiag(n-1,a,b,c,r,y2,WITH_PIVOTING);
    }
    else {
      /* y2 needed at both ends */
      tridiag(n,a,b,c,r,y2,WITH_PIVOTING);
    }
  }

  for (j = 0; j < n; j++) {
    (table+j)->z = y2[j];
  }

  /* 
   * Free allocated memory.
   */
  free_fvector(y2,0,n-1,dbmsname);
  free_fvector(r, 0,n-1,dbmsname);
  free_fvector(c, 0,n-1,dbmsname);
  free_fvector(b, 0,n-1,dbmsname);
  free_fvector(a, 0,n-1,dbmsname);
  
  return;
}

/*======================= end of spline() ===================================*/

/*======================= splint() ==========================================*/

/*  
 *  Evaluates cubic-spline interpolations.
 *  This version assumes you have already found the correct position
 *  in the tables, unlike the Numerical Recipes version. 
 *  The function find_place_in_table() may be used to find the position.
 *
 *  NOTE: We stripe the data into one array with float_triplet to get
 *        a cache-aware memory layout.
 */

inline FLOAT splint(register FLOAT          xx, 
                    register float_triplet *table,
                    register FLOAT          dx)
{
  register FLOAT  
    a, 
    b,
    ans;

  a = ( (table+1)->x - xx       )/dx;
  b = (     xx       - table->x )/dx;

  ans = a*table->y + b*(table+1)->y +( (a*a*a-a)*(table  )->z
                                      +(b*b*b-b)*(table+1)->z )*dx*dx/6;
    
  return ans;
}

/*======================= end of splint() ===================================*/

/*======================= linint() ==========================================*/

/*
 * Evaluates linear interpolation using same arguments as splint(),
 * to make it easy to switch between the two.
 */

FLOAT linint(register FLOAT          xx,
                    register float_triplet *table,
                    register FLOAT          dx)
{
  FLOAT
    ans;

  ans = table->y + ((table+1)->y - table->y)*(xx - table->x)/dx;
  
  return ans;
}

/*======================= end of linint() ===================================*/

/*=========================== spline_pchip() ================================*/

/*
 *  Purpose:
 *
 *    Sets derivatives for a piecewise cubic Hermite interpolant.
 *
 *  Discussion:
 *
 *    This routine computes what would normally be called a Hermite 
 *    interpolant.  However, the user is only required to supply function
 *    values, not derivative values as well.  This routine computes
 *    "suitable" derivative values, so that the resulting Hermite interpolant
 *    has desirable shape and monotonicity properties.
 *
 *    The interpolant will have an extremum at each point where
 *    monotonicity switches direction.
 *
 *    The resulting piecewise cubic Hermite function may be evaluated
 *    by splint_pchip().
 *
 *    This routine was originally called "PCHIM".
 *
 *    The acronym PCHIP means Piecewise Cubic Hermite Interpolating Polynomial.
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Modified for use in the EPIC model by Tim Dowling.
 *
 *  References:
 *
 *    Fred Fritsch, Ralph Carlson,
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *    Fred Fritsch, Judy Butland,
 *    A Method for Constructing Local Monotone Piecewise 
 *    Cubic Interpolants,
 *    SIAM Journal on Scientific and Statistical Computing,
 *    Volume 5, Number 2, 1984, pages 300-304.
 *
 *  Parameters:
 *
 *    Input, int n, the number of data points; n must be at least 2.
 *
 *    Input, 
 *      float_triplet table[n]:
 *      table[n].x are the strictly increasing independent variable values.
 *
 *      table[n].y are the dependent variable values to be interpolated. 
 *      This routine is designed for monotonic data, but it will work for
 *      any data. It will force extrema at points where monotonicity switches
 *      direction.
 *
 *    Output, 
 *      table[n].z are the derivative values at the data points.  If the
 *      data are monotonic, these values will determine a monotone cubic
 *      Hermite function.
 */

void spline_pchip(int            n,
                  float_triplet *table)
{
  int
    i,nless1;
  register FLOAT
    del1,del2,
    dmax,dmin,
    drat1,drat2,
    h1,h2,hsum,hsumt3,
    w1,w2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="spline_pchip";

  /*
   * Check the arguments.
   */
  if ( n < 2 ) {
    fprintf(stderr,"**error:%s: number of data points, n=%d, is less than 2.\n",
                    dbmsname,n);
    exit (-1);
  }

  for ( i = 1; i < n; i++ ) {
    if ( table[i].x <= table[i-1].x ) {
      fprintf(stderr,"**error:%s: x array is not strictly increasing, table[%d].x=%g <= table[%d].x=%g\n",
                      dbmsname,i,table[i].x,i-1,table[i-1].x);
      exit (-3);
    }
  }

  nless1 = n - 1;
  h1     =   table[1].x - table[0].x;
  del1   = ( table[1].y - table[0].y ) / h1;

  /*
   * Special case n = 2, use linear interpolation.
   */
  if ( n == 2 ) {
    table[0  ].z = del1;
    table[n-1].z = del1;
    return;
  }

  /*
   *  Normal case, 3 <= n.
   */
  h2   =   table[2].x - table[1].x;
  del2 = ( table[2].y - table[1].y ) / h2;

  /*
   * Set table[0].z via non-centered three point formula, adjusted to be shape preserving.
   */
  hsum       = h1 + h2;
  w1         = ( h1 + hsum ) / hsum;
  w2         = -h1 / hsum;
  table[0].z = w1 * del1 + w2 * del2;

  if ( pchst( table[0].z, del1 ) <= 0.0 ) {
    table[0].z = 0.0;
  }
  else if ( pchst( del1, del2 ) < 0.0 ) {
    /*
     *  Need to do this check only if monotonicity switches.
     */
    dmax = 3.0 * del1;

    if ( fabs( dmax ) < fabs( table[0].z ) ) {
      table[0].z = dmax;
    }
  }

  /*
   *  Loop through interior points.
   */
  for ( i = 2; i <= nless1; i++ ) {
    if ( 2 < i ) {
      h1   = h2;
      h2   = table[i].x - table[i-1].x;
      hsum = h1 + h2;
      del1 = del2;
      del2 = ( table[i].y - table[i-1].y ) / h2;
    }

    /*
     *  Set table[i-1].z = 0 unless data are strictly monotonic.
     */
    table[i-1].z = 0.0;
    if ( pchst( del1, del2 ) > 0.0 ) {
     /*
      *  Use Brodlie modification of Butland formula.
      */
      hsumt3       = 3.0 * hsum;
      w1           = ( hsum + h1 ) / hsumt3;
      w2           = ( hsum + h2 ) / hsumt3;
      dmax         = MAX ( fabs ( del1 ), fabs ( del2 ) );
      dmin         = MIN ( fabs ( del1 ), fabs ( del2 ) );
      drat1        = del1 / dmax;
      drat2        = del2 / dmax;
      table[i-1].z = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }

  /*
   *  Set table[n-1].z via non-centered three point formula, adjusted to be
   *  shape preserving.
   */
  w1           = -h2 / hsum;
  w2           = ( h2 + hsum ) / hsum;
  table[n-1].z = w1 * del1 + w2 * del2;

  if ( pchst( table[n-1].z, del2 ) <= 0.0 ) {
    table[n-1].z = 0.0;
  }
  else if ( pchst( del1, del2 ) < 0.0 ) {
    /*
     *  Need to do this check only if monotonicity switches.
     */
    dmax = 3.0 * del2;

    if ( fabs( dmax ) < fabs( table[n-1].z ) ) {
      table[n-1].z = dmax;
    }
  }

  return;
}

/*========================= end of spline_pchip() ==============================*/

/*========================= splint_pchip() =====================================*/

/*
 *  Purpose:
 *
 *    Evaluates a cubic polynomial given in Hermite form.
 *
 *  Discussion:
 *
 *    This routine evaluates a cubic polynomial given in Hermite form, and
 *    is designed to work with spline_pchip().
 *
 *    The cubic polynomial is determined by function values
 *    table[0].y, table[1].y and derivatives table[0].z, table[1].z
 *    on the interval [table[0].x,table[1].x].
 *
 *    This routine was originally called "CHFEV".
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Modified for use in the EPIC model by Tim Dowling.
 *
 *  References:
 *
 *    Fred Fritsch, Ralph Carlson, 
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *    David Kahaner, Cleve Moler, Steven Nash,
 *    Numerical Methods and Software,
 *    Prentice Hall, 1989,
 *    ISBN: 0-13-627258-4,
 *    LC: TA345.K34.
 *
 *  Parameters:
 *
 *    Input, FLOAT xx, the position of the interpolation
 *
 *    Input, float_triplet *table:
 *       The use of the data type float_triplet facilitates
 *       cache-aware memory allocation. 
 *
 *       NOTE: This function does not search the table,
 *             but assumes xx is between the given inputs as
 *             returned by find_place_in_table().
 *
 *       table[0].x and table[1].x are the endpoints of the interval of
 *       definition of the cubic; they must be distinct.
 *
 *       table[0].y and table[1].y are the corresponding values of the function.
 *
 *       table[0].z and table[1].z are the corresponding derivative values.
 *
 *    Input, h = dx, the x interval, as returned by find_place_in_table().
 *    This is included here to make it easy to switch with splint().
 *
 *    Returns the value of the cubic function at the point xx.
 */

FLOAT splint_pchip(FLOAT          xx,
                          float_triplet *table,
                          FLOAT          h)
{
  register FLOAT 
    c2,c3,
    del1,del2,delta,
    x,f;

  /*
   *  Compute cubic coefficients expanded about table[0].x.
   */
  delta =  ( table[1].y - table[0].y ) / h;
  del1  =  ( table[0].z - delta      ) / h;
  del2  =  ( table[1].z - delta      ) / h;
  c3    =  ( del1 + del2             ) / h;
  c2    = -( del1 + del1 + del2 );

  /*
   *  Evaluation.
   */
  x = xx - table[0].x;
  f = table[0].y + x * ( table[0].z + x * ( c2 + x * c3 ) );

  return f;
}

/*========================= end of splint_pchip() ==============================*/

/*========================= pchst() ============================================*/

/*
 *  Purpose:
 *
 *    Sign-testing routine.
 *
 *  Discussion:
 *
 *    This routine essentially computes the sign of arg1 * arg2.
 *
 *    The object is to do this without multiplying arg1 * arg2, to avoid
 *    possible over/underflow problems.
 *
 *  Modified:
 *
 *    25 August 2007
 *
 *  Author:
 *
 *    Fred Fritsch,
 *    Mathematics and Statistics Division,
 *    Lawrence Livermore National Laboratory.
 *
 *    C translation by John Burkardt.
 *
 *    Adapted to the EPIC model by Tim Dowling.
 *
 *  Reference:
 *
 *    Fred Fritsch, Ralph Carlson, 
 *    Monotone Piecewise Cubic Interpolation,
 *    SIAM Journal on Numerical Analysis,
 *    Volume 17, Number 2, April 1980, pages 238-246.
 *
 *  Parameters:
 *
 *    Input, FLOAT arg1, arg2, two values to check.
 *
 *    Output, FLOAT pchst,
 *    -1.0, if arg1 and arg2 are of opposite sign.
 *     0.0, if either argument is zero.
 *    +1.0, if arg1 and arg2 are of the same sign.
 */

extern inline FLOAT pchst(FLOAT arg1,
                          FLOAT arg2)
{
  FLOAT
    value;

  if ( arg1 == 0.0 ) {
    value = 0.0;
  }
  else if ( arg1 < 0.0 ) {
    if ( arg2 < 0.0 ) {
      value = 1.0;
    }
    else if ( arg2 == 0.0 ) {
      value = 0.0;
    }
    else if ( 0.0 < arg2 ) {
      value = -1.0;
    }
  }
  else if ( 0.0 < arg1 ) {
    if ( arg2 < 0.0 ) {
      value = -1.0;
    }
    else if ( arg2 == 0.0 ) {
      value = 0.0;
    }
    else if ( 0.0 < arg2 ) {
      value = 1.0;
    }
  }

  return value;
}

/*======================= end of pchst() ====================================*/

/*======================= lagrange_interp() =================================*/
/*
 * Returns the lagrangian interpolation of f(x).
 *
 * Input: "order" is the order of polynomial interpolation 
 *        (0=copy,1=linear,2=quadradic,3=cubic,etc.).
 *
 *        f[order+2] where f[0] is empty,
 *                         f[1] = f(x1),
 *                         f[2] = f(x2),
 *                         f[3] = f(x3), etc.
 *
 *        x[order+2] where x[0] is the point of interpolation, i.e., "x" in f(x),
 *                          x[1] = x1,
 *                          x[2] = x2,
 *                          x[3] = x3, etc.
 */

extern inline FLOAT lagrange_interp(register FLOAT *f,
                             register FLOAT *x,
                             register int order)
{
  register int
    jj, il;
  FLOAT 
    lp[ order+2 ];

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lagrange_interp";

  if ( order < 0 ) {
    fprintf(stderr,"**error:%s, order=%d is not in valid (must be >= 0)\n",dbmsname,order);
    exit(1);
  }

  for (jj=1; jj<=order+1; jj++) {
    lp[jj] = 1.;
    for (il=1; il<=order+1; il++) {
      if (il != jj)
      lp[jj] *= (x[0]-x[il])/(x[jj]-x[il]);
    }
  }

  f[0] = 0.0;
  for (jj=1; jj<=order+1; jj++) {
    f[0] += f[jj]*lp[jj];
  }

  return f[0];
}

/*==================== end of lagrange_interp() =============================*/

/*=============================== gamma_nr() ================================*/

/*
 * CJP 06/2003 *A*  
 * Based on "Numerical Recipes in C", 2nd Ed., p 213-214, Cambridge.
 */

FLOAT gamma_nr(FLOAT xx_input)
{
  register int
    j;
  register double 
    x,xx,
    y,tmp,ser,loggamma;
  static double 
    cof[6]={76.18009172947146,     -86.505320032941677,
            24.01409824083091,      -1.231739572450155,
             0.1208650973866179e-2  -0.5395239384953e-5};
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="gamma_nr";
    
  if (xx_input <= 0.) {
    fprintf(stderr,"%s: Negative or zero argument: xx=%e\n",dbmsname,xx_input);
  }  
  /*
   * Use double precision for internal calculations.
   */
  xx   = (double)xx_input;
  y    = x = xx;
  tmp  = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser  = 1.000000000190015;
  for  (j = 0; j <= 5; j++) {
    ser += cof[j]/++y;
  }
  loggamma = -tmp+log(2.5066282746310005*ser/x);

  return (FLOAT)exp(loggamma);
}

/*============================end of gamma_nr() =============================*/

/*======================= sech2() ===========================================*/

FLOAT sech2(FLOAT xx)
     /*
      *  Evaluates the square of the hyperbolic secant
      */
{
  FLOAT 
    a;

  a = 1./cosh(xx);
  return a*a;
}

/*======================= end of sech2() =====================================*/

/*======================= normed_legendre() ==================================*/

/*
 * Return the value of the normalized associated Legendre polynomial.
 * Used in spherical harmonics.
 * Adapted from Numerical Recipes in C, 2nd Ed., p. 254.
 */

FLOAT normed_legendre(int l,
                      int m,
                      FLOAT x)
{
  int
    i,ll;
  /* 
   * NOTE: The internal values need to be double precision
   *       to prevent roundoff problems.
   */
  double 
    factor1,factor2,
    somx2,
    pll,pmm,pmmp1;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="normed_legendre";

  /*
   * Check validity of input arguments.
   */
  if (m < 0 || m > l || fabs(x) > 1.) {
    fprintf(stderr,"**error:%s, bad argument l,m,x = %d %d %f\n",dbmsname,l,m,x);
    exit(1);
  }

  /*
   * Normalization convention is such that the integral of squared
   * spherical harmonics over the sphere yields 4.*M_PI, and there is
   * a 180 deg shift in longitude ( no (-1)^m factor) compared to the 
   * Numerical Recipes routine plgndr().
   */
  if (m == 0) {
    pmm = sqrt((2.*(FLOAT)l+1.));
  }
  else {
    pmm = sqrt(2.*(2.*(FLOAT)l+1.));
  }

  if (m > 0) {
    somx2   = sqrt((1.-x)*(1.+x));
    factor1 = 1.;
    factor2 = 1.+l-m;
    for (i = 1; i <= m; i++) {
      pmm     *= factor1*somx2/sqrt(factor2*(factor2+1));
      factor1 += 2.;
      factor2 += 2.;
    }
  }
  if (l == m) {
    return (FLOAT)pmm;
  }
  else {
    pmmp1 = x*(2*m+1)*pmm;
    if (l == (m+1)) {
      return (FLOAT)pmmp1;
    }
    else {
      for (ll = m+2; ll <= l; ll++) {
        pll   = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm   = pmmp1;
        pmmp1 = pll;
      }
      return (FLOAT)pll;
    }
  }
}

/*======================= end of normed_legendre() ===========================*/

/*======================= machine_epsilon() ==================================*/

/*
 * Calculate machine's smallest meaningful floating-point number.
 * From Dennis and Schnabel (1996), Algorithm A1.3.1.
 */

FLOAT machine_epsilon(void)
{
  static int
    initialized = FALSE;
  static FLOAT
    eps;

  if (!initialized) {
    eps = 1.;
    while (1.+eps != 1.) {
      eps *= .5;
    }
    eps *= 2.;

    initialized = TRUE;
  }

  return eps;
}

/*======================= end of machine_epsilon() ===========================*/

/*======================= find_root() ========================================*/

/*
 * Adapted from zridder(), Numerical Recipes in C, 2nd ed., p.358.
 * Returns 0 if root is found,
 *        -1 if fabs(func(x1)) <  fabs(func(x2)) and zero is not bracketed,
 *         1 if fabs(func(x2)) <= fabs(func(x1)) and zero is not bracketed.
 */

#undef  MAX_IT
#define MAX_IT 100

#undef  UNLIKELY_VAL
#define UNLIKELY_VAL -1.11111e+30

int find_root(FLOAT  x1,
              FLOAT  x2,
              FLOAT  xacc,
              FLOAT *x_root,
              FLOAT  (*func)(FLOAT))
{
  register int
    iter,
    compare;
  register FLOAT
    fh,fl,fm,fnew,
    s,xh,xl,xm,xnew;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="find_root";

  fl = (*func)(x1);
  fh = (*func)(x2);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl      = x1;
    xh      = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5*(xl+xh);
      fm = func(xm);
      s  = sqrt(fm*fm-fl*fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm+(xm-xl)*((fl > fh ? 1. : -1.)*fm/s);

      if (fabs(xnew-*x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew    = func(*x_root);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      }
      else {
        fprintf(stderr,"**error:%s, should never get here\n",dbmsname);
        exit(1);
      }
      if (fabs(xh-xl) <= xacc) {
        return 0;
      }
    }
    fprintf(stderr,"**error:%s, exceeded MAX_IT = %d, current root calc = %e\n", dbmsname, MAX_IT, *x_root);
    exit(1);
  }
  else {
    if (fl == 0.) {
      *x_root = x1;
      return 0;
    }
    if (fh == 0.) {
      *x_root = x2;
      return 0;
    }

    compare = fcmp(fabs(fl),fabs(fh));
    if (compare < 0) {
      return -1;
    }
    else {
      return 1;
    }
  }

  /* Should never get here. */
  fprintf(stderr,"**error:%s, should never reach here\n",dbmsname);
  exit(1);
}

/*======================= end of find_root() =================================*/

/*======================= broyden_root() =====================================*/

/*
 * Finds a vector root, x[], of the vector function vecfunc().
 * This globally convergent algorithm is adapted from 
 * Numerical Recipes in C, 2nd ed., p. 389-392, which is based on
 * Dennis and Schnabel (1996).
 *
 * Call with an initial guess for x[].
 * Assumes zero-based indexing.
 *
 * The function global_step() performs the globally-convergent step. 
 * The argument step_type specifies the type of step taken. Currently, the
 * valid choices are DS_LINE_STEP and DS_DOGLEG_STEP, the latter being more
 * sophisticated and reliable (the "DS" stands for Dennis and Schnabel 1996).
 *
 * Returns 0 on normal execution and an error code if
 * the routine has failed, has converged to a local minimum, or can make
 * no further progress, in which case one should retry with a
 * different initial guess.
 */

int broyden_root(int    n,
                 FLOAT *x,
                 void  (*vecfunc)(int,FLOAT *,FLOAT *),
                 FLOAT  tol_f,
                 int    max_it)
{
  int
    k,j,i,
    restart,singular,skip,
    num_bytes,
    old_max_taken,
    it              = 0,
    max_taken       = FALSE,
    count_max_taken = 0,
    status          = DS_X_ACCEPTED;
  static int
    nold = 0;
  FLOAT
    sum,denom,
    f,f_old,
    max_step,
    test,h,
    tmp,
    delta = -1;
  static FLOAT
    *c,*d,
    *x_old,
    *fvec,*fvec2,*fvec_old,
    *g,*sn,*s,*t,*w,
    *qt,*r;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="broyden_root";

  if (n > nold) {
    if (nold != 0) {
      /* 
       * Free previously allocated memory: 
       */
      free_fvector(r,       0,nold*nold-1,dbmsname);
      free_fvector(qt,      0,nold*nold-1,dbmsname);
      free_fvector(w,       0,nold-1,dbmsname);
      free_fvector(t,       0,nold-1,dbmsname);
      free_fvector(s,       0,nold-1,dbmsname);
      free_fvector(sn,      0,nold-1,dbmsname);
      free_fvector(g,       0,nold-1,dbmsname);
      free_fvector(fvec_old,0,nold-1,dbmsname);
      free_fvector(fvec2,   0,nold-1,dbmsname);
      free_fvector(fvec,    0,nold-1,dbmsname);
      free_fvector(x_old,   0,nold-1,dbmsname);
      free_fvector(d,       0,nold-1,dbmsname);
      free_fvector(c,       0,nold-1,dbmsname);
    }
    /*
     * Allocate memory: 
     */
    c        = fvector(0,n-1,dbmsname);
    d        = fvector(0,n-1,dbmsname);
    x_old    = fvector(0,n-1,dbmsname);
    fvec     = fvector(0,n-1,dbmsname);
    fvec2    = fvector(0,n-1,dbmsname);
    fvec_old = fvector(0,n-1,dbmsname);
    g        = fvector(0,n-1,dbmsname);
    sn       = fvector(0,n-1,dbmsname);
    s        = fvector(0,n-1,dbmsname);
    t        = fvector(0,n-1,dbmsname);
    w        = fvector(0,n-1,dbmsname);
    qt       = fvector(0,n*n-1,dbmsname);
    r        = fvector(0,n*n-1,dbmsname);

    nold = n;
  }
  else {
    /* Clear working memory: */
    num_bytes = n*sizeof(FLOAT);
    memset(c,       0,num_bytes);
    memset(d,       0,num_bytes);
    memset(x_old,   0,num_bytes);
    memset(fvec,    0,num_bytes);
    memset(fvec2,   0,num_bytes);
    memset(fvec_old,0,num_bytes);
    memset(g,       0,num_bytes);
    memset(sn,      0,num_bytes);
    memset(s,       0,num_bytes);
    memset(t,       0,num_bytes);
    memset(w,       0,num_bytes);
    num_bytes *= n;
    memset(qt,      0,num_bytes);
    memset(r,       0,num_bytes);
  }

  /*
   * Calculate fvec[].
   */
  (*vecfunc)(n,x,fvec);

  f = 0.;
  for (i = 0; i < n; i++) {
    f += fvec[i]*fvec[i];
  }
  f *= 0.5;

  /*
   * Test if initial guess is a root.
   * NOTE: We do not compare to the more stringent 0.01*tol_f used in 
   *       NR and DS96 at iteration zero, so that our solutions are 
   *       consistent.
   */
  test = 0.;
  for (i = 0; i < n; i++) {
    test = MAX(test,fabs(fvec[i]));
  }
  if (test < tol_f) {
    /* initial x[] is a root. */
    return status;
  }

  /*
   * Calculate max_step for globally convergent step.
   */
  sum = 0.;
  for (i = 0; i < n; i++) {
    sum += x[i]*x[i];
  }
  max_step = DS_MAX_STEP*MAX(sqrt(sum),(FLOAT)n);

  /*
   * Main iteration loop.
   */
  restart = TRUE;
  for (it = 0; it < max_it; it++) {
    if (restart == TRUE) {
      /*
       * Compute forward-difference approximation to Jacobian.
       */
      for (j = 0; j < n; j++) {
        tmp  = x[j];
        h    = DS_SQRT_EPS*fabs(tmp);
        if (h == 0.) {
          h = DS_SQRT_EPS;
        }
        x[j] = tmp+h;

        (*vecfunc)(n,x,fvec2);

        h    = x[j]-tmp;
        x[j] = tmp;
        for (i = 0; i < n; i++) {
          R(i,j) = (fvec2[i]-fvec[i])/h;
        }
      }

      /*
       * Calculate QR decomposition.
       */
      singular = qr_decompose(n,r,c,d);
      if (singular) {
        status = DS_SINGULAR_JACOBIAN;
        return status;
      }


      /* Compute transpose, QT. */
      memset(qt,0,n*n*sizeof(FLOAT));
      for (i = 0; i < n; i++) {
        QT(i,i) = 1.;
      }
      for (k = 0; k < n; k++) {
        if (c[k]) {
          for (j = 0; j < n; j++) {
            sum = 0.;
            for (i = k; i < n; i++) {
              sum += R(i,k)*QT(i,j);
            }
            sum /= c[k];
            for (i = k; i < n; i++) {
              QT(i,j) -= sum*R(i,k);
            }
          }
        }
      }
      /* Form R explicitly. */
      for (i = 0; i < n; i++) {
        R(i,i) = d[i];
        for (j = 0; j < i; j++) {
          R(i,j) = 0.;
        }
      }
    }
    else {
      for (i = 0; i < n; i++) {
        s[i] = x[i]-x_old[i];
      }
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = i; j < n; j++) {
          sum += R(i,j)*s[j];
        }
        t[i] = sum;
      }
      skip = TRUE;
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = 0; j < n; j++) {
          sum += QT(j,i)*t[j];
        }
        w[i] = fvec[i]-fvec_old[i]-sum;
        if (fabs(w[i]) >= DS_EPS*(fabs(fvec[i])+fabs(fvec_old[i]))) {
          skip = FALSE;
        }
        else {
          w[i] = 0.;
        }
      }
      if (skip == FALSE) {
        for (i = 0; i < n; i++) {
          sum = 0.;
          for (j = 0; j < n; j++) {
            sum += QT(i,j)*w[j];
          }
          t[i] = sum;
        }

        denom = 0.;
        for (i = 0; i < n; i++) {
          denom += s[i]*s[i];
        }
        /* Store s/(s.s) in s. */
        for (i = 0; i < n; i++) {
          s[i] /= denom;
        }

        /* 
         * Update r and qt. 
         */
        qr_update(n,r,qt,t,s);

        for (i = 0; i < n; i++) {
          if (R(i,i) == 0.) {
            fprintf(stderr,"**error:%s, R(%d,%d) singular\n",dbmsname,i,i);
            exit(1);
          }
          d[i] = R(i,i);
        }
      }
    }

    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i,j)*fvec[j];
      }
      g[i] = sum;
    }
    for (i = n-1; i >= 0; i--) {
      sum = 0.;
      for (j = 0; j <= i; j++) {
        sum += R(j,i)*g[j];
      }
      g[i] = sum;
    }

    /*
     * Store old x,fvec,f.
     */
    for (i = 0; i < n; i++) {
      /*
       * Screen for NaN: 
       */
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        fprintf(stderr,"**error:%s, x[%d]=%g fvec[%d]=%g\n",dbmsname,i,x[i],i,fvec[i]);
        exit(1);
      }
      x_old[   i] = x[   i];
      fvec_old[i] = fvec[i];
    }
    f_old = f;

    /*
     * Compute right-hand side of linear equations, sn[].
     */
    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i,j)*fvec[j];
      }
      sn[i] = -sum;
    }

    /*
     * Solve R.x = sn.
     * See Numerical Recipes in C, 2nd ed., p. 100.
     */
    sn[n-1] /= d[n-1];
    for (i = n-2; i >= 0; i--) {
      sum = 0.;
      for (j = i+1; j < n; j++) {
        sum += R(i,j)*sn[j];
      }
      sn[i] = (sn[i]-sum)/d[i];
    }

    /*
     * Calculate new x,f,fvec[].
     */
    old_max_taken = max_taken;
    max_taken = global_step(n,x_old,f_old,g,r,sn,max_step,&delta,
                            DS_DOGLEG_STEP,&status,x,&f,fvec,vecfunc);

    /*
     * Screen for NaN: 
     */
    for (i = 0; i < n; i++) {
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        fprintf(stderr,"**error:%s, after global_step(): x[%d]=%g fvec[%d]=%g\n",dbmsname,i,x[i],i,fvec[i]);
        exit(1);
      }
    }

    /*
     * Screen for repeated maximum steps.
     */
    if (max_taken == TRUE) {
      if (old_max_taken == TRUE) {
        count_max_taken++;
      }
      else {
        count_max_taken = 1;
      }
    }
    else {
      count_max_taken = 0;
    }

    /*
     * Test for convergence.
     */
    test = 0.;
    for (i = 0; i < n; i++) {
      test = MAX(test,fabs(fvec[i]));
    }
    if (test < tol_f) {
      status = DS_X_ACCEPTED;
      return status;
    }

    if (count_max_taken >= 5) {
      status = DS_MAX_TAKEN_5;
      return status;
    }

    if (status == DS_X_NO_PROGRESS) {
      if (restart == TRUE) {
        return status;
      }
      else {
        test  = 0.;
        denom = MAX(f,0.5*(FLOAT)n);
        for (i = 0; i < n; i++) {
          tmp  = fabs(g[i])*MAX(fabs(x[i]),1.)/denom;
          test = MAX(test,tmp);
        }
        if (test < DS_TOL_MIN) {
          return status;
        }
        else {
          /*
           * Try reinitializing the Jacobian.
           */
          restart = TRUE;
        }
      }
    }
    else {
      restart = FALSE;
      test    = 0.;
      for (i = 0; i < n; i++) {
        tmp  = (fabs(x[i]-x_old[i]))/MAX(fabs(x[i]),1.);
        test = MAX(test,tmp);
        if (test < DS_TOL_X) {
          /* 
           * Convergence. 
           */
          return status;
        }
      }
    }
  }

  status = DS_MAX_IT_EXCEEDED;
  return status;
}

/*======================= end of broyden_root() ==============================*/

/*======================= global_step() ======================================*/


/*
 * Take a globally convergent step towards a vector root.
 * Returns max_taken.
 *
 * The function line_search() backtracks along the quasi-Newton direction.
 *
 * The function dogleg_driver() uses the model trust-region approach, where
 * delta is the radius of the trust region. A step is taken 
 * in the steepest descent direction for small delta, in the quasi_Newton 
 * direction for large delta, and on the connecting line segment for  
 * intermediate delta.
 * See Dennis and Schnabel (1996, DS96).
 *
 * The value of step_type can be DS_LINE_STEP, DS_HOOK_STEP, or DS_DOGLEG_STEP.
 *
 * NOTE: Unlike in DS96, here R is R of QR, not the transpose of R. 
 */

int global_step(int    n,
                FLOAT *x_old,
                FLOAT  f_old,
                FLOAT *g,
                FLOAT *r,
                FLOAT *sn,
                FLOAT  max_step,
                FLOAT *delta,
                int    step_type,
                int   *status,
                FLOAT *x,
                FLOAT *f,
                FLOAT *fvec,
                void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    max_taken = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="global_step";

  if (step_type == DS_LINE_STEP) {
    max_taken = line_search(n,x_old,f_old,g,sn,max_step,
                            status,x,f,fvec,vecfunc);
  }
  else if (step_type == DS_HOOK_STEP) {
    fprintf(stderr,"**error:%s, DS_HOOK_STEP not yet implemented\n",dbmsname);
    exit(1);
  }
  else if (step_type == DS_DOGLEG_STEP) {
    max_taken = dogleg_driver(n,x_old,f_old,g,r,sn,max_step,delta,
                              status,x,f,fvec,vecfunc);
  }
  else {
    fprintf(stderr,"**error:%s, unrecognized step_type=%d\n",dbmsname,step_type);
    exit(1);
  }

  return max_taken;
}

/*======================= end of global_step() ===============================*/

/*======================= line_search() ======================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 385-386, which is
 * an implementation of Algorithm A6.3.1 of Dennis and Schnabel (1996).
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

int line_search(int    n,
                FLOAT *x_old,
                FLOAT  f_old,
                FLOAT *g,
                FLOAT *sn,
                FLOAT  max_step,
                int   *status,
                FLOAT *x,
                FLOAT *f,
                FLOAT *fvec,
                void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,
    max_taken = FALSE;
  FLOAT
    a,b,
    lambda,lambda_prev,lambda_min,
    disc,f_prev,
    rhs1,rhs2,
    initial_slope,
    newt_length,
    rel_step_length,
    tmp,
    tmp_lambda; 
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="line_search";

  *status = DS_X_ACCEPTED;

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i]*sn[i];
  }
  newt_length = sqrt(tmp);

  if (newt_length > max_step) {
    tmp = max_step/newt_length;
    for (i = 0; i < n; i++) {
      sn[i] *= tmp;
    }
  }

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i]*sn[i];
  }

  rel_step_length = 0.;
  for (i = 0; i < n; i++) {
    rel_step_length = MAX( rel_step_length,
                           fabs(sn[i])/MAX(fabs(x_old[i]),1.) );
  }

  lambda_min = DS_TOL_X/rel_step_length;
  lambda     = 1.;

  /*
   * Iteration loop.
   */
  while (TRUE) {
    for (i = 0; i < n; i++) {
      x[i] = x_old[i]+lambda*sn[i];
    }

    /*
     * Calculate fvec[].
     */
    (*vecfunc)(n,x,fvec);

    *f = 0.;
    for (i = 0; i < n; i++) {
      *f += fvec[i]*fvec[i];
    }
    *f *= .5;

    if (lambda < lambda_min) {
      /* 
       * Convergence on dx. 
       * For zero finding, calling program should verify the convergence.
       */
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
      *status = DS_X_NO_PROGRESS;
      return max_taken;
    }
    else if (*f <= f_old+DS_ALPHA*lambda*initial_slope) {
      /* Sufficient function decrease. */
      if (lambda == 1. && newt_length > 0.99*max_step) {
        max_taken = TRUE;
      }
      return max_taken;
    }
    else {
      if (lambda == 1.) {
        /* First backtrack uses a quadratic model. */
        tmp_lambda = -initial_slope/(2.*(*f-f_old-initial_slope));
      }
      else {
        /* Subsequent backtracks use a cubic model. */
        rhs1 =     *f-initial_slope*lambda     -f_old;
        rhs2 = f_prev-initial_slope*lambda_prev-f_old;
        a    = (rhs1/(lambda*lambda)-rhs2/(lambda_prev*lambda_prev))
               /(lambda-lambda_prev);
        b    = ( -lambda_prev*rhs1/(lambda*lambda)
                 +lambda*rhs2/(lambda_prev*lambda_prev) )/(lambda-lambda_prev);
        if (a == 0.) {
          tmp_lambda = -initial_slope/(2.*b);
        }
        else {
          disc = b*b-3.*a*initial_slope;
          if (disc < 0.) {
            fprintf(stderr,"**error:%s,roundoff problem\n",dbmsname);
            exit(1);
          }
          else {
            tmp_lambda = (-b+sqrt(disc))/(3.*a);
          }
        }
        if (tmp_lambda > .5*lambda) {
          tmp_lambda = .5*lambda;
        }
      }
    }
    lambda_prev = lambda;
    f_prev      = *f;
    lambda      = MAX(tmp_lambda,.1*lambda);
  }

  /* Never get here. */
}

/*======================= end of line_search() ===============================*/

/*======================= dogleg_driver() ====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.3.
 * Returns max_taken;
 * Assumes zero-based indexing.
 */
int dogleg_driver(int    n,
                  FLOAT *x_old,
                  FLOAT  f_old,
                  FLOAT *g,
                  FLOAT *r,
                  FLOAT *sn,
                  FLOAT  max_step,
                  FLOAT *delta,
                  int   *status,
                  FLOAT *x,
                  FLOAT *f,
                  FLOAT *fvec,
                  void   (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,
    max_taken,
    newt_taken,
    first_dog = TRUE;
  FLOAT
    newt_length,
    f_prev,
    tmp,
   *s,
   *s_hat,
   *nu_hat,
   *x_prev;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dogleg_driver";

  /*
   * Allocate memory.
   */
  s      = fvector(0,n-1,dbmsname);
  s_hat  = fvector(0,n-1,dbmsname);
  nu_hat = fvector(0,n-1,dbmsname);
  x_prev = fvector(0,n-1,dbmsname);

  *status = DS_INITIAL;

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i]*sn[i];
  }
  newt_length = sqrt(tmp);

  while (*status >= DS_REDUCE_DELTA) {
    /*
     * Find new step.
     */
    newt_taken = dogleg_step(n,g,r,sn,newt_length,max_step,
                             delta,&first_dog,s_hat,nu_hat,s);
    /*
     * Check new point and update trust region.
     */
    max_taken = trust_region(n,x_old,f_old,g,s,newt_taken,max_step,
                             DS_DOGLEG_STEP,r,delta,status,x_prev,&f_prev,
                             x,f,fvec,vecfunc);
  }

  /*
   * Free allocated memory.
   */
  free_fvector(x_prev,0,n-1,dbmsname);
  free_fvector(nu_hat,0,n-1,dbmsname);
  free_fvector(s_hat, 0,n-1,dbmsname);
  free_fvector(s,     0,n-1,dbmsname);

  return max_taken;
}

/*======================= end of dogleg_driver() =============================*/

/*======================= dogleg_step() ======================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.4.
 * Returns newt_taken.
 * Assumes zero-based indexing.
 */
int dogleg_step(int    n,
                FLOAT *g,
                FLOAT *r,
                FLOAT *sn,
                FLOAT  newt_length,
                FLOAT  max_step,
                FLOAT *delta,
                int   *first_dog,
                FLOAT *s_hat,
                FLOAT *nu_hat,
                FLOAT *s)
{
  int
    i,j,
    newt_taken;
  static int
    eta_warned = FALSE;
  FLOAT
    alpha,beta,al_be,
    eta,
    lambda,
    tmp,tmp_nu,tmp_cauchy;
  static FLOAT
    cauchy_length;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dogleg_step";

  if (newt_length <= *delta) {
    /*
     * s is Newton step.
     */
    newt_taken = TRUE;
    for (i = 0; i < n; i++) {
      s[i] = sn[i];
    }
    *delta = newt_length;
  }
  else {
    /*
     * Newton step is too long, find s on double-dogleg curve.
     */
    newt_taken = FALSE;
    if (*first_dog == TRUE) {
      /*
       * Calculate double-dogleg curve.
       */
      *first_dog = FALSE;

      alpha      = 0.;
      for (i = 0; i < n; i++) {
        alpha += g[i]*g[i];
      }

      beta = 0.;
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i,j)*g[j];
        }
        beta += tmp*tmp;
      }
      al_be         = alpha/beta;
      cauchy_length = al_be*sqrt(alpha);

      tmp = 0.;
      for (i = 0; i < n; i++) {
        s_hat[i] = -g[i]*al_be;
        tmp     +=  g[i]*sn[i];
      }
      eta = 0.2+0.8*alpha*al_be/fabs(tmp);

      /*
       * Check range of eta, which should be [0,1]:
       */
      if (eta > 1.) {
        /* 
         * If eta > 1.0 print a one-time warning. 
         * This behavior is known to be caused by optimization (-O) 
         * for gcc (at least versions 2.7.2.3.f.1 and 2.96).
         */
        if (!eta_warned) {
          fprintf(stderr,"eta=%g > 1, setting eta=1. This warning will not be repeated.\n"
                          "Check for optimization error (-O vs -g).",eta);
          eta_warned = TRUE;
        }
        eta = 1.;
      }

      for (i = 0; i < n; i++) {
        nu_hat[i] = eta*sn[i]-s_hat[i];
      }

      if (*delta == -1) {
        /*
         * First iteration, and no initial trust region was
         * provided by the user.
         */
        *delta = MIN(cauchy_length,max_step);
      }
    }

    if (eta*newt_length <= *delta) {
      /*
       * Take partial step in Newton direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta)/newt_length)*sn[i];
      }
    }
    else if ((cauchy_length) >= (*delta)) {
      /*
       * Take step in steepest descent direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta)/(cauchy_length))*s_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          fprintf(stderr,"**error:%s, s[%d]=%g,*delta=%g,cauchy_length=%g,s_hat=%g",
                           dbmsname,i,s[i],*delta,cauchy_length,s_hat[i]);
          exit(1);
        }
      }
    }
    else {
      /*
       * Take convex-combination step.
       */
      tmp    = 0.;
      tmp_nu = 0.;
      for (i = 0; i < n; i++) {
        tmp    += nu_hat[i]*s_hat[i];
        tmp_nu += nu_hat[i]*nu_hat[i];
      }
      tmp_cauchy = cauchy_length*cauchy_length-(*delta)*(*delta);
      lambda     = (-tmp+sqrt(tmp*tmp-tmp_nu*tmp_cauchy))/tmp_nu;
      for (i = 0; i < n; i++) {
        s[i] = s_hat[i]+lambda*nu_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          fprintf(stderr,"**error:%s, s[%d]=%g,lambda=%g,nu_hat=%g,s_hat=%g",
                           dbmsname,i,s[i],lambda,nu_hat[i],s_hat[i]);
          exit(1);
        }
      }
    }
  }


  return newt_taken;
}


/*======================= end of dogleg_step() ===============================*/

/*======================= trust_region() =====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.5.
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

int trust_region(int    n,
                 FLOAT *x_old,
                 FLOAT  f_old,
                 FLOAT *g,
                 FLOAT *s,
                 int    newt_taken,
                 FLOAT  max_step,
                 int    step_type,
                 FLOAT *r,
                 FLOAT *delta,
                 int   *status,
                 FLOAT *x_prev,
                 FLOAT *f_prev,
                 FLOAT *x,
                 FLOAT *f,
                 FLOAT *fvec,
                 void  (*vecfunc)(int,FLOAT *,FLOAT *))
{
  int
    i,j,
    max_taken = FALSE;
  FLOAT
    initial_slope,
    step_length,
    rel_step_length,
    delta_tmp,
    tmp;
  FLOAT
     df,
     df_tol,
     df_pred;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="trust_region";

  /*
   * Screen for NaN.
   */
  for (i = 0; i < n; i++) {
    if (!isfinite(s[i])) {
      fprintf(stderr,"**error:%s,s[%d]=%g\n",dbmsname,i,s[i]);
      exit(1);
    }
  }

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += s[i]*s[i];
  }
  step_length = sqrt(tmp);

  /* Take step. */
  for (i = 0; i < n; i++) {
    x[i] = x_old[i]+s[i];
  }

  /* Calculate fvec[]. */
  (*vecfunc)(n,x,fvec);

  /* Compute f. */
  *f = 0.;
  for (i = 0; i < n; i++) {
    *f += fvec[i]*fvec[i];
  }
  *f *= .5;

  df = *f-f_old;

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i]*s[i];
  }

  if (*status != DS_INCREASE_DELTA) {
    *f_prev = 0.;
  }

  df_tol = DS_ALPHA*initial_slope;

  if (*status == DS_INCREASE_DELTA && (*f >= *f_prev || df > df_tol)) {
    /*
     * Retreat.
     */
    *status = DS_X_ACCEPTED;
    for (i = 0; i < n; i++) {
      x[i] = x_prev[i];
    }
    *f      = *f_prev;
    *delta *= .5;
  }
  else if (df >= df_tol) {
    /*
     * The value of f is too large.
     */
    rel_step_length = 0.;
    for (i = 0; i < n; i++) {
      rel_step_length = MAX( rel_step_length,
                             fabs(s[i])/MAX(fabs(x[i]),1.) );
    }

    if (rel_step_length < DS_TOL_X) {
      /*
       * The step is too small.
       */
      *status = DS_X_NO_PROGRESS;
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
    }
    else {
      /*
       * Reduce delta.
       */
      *status   = DS_REDUCE_DELTA;
      delta_tmp = (-initial_slope*step_length)/(2.*(df-initial_slope));
      if (delta_tmp < (*delta)*0.1) {
        *delta *= 0.1;
      }
      else if (delta_tmp > (*delta)*0.5) {
        *delta *= 0.5;
      }
      else {
        *delta = delta_tmp;
      }
    }
  }
  else {
    /*
     * The value of f is sufficiently small.
     */
    df_pred = initial_slope;

    if (step_type == DS_HOOK_STEP) {
      fprintf(stderr,"**error:%s, DS_HOOK_STEP not yet implemented\n",dbmsname);
      exit(1);
    }
    else if (step_type == DS_DOGLEG_STEP) {
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i,j)*s[j];
        }
        df_pred += tmp*tmp*.5;
      }
    }
    else {
      fprintf(stderr,"**error:%s, unrecognized step_type=%d\n",dbmsname,step_type);
      exit(1);
    }

    if ( ((*status) != DS_REDUCE_DELTA && fabs(df_pred-df) <= 0.1*fabs(df)) ||
         (df <= initial_slope && newt_taken == FALSE && (*delta) <= 0.99*max_step) ) {
      /*
       * Double delta.
       */
      *status = DS_INCREASE_DELTA;
      for (i = 0; i < n; i++) {
        x_prev[i] = x[i];
      }
      *f_prev = *f;
      *delta  = MIN((*delta)*2.,max_step);
    }
    else {
      /*
       * Accept x, choose delta for next iteration.
       */
      *status = DS_X_ACCEPTED;
      if (step_length > 0.99*max_step) {
        max_taken = TRUE;
      }
      if (df >= 0.1*df_pred) {
        /*
         * Decrease delta.
         */
        *delta *= 0.5;
      }
      else if (df <= 0.75*df_pred) {
        /*
         * Increase delta.
         */
        *delta = MIN((*delta)*2.,max_step);
      }
      else {
        /*
         * Leave delta unchanged.
         */
        ;
      }
    }
  }

  return max_taken;
}

/*======================= end of trust_region() ==============================*/

/*======================= qr_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p.99.
 * Returns FALSE if normal and TRUE if singular.
 * Assumes zero-based indexing.
 */

int qr_decompose(int    n,
                 FLOAT *r,
                 FLOAT *c,
                 FLOAT *d)
{
  int
    i,j,k,
    singular;
  FLOAT
    sigma,sum,tau,scale;

  singular = FALSE;

  for (k = 0; k < n-1; k++) {
    /* 
     * Put scale=0. inside k loop as in Dennis & Schnabel (1996), 
     * Algorithm A3.2.1, p.305, rather than outside the loop as in 
     * Numerical Recipes' qrdcmp(), p. 99.
     */
    scale = 0.;
    for (i = k; i < n; i++) {
      scale = MAX(scale,fabs(R(i,k)));
    }
    if (scale == 0.) {
      /*
       * Singular case.
       */
      c[k]     = 0.;
      d[k]     = 0.;
      singular = TRUE;
    }
    else {
      for (i = k; i < n; i++) {
        R(i,k) /= scale;
      }
      sum = 0.;
      for (i = k; i < n; i++) {
        sum += R(i,k)*R(i,k);
      }
      sigma = NR_SIGN(sqrt(sum),R(k,k)); 
      R(k,k) += sigma;
      c[k]    = sigma*R(k,k);
      d[k]    = -scale*sigma;
      for (j = k+1; j < n; j++) {
        sum = 0.;
        for (i = k; i < n; i++) {
          sum += R(i,k)*R(i,j);
        }
        tau = sum/c[k];
        for (i = k; i < n; i++) {
          R(i,j) -= tau*R(i,k);
        }
      }
    }
  }
  d[n-1] = R(n-1,n-1);

  if (d[n-1] == 0.) {
    singular = TRUE;
  }

  return singular;
}

/*======================= end of qr_decompose() ==============================*/

/*======================= qr_update() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_update(int    n,
               FLOAT *r,
               FLOAT *qt,
               FLOAT *u,
               FLOAT *v)
{
  int
    i,j,k;
  FLOAT
    tmp;

  /* Find largest k such that u[k] != 0. */
  for (k = n-1; k > 0; k--) {
    if (u[k]) {
      break;
    }
  }

  for (i = k-1; i >= 0; i--) {
    qr_rotate(n,r,qt,i,u[i],-u[i+1]);
    if (u[i] == 0.) {
      u[i] = fabs(u[i+1]);
    }
    else if (fabs(u[i]) > fabs(u[i+1])) {
      tmp  = u[i+1]/u[i];
      u[i] = fabs(u[i])*sqrt(1.+tmp*tmp);
    }
    else {
      tmp  = u[i]/u[i+1];
      u[i] = fabs(u[i+1])*sqrt(1.+tmp*tmp);
    }
  }

  for (j = 0; j < n; j++) {
    R(0,j) += u[0]*v[j];
  }
  for (i = 0; i < k; i++) {
    qr_rotate(n,r,qt,i,R(i,i),-R(i+1,i));
  }

  return;
}

/*======================= end of qr_update() =================================*/

/*======================= qr_rotate() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_rotate(int    n,
               FLOAT *r,
               FLOAT *qt,
               int    i,
               FLOAT  a,
               FLOAT  b)
{
  int
    j;
  FLOAT
    c,factor,
    s,w,y;

  if (a == 0.) {
    c = 0.;
    s = b > 0. ? 1. : -1.;
  }
  else if (fabs(a) > fabs(b)) {
    factor = b/a;
    c      = NR_SIGN(1./sqrt(1.+factor*factor),a);
    s      = factor*c;
  }
  else {
    factor = a/b;
    s      = NR_SIGN(1./sqrt(1.+factor*factor),b);
    c      = factor*s;
  }

  for (j = 0; j < n; j++) {
    y        = R(i,  j);
    w        = R(i+1,j);
    R(i,  j) = c*y-s*w;
    R(i+1,j) = s*y+c*w;
  }

  for (j = 0; j < n; j++) {
    y         = QT(i,  j);
    w         = QT(i+1,j);
    QT(  i,j) = c*y-s*w;
    QT(i+1,j) = s*y+c*w;
  }

  return;
}

/*======================= end of qr_rotate() =================================*/

/*======================= lu_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, pp. 46-47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_decompose(int   n,
                 FLOAT *a,
                 int   *index,
                 FLOAT *d)
{
  int
    i,imax,j,k;
  static int
    n_max = 0;
  FLOAT
    tiny = 1.e-20,
    big,dum,sum,temp;
  static FLOAT
   *vv;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_decompose";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    vv    = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(vv,0,n_max-1,dbmsname);
    n_max = n;
    vv    = fvector(0,n_max-1,dbmsname);
  }

  *d = 1.;
  for (i = 0; i < n; i++) {
    big = 0.;
    for (j = 0; j < n; j++) {
      temp = fabs(A(i,j));
      if (temp > big) {
        big = temp;
      }
    }
    if (big == 0.) {
      fprintf(stderr,"**error:%s, matrix singular\n",dbmsname);
      exit(1);
    }
    vv[i] = 1./big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = A(i,j);
      for (k = 0; k < i; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
    }
    big = 0.;
    for (i = j; i < n; i++) {
      sum = A(i,j);
      for (k = 0; k < j; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
 
      dum = vv[i]*fabs(sum);
      if (dum >= big) {
        big  = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum       = A(imax,k);
        A(imax,k) = A(j,   k);
        A(j,   k) = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    index[j] = imax;

    if (A(j,j) == 0.) {
      A(j,j) = tiny;
      fprintf(stderr,"Warning, %s(), matrix is singular\n",dbmsname);
    }

    if (j != n-1) {
      dum = 1./A(j,j);
      for (i = j+1; i < n; i++) {
        A(i,j) *= dum;
      }
    }
  }
 
  return;
}

/*======================= end of lu_decompose() ==============================*/

/*======================= lu_backsub() =======================================*/

/*
 * Adapted from Numerical Recipes in C, p. 47.
 * Assumes zero-based indexing, with (i,j) referring to (row,column).
 */

#undef  A
#define A(i,j) a[(j)+n*(i)]

void lu_backsub(int    n,
                FLOAT *a,
                int   *index,
                FLOAT *b)
{
  int
    ii = -1,
    i,ip,j; 
  FLOAT
    sum;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_backsub";

  for (i = 0; i < n; i++) {
    ip    = index[i];
    sum   = b[ip];
    b[ip] = b[i];
    if (ii != -1) {
      for (j = ii; j <= i-1; j++) {
        sum -= A(i,j)*b[j];
      }
    }
    else if (sum != 0.) {
      ii = i;
    }
    b[i] = sum;
  } 

  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i+1; j < n; j++) {
      sum -= A(i,j)*b[j];
    }
    b[i] = sum/A(i,i);
  }

  return;
}

/*======================= end of lu_backsub() ================================*/

/*======================= lu_improve() =======================================*/

/*
 * Adapted from Numerical Recipes in C, p. 56.
 * Assumes zero-based indexing.
 *
 * NOTE: May not have much effect if lu_decompose(), lu_backsub()
 *       are in double precision.
 */

void lu_improve(int    n,
                FLOAT *a,
                FLOAT *alu,
                int   *index,
                FLOAT *b,
                FLOAT *x)
{
  int
    j,i;
  static int
    n_max = 0;
  static FLOAT
   *r;
  double
    sdp;  /* NOTE: sdp must be double precision. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lu_improve";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(r,0,n_max-1,dbmsname);
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }

  for (i = 0 ; i < n; i++) {
    sdp = -b[i];
    for (j = 0; j < n; j++) {
      sdp += A(i,j)*x[j];
    }
    r[i] = sdp;
  }
  lu_backsub(n,alu,index,r);
  for (i = 0; i < n; i++) {
    x[i] -= r[i];
  }

  return;
}

/*======================= end of lu_improve() ================================*/

/*======================= find_place_in_table() ==============================*/

/*
 * Find place in table using bisection.
 * Adapted from Numerical Recipes in C, p. 117.
 *
 * NOTE: Finds place relative to the x component of the float_triplet table.
 */
int find_place_in_table(int            n,
                        float_triplet *table,
                        FLOAT          x,
                        FLOAT         *dx)
{
  register int
    il,im,iu,
    ascend,
    cmplo,cmphi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="find_place_in_table";

  il     = -1;
  iu     =  n;
  ascend = (table[n-1].x >= table[0].x);
  cmplo  = fcmp(x,table[0  ].x);
  cmphi  = fcmp(x,table[n-1].x);

  /*
   * Check x against table endpoints.
   */
  if (ascend) {
    if (cmplo <= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi >= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }
  else {
    if (cmplo >= 0) {
      il  = 0;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
    if (cmphi <= 0) {
      il  = n-2;
      *dx = table[il+1].x-table[il].x;
      return il;
    }
  }

  /*
   * Use bisection to search table.
   */
  while (iu-il > 1) {
    im = (iu+il) >> 1;
    if (fcmp(x,table[im].x) >= 0 == ascend) {
      il = im;
    }
    else {
      iu = im;
    }
  }

  *dx = table[il+1].x-table[il].x;
  return il;
}

/*======================= end of find_place_in_table() =======================*/

/*======================= tridiag() =========================================*/

/*
 * Solve a tridiagnonal matrix system.
 * Adapted from Numerical Recipes in C, 2nd ed., pp. 51-54.
 * If pivot_type = WITH_PIVOTING, use band_decomp() and band_back_sub().
 * Assumes zero-based indexing.
 */

#undef  AA
#define AA(i,j) aa[(m1+m2+1)*(i)+(j)]
#undef  AAORIG
#define AAORIG(i,j) aaorig[(m1+m2+1)*(i)+(j)]
#undef  AAL
#define AAL(i,j) aal[m1*(i)+(j)]

void tridiag(int    n,
             FLOAT *a,
             FLOAT *b,
             FLOAT *c,
             FLOAT *r,
             FLOAT *u,
             int    pivot_type)
{
  int
    j,
    m1,m2,mm,
    *index;
  FLOAT
    bet,
    *gam,
    *aa,
    *aaorig,
    *aal,
     d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="tridiag";

  /*
   * Check validity of n.
   */
  if (n <= 0) {
    fprintf(stderr,"**error:%s, n = %d\n",dbmsname,n);
    exit(1);
  }

  if (pivot_type == WITHOUT_PIVOTING) {
    /* Allocate memory. */
    gam = fvector(0,n-1,dbmsname);

    if (b[0] == 0.0) {
      fprintf(stderr,"**error:%s,b[0] = 0\n"
                     "Rewrite equations as a set of order n-1, with u[1]\n"
                     "trivially eliminated\n",dbmsname);
      exit(1);
    }
    bet  = b[0];
    u[0] = r[0]/bet;
    for (j = 1; j < n; j++) {
      gam[j] = c[j-1]/bet;
      bet    = b[j]-a[j]*gam[j];
      if (bet == 0.) {
        /* 
         * Encountered a zero pivot. 
         * Try again using pivot_type = WITH_PIVOTING.
         */
        fprintf(stderr,"Warning: tridiag(): retrying with pivoting.\n");
        /* Free allocated memory. */
        free_fvector(gam,0,n-1,dbmsname);
        tridiag(n,a,b,c,r,u,WITH_PIVOTING);
        return;
      }
      u[j] =(r[j]-a[j]*u[j-1])/bet;
    }

    /* Backsubstitution: */
    for (j = n-2; j >= 0; j--) {
      u[j] -= gam[j+1]*u[j+1];
    }
    /* Free allocated memory. */
    free_fvector(gam,0,n-1,dbmsname);
    return;
  }
  else if (pivot_type == WITH_PIVOTING) {
    /*
     * Use band_decomp() and band_back_sub().
     */
    m1 = 1;
    m2 = 1;
    mm = m1+m2+1;
    /*
     * Allocate memory.
     */
    aa     = fvector(0,n*(m1+m2+1)-1,dbmsname);
    aaorig = fvector(0,n*(m1+m2+1)-1,dbmsname);
    aal    = fvector(0,n*m1-1,dbmsname);
    index  = ivector(0,n-1,dbmsname);

    /*
     * Load matrix AA and keep copy AAORIG.
     */
    for (j = 0; j < n; j++) {
      AA(j,m1+1) = AAORIG(j,m1+1) = c[j];
      AA(j,m1  ) = AAORIG(j,m1  ) = b[j];
      AA(j,m1-1) = AAORIG(j,m1-1) = a[j];
    }
    
    band_decomp(n,m1,m2,aa,aal,index,&d);

    /* 
     * Since tridiag() does not overwrite the input rhs vector, r,
     * with the answer, u, but band_back_sub() does, copy r into u
     * before calling band_back_sub().
     */
    for (j = 0; j < n; j++) {
      u[j] = r[j];
    }

    band_back_sub(n,m1,m2,aa,aal,index,u);

    /*
     *  Reduce roundoff errors with call to band_improve().
     */
    band_improve(n,m1,m2,aaorig,aa,aal,index,r,u);

    /*
     * Free allocated memory.
     */
    free_fvector(aa,    0,n*(m1+m2+1)-1,dbmsname);
    free_fvector(aaorig,0,n*(m1+m2+1)-1,dbmsname);
    free_fvector(aal,   0,n*m1-1,       dbmsname);
    free_ivector(index, 0,n-1,          dbmsname);

    return;
  }
  else {
    fprintf(stderr,"**error:%s, unrecognized pivot_type=%d\n",dbmsname,pivot_type);
    exit(1);
  }
}

/*======================= end of tridiag() ===================================*/

/*======================= band_decomp() ======================================*/

/*
 * Decompose a banded matrix.
 * Adapted from Numerical Recipes in C, 2nd ed., p. 53.
 * The input matrix must be stored in compact form, as described on p. 52.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}
#undef  TINY
#define TINY 1.e-20

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_decomp(int     n,
                 int     m1,
                 int     m2,
                 FLOAT  *a,
                 FLOAT  *al,
                 int    *index,
                 FLOAT  *d)
{
  int
    i,j,k,l,
    mm;
  FLOAT
    tmp;
  static int
    warned = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_decomp";

  mm = m1+m2+1;
  l  = m1;
  for (i = 0; i < m1; i++) {
    for (j = m1-i; j < mm; j++) {
      A(i,j-l) = A(i,j);
    }
    l--;
    for (j = mm-l-1; j < mm; j++) {
      A(i,j) = 0.;
    }
  }
  *d = 1.;
  l  = m1;
  for (k = 0; k < n; k++) {
    tmp = A(k,0);
    i   = k;
    if (l < n) {
      l++;
    }
    for (j = k+1; j < l; j++) {
      if (fabs(A(j,0)) > fabs(tmp)) {
        tmp = A(j,0);
        i   = j;
      }
    }
    index[k] = i;
    if (tmp == 0.) {
      /*
       * Matrix is algorithmically singular.
       */
      A(k,0) = TINY;
      if (!warned) {
        fprintf(stderr,"**warning: %s, matrix is algorithmically singular (future warnings suppressed)\n",dbmsname);
        warned = TRUE;
      }
    }
    if (i != k) {
      *d = -(*d);
      for (j = 0; j < mm; j++) {
        SWAP(A(k,j),A(i,j))
      }
    }
    for (i = k+1; i < l; i++) {
      tmp          = A(i,0)/A(k,0);
      AL(k,i-k-1) = tmp;
      for (j = 1; j < mm; j++) {
        A(i,j-1) = A(i,j)-tmp*A(k,j);
      }
      A(i,mm-1) = 0.;
    }
  }

  return;
}

/*======================= end of band_decomp() ===============================*/

/*======================= band_back_sub() ====================================*/

/*
 * Back substitute for a banded matrix using the output from band_decomp().
 * Adapted from Numerical Recipes in C, 2nd ed., p. 54.
 * Assumes zero-based indexing.
 */

#undef  SWAP
#define SWAP(a,b) {tmp=(a);(a)=(b);(b)=tmp;}

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void band_back_sub(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *al,
                   int    *index,
                   FLOAT  *b)
{
  int
    i,k,l,
    mm;
  FLOAT
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_back_sub";

  mm = m1+m2+1;
  l  = m1;
  for (k = 0; k < n; k++) {
    i = index[k];
    if (i != k) {
      SWAP(b[k],b[i]);
    }
    if (l < n) {
      l++;
    }
    for (i = k+1; i < l; i++) {
      b[i] -= AL(k,i-k-1)*b[k];
    }
  }
  l = 1;
  for (i = n-1; i >= 0; i--) {
    tmp = b[i];
    for (k = 1; k < l; k++) {
      tmp -= A(i,k)*b[k+i];
    }
    b[i] = tmp/A(i,0);
    if (l < mm) {
      l++;
    }
  }

  return;
}

/*======================= end of band_back_sub() =============================*/

/*======================= band_multiply() ====================================*/
/*
 *  Compute matrix muliplication b = A.x, where A is in the compact-storage
 *  form of a band-diagonal matrix. 
 *  Based on Numerical Recipes in C, banmul(), p. 52.
 *  Assumes zero-based indexing.
 */

#undef  A
#define A(i,j) a[(m1+m2+1)*(i)+(j)]

void band_multiply(int     n,
                   int     m1,
                   int     m2,
                   FLOAT  *a,
                   FLOAT  *x,
                   FLOAT  *b)
{
  int
    i,j,k,tmploop;

  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    b[i]    = 0.;
    for (j = IMAX(0,-k); j < tmploop; j++) {
      b[i] += A(i,j)*x[j+k];
    }
  }

  return;
}
                   
/*======================= end of band_multiply() =============================*/

/*======================= band_improve() =====================================*/
/*
 * Based on Numerical Recipes in C, Secion 2.5, Iterative Improvement of
 * a Solution to Linear Equations.
 * This is for the band-diagnonal matrix case, and is analogous to lu_improve().
 * AORIG is the original matrix in compact form, whereas A and AL are the
 * matrices returned from band_decomp().
 * NOTE: The functionality of band_multiply() is echoed here because of the
 *       requirement of double precision.
 * Assumes zero-based indexing.
 */

#undef  AORIG
#define AORIG(i,j) aorig[(m1+m2+1)*(i)+(j)]

void band_improve(int     n,
                  int     m1,
                  int     m2,
                  FLOAT  *aorig,
                  FLOAT  *a,
                  FLOAT  *al,
                  int    *index,
                  FLOAT  *b,
                  FLOAT  *x)
{
  int
    k,j,i,tmploop;
  static int
    n_max = 0;
  static FLOAT
   *r;
  double
    sdp;  /* NOTE: sdp must be double precision. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="band_improve";

  /*
   * Allocate memory.
   */
  if (n_max == 0) {
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }
  else if (n > n_max) {
    free_fvector(r,0,n_max-1,dbmsname);
    n_max = n;
    r     = fvector(0,n_max-1,dbmsname);
  }

  /*
   * The band-diagonal indexing is as in band_multiply().
   */
  for (i = 0; i < n; i++) {
    k       = i-m1;
    tmploop = IMIN(m1+m2+1,n-k);
    sdp     = -b[i];
    for (j = IMAX(0,-k); j < tmploop; j++) {
      sdp += AORIG(i,j)*x[j+k];
    }
    r[i] = sdp;
  }

  band_back_sub(n,m1,m2,a,al,index,r);

  for (i = 0; i < n; i++) {
    x[i] -= r[i];
  }

  return;
}

/*======================= end of band_improve() ==============================*/

/*======================= poly_interp() ======================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., pp. 109-110.
 * Assumes zero-offset arrays.
 */

#undef  N_MAX
#define N_MAX 10

FLOAT poly_interp(int     n,
                  FLOAT  *xa,
                  FLOAT  *ya,
                  FLOAT   x,
                  FLOAT  *dy)
{
  int
    i,m,ns;
  static int
    initialized=0;
  FLOAT
    den,dif,dift,
    ho,hp,w,
    y;
  static FLOAT
    *c,
    *d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="poly_interp";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    c = fvector(0,N_MAX-1,dbmsname);
    d = fvector(0,N_MAX-1,dbmsname);
    initialized = 1;
  }

  dif = FLOAT_MAX;
  for (i = 0; i < n; i++) {
    dift = fabs(x-xa[i]);
    if (dift < dif) {
      ns  = i;
      dif = dift;
    }
  }
  y = ya[ns--];

  if (dif == 0.) {
    *dy = 0.;
  }
  else {

    if (n > N_MAX) {
      /*
       * Some systems have trouble with too many calls to calloc(),
       * so we allocate memory once for this function, assuming
       * n <= N_MAX.  Screen for this limit.
       */
      fprintf(stderr,"**error:%s, n=%d > N_MAX=%d, increase N_MAX.\n",dbmsname,n,N_MAX);
      exit(1);
    }

    for (i = 0; i < n; i++) {
      c[i] = ya[i];
      d[i] = ya[i];
    }

    for (m = 0; m <= n-2; m++) {
      for (i = 0; i <= n-2-m; i++) {
        ho  = xa[i    ]-x;
        hp  = xa[i+m+1]-x;
        w   = c[i+1]-d[i];
        den = ho-hp;
        if (den == 0.) {
          fprintf(stderr,"**error:%s, xa[%d]=%g, xa[%d]=%g\n",dbmsname,i,xa[i],i+m+1,xa[i+m+1]);
          exit(1);
        }
        den  = w/den;
        d[i] = hp*den;
        c[i] = ho*den;
      }
      *dy  = (2*ns < n-m-3) ? c[ns+1] : d[ns--];
      y   += *dy;
    }
  }

  return y;
}

/*======================= end of poly_interp() ===============================*/

/*======================= nth_trapezoidal() ==================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., p. 137.
 * Usage note: Call with successively increasing n, starting at n = 1.
 */

FLOAT nth_trapezoidal(int   n,
                      FLOAT (*func)(FLOAT),
                      FLOAT a,
                      FLOAT b)
{
  register int
    it,j;
  register FLOAT
    x,tnm,sum,dx;
  static FLOAT
    s;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="nth_trapezoidal";

  if (n < 1) {
    fprintf(stderr,"**error:%s, n=%d\n",dbmsname,n);
    exit(1);
  }
  else if (n == 1) {
    s = 0.5*(b-a)*(func(a)+func(b));
    return s;
  }
  else {
    it  = (1 << (n-2));
    tnm = (FLOAT)it;
    dx  = (b-a)/tnm;
    x   = a+0.5*dx;
    sum = 0.;
    for (j = 0; j < it; j++) {
      sum += func(x);
      x   += dx;
    }
    s = 0.5*(s+(b-a)*sum/tnm);
    return s;
  }

  /* Never get here. */
}

/*======================= end of nth_trapezoidal() ===========================*/

/*======================= romberg_integral() =================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd Ed., p. 140.
 * The argument tol is the fractional accuracy (tolerance) desired.
 */
#undef  MAX_IT
#define MAX_IT  30
#undef  NUM_PTS
#define NUM_PTS 5

FLOAT romberg_integral(FLOAT (*func)(FLOAT),
                       FLOAT a,
                       FLOAT b,
                       FLOAT tol)
{
  int
    it,
    n;
  FLOAT
    ss,
    dss,
    s[MAX_IT+1],
    h[MAX_IT+1];

  if (a == b) {
    return 0.;
  }

  h[0]   = 1.;
  for (it = 0; it < MAX_IT; it++) {
    n = it+1;
    s[it] = nth_trapezoidal(n,func,a,b);
    if (n >= NUM_PTS) {
      ss = poly_interp(NUM_PTS,h+(n-NUM_PTS),s+(n-NUM_PTS),0.,&dss);
      if (fabs(ss) < tol && fabs(dss) < tol) {
        return ss;
      }
      else if (fabs(dss) <= tol*fabs(ss)) {
        return ss;
      }
    }
    s[n] = s[it];
    h[n] = 0.25*h[it];
  }

  fprintf(stderr,"Warning: romberg_integral(): reached it=%d, a,b,tol=%g %g %g ans=%g\n",
                  MAX_IT,a,b,tol,ss);
  return ss;
}

/*======================= end of romberg_integral() ==========================*/

/*======================= crank_nicolson() ===================================*/

/*
 * Apply diffusion to the 1D input vector, A[j], for one timestep using the Crank-Nicolson scheme.
 *
 * The equation integrated in time is:
 *
 *      dA     1   d   /     dA  \
 *      --- = --- --- (  mu ----  )
 *      dt    rho  dz  \     dz  /
 *
 * The algorithm does not assume a regularly spaced grid, but
 * it does assume that the grid position, z[j], does not change during the timestep.
 * See Blottner (1980, Computer and Fluids 8, 421-434) for a description of how the variable grid
 * is handled.
 *
 * For input A[j], j = 0 and j = n+1 are boundary points.  Set to 1.e+20
 * to signal a no-flux boundary, otherwise set to the boundary value.
 *
 * One may specify the case rho[j] = 1.0 for all j by inputting NULL for rho.
 *
 * The input vector mu[j] may vary spatially, but is assumed to not change during the timestep.
 * The position of mu[j] is centered between A[j] and A[j+1].
 *
 * Input vectors z, A and D are assumed to be continguous 1D vectors with memory
 * allocations in the range [0,n+1]. The grid positions z are expected to be monotonically
 * increasing, and an error is generated if this is not the case.
 * The result is returned in ANS[j] for j = 1 to n, such that 
 * ANS[0] and ANS[n+1] are not altered. Memory for ANS should not overlap A.
 *
 * NOTE: Assumes no domain decomposition (ie. not MPI ready).
 */

void crank_nicolson(int    n,
                    FLOAT  dt,
                    FLOAT *z,
                    FLOAT *A,
                    FLOAT *mu,
                    FLOAT *rho,
                    FLOAT *ANS)
{
  int
    j,jstart,jend;
  FLOAT
    coeff1,coeff3,factor;
  static int
    nold = -1;
  static FLOAT
   *a,
   *b,
   *c,
   *r;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="crank_nicolson";

  /*
   * Screen for invalid n.
   */
  if (n < 1) {
    fprintf(stderr,"**error: %s: called with n=%d < 1\n",dbmsname,n);
    exit(1);
  }

  /* 
   * Screen for non-monotonically increasing z.
   */
  if (A[0] > .99e+20) {
    /* z[1]-z[0] is not used for no-flux b.c. */
    jstart = 2;
  }
  else {
    jstart = 1;
  }
  if (A[n+1] > .99e+20) {
    /* z[n+1]-z[n] is not used for no-flux b.c. */
    jend = n;
  }
  else {
    jend = n+1;
  }
  for (j = jstart; j <= jend; j++) {
    if (z[j] <= z[j-1]) {
      fprintf(stderr,"**error: %s: z[%d]=%g <= z[%d]=%g\n",dbmsname,j,z[j],j-1,z[j-1]);
      exit(1);
    }
  }

  /* Allocate memory for local vectors. */
  if (nold == -1) {
    a    = fvector(0,n-1,dbmsname);
    b    = fvector(0,n-1,dbmsname);
    c    = fvector(0,n-1,dbmsname);
    r    = fvector(0,n-1,dbmsname);
    nold = n;
  }
  else if (nold < n) {
    free_fvector(a,0,nold-1,dbmsname);
    free_fvector(b,0,nold-1,dbmsname);
    free_fvector(c,0,nold-1,dbmsname);
    free_fvector(r,0,nold-1,dbmsname);
    a    = fvector(0,n-1,dbmsname);
    b    = fvector(0,n-1,dbmsname);
    c    = fvector(0,n-1,dbmsname);
    r    = fvector(0,n-1,dbmsname);
    nold = n;
  }
  
  if (A[0] > .99e+20) {
    /* No-flux boundary condition. */
    coeff3 = 0.;
  }
  else {
    coeff3 = mu[0]/(z[1]-z[0]);
  }
  for (j = 1; j <= n; j++) {
    coeff1 = coeff3;
    if (j == 1) {
      coeff3 = mu[j]/(z[j+1]-z[j]);
      if (rho) {
        factor = dt/(rho[j]*(z[j+1]-z[j-1]));
      }
      else {
        factor = dt/(z[j+1]-z[j-1]);
      }
      b[j-1] = 1.+(coeff1+coeff3)*factor;
      c[j-1] = -coeff3*factor;
      r[j-1] = 2.*A[j-1]*coeff1*factor
                 +A[j  ]*(1.-(coeff1+coeff3)*factor)
                 +A[j+1]*coeff3*factor;
    }
    else if (j == n) {
      if (A[n+1] > .99e+20) {
        /* No-flux boundary condition. */
        coeff3 = 0.;
      }
      else {
        coeff3 = mu[j]/(z[j+1]-z[j]);
      }
      if (rho) {
        factor = dt/(rho[j]*(z[j+1]-z[j-1]));
      }
      else {
        factor = dt/(z[j+1]-z[j-1]);
      }
      a[j-1] = -coeff1*factor;
      b[j-1] = 1.+(coeff1+coeff3)*factor;
      r[j-1] =    A[j-1]*coeff1*factor
                 +A[j  ]*(1.-(coeff1+coeff3)*factor)
              +2.*A[j+1]*coeff3*factor;
    }
    else {
      coeff3 = mu[j]/(z[j+1]-z[j]);
      if (rho) {
        factor = dt/(rho[j]*(z[j+1]-z[j-1]));
      }
      else {
        factor = dt/(z[j+1]-z[j-1]);
      }
      a[j-1] = -coeff1*factor;
      b[j-1] = 1.+(coeff1+coeff3)*factor;
      c[j-1] = -coeff3*factor;
      r[j-1] = A[j-1]*coeff1*factor
              +A[j  ]*(1.-(coeff1+coeff3)*factor)
              +A[j+1]*coeff3*factor;
    }
  }

  tridiag(n,a,b,c,r,ANS+1,WITH_PIVOTING);

  return;
}

/*======================= end of crank_nicolson() ============================*/

/*======================= hqr() ==============================================*/

/* 
 * From Numerical Recipes in C, 2nd Ed., p. 491-492.
 * Finds eigenvalues of a real Hessenberg matrix. 
 * Assumes zero-based indexing.
 */

void hqr(FLOAT *a_eig, 
         int    n, 
         FLOAT *wr, 
         FLOAT *wi)
{
  int 
    nn,m,l,k,j,
    klen,its,i,mmin;
  FLOAT 
    z,y,x,w,v,u,
    t,s,r,q,p,
    anorm;

  /* Need klen for A_EIG() macro. */
  klen  = n;

  anorm = fabs(A_EIG(0,0));
  for (i = 1; i < n; i++) {
    for (j= i-1; j < n; j++) {
      anorm += fabs(A_EIG(i,j));
    }
  }

  nn = n-1;
  t  = 0.;
  while (nn >= 0) {
    its = 0;
    do {
      for (l = nn; l >= 1; l--) {
	s = fabs(A_EIG(l-1,l-1))+fabs(A_EIG(l,l));
	if (s == 0.0) {
          s = anorm;
        }
	if ((fabs(A_EIG(l,l-1))+s) == s) {
          break;
        }
      }
      x = A_EIG(nn,nn);
      if (l == nn) {
	wr[nn  ] = x+t;
	wi[nn--] = 0.0;
      } 
      else {
	y = A_EIG(nn-1,nn-1);
	w = A_EIG(nn,nn-1)*A_EIG(nn-1,nn);
	if (l == nn-1) {
	  p = 0.5*(y-x);
	  q = p*p+w;
	  z = sqrt(fabs(q));
	  x += t;
	  if (q >= 0.) {
	    z        = p+NR_SIGN(z,p);
	    wr[nn-1] = wr[nn]=x+z;
	    if (z) {
              wr[nn] = x-w/z;
            }
	    wi[nn-1] = wi[nn] = 0.0;
	  } 
          else {
	    wr[nn-1] =   wr[nn] = x+p;
	    wi[nn-1] = -(wi[nn] = z);
	  }
	  nn -= 2;
	} 
        else {
	  if (its == 30){
	    fprintf(stderr, "Warning: hqr(): too many iterations\n");
	    return;
	  }
	  if (its == 10 || its == 20) {
	    t += x;
	    for (i = 0; i < nn; i++) {
              A_EIG(i,i) -= x;
            }
	    s = fabs(A_EIG(nn,nn-1))+fabs(A_EIG(nn-1,nn-2));
	    y = x = 0.75*s;
	    w = -0.4375*s*s;
	  }
	  its++;
	  for (m = nn-2; m >= 0; m--) {
	    z = A_EIG(m,m);
	    r = x-z;
	    s = y-z;
	    p = (r*s-w)/A_EIG(m+1,m)+A_EIG(m,m+1);
	    q = A_EIG(m+1,m+1)-z-r-s;
	    r = A_EIG(m+2,m+1);
	    s = fabs(p)+fabs(q)+fabs(r);
	    p /= s;
	    q /= s;
	    r /= s;
	    if (m == l) {
              break;
            }
	    u = fabs(A_EIG(m,m-1))*(fabs(q)+fabs(r));
	    v = fabs(p)*(fabs(A_EIG(m-1,m-1))+fabs(z)+fabs(A_EIG(m+1,m+1)));
	    if ((u+v) == v) {
              break;
            }
	  }
	  for (i = m+2; i <= nn; i++) {
	    A_EIG(i,i-2) = 0.;
	    if (i != (m+2)) {
              A_EIG(i,i-3) = 0.;
            }
	  }
	  for (k = m; k <= nn-1; k++) {
	    if (k != m) {
	      p = A_EIG(k,  k-1);
	      q = A_EIG(k+1,k-1);
	      r = 0.0;
	      if (k != nn-1) {
                r = A_EIG(k+2,k-1);
              }
	      if ((x = fabs(p)+fabs(q)+fabs(r)) != 0.) {
		p /= x;
		q /= x;
		r /= x;
	      }
	    }
	    if ((s = NR_SIGN(sqrt(p*p+q*q+r*r),p)) != 0.) {
	      if (k == m) {
		if (l != m) {
		  A_EIG(k,k-1) = -A_EIG(k,k-1);
                }
	      } 
              else {
		A_EIG(k,k-1) = -s*x;
              }
	      p += s;
	      x  = p/s;
	      y  = q/s;
	      z  = r/s;
	      q /= p;
	      r /= p;
	      for (j = k; j <= nn; j++) {
		p = A_EIG(k,j)+q*A_EIG(k+1,j);
		if (k != nn-1) {
		  p            += r*A_EIG(k+2,j);
		  A_EIG(k+2,j) -= p*z;
		}
		A_EIG(k+1,j) -= p*y;
		A_EIG(k,  j) -= p*x;
	      }
	      mmin = (nn < k+3) ? nn : k+3;
	      for (i = l; i <= mmin; i++) {
		p = x*A_EIG(i,k)+y*A_EIG(i,k+1);
		if (k != nn-1) {
		  p            += z*A_EIG(i,k+2);
		  A_EIG(i,k+2) -= p*r;
		}
		A_EIG(i,k+1) -= p*q;
		A_EIG(i,k  ) -= p;
	      }
	    }
	  }
	}
      }
    } while (l < nn-1);
  }
}

/*======================= end of hqr() =======================================*/

/*======================= quicksort() ========================================*/
/* 
 * From K & R 2nd ed., p. 87: 
 */
void quicksort(FLOAT *mag, 
               int    left, 
               int    right)
{
  int 
    i,last;
  void 
    swap(FLOAT *mag,int i,int j);

  if (left >= right) {
    return;
  }

  swap(mag,left,(left+right)/2);
  last = left;
  for (i = left+1; i <= right; i++) {
    if (mag[i] > mag[left]) {
      swap(mag, ++last, i);
    }
  }
  swap(mag,left,last);

  quicksort(mag,left,  last-1);
  quicksort(mag,last+1,right);

  return;
}

/*====================== end of quicksort() =================================*/

/*====================== swap() =============================================*/

void swap(FLOAT *mag,
          int    i,
          int    j)
{
  FLOAT 
    temp;

  temp   = mag[i];
  mag[i] = mag[j];
  mag[j] = temp;

  return;
}

/*======================= end of swap() ======================================*/

/*======================= four1() ============================================*/

      /*
       *  FFT routine.
       *  Numerical Recipes in C, 2nd ed, p. 507.
       *  Assumes data length is a power of two.
       */

#undef  SWAP
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(FLOAT         data[], 
           unsigned long nn, 
           int           isign)
{
  unsigned long 
    n,mmax,m,j,istep,i;
  FLOAT 
    tempr,tempi;
  FLOAT 
    wtemp,wr,wpr,wpi,wi,theta;

  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i+= 2) {
    if (j > i) {
       SWAP(data[j  ],data[i  ]);
       SWAP(data[j+1],data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(2*M_PI/mmax);
    wtemp = sin(0.5*theta);
    wpr   = -2.0*wtemp*wtemp;
    wpi   = sin(theta);
    wr    = 1.0;
    wi    = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr      = wr*data[j]-wi*data[j+1];
        tempi      = wr*data[j+1]+wi*data[j];
        data[j  ]  = data[i  ]-tempr;
        data[j+1]  = data[i+1] - tempi;
        data[i  ] += tempr;
        data[i+1] += tempi;
      }
      wr = (wtemp = wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}
      
/*======================= end of four1() ====================================*/

/*======================= realft() ==========================================*/

      /*
       *  Real FFT routine.
       *  Numerical Recipes in C, 2nd ed, p. 513.
       *  Assumes data length is a power of 2.
       *  Result of inverse must be multiplied by 2/n.
       */
void realft(FLOAT data[], unsigned long n, int isign)
{
  void 
    four1(FLOAT data[], unsigned long nn, int isign);
  unsigned long 
    i, i1, i2, i3, i4, np3;
  FLOAT 
    c1=0.5,c2,h1r,h1i,h2r,h2i;
  FLOAT 
    wr,wi,wpr,wpi,wtemp,theta;

  theta = M_PI/(FLOAT) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, n>>1, 1);
  } 
  else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5*theta);
  wpr   = -2.0*wtemp*wtemp;
  wpi   = sin(theta);
  wr    = 1.0+wpr;
  wi    = wpi;
  np3   = n+3;
  for (i = 2; i <= (n>>2); i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r= -c2*(data[i2]+data[i4]);
    h2i = c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4]= -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r - data[2];
  }
  else {
    data[1] = c1*((h1r=data[1])+data[2]);
    data[2] = c1*(h1r-data[2]);
    four1(data, n>>1, -1);
  }
}

/*======================= end of realft() ====================================*/

/*======================= c_num() ============================================*/

complex c_num(FLOAT x, FLOAT y)
{
  complex
    ans;

  ans.x = x;
  ans.y = y;

  return ans;
}

/*======================= end of c_num() =====================================*/

/*======================= c_mult() ===========================================*/

complex c_mult(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)*(z2.x)-(z1.y)*(z2.y);
  ans.y = (z1.x)*(z2.y)+(z1.y)*(z2.x);

  return ans;
}
/*======================= end of c_mult() ====================================*/

/*======================= c_add() ============================================*/

complex c_add(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)+(z2.x);
  ans.y = (z1.y)+(z2.y);

  return ans;
}

/*======================= end of c_add() =====================================*/

/*======================= c_sub() ============================================*/

complex c_sub(complex z1,complex z2)
{
  complex
    ans;

  ans.x = (z1.x)-(z2.x);
  ans.y = (z1.y)-(z2.y);

  return ans;
}

/*======================= end of c_sub() =====================================*/

/*======================= c_exp() ============================================*/

complex c_exp(complex z)
{
  complex 
    ans;

  ans.x = exp(z.x)*cos(z.y);
  ans.y = exp(z.x)*sin(z.y);

  return ans;
}

/*======================= end of c_exp() =====================================*/

/*======================= c_abs() ============================================*/

/*
 * NOTE: For LINUX with -D_BSD_SOURCE, cabs() is defined, such that a 
 * type-mismatch error occurs if we call this function cabs().
 */
 
FLOAT c_abs(complex z)
{

  return sqrt((z.x)*(z.x)+(z.y)*(z.y));
}

/*======================= end of c_abs() =====================================*/

/*======================= c_real() ===========================================*/

FLOAT c_real(complex z) 
{
  return (z.x);
}

/*======================= end of c_real() ====================================*/

/*======================= c_imag() ===========================================*/

FLOAT c_imag(complex z)
{
  return (z.y);
}

/*======================= end of c_imag() ====================================*/

/*======================= fcmp() =============================================*/

/*
 * Derived from fcmp(), version 1.2.2, 
 * Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * <mailto:Ted.Belding@umich.edu>
 * <http://fcmp.sourceforge.net>
 *
 * The major modification we have made is to remove the "epsilon" argument
 * and set epsilon inside the fcmp() function.
 *
 * Description:
 *   It is generally not wise to compare two floating-point values for
 *   exact equality, for example using the C == operator.  The function
 *   fcmp() implements Knuth's suggestions for safer floating-point
 *   comparison operators, from:
 *   Knuth, D. E. (1998). The Art of Computer Programming.
 *   Volume 2: Seminumerical Algorithms. 3rd ed. Addison-Wesley.
 *   Section 4.2.2, p. 233. ISBN 0-201-89684-2.
 *
 * Input parameters:
 *   x1, x2: numbers to be compared
 *
 * This routine may be used for both single and double precision.
 *
 * Returns:
 *   -1 if x1 < x2
 *    0 if x1 == x2
 *    1 if x1 > x2		
 */

int fcmp(double x1, double x2) {
  int 
    exponent;
  double
    delta,
    difference;
#if EPIC_PRECISION == DOUBLE_PRECISION
  const double
    epsilon = DBL_EPSILON;
#else
  const double
    epsilon = FLT_EPSILON;
#endif
  
  /* 
   * Get exponent(max(fabs(x1),fabs(x2))) and store it in exponent. 
   *
   * If neither x1 nor x2 is 0,
   * this is equivalent to max(exponent(x1),exponent(x2)).
   *
   * If either x1 or x2 is 0, its exponent returned by frexp would be 0,
   * which is much larger than the exponents of numbers close to 0 in
   * magnitude. But the exponent of 0 should be less than any number
   * whose magnitude is greater than 0.
   *
   * So we only want to set exponent to 0 if both x1 and x2 are 0. 
   * Hence, the following works for all x1 and x2. 
   */
  frexp(fabs(x1) > fabs(x2) ? x1 : x2,&exponent);

  /* 
   * Do the comparison.
   *
   * delta = epsilon*pow(2,exponent)
   *
   * Form a neighborhood around x2 of size delta in either direction.
   * If x1 is within this delta neighborhood of x2, x1 == x2.
   * Otherwise x1 > x2 or x1 < x2, depending on which side of
   * the neighborhood x1 is on.
   */
  delta      = ldexp(epsilon,exponent); 
  difference = x1-x2;

  if (difference > delta) {
    /* x1 > x2 */
    return 1;
  }
  else if (difference < -delta) {
    /* x1 < x2 */
    return -1;
  }
  else  {
    /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
  }
}

/*======================= end of fcmp() ======================================*/

/*======================= least_squares() ====================================*/

/*
 * Code for fitting data to a straight line using the least squares technique.
 *   Csaba J. Palotai *A*  4/14/2004
 * Adapted from "Numerical recipes in C," p665.
 *
 * Given a set of data points x[0..n-1],y[0..n-1] the code fits them to a 
 * straight line y= a+bx by minimizing X^2. Returned are a,b.
 *
 * Commented out is code that can be used to calculate the uncertainties 
 * siga and sigb, and the chi-square, chi2.
 */

void least_squares(FLOAT *x,
                   FLOAT *y,
                   int    n,
                   FLOAT *a)

{
  
  register int
    i;
  register FLOAT 
    siga,sigb,chi2,q,t,
    sxoss,ss,sigdat,
    sx  = 0.0,
    sy  = 0.0,
    st2 = 0.0;

  a[1] = 0.0;

  for (i = 0; i < n; i++) {
    sx += x[i];
    sy += y[i];
  }
  
  ss    = n;
  sxoss = sx/ss;
  
  for (i = 0; i < n; i++) {
    t     = x[i]-sxoss;
    st2  += t*t;
    a[1] += t*y[i];
  }
  
  a[1] /= st2;
  a[0]  = (sy-sx*a[1])/ss;

  /*  
   * Code for chi2 and uncertainties for a and b.
   *
  siga = sqrt((1.+sx*sx/(ss*st2))/ss);
  sigb = sqrt(1./st2);
  chi2 = 0.0;
  q    = 1.0;
  for (i = 0; i < n; i++) {
    t      = (y[i]-(*a)-(*b)*x[i]);
    chi2  += t*t;
    sigdat = sqrt((chi2)/(n-2));
    siga  *= sigdat;
    sigb  *= sigdat;
  }
   *
   *
   */

  return;
}  

/*======================= end of least_squares() =============================*/

/*======================= savitzky_golay() ===================================*/

/*
 * Calculate Savitzky-Golay weighting coefficients to compute a smooth value
 * or smooth derivative of noisy data.
 *
 * Based on Numerical Recipes in C, p. 652.
 * Note the wrap-around ordering of the coefficients returned in c[].
 *
 * Assumes zero-based arrays.
 */

#undef  A
#define A(i,j) a[j+(m+1)*i]
 
void savitzky_golay(FLOAT *c,
                    int    np,
                    int    nl,
                    int    nr,
                    int    ld,
                    int    m)
{
  int
    imj,ipj,j,k,kk,mm;
  int
   *index;
  FLOAT
    d,fac,sum;
  FLOAT
    *a,
    *b;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="savitsky_golay";

  /*
   * Screen for inconsistent arguments.
   */
  if (np < nl+nr+1) {
    fprintf(stderr,"**error: %s: np=%d < nl+nr+1=%d",np,nl+nr+1);
    exit(1);
  }
  else if (nl < 0) {
    fprintf(stderr,"**error: %s: nl=%d < 0",dbmsname,nl);
    exit(1);
  }
  else if (nr < 0) {
    fprintf(stderr,"**error: %s: nr=%d < 0",dbmsname,nr);
    exit(1);
  }
  else if (ld > m) {
    fprintf(stderr,"**error: %s: ld=%d > m=%d",dbmsname,ld,m);
    exit(1);
  }
  else if (nl+nr < m) {
    fprintf(stderr,"**error: %s: nl+nr=%d < m=%d",dbmsname,nl+nr,m);
    exit(1);
  }

  /* Allocate memory */
  index = ivector(0,m,dbmsname);
  a     = fvector(0,(m+1)*(m+1)-1,dbmsname);
  b     = fvector(0,m,dbmsname);

  for (ipj = 0; ipj <= (m << 1); ipj++) {
    sum = (ipj ? 0. : 1.);
    for (k = 1; k <= nr; k++) {
      sum += pow((double)k,(double)ipj);
    }
    for (k = 1; k <= nl; k++) {
      sum += pow((double)-k,(double)ipj);
    }
    mm = IMIN(ipj,2*m-ipj);
    for (imj = -mm; imj <= mm; imj+=2) {
      A((ipj+imj)/2,(ipj-imj)/2) = sum;
    }
  }

  lu_decompose(m+1,a,index,&d);

  for (j = 0; j < m+1; j++) {
    b[j] = 0.;
  }
  b[ld] = 1.0;

  lu_backsub(m+1,a,index,b);

  for (kk = 0; kk < np; kk++) {
    c[kk] = 0.;
  }

  for (k = -nl; k <= nr; k++) {
    sum = b[0];
    fac = 1.;
    for (mm = 0; mm < m; mm++) {
      sum += b[mm+1]*(fac *= k);
    }
    kk = ((np-k)%np);
    c[kk] = sum;
  }

  /* Free allocated memory. */
  free_ivector(index,0,m,dbmsname);
  free_fvector(a,0,(m+1)*(m+1)-1,dbmsname);
  free_fvector(b,0,m,dbmsname);

  return;
}

/*======================= end of savitzky_golay() ============================*/

/*======================= random_number() ====================================*/

/*
 * Based on ran1() in Numerical Recipes in C, p. 280.
 * 
 * Returns a uniform deviate between 0.0 and 1.0 (exclusive of the endpoint values).
 *
 * Call with idum a negative number to initialize; thereafter, do not alter
 * idum between successive deviates in a sequence.
 */

#undef  IA
#define IA   16807
#undef  IM
#define IM   2147483647
#undef  AM
#define AM   (1./IM)
#undef  IQ
#define IQ   127773
#undef  IR
#define IR   2836
#undef  NTAB
#define NTAB 32
#undef  NDIV
#define NDIV (1+(IM-1)/NTAB)
#undef  EPS
#define EPS   1.2e-7
#undef  RNMX
#define RNMX (1.-EPS)

FLOAT random_number(long *idum)
{
  int
    j;
  long
    k;
  static long
    iy=0,
    iv[NTAB];
  FLOAT
    temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k     = (*idum)/IQ;
      *idum = IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j     = iy/NDIV;
  iy    = iv[j];
  iv[j] = *idum;
  if ((temp = AM*iy) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}

/*======================= end of random_number() =============================*/

#undef DEG
#define DEG (M_PI/180.)

/*======================= lat_centric_to_graphic() ===========================*/

/*
 * lat: [Deg]
 * rerp: (equatorial radius)/(polar radius)
 */
EPIC_FLOAT lat_centric_to_graphic(EPIC_FLOAT lat,
                                  EPIC_FLOAT rerp)
{
  return (fabs(lat) == 90.) ? lat : atan(rerp*rerp*tan(lat*DEG))/DEG;
}

/*======================= end of lat_centric_to_graphic() ====================*/

/*======================= lat_graphic_to_centric() ===========================*/

/*
 * lat: [Deg]
 * rerp: (equatorial radius)/(polar radius)
 */
EPIC_FLOAT lat_graphic_to_centric(EPIC_FLOAT lat,
                                  EPIC_FLOAT rerp)
{
  return (fabs(lat) == 90.) ? lat : atan(tan(lat*DEG)/(rerp*rerp))/DEG;
}

/*======================= end of lat_graphic_to_centric() ====================*/


/* * * * * * * * * * *  end of epic_funcs_util.c  * * * * * * * * * * * * * * */

