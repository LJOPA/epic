/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998 Joseph Matarese                              *
 *                                                                 *
 * This program is free software; you can redistribute it and/or   *
 * modify it under the terms of the GNU General Public License     *
 * as published by the Free Software Foundation; either version 2  *
 * of the License, or (at your option) any later version.          *
 * A copy of this License is in the file:                          *
 *   $EPIC_PATH/License.txt                                        *
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

/* map_var.c */
/* Custom Fortran interface file */
#include "mpiimpl.h"

#ifdef POINTER_64_BITS
extern void *MPIR_ToPointer();
extern int MPIR_FromPointer();
extern void MPIR_RmPointer();
#else
#define MPIR_ToPointer(a) (a)
#define MPIR_FromPointer(a) (int)(a)
#define MPIR_RmPointer(a)
#endif

#ifdef MPI_BUILD_PROFILING
#ifdef FORTRANCAPS
#define mpg_cart_varoffset_ PMPG_CART_VAROFFSET
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_varoffset_ pmpg_cart_varoffset__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_varoffset_ pmpg_cart_varoffset
#else
#define mpg_cart_varoffset_ pmpg_cart_varoffset_
#endif
#else
#ifdef FORTRANCAPS
#define mpg_cart_varoffset_ MPG_CART_VAROFFSET
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_varoffset_ mpg_cart_varoffset__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_varoffset_ mpg_cart_varoffset
#endif
#endif

#include "qa.h"

 void mpg_cart_varoffset_ ( ndims, start, end, pad,
			    ng, offset, ierr )
int              *ndims;
int              *start;
int              *end;
int              *pad;
int              *ng;
int              *offset;
int              *ierr;
{
  int   size, i;

  *offset = 1;
  for (size=1, i=0; i<*ndims; i++) {
    *offset += pad[i]*size;
    size *= (end[i] - start[i] + 1 + 2*pad[i]);
  }

  debug ('m', "space available = %d; ", *ng);
  debug ('m', "space needed = %d\n", size);

  if (size > *ng) *ierr = 1;
  else *ierr = 0;

}

