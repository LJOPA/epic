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

/* map_decomp.c */
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
#define mpg_cart_decomp_ PMPG_CART_DECOMP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_decomp_ pmpg_cart_decomp__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_decomp_ pmpg_cart_decomp
#else
#define mpg_cart_decomp_ pmpg_cart_decomp_
#endif
#else
#ifdef FORTRANCAPS
#define mpg_cart_decomp_ MPG_CART_DECOMP
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_decomp_ mpg_cart_decomp__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_decomp_ mpg_cart_decomp
#endif
#endif

 void mpg_cart_decomp_ ( comm_old, ndims, length, dims, periods, pad,
			start, end, stride, comm_cart, ierr )
MPI_Comm          comm_old;
int              *ndims;
int              *length;
int              *dims;
int              *periods;
int              *pad;
int              *start;
int              *end;
int              *stride;
MPI_Comm         *comm_cart;
int              *ierr;
{
  int lperiods[20], i;
  MPI_Comm lcomm_old = (MPI_Comm) MPIR_ToPointer( *((int *)comm_old) );
  MPI_Comm lcomm_cart;
  
  if (*ndims > 20) {
    *ierr = MPIR_ERROR( lcomm_old, MPI_ERR_LIMIT, "Too many dimensions" );
    return;
  }
  for (i=0; i<*ndims; i++) 
    lperiods[i] = MPIR_FROM_FLOG(periods[i]);

  MPG_Cart_decomp( lcomm_old, *ndims, length, dims,
		  lperiods, pad, start, end,
		  stride, &lcomm_cart);
  *ierr = 0;

  if (*ierr == 0) 
    *(int *)comm_cart = MPIR_FromPointer( lcomm_cart );

}
