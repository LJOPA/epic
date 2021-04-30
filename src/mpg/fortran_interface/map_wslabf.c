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

/* map_wslab.c */
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
#define mpg_cart_wslab_ PMPG_CART_WSLAB
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_wslab_ pmpg_cart_wslab__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_wslab_ pmpg_cart_wslab
#else
#define mpg_cart_wslab_ pmpg_cart_wslab_
#endif
#else
#ifdef FORTRANCAPS
#define mpg_cart_wslab_ MPG_CART_WSLAB
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_wslab_ mpg_cart_wslab__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_wslab_ mpg_cart_wslab
#endif
#endif

void mpg_cart_wslab_ ( file, comm_cart, ndims, slab_hi, slab_lo, length,
		      pad, datatype, buffer, ierr )
const char       *file;
MPI_Comm          comm_cart;
int              *ndims;
int              *slab_hi;
int              *slab_lo;
int              *length;
int              *pad;
MPI_Datatype      datatype;
void             *buffer;
int              *ierr;
{
  FILE  *fp;
  int    eflag, myid;
  MPI_Comm  lcomm_cart = (MPI_Comm) MPIR_ToPointer( *((int *)comm_cart) );

  MPI_Comm_rank (lcomm_cart, &myid);

  eflag = 0;
  if (myid == 0 && (fp=fopen(file,"w")) == NULL) {
    fprintf (stderr, "cannot open file %s\n", file);
    eflag = 1;
  }
  MPI_Bcast (&eflag, 1, MPI_INT, 0, lcomm_cart);
  if (eflag) MPI_Abort (MPI_COMM_WORLD, 0);
  
  MPG_Cart_wslab( fp, lcomm_cart, *ndims, slab_hi, slab_lo, length, pad,
		 (MPI_Datatype)MPIR_ToPointer( *(int*)(datatype) ),
		 MPIR_F_PTR(buffer));

  eflag = 0;
  if (myid == 0 && fclose(fp) == EOF) {
    fprintf (stderr, "cannot close file %s\n", file);
    eflag = 1;
  }
  MPI_Bcast (&eflag, 1, MPI_INT, 0, lcomm_cart);
  if (eflag) MPI_Abort (MPI_COMM_WORLD, 0);

  *ierr = 0;

}
