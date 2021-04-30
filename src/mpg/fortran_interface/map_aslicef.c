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

/* map_aslice.c */
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
#define mpg_cart_aslice_ PMPG_CART_ASLICE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_aslice_ pmpg_cart_aslice__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_aslice_ pmpg_cart_aslice
#else
#define mpg_cart_aslice_ pmpg_cart_aslice_
#endif
#else
#ifdef FORTRANCAPS
#define mpg_cart_aslice_ MPG_CART_ASLICE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define mpg_cart_aslice_ mpg_cart_aslice__
#elif !defined(FORTRANUNDERSCORE)
#define mpg_cart_aslice_ mpg_cart_aslice
#endif
#endif

 void mpg_cart_aslice_ ( file, comm_cart, ndims, slice, length, pad, datatype,
                      buffer, ierr )
const char       *file;
MPI_Comm          comm_cart;
int              *ndims;
int              *slice;
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
  if (myid == 0 && (fp=fopen(file,"a")) == NULL) {
    fprintf (stderr, "cannot open file %s\n", file);
    eflag = 1;
  }
  MPI_Bcast (&eflag, 1, MPI_INT, 0, lcomm_cart);
  if (eflag) MPI_Abort (MPI_COMM_WORLD, 0);
  
  MPG_Cart_wslice( fp, lcomm_cart, *ndims, slice, length, pad,
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
