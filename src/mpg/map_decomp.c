/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998 Joseph Matarese                              *
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

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "qa.h"

#define TRUE 1

/* Function prototypes: */
void balance (int numproc, 
              int dim, 
              int *dims, 
              int *length);
void factor (int numproc, 
             int *nf, 
             int **factors);
int  comp (const void *i, 
           const void *j);

/*=================== MPG_Cart_decomp() ===========================*/

void  MPG_Cart_decomp(MPI_Comm  comm,
                      int       dim,
                      int      *length,
                      int      *dims, 
		      int      *periods,
                      int      *pad,
                      int      *start,
                      int      *end,
		      int      *stride,
                      MPI_Comm *comm_cart) 
{
  int  
    i,numproc,
    *coords;

  /* Allocate memory for coords: */
  if ((coords = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr,"Cannot allocate space\n");
    MPI_Abort(MPI_COMM_WORLD,0);
  }

  MPI_Comm_size(comm,&numproc);

  balance(numproc,dim,dims,length);

  for (i = 0; i < dim; i++) {
    debug ('v',"i = %d; ",i);
    debug ('v',"dims[i] = %d\n",dims[i]);
  }

  MPI_Dims_create(numproc,dim,dims);
  MPI_Cart_create(comm,dim,dims,periods,TRUE,comm_cart);
  MPI_Cart_get(*comm_cart,dim,dims,periods,coords);

  for (i = 0; i < dim; i++) {
    MPE_Decomp1d(length[i],dims[i],coords[i],&start[i],&end[i]);
    stride[i] = end[i]-start[i]+1+2*pad[i];
    start[i]--;  
    end[  i]--;

    debug ('v',"i = %d; ",               i);
    debug ('v',"dims[i] = %d; ",    dims[i]);
    debug ('v',"start[i] = %d; ",  start[i]);
    debug ('v',"end[i] = %d; ",      end[i]);
    debug ('v',"stride[i] = %d\n",stride[i]);
  }

  (void)free ((char *)coords);

  return;
}

/*=================== end of MPG_Cart_decomp() ====================*/

/*=================== balance() ===================================*/

void balance (int numproc, 
              int dim, 
              int *dims, 
              int *length) 
{
  int  
    i,j,nf,
    max_dim,max_length,numleft,
    *dims_copy,*length_copy,*factors;

  /* Allocate memory for dims_copy,length_copy: */
  if ((dims_copy   = (int *)malloc(dim*sizeof(int))) == NULL ||
      (length_copy = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr, "mpg:map_decomp:cannot allocate memory\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

/*---------------------------------------------------------------------------*
 * if dim[i] == 0, leave it alone
 * if dim[i] <  0, "fill it in" by rotating bisections (greatest dimlen)
 *---------------------------------------------------------------------------*/

  for (i = 0; i < dim; i++) {
    if (dims[i] == 0) {
      /* Free allocated memory: */
      (void)free((char *)dims_copy);
      (void)free((char *)length_copy);
      return;
    }
  }

  numleft = numproc;
  for (i = 0; i < dim; i++) {
    dims_copy[  i] = dims[  i];
    length_copy[i] = length[i];
    if (dims[i] < 0) {
      dims[i] = 1;
    }
    numleft /= dims[i];
    if (numleft < 1) {
      fprintf(stderr,"mpg:map_decomp:cannot subdivide grid further\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  }

  factor(numleft,&nf,&factors);

  for (j = 0; j < nf; j++) {
    max_dim    = -1;  
    max_length =  0;
    for (i = 0; i < dim; i++) {
      if (dims_copy[i] < 0 && max_length <= length_copy[i]) {
        max_dim    = i;
	max_length = length_copy[i];
      }
    }
    if (max_dim != -1) {
      length_copy[max_dim] /= factors[j]; 
      if (length_copy[max_dim] < 1) {
        fprintf (stderr,"mpg:map_decomp:cannot subdivide grid further\n");
        MPI_Abort(MPI_COMM_WORLD,0);
      }
      dims[max_dim] *= factors[j];
    }
  }

  /* Free allocated memory: */
  (void)free((char *)dims_copy);
  (void)free((char *)length_copy);

  return;
}

/*=================== end of balance() ============================*/

/*=================== factor() ====================================*/

void  factor (int num, 
              int *nf, 
              int **factors) {
  int  
    tmp,i,divisor;

  tmp = num;
  for (i = 0; tmp > 2; i++) {
    for (divisor = 2; divisor < tmp/2 && tmp%divisor; divisor++);
    if (divisor > tmp/2) {
      break;
    }
    else {
      tmp /= divisor;
    }
  }

  *nf = i+1;
  /* Allocate memory for factors: */
  if ((*factors = (int *)malloc((*nf)*sizeof(int))) == NULL) {
    fprintf (stderr, "Cannot allocate space\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

  tmp = num;
  for (i = 0; tmp > 2; i++) {
    for (divisor=2; divisor < tmp/2 && tmp%divisor; divisor++);
    if (divisor > tmp/2) {
      break;
    }
    else {
      tmp /= divisor;
      (*factors)[i] = divisor;
    }
  }
  (*factors)[i] = tmp;
  
  qsort (*factors,*nf,sizeof(int),comp);

  return;
}

/*=================== end of factor() =============================*/

/*=================== comp() ======================================*/

int comp(const void *i,
         const void *j) 
{
  return (*(int *)j - *(int *)i);
}

/*=================== end of comp () ==============================*/

