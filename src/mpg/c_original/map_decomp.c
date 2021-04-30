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

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "qa.h"

#pragma weak MPG_Cart_decomp = _MPG_Cart_decomp
#pragma weak MPG_CART_DECOMP = _MPG_CART_DECOMP

#define TRUE 1

void balance (int numproc, int dim, int *dims, int *length);
void factor (int numproc, int *nf, int **factors);
int  comp (const int *i, const int *j);

void  _MPG_Cart_decomp (MPI_Comm comm, int dim, int *length, int *dims, 
		       int *periods, int *pad, int *start, int *end,
		       int *stride, MPI_Comm *comm_cart) {
  int  i, *coords, numproc;

  if ((coords = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr, "Cannot allocate space\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

  MPI_Comm_size (comm, &numproc);

  balance (numproc, dim, dims, length);
  for (i=0; i<dim; i++) {
    debug ('v', "i = %d; ", i);
    debug ('v', "dims[i] = %d\n", dims[i]);
  }

  MPI_Dims_create (numproc, dim, dims);
  MPI_Cart_create (comm, dim, dims, periods, TRUE, comm_cart);
  MPI_Cart_get (*comm_cart, dim, dims, periods, coords);
  for (i=0; i<dim; i++) {
    MPE_Decomp1d (length[i], dims[i], coords[i], &start[i], &end[i]);
    stride[i] = end[i] - start[i] + 1 + 2*pad[i];
    start[i]--;  end[i]--;

    debug ('v', "i = %d; ", i);
    debug ('v', "dims[i] = %d; ", dims[i]);
    debug ('v', "start[i] = %d; ", start[i]);
    debug ('v', "end[i] = %d; ", end[i]);
    debug ('v', "stride[i] = %d\n", stride[i]);
  }

  (void)free ((char *)coords);

  return;
}

void  balance (int numproc, int dim, int *dims, int *length) {
  int  *tmp, i, j, nf, *factors, max_length, max_dim, rule1, rule2;
  int  *dims_copy, numleft;

  if ((tmp = (int *)malloc(dim*sizeof(int))) == NULL ||
      (dims_copy = (int *)malloc(dim*sizeof(int))) == NULL) {
    fprintf (stderr, "Cannot allocate space\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

/*---------------------------------------------------------------------------*
 * if dim[i] >= 0, leave it alone...
 * if dim[i] < 0, "fill it in" according to rule:
 *   rotate bisections (greatest dimlen)
 *---------------------------------------------------------------------------*/

  for (rule1=rule2=0, numleft=numproc, i=0; i<dim; i++) {
    if (dims[i] == 0) {
      if (rule2) {
	fprintf (stderr,
		 "Cannot mix subdivision rules\n");
	MPI_Abort (MPI_COMM_WORLD, 0);
      }
      (void)free ((char *)tmp);
      (void)free ((char *)dims_copy);
      return;
    }
    if (dims[i] == -1) {
      if (rule1) {
	fprintf (stderr,
		 "Cannot mix subdivision rules\n");
	MPI_Abort (MPI_COMM_WORLD, 0);
      }
      rule2 = 1;
    }
    if (dims[i] > 0) numleft /= dims[i];
    if (numleft == 0) {
      fprintf (stderr,
	       "Cannot subdivide grid further; use more processors\n");
      MPI_Abort (MPI_COMM_WORLD, 0);
    }
  }

  for (i=0; i<dim; i++) {
    dims_copy[i] = dims[i];
    if (dims[i] == -1) { dims[i] = 1;  tmp[i] = length[i]; }
  }
  factor (numleft, &nf, &factors);
  for (j=0; j<nf; j++) {
    for (max_length=0, i=0; i<dim; i++) {
      if (dims_copy[i] == -1 && max_length <= tmp[i]) {
	max_length = tmp[i];
	max_dim = i;
      }
    }
    tmp[max_dim] /= factors[j];
    if (tmp[max_dim] < 1) {
      fprintf (stderr,
	       "Cannot subdivide grid further; use fewer processors\n");
      MPI_Abort (MPI_COMM_WORLD, 0);
    }
    dims[max_dim] *= factors[j];
  }

  (void)free ((char *)tmp);
  (void)free ((char *)dims_copy);

  return;
}

void  factor (int num, int *nf, int **factors) {
  int  tmp, i, divisor;

  for (tmp=num, i=0; tmp>2; i++) {
    for (divisor=2; divisor<tmp/2 && tmp%divisor; divisor++);
    if (divisor > tmp/2) break;
    else tmp /= divisor;
  }

  *nf = i+1;
  if ((*factors = (int *)malloc((*nf)*sizeof(int))) == NULL) {
    fprintf (stderr, "Cannot allocate space\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

  for (tmp=num, i=0; tmp>2; i++) {
    for (divisor=2; divisor<tmp/2 && tmp%divisor; divisor++);
    if (divisor > tmp/2) break;
    else {
      tmp /= divisor;
      (*factors)[i] = divisor;
    }
  }
  (*factors)[i] = tmp;
  
  qsort (*factors, *nf, sizeof(int), comp);

  return;
}

int comp(const int *i,const int *j) {
  return (*j - *i);
}

