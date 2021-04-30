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

/*---------------------------------------------------------------------------*
 * Logic:
 *
 * 1. Node 0 reads row.
 * 2. Node 0 broadcasts row.
 * 3. Node i!=0 copy portion of row to local array, if applicable.
 * 
 *
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "qa.h"

#define M_NEED_ROW 240

extern void  index_row (int row, int *length, int *start, int dim, int *index);
extern void  offset_index (int *end, int *start, int *pad, int *index, int dim,
			   int *offset);
extern void  map_error (const char *s);

void  MPG_Cart_wslab (
     FILE  *stream,                        /* stream from which to read data */
     MPI_Comm  comm_cart,                 /* communicator for cartesian grid */
     int  ndim,                              /* dimensions of Cartesian grid */
     int  *slab_hi,                       /* index of maximum value to write */
     int  *slab_lo,                       /* index of minimum value to write */
     int  *length,                     /* number of elements along each axis */
     int  *pad,
     MPI_Datatype  etype,           /* data type of grid elements to be read */
     void  *local                            /* local array for storing data */
		    )
{
  int                           inode, nnode;       /* ith node, total nodes */
  int                           idim;                       /* ith dimension */
  MPI_Aint                      tsize;              /* stores data type size */
  int                           esize;       /* stores data type size as int */
  int                           irow, nrow;  /* ith row, number of grid rows */
  int                           row_node, row_nodes, tag;
  int                           ne;        /* number of elements in grid row */
  MPI_Status                    stat;
  char                         *buffer;                /* buffer to hold row */
  int                          *periods, *dims, *coords;       /* dummy ptrs */
  int                          *index;      /* index of starting row element */
  int                          *start, *end;        /* start and end indices */
  int                           my_row;                    /* row on my node */
  int                           my_slice, slice_offset, slice_row, count,
                                orow=0;
  int                           offset, size, i;
  double                        t1, t2;

  for (my_slice=0, idim=0; idim<ndim; idim++)
    if (slab_hi[idim] >= slab_lo[idim] &&
	slab_lo[idim] >= 0 && slab_hi[idim] < length[idim])
      { my_slice = 1;  break; }
  if (!my_slice)
    map_error ("no slice selected");

  MPI_Comm_rank (comm_cart, &inode);
  MPI_Comm_size (comm_cart, &nnode);

  MPI_Type_size (etype, &tsize);
  esize = (int)tsize;

  if ((start = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (end = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (index = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (dims = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (periods = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (coords = (int *)malloc(ndim*sizeof(int))) == NULL)
      map_error ("Cannot allocate space");

  ne = slab_hi[0] - slab_lo[0] + 1;
  MPI_Cart_get (comm_cart, ndim, dims, periods, coords);
  for (nrow=1, idim=0; idim<ndim; idim++) {
    if (idim) nrow *= length[idim];
    MPE_Decomp1d (length[idim], dims[idim], coords[idim],
		  &start[idim], &end[idim]);
    start[idim]--;  end[idim]--;
  }

  row_nodes = dims[0];  row_node = coords[0];
/*
  fprintf (stderr, "node %d has row_node %d of %d\n",
	   inode, row_node, row_nodes);
 */

/*- Allocate row buffer -----------------------------------------------------*/

  if ((buffer = malloc(ne*esize)) == NULL)
    map_error ("Cannot allocate space");
  
  MPI_Barrier (comm_cart);
  t1 = MPI_Wtime();

  for (irow=0; irow<nrow; irow++) {
    
/*- Node i determines whether it has row ------------------------------------*/
/*- Node i determines whether row intersects slice --------------------------*/
/*- Node 0 determines whether node i has slice ------------------------------*/

    index_row (irow, length, start, ndim, index);
  
    for (my_row=1, idim=1; idim<ndim; idim++)
      if (index[idim] < start[idim] || index[idim] > end[idim]) {
	my_row = 0;  break;
      }
  
    for (slice_row=1, idim=1; idim<ndim; idim++)
      if (index[idim] < slab_lo[idim] || index[idim] > slab_hi[idim]) {
	slice_row = 0;  break;
      }
  
    if (slice_row &&
	(slab_hi[0] >= start[0] && slab_lo[0] <= end[0]))
      my_slice = my_row;
    else my_slice = 0;

    if (slice_row) {

      MPI_Barrier (comm_cart);

/*- Node i identifies portion of local memory belonging to row. -------------*/

      count = offset = slice_offset = 0;
      if (my_slice) {
	/* debug ('i', "node %d has slice\n", inode); */
	offset_index (end, start, pad, index, ndim, &offset);
	slice_offset = slab_lo[0]-start[0];
	if (slice_offset < 0) slice_offset = 0;
	if (slab_hi[0] > end[0]) count = end[0]-(start[0]+slice_offset)+1;
	else count = slab_hi[0]-(start[0]+slice_offset)+1;
      }

/*- On node !=0 if my_row and my_slice---------------------------------------*/

      if (inode) {
	if (my_slice) {
	  /*- Send portion size to node 0 with message tag = row position ---*/
	  tag = row_node;
	  MPI_Send (&count, 1, MPI_INT, 0, tag, comm_cart);
	  /*- Node !=0 sends portion to node 0 w/tag = nnodes + row position */
	  tag += row_nodes;
	  MPI_Send ((char *)local + (offset+slice_offset)*esize, count, etype,
		    0, tag, comm_cart);
	}
      }

/*- On node 0, for each of nnodes messages containing portion sizes ---------*/

      if (inode == 0) {
	if (my_slice) {
	  if (fwrite ((char *)local + (offset+slice_offset)*esize, esize,
		      count, stream) != count)
	    map_error ("Cannot write expected number of elements");
	}
	for (i=my_row; i<row_nodes; i++) {
	  /*- Read portion size message -------------------------------------*/
	  tag = i;
	  MPI_Recv (&count, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm_cart,
		    &stat);
	  /*- Read portion message ------------------------------------------*/
	  tag += row_nodes;
	  MPI_Recv (buffer, count, etype, MPI_ANY_SOURCE, tag, comm_cart,
		    &stat);
	  /*- Write portion to file -----------------------------------------*/
	  if (fwrite (buffer, esize, count, stream) != count)
	    map_error ("Cannot write expected number of elements");
	}
      }

      if (inode == 0) {
	debug ('i', "%d of ", ++orow);
	debug ('i', "%d rows written\r", nrow);
      }
    }
  }

  if (inode == 0) {
    debug ('i', "%s", "\n");
  }

  MPI_Barrier (comm_cart);
  t2 = MPI_Wtime();

  if (inode == 0) {
    debug ('t', "MPI_Cart_write:  elapsed time = %g\n", t2-t1);
  }

  (void)free ((char *)buffer);

  (void)free ((char *)start);
  (void)free ((char *)end);
  (void)free ((char *)index);
  (void)free ((char *)dims);
  (void)free ((char *)periods);
  (void)free ((char *)coords);

  return;
}
