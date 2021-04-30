/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *				MPICH COPYRIGHT                                    *
 *                                                                                 *
 * The following is a notice of limited availability of the code, and disclaimer   *
 * which must be included in the prologue of the code and in all source listings   *
 * of the code.                                                                    *
 *                                                                                 *
 * Copyright Notice                                                                *
 *  + 1993 University of Chicago                                                   *
 *  + 1993 Mississippi State University                                            *
 *                                                                                 *
 * Permission is hereby granted to use, reproduce, prepare derivative works, and   *
 * to redistribute to others.  This software was authored by:                      *
 *                                                                                 *
 * Argonne National Laboratory Group                                               *
 *   W. Gropp: (630) 252-4318; FAX: (708) 252-7852; e-mail: gropp@mcs.anl.gov      *
 *   E. Lusk: (630) 252-7852; FAX: (708) 252-7852; e-mail: lusk@mcs.anl.gov        *
 *   Mathematics and Computer Science Division                                     *
 *   Argonne National Laboratory, Argonne IL 60439                                 *
 *                                                                                 *
 * Mississippi State Group                                                         *
 *   N. Doss:  (601) 325-2565; FAX: (601) 325-7692; e-mail: doss@erc.msstate.edu   *
 *   A. Skjellum:(601) 325-8435; FAX: (601) 325-8997; e-mail: tony@erc.msstate.edu *
 *   Mississippi State University, Computer Science Department &                   *
 *   NSF Engineering Research Center for Computational Field Simulation            *
 *   P.O. Box 6176, Mississippi State MS 39762                                     *
 *                                                                                 *
 *			      GOVERNMENT LICENSE                                   *
 *                                                                                 *
 * Portions of this material resulted from work developed under a U.S.             *
 * Government Contract and are subject to the following license: the Government    *
 * is granted for itself and others acting on its behalf a paid-up, nonexclusive,  *
 * irrevocable worldwide license in this computer software to reproduce, prepare   *
 * derivative works, and perform publicly and display publicly.                    *
 *                                                                                 *
 * 				  DISCLAIMER                                       *
 *                                                                                 *
 * This computer code material was prepared, in part, as an account of work        *
 * sponsored by an agency of the United States Government.  Neither the United     *
 * States, nor the University of Chicago, nor Mississippi State University, nor    *
 * any of their employees, makes any warranty express or implied, or assumes any   *
 * legal liability or responsibility for the accuracy, completeness, or            *
 * usefulness of any information, apparatus, product, or process disclosed, or     *
 * represents that its use would not infringe privately owned rights.              *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * MPE_Decomp1d - Compute a balanced decomposition of a 1-D array.
 *
 * This file contains a routine for producing a decomposition of a 1-d array
 * when given a number of processors.  It may be used in "direct" product
 * decomposition.  The values returned assume a "global" domain in [1:n]
 *
 *  Input Parameters:
 *       n - Length of the array
 *    size - Number of processors in decomposition
 *    rank - Rank of this processor in the decomposition (0 <= rank < size)
 *
 *  Output Parameters:
 *     s,e - Array indices are s:e, with the original array considered as 1:n.  
 */

#include <stdio.h>
#include "mpi.h"

int MPE_Decomp1d(int n, 
                 int size, 
                 int rank, 
                 int *s, 
                 int *e)
{
  int 
    nlocal,
    deficit;
 
  nlocal  = n/size;
  *s      = rank*nlocal+1;
  deficit = n%size;
  *s      = *s+((rank < deficit) ? rank : deficit);
  if (rank < deficit) {
    nlocal++;
  }
  *e = *s+nlocal-1;
  if (*e > n || rank == size-1) {
    *e = n;
  }

  return MPI_SUCCESS;
}
