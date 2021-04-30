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


T. Dowling, 1/15/98-------------------------------------------------------
The files in this directory are custom fortran interface functions.
They include mpiimpl.h, which causes the following problems:
 i) mpiimpl.h is not copied to $(MPIR_HOME)/include when mpi is installed.
    Thus, the user has to hand copy the file from mpich/include.
    
ii) Several of the include files called by mpiimpl.h have the same problem
    as i). Others appear to be missing.
    
I moved these *f.c files to this fortran_interface directory, and
commented-out their compilation lines in the makefile.  We give up this
fortran functionality in order to get the code to compile.
---------------------------------------------------------------------------
