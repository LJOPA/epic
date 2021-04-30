/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling                      *
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

#include "epic_datatypes.h"

#ifndef EPIC_SUBGRID_H
#define EPIC_SUBGRID_H

/* * * * * * * * * * * * * epic_subgrid.h  * * * * * * * * * * * * * 
 *                                                                 *
 *       Vimal Kumar Parimi and Raymond P. LeBeau                  *
 *                                                                 *
 *       Header file for epic_subgrid.c                            *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Defines.
 */
#define NU_TURB_EPSILON (1.e-3)
#define NU_TURB_MIN (1.e-3*planet->kinvisc)

/*
 * Shift macros.
 */
#define D_WALL(k,j,i)   d_wall[   i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define DIFF_COEF(j,i)  diff_coef[i+(j)*Iadim-Shift2d]
#define TAU11(j,i)      tau11[    i+(j)*Iadim-Shift2d]
#define TAU22(j,i)      tau22[    i+(j)*Iadim-Shift2d]
#define TAU33(j,i)      tau33[    i+(j)*Iadim-Shift2d]
#define TAU12(j,i)      tau12[    i+(j)*Iadim-Shift2d]
#define DIE(j,i)        die[     i+(j)*Iadim-Shift2d]
#define KIE(j,i)        kie[     i+(j)*Iadim-Shift2d]
#define TAU_WALL(j,i)   tau_wall[i+(j)*Iadim-Shift2d]

/*
 * Function prototypes.
 */
void set_max_nu(double *max_nu);

void init_viscosity(planetspec  *planet);

void scalar_horizontal_diffusion(planetspec  *planet,
                                 EPIC_FLOAT **Buff2D);

void scalar_horizontal_subgrid(planetspec  *planet,
                               EPIC_FLOAT **Buff2D);

void zonal_hyperviscosity(planetspec  *planet,
                          int          is,
                          EPIC_FLOAT  *variable,
                          int          dim,
                          EPIC_FLOAT **Buff2D);

void meridional_hyperviscosity(planetspec  *planet,
                               int          is,
                               EPIC_FLOAT  *variable,
                               int          dim,
                               EPIC_FLOAT **Buff2D);

void laplacian_h(planetspec *planet,
                 EPIC_FLOAT *hh,
                 EPIC_FLOAT *diff_coeff,
                 EPIC_FLOAT *lph,
                 EPIC_FLOAT *buff1,
                 EPIC_FLOAT *buff2);

void d2dy2(planetspec *planet,
           int         index,
           EPIC_FLOAT *hh,
           EPIC_FLOAT  diff_coeff,
           EPIC_FLOAT *d2dy2);

void scalar_vertical_subgrid(planetspec  *planet,
			     EPIC_FLOAT **Buff2D);

void scalar_vertical_diffusion(planetspec  *planet,
			       EPIC_FLOAT **Buff2D);

void uv_horizontal_subgrid(planetspec  *planet,
                           EPIC_FLOAT **Buff2D);

void divergence_damping(planetspec *planet);

void uv_horizontal_diffusion(planetspec  *planet,
                             EPIC_FLOAT **Buff2D);

void uv_vertical_subgrid(planetspec  *planet,
                         EPIC_FLOAT **Buff2D);

void uv_vertical_diffusion(planetspec  *planet,
			   EPIC_FLOAT **Buff2D);

void make_arrays_subgrid(void);

void free_arrays_subgrid(void);

void init_subgrid(planetspec *planet);

void set_diffusion_coef(planetspec *planet);

inline void source_sink_turb(planetspec  *planet,
	                     EPIC_FLOAT **Buff2D);

void source_sink_SA(planetspec  *planet,
		    EPIC_FLOAT **Buff2D);

void dwall_SA(planetspec *planet,
              EPIC_FLOAT *d_wall);

void fp_init_prof(planetspec *planet);
      
EPIC_FLOAT delta_SA(planetspec *planet, 
                    int         K, 
                    int         J, 
                    int         I);

void set_bc_nu_turb(planetspec *planet);

EPIC_FLOAT tau_surface(planetspec  *planet,
                       int          index,
		       EPIC_FLOAT  *tau_wall,
                       EPIC_FLOAT  *buffji); 

EPIC_FLOAT law_of_the_wall(planetspec *planet,
                           int         K,
                           int         J,
                           int         I,
                           int         index,
			   EPIC_FLOAT  t_vis,
			   EPIC_FLOAT  u_tan);

EPIC_FLOAT func_utau(EPIC_FLOAT u_tau,
                     EPIC_FLOAT u_tan,
                     EPIC_FLOAT dwall);

EPIC_FLOAT invert_fv1(planetspec *planet,
                      EPIC_FLOAT t_vis);

/* * * * * * * * * * * * * * * end of epic_subgrid.h * * * * * * * */
#endif
