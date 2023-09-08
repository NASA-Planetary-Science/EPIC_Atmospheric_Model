/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2018 Timothy E. Dowling                      *
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
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
 * Boston, MA 02110-1301, USA.                                     *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "epic_datatypes.h"

#ifndef EPIC_SUBGRID_H
#define EPIC_SUBGRID_H

/* * * * * * * * * * * * * epic_subgrid.h  * * * * * * * * * * * * * 
 *                                                                 *
 *       TE Dowling, VK Parimi, RP LeBeau,                         *
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
#define DIE(j,i)        die[      i+(j)*Iadim-Shift2d]
#define KIE(j,i)        kie[      i+(j)*Iadim-Shift2d]
#define TAU_WALL(j,i)   tau_wall[ i+(j)*Iadim-Shift2d]

#define LPHH(j,i)  lphh[ i+(j)*Iadim-Shift2d]
#define LPUU(j,i)  lpuu[ i+(j)*Iadim-Shift2d]
#define LPVV(j,i)  lpvv[ i+(j)*Iadim-Shift2d]

#define BUFF1(j,i) buff1[i+(j)*Iadim-Shift2d]
#define BUFF2(j,i) buff2[i+(j)*Iadim-Shift2d]

#define ZE(j,i) ze[i+(j)*Iadim-Shift2d]
#define DI(j,i) di[i+(j)*Iadim-Shift2d]

#define COEF(itmp,j,i) coef[itmp][i+(j)*Iadim-Shift2d]

/*
 * Function prototypes.
 */
EPIC_FLOAT max_nu_nondim(int order);

void set_max_nu(double *max_nu_horizontal);

void set_hyperviscosity(void);

void scalar_horizontal_subgrid(planetspec  *planet,
                               EPIC_FLOAT **Buff2D);

void scalar_horizontal_diffusion(planetspec  *planet,
                                 EPIC_FLOAT **Buff2D);

void scalar_hyperviscosity(int          nu_order,
                           double       nu_hyper,
                           EPIC_FLOAT **Buff2D,
                           int          kstart,
                           int          kend,
                           EPIC_FLOAT  *h);

void adiabatic_adjustment(planetspec  *planet,
                          EPIC_FLOAT **Buff2D);

void laplacian_h(int         kk,
                 EPIC_FLOAT *hh,
                 EPIC_FLOAT *diff_coeff,
                 EPIC_FLOAT *lph,
                 EPIC_FLOAT *buff1,
                 EPIC_FLOAT *buff2);

void scalar_vertical_subgrid(planetspec  *planet,
			     EPIC_FLOAT **Buff2D);

void scalar_vertical_diffusion(planetspec  *planet,
			       EPIC_FLOAT **Buff2D);

void uv_horizontal_subgrid(planetspec  *planet,
                           EPIC_FLOAT **Buff2D);

void uv_horizontal_diffusion(planetspec  *planet,
                             EPIC_FLOAT **Buff2D);

void divergence_damping(planetspec  *planet,
                        int          K,
                        EPIC_FLOAT   nudiv_nondim,
                        EPIC_FLOAT **Buff2D);

void uv_hyperviscosity(int          nu_order,
                       double       nu_hyper,
                       EPIC_FLOAT **Buff2D);

void laplacian_uv(int         K,
                  EPIC_FLOAT *uu,
                  EPIC_FLOAT *vv,
                  EPIC_FLOAT  viscosity,
                  EPIC_FLOAT *lpuu,
                  EPIC_FLOAT *lpvv,
                  EPIC_FLOAT *buff1,
                  EPIC_FLOAT *buff2);

void uv_vertical_subgrid(planetspec  *planet,
                         EPIC_FLOAT **Buff2D);

void uv_vertical_diffusion(planetspec  *planet,
			   EPIC_FLOAT **Buff2D);

void make_arrays_subgrid(void);

void free_arrays_subgrid(void);

void init_subgrid(planetspec *planet);

void set_diffusion_coef(planetspec *planet);

void source_sink_turb(planetspec  *planet,
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

void tau_surface(planetspec  *planet,
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
