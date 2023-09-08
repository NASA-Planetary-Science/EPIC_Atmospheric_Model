
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2002-2011 Csaba J. Palotai                        *
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

#ifndef EPIC_MICRO_H
#define EPIC_MICRO_H
/* * * * * * * * * * * * * epic_microphysics.h * * * * * * * * * * * * * * * 
 *                                                                         *
 *                        Csaba J. Palotai  *A*                            *
 *                                                                         *
 *              Header file for the EPIC microphysics model                *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * *
 *                                             *
 *       Structures, Constants and Macros      *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * NOTE: This file does not have the EPIC model definitions like EPIC_FLOAT
 *       available to it (they are defined "downstream"). For cloud microphysics functions
 *       that require these definitions, their prototypes should go in epic.h.
 */

/*
 * CJP 3/3/05
 * introducing boolean type for faster evaluation of frequent fcmp() terms
 */
typedef int boolean;

/*
 * defining constants and macros for the microphysical processes
 */

#define EPSILON machine_epsilon()

typedef struct{
  double
    x,y,z;
} triplet;

/*
 * CJP 10/1/2003
 * Structure for microphysical properties, i.e. threshold values, densities
 * We can specify different or same values for each species.             
 * The values are being assigned in set_species_thermo_data() 
 *
 */
  
/* 
 * Microphysical properties structure, CJP 10/15/2003
 * Modified for multi-planet applications, TD 4/30/09
 *
 * NOTE: The following have been removed from microphysics_spec, because the scheme currently does not use them:
 *   n_c,     Droplet concentration [m-3],  H04 p109; typically 3.e+8
 *   r_cr,    Critical mean droplet radius,  H04 p109; typically 8.e-6
 *   m_i0,    Initial mass of a new ice crystal, D89 p3103; typically 1.e-12
 */

typedef struct {
  double
    rain_density,
    snow_density,
    rain_threshold,            /*  Threshold specific mass value for autoconversion LIQUID-->RAIN [kg/kg] */
    snow_threshold,            /*  Threshold specific mass value for autoconversion SOLID -->SNOW [kg/kg] */
    powerlaw_x_i,              /* In powerlaw V_i(D)  = x*(D)^y  for ice crystals */
    powerlaw_y_i,
    powerlaw_x_s,              /* In powerlaw V_s(D)  = x*(D)^y  for snow */
    powerlaw_y_s,
    powerlaw_x_r,              /* In powerlaw V_r(D)  = x*(D)^y  for rain */
    powerlaw_y_r,
    powerlaw_a,                /* In powerlaw     V   = a*(RHOQ)^b        */ 
    powerlaw_b,
    powerlaw_c,                /* In powerlaw     N_I = c*(RHOQ)^d        */ 
    powerlaw_d,
    powerlaw_m,                /* In powerlaw     M   = alpha*(D)^beta    */
    powerlaw_n,
    n_0r,                      /* Intercept value in raindrop size distribution [m-4], FRR96 p526, D89 p3103 */
    e_r,                       /* Accretion efficiency for rain,  D89 p3104 */
    d_icrit,                   /* Maximum diameter of cloud ice crystal H04 p108 */
    f1r,                       /* Coefficient in the ventillation factor, D89 p3103 */
    f2r,                       /* Coefficient in the ventillation factor, D89 p3103 */
    f1s,                       /* Coefficient in the ventillation factor, D89 p3103 */
    f2s,                       /* Coefficient in the ventillation factor, D89 p3103 */
    p_exp_liq,                 /* Exponent in the (p0/p)^x term for grub:cloud liquid? */
    p_exp_ice,                 /* Exponent in the (p0/p)^x term for grub:cloud ice? */
    p_exp_snow,                /* Exponent in the (p0/p)^x term for snow */
    alpha_raut,                /* Rate coeff. for autoconversion of cloud water [sec-1]. Fowler et al (1996) and Kessler (1969) */
    q_liq_0,                   /* Threshold value for rain autoconversion [kg kg-1]. Fowler et al. (1996) */
    t_00;                      /* Supercooled liquid threshold temperature */
} microphysics_spec;


/* * * * * * * * * * * * * * * * * *
 *                                 *
 * Function prototypes.            *
 *                                 *
 * * * * * * * * * * * * * * * * * */

/*
 * Functions for the microphysical processes
 */
void instantaneous_processes(int   is,
                             int   K,
                             int   J,
                             int   I,
			     int   non_precip,
			     int   no_heating);

void finite_rate_processes(int     is,
                           int     K,
                           int     J,
	                   int     I,
                           double  rh_ji,
			   int     non_precip,
			   int     no_heating);

double terminal_velocity(int    is,
                         int    ip,
                         double pressure,
                         double temperature,
                         double precip_density);

void set_species_thermo_data(void);

void set_microphysics_params(int planet_index);

void restore_mass_min(int  is,
                      int  K);

double dynvisc(  char *globe,
               double  temp);

double conductivity(  char *globe,
                    double  temp);

/* * * * * * * * * * end of epic_microphysics.h * * * * * * * * * * * * * */ 
#endif
