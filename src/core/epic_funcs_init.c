/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2023 Timothy E. Dowling                      *
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

/* * * * * * * * * epic_funcs_init.c * * * * * * * * * * * * * * * * * * * * * 
 *                                                                           *
 *       Include here functions used to initialize or change the model       *
 *       that are not called while the model is running.  These functions    *
 *       do not need to be MPI ready.                                        *
 *                                                                           *
 *       This file contains the following functions:                         *
 *                                                                           *
 *           set_re_rp()                                                     *
 *           init_phi_surface()                                              *
 *           init_with_u()                                                   *
 *           init_with_ref()                                                 *
 *           init_fpara_as_fpe()                                             *
 *           init_species()                                                  *
 *           init_vmr_via_deep_value()                                       *
 *           init_vmr_via_obs()                                              *
 *           setup_mu_p(), mu_p()                                            *
 *           t_yp(),fp_yp()                                                  *
 *           p_phi()                                                         *
 *                                                                           *
 *           For processing non-EPIC input:                                  *
 *             openmars, emars, weizmann                                     *
 *               *_make_arrays()                                             *
 *               *_free_arrays()                                             *
 *               *_var_read()                                                *
 *               *_epic_nc()                                                 *
 *               *_conversion()                                              *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>
#include <sys/wait.h>

/*====================== set_re_rp() =========================================*/

void set_re_rp(planetspec *planet)
{

  /*
   * T. Dowling 3/25/2022
   *
   * This version of set_re_rp() fixes re and rp to be constants, not functions of the
   * vertical index, K.
   *
   * We reverted to this because our implementation of variable re and rp
   * appears to introduce biases that contribute to numerical instability.
   * In other words, it was not self-consistent throughout the model.
   */
  register int
    K;

  for (K = 0; K <= KHI; K++) {
    grid.re[K] = planet->re;
    grid.rp[K] = planet->rp;
  }

  return;
}

/*====================== end of set_re_rp() ==================================*/

// /*====================== set_re_rp() =========================================*/
// 
// /*
//  * Determine the equatorial and polar radii as a function of vertical index K.
//  * Assume a uniformly rotating fluid body with no relative winds.  The values
//  * for grid.re[K] and grid.rp[K] are included in the model input/output, hence
//  * this function needs to be called only in epic_initial.c.
//  * 
//  * NOTE: Oblateness of the geopotential increases with altitude.
//  */
// 
// #undef  MAX_ITER
// #define MAX_ITER 10
// #undef  TOL
// #define TOL      (1.e-9)
// 
// /*
//  * Need to smooth grid.re[K], grid.rp[K].
//  * Use a Savitzky-Golay filter of order m, width np.
//  */
// #undef  SAVITZKY_GOLAY_M
// #define SAVITZKY_GOLAY_M   2
// #undef  SAVITZKY_GOLAY_NP
// #define SAVITZKY_GOLAY_NP (2*(SAVITZKY_GOLAY_M+3)+1)
// 
// void set_re_rp(planetspec *planet)
// {
//   register int
//     K,K0,kk,iter,
//     kay,j,i;
//   int
//     nr,nl,
//     nc = SAVITZKY_GOLAY_NP/2;
//   EPIC_FLOAT
//     old_diff,diff,
//     GM,GM_mid,
//     spin_factor,area,a,c,
//     ge,ge_mid,gp,gp_mid,g_avg,m_layer,dp,
//     change,old_ge,prev_gavg,prev_ge,prev_GM,
//     dy;
//   EPIC_FLOAT
//     c_sav[SAVITZKY_GOLAY_NP],
//     smooth_re[grid.nk+1],
//     smooth_rp[grid.nk+1],
//     extended_re[grid.nk+1+2*nc],
//     extended_rp[grid.nk+1+2*nc],
//     xleft[ SAVITZKY_GOLAY_M+1],
//     xright[SAVITZKY_GOLAY_M+1],
//    *padded_re,
//    *padded_rp;
//   /* 
//    * The following are part of DEBUG_MILESTONE(.) statements: 
//    */
//   int
//     idbms=0;
//   static char
//     dbmsname[]="set_re_rp";
// 
//   /*
//    * Shift padded_re, padded_rp.
//    */
//   padded_re = extended_re+nc;
//   padded_rp = extended_rp+nc;
// 
//   /*
//    * Assume grid.p_ref[kk] and grid.rho_ref[kk] refer to a representative sounding at the equator.
//    * Assume planet->re is the equatorial radius for the layer K = K0 where grid.p_ref[2*K+1] is 
//    * closest to 1000 hPa.
//    */
//   old_diff = FLOAT_MAX, K0 = INT_MAX;
//   for (K = 0; K <= KHI; K++) {
//     kk   = 2*K+1;
//     diff = fabs(grid.p_ref[kk]-1000.*100.);
//     if (diff < old_diff) {
//       K0       = K;
//       old_diff = diff;
//     }
//   }
//   grid.re[K0] = a = planet->re;
//   /*
//    * To get started, estimate GM = planet->GM, c = planet->rp;
//    */
//   c            = planet->rp;
//   GM           = planet->GM;
//   spin_factor  = (planet->omega_sidereal*a)*(planet->omega_sidereal*a)*(a/GM);
//   ge           = (1.+1.5*planet->J2-spin_factor)*GM/(a*a);
//   gp           = (1.-3.*(a/c)*(a/c)*planet->J2)*GM/(c*c);
//   area         = surface_area_oblate(a,c);
//   g_avg        = .5*(ge+gp);
//   m_layer      = (grid.p_ref[2*K0+1]-0.)*area/g_avg;
//   /* 
//    * Update estimates of GM, c
//    */
//   GM          -= G_NEWTON*m_layer;
//   grid.rp[K0]  = c = polar_radius(a,planet->J2,GM,planet->omega_sidereal);
// 
//   /*
//    * Need middle values for continuity with downward integration below.
//    */
//   GM_mid = GM;
//   ge_mid = ge;
//   gp_mid = gp;
// 
//   /*
//    * Integrate the hydrostatic balance equation upwards from K0.
//    */
//   for (K = K0-1; K >= 0; K--) {
//     dp        = (grid.p_ref[2*K+1]-grid.p_ref[2*(K+1)+1]);
//     prev_ge   = ge;
//     prev_gavg = .5*(ge+gp);
//     prev_GM   = GM;
//     for (iter = 0; iter < MAX_ITER; iter++) {
//       /*
//        * Iterate to improve self-consistency
//        */
//       old_ge     = ge;
//       grid.re[K] = a = grid.re[K+1]-dp/(.5*(ge+prev_ge)*grid.rho_ref[2*(K+1)]);
//       grid.rp[K] = c = polar_radius(a,planet->J2,GM,planet->omega_sidereal);
//       /*
//        * Update estimate of GM, ge, gp.
//        */
//       area         = .5*(surface_area_oblate(a,c)
//                         +surface_area_oblate(grid.re[K+1],grid.rp[K+1]));
//       g_avg        = .5*(.5*(ge+gp)+prev_gavg);
//       m_layer      = -dp*area/g_avg;
//       GM           = prev_GM+G_NEWTON*m_layer;
//       spin_factor  = SQR(planet->omega_sidereal*a)*(a/GM);
//       ge           = (1.+1.5*planet->J2-spin_factor)*GM/(a*a);
//       gp           = (1.-3.*(a/c)*(a/c)*planet->J2)*GM/(c*c);
//       change       = (ge-old_ge)/old_ge;
//       if (fabs(change) <= TOL) {
//         break;
//       }
//     }
//   }
// 
//   /*
//    * Integrate the hydrostatic balance equation downwards from K0.
//    */
//   GM = GM_mid;
//   ge = ge_mid;
//   gp = gp_mid;
//   for (K = K0+1; K <= KHI; K++) {
//     dp        = (grid.p_ref[2*(K-1)+1]-grid.p_ref[2*K+1]);
//     prev_ge   = ge;
//     prev_gavg = .5*(ge+gp);
//     prev_GM   = GM;
//     for (iter = 0; iter < MAX_ITER; iter++) {
//       /*
//        * Iterate to improve self-consistency
//        */
//       old_ge     = ge;
//       grid.re[K] = a = grid.re[K-1]+dp/(.5*(ge+prev_ge)*grid.rho_ref[2*K]);
//       grid.rp[K] = c = polar_radius(a,planet->J2,GM,planet->omega_sidereal);
//       /*
//        * Update estimate of GM, ge, gp.
//        */
//       area         = .5*(surface_area_oblate(a,c)
//                         +surface_area_oblate(grid.re[K-1],grid.rp[K-1]));
//       g_avg        = .5*(.5*(ge+gp)+prev_gavg);
//       m_layer      = -dp*area/g_avg;
//       GM           = prev_GM-G_NEWTON*m_layer;
//       spin_factor  = SQR(planet->omega_sidereal*a)*(a/GM);
//       ge           = (1.+1.5*planet->J2-spin_factor)*GM/(a*a);
//       gp           = (1.-3.*(a/c)*(a/c)*planet->J2)*GM/(c*c);
//       change       = (ge-old_ge)/old_ge;
//       if (fabs(change) <= TOL) {
//         break;
//       }
//     }
//   }
// 
//   if (grid.nk >= SAVITZKY_GOLAY_NP-1) {
//     /*
//      * Smooth grid.re[K] and grid.rp[K].
//      * See Savitzky-Golay parameters defined at top of function.
//      *
//      * NOTE: Near the endpoints, we tried off-centering the filter in a continuous manner,
//      *       but this generated glitches in the second derivative of re[K] and rp[K].
//      *       Instead, we pad the endpoints with a polynomial fit of order SAVITZKY_GOLAY_M,
//      *       which allows us to run a centered filter across the array from K = 0 to grid.nk,
//      *       resulting in a smooth second-derivative for grid.re[K] and grid.rp[K].
//      *
//      * Extend ends of re[K], rp[K] with a polynomial fit of order SAVITZKY_GOLAY_M,
//      * fitted to the respective end points.
//      */
//     for (kay = 0; kay <= SAVITZKY_GOLAY_M; kay++) {
//       xleft[ kay] = (EPIC_FLOAT)kay;
//       xright[kay] = (EPIC_FLOAT)(kay+grid.nk-SAVITZKY_GOLAY_M);
//     }
//     for (kay = -nc; kay < 0; kay++) {
//       padded_re[kay]         = poly_interp(SAVITZKY_GOLAY_M+1,xleft, grid.re,                         (EPIC_FLOAT)kay,          &dy);
//       padded_rp[kay]         = poly_interp(SAVITZKY_GOLAY_M+1,xleft, grid.rp,                         (EPIC_FLOAT)kay,          &dy);
//       padded_re[grid.nk-kay] = poly_interp(SAVITZKY_GOLAY_M+1,xright,grid.re+grid.nk-SAVITZKY_GOLAY_M,(EPIC_FLOAT)(grid.nk-kay),&dy);
//       padded_rp[grid.nk-kay] = poly_interp(SAVITZKY_GOLAY_M+1,xright,grid.rp+grid.nk-SAVITZKY_GOLAY_M,(EPIC_FLOAT)(grid.nk-kay),&dy);
//     }
//     for (K = 0; K <= grid.nk; K++) {
//       padded_re[K] = grid.re[K];
//       padded_rp[K] = grid.rp[K];
//     }
// 
//     /* Zero smooth arrays. */
//     memset(smooth_re,0,(grid.nk+1)*sizeof(EPIC_FLOAT));
//     memset(smooth_rp,0,(grid.nk+1)*sizeof(EPIC_FLOAT));
// 
//     /* Centered filter */
//     nl = nr = nc;
//     savitzky_golay(c_sav,SAVITZKY_GOLAY_NP,nl,nr,0,SAVITZKY_GOLAY_M);
//     for (K = 0; K <= KHI; K++) {
//       for (kay = K-nl,i = 0; kay <= K+nr; kay++,i++) {
//         /* c_sav[] uses a wrap-around index */
//         j             = (SAVITZKY_GOLAY_NP+nl-i)%SAVITZKY_GOLAY_NP;
//         smooth_re[K] += padded_re[kay]*c_sav[j];
//         smooth_rp[K] += padded_rp[kay]*c_sav[j];
//       }
//     }
// 
//     for (K = 0; K <= KHI; K++) {
//       grid.re[K] = smooth_re[K];
//       grid.rp[K] = smooth_rp[K];
//     }
//   }
// 
//   return;
// }
// 
// /*====================== end of set_re_rp() ==================================*/

/*====================== init_phi_surface() ==================================*/

/*
 * Set surface geopotential, phi, for terrestrial planets, 
 * using spherical harmonic gravity and surface data.
 * The level phi = 0 corresponds to an appropriate reference radius.
 */

void init_phi_surface(planetspec       *planet,
                      init_defaultspec *def)
{
  int
    K,J,I,
    jlocal,ilocal,
    l,m,
    max_l = -1, 
    i;
  int
    num_files,nc_err,nc_id,
    nc_grid_ni,
    nc_grid_nj,
    file_match,
    node,num_nodes,
    jdim,idim,
    start[TWODIM],
    end[TWODIM];
  size_t
    index[1];
  EPIC_FLOAT
    nc_grid_globe_lonbot,
    nc_grid_globe_lontop,
    nc_grid_globe_latbot,
    nc_grid_globe_lattop;
  char
    nc_grid_geometry[GEOM_STR],
    nc_planet_name[N_STR],  /* NOTE: If [16], on Darwin next variable sometimes overwritten */
    phi_surface_nc[N_STR];
  EPIC_FLOAT
    lonr,latr,sin_lat,cos_lat,
    sum,r,phi,
    rerp,
    r0 = 0.,
    two_l,lplusm,lminusm,
    ge,gp,
    spin_factor,sinlat2,coslat2;
  EPIC_FLOAT
   ***p_lm,
    **cos_mlon,
    **sin_mlon,
    **c_lm_r,     /* r for radius */
    **s_lm_r;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *gravity;            /* Some planets need an estimate of surface gravity vs latitude */
  char
    data_file_name[N_STR],
    header[N_STR];
  nc_type
    nc_float_type;
  FILE
    *infile_r;
  struct dirent
    **namelist;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_phi_surface";

#if EPIC_PRECISION == DOUBLE_PRECISION
   nc_float_type  = NC_DOUBLE;
#else
   nc_float_type  = NC_FLOAT;
#endif

  if (!initialized) {
    if (strcmp(planet->type,"terrestrial") != 0) {
      sprintf(Message,"planet->type=%s not implemented",planet->type);
      epic_error(dbmsname,Message);
    }

    /* Allocate memory */
    gravity = fvector(0,JHI-JLO,dbmsname);

    /*
     * Set up surface gravity vs latitude, for cases that need it.
     *
     * NOTE: grid.g[kk][jj] not available.
     */
    spin_factor = (planet->omega_sidereal*planet->re)*(planet->omega_sidereal*planet->re)*(planet->re/planet->GM);
    ge          = (1.+1.5*planet->J2-spin_factor)*planet->GM/(planet->re*planet->re);
    gp          = (1.-3.*(planet->re/planet->rp)*(planet->re/planet->rp)*planet->J2)*planet->GM/(planet->rp*planet->rp);
    for (J = JLO; J <= JHI; J++) {
      sinlat2        = sin(grid.lat[2*J+1]*DEG);
      sinlat2       *= sinlat2;
      coslat2        = cos(grid.lat[2*J+1]*DEG);
      coslat2       *= coslat2;
      gravity[J-JLO] = (planet->re*ge*coslat2+planet->rp*gp*sinlat2)/sqrt(planet->re*planet->re*coslat2+planet->rp*planet->rp*sinlat2);
    }

    initialized = TRUE;
  }

  /*
   * Handle special cases first.
   */
  if (strcmp(planet->name,"Held_Suarez") == 0) {
    /*
     * Held-Suarez test case has no topography.
     * Set PHI_SURFACE(J,I) = 0. and return.
     */
    memset(var.phi_surface.value,0,sizeof(EPIC_FLOAT)*Nelem2d);
    return;
  }
  else if (strcmp(planet->name,"Goldmine") == 0) {
    /*
     * Goldmine test case has a user defined idealized topography.
     * Insert a Gaussian hill centered at 0-degrees longitude and
     * centered in latitude on the channel chosen via initial program.
     *
     * NOTE: grid.rln[kk][jj], grid.rlt[kk][jj] not available.
     */
    EPIC_FLOAT
      x,y,r2,a,height;

    grid.aux_fa = height = 1000.*input_float("Height of mountain [km]\n",8.);
    if (height == 0.){
      memset(var.phi_surface.value,0,sizeof(EPIC_FLOAT)*Nelem2d);
      return;
    }
    else {
      grid.aux_fb = a = 1000.*input_float("Gaussian width [km]\n",80.);
      for (J = JLO; J <= JHI; J++){
        y = (grid.lat[2*J+1]-grid.lat[2*(grid.nj/2)-1])*DEG*planet->re;
        for (I = ILO; I <= IHI; I++){
          x = (grid.lon[2*I+1]-0.)*DEG*planet->re*cos(grid.lat[2*J+1]*DEG);
          r2 = x*x + y*y;
          PHI_SURFACE(J,I)=gravity[J-JLO]*height*exp(-r2/(a*a));
        }
      }
    }
    /* Need to apply bc_lateral() here */
    bc_lateral(var.phi_surface.value,TWODIM);
    
    return;
  }

  /* 
   * The flag file_match indicates whether an appropriate phi_surface.nc file already exists. 
   */
  file_match = FALSE;

  if (IAMNODE == NODE0) {
    /*
     * Search to see if appropriate phi_surface.nc data file
     * already exists.
     */

    /*
     * On sun4, rs6000, and sp2 platforms, we find that calling the function alphasort 
     * in the last argument of scandir() generates errors, even though it is
     * described in their man pages.  Since we don't need the directory
     * entries to be alphabetized, we fall back to using the NULL argument.
     */
    num_files = scandir(".",&namelist,is_phi_surface_file,NULL);

    if (num_files < 0) {
      /*
       * Problem encountered reading the files.
       */
      perror("init_phi_surface():scandir");
      exit(1);
    }
    else if (num_files > 0) {
      /*
       * There are a total of num_files phi_surface netCDF data files.
       * Search them to see if one matches what is needed.
       */
      for (l = 0; l < num_files; l++) {
        /*
         * Open file.
         */
        sprintf(phi_surface_nc,"./%s",namelist[l]->d_name);
        nc_err = nc_open(phi_surface_nc,NC_NOWRITE,&nc_id);
        if (nc_err != NC_NOERR) {
          sprintf(Message,"%s, %s",nc_strerror(nc_err),phi_surface_nc);
          epic_error(dbmsname,Message);
        }

        /*
         * Check to see if this file matches what is needed.
         */

        /* name */
        nc_err = nc_get_att_text(nc_id,NC_GLOBAL,"planet_name",nc_planet_name);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (strcmp(nc_planet_name,planet->name) != 0) {
          nc_close(nc_id); 
          continue;
        }

        /* geometry */
        nc_err = nc_get_att_text(nc_id,NC_GLOBAL,"grid_geometry",nc_grid_geometry);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (strcmp(nc_grid_geometry,grid.geometry) != 0) {
          nc_close(nc_id); 
          continue;
        }

        /* ni */
        nc_err = nc_get_att_int(nc_id,NC_GLOBAL,"grid_ni",&nc_grid_ni);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_ni != grid.ni) {
          nc_close(nc_id); 
          continue;
        }

        /* nj */
        nc_err = nc_get_att_int(nc_id,NC_GLOBAL,"grid_nj",&nc_grid_nj);
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_nj != grid.nj) {
          nc_close(nc_id); 
          continue;
        }

        /* latbot */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_latbot",&nc_grid_globe_latbot);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_latbot",&nc_grid_globe_latbot);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_latbot != grid.globe_latbot) {
          nc_close(nc_id); 
          continue;
        }

        /* lattop */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lattop",&nc_grid_globe_lattop);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lattop",&nc_grid_globe_lattop);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lattop != grid.globe_lattop) {
          nc_close(nc_id); 
          continue;
        }

        /* lonbot */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lonbot",&nc_grid_globe_lonbot);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lonbot",&nc_grid_globe_lonbot);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lonbot != grid.globe_lonbot) {
          nc_close(nc_id); 
          continue;
        }

        /* lontop */
#if EPIC_PRECISION == DOUBLE_PRECISION
        nc_err = nc_get_att_double(nc_id,NC_GLOBAL,"grid_globe_lontop",&nc_grid_globe_lontop);
#else
        nc_err = nc_get_att_float(nc_id,NC_GLOBAL,"grid_globe_lontop",&nc_grid_globe_lontop);
#endif
        if (nc_err != NC_NOERR) {
          fprintf(stderr,"Warning: init_phi_surface(): %s\n", nc_strerror(nc_err));
          continue;
        }
        if (nc_grid_globe_lontop != grid.globe_lontop) {
          nc_close(nc_id); 
          continue;
        }

        /* 
         * Everything matches. 
         */
        file_match = TRUE;
        nc_close(nc_id);
        break;
      }
    }
  }
  
  if (file_match == TRUE) {
    /*
     * Input PHI_SURFACE(J,I) from file.
     */
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\nReading PHI_SURFACE(J,I) from %s \n\n",
                     phi_surface_nc);
      /* Open file. */
      nc_err = nc_open(phi_surface_nc,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        fprintf(stderr,"Cannot find input file %s \n",phi_surface_nc);
        exit(1);
      }
    }
    num_nodes = setup_read_array();
    /*
     * Load start and end vectors.
     */
    start[0]  = 1;
    end[  0]  = grid.ni;

    for (node = 0; node < num_nodes; node++) {
      get_jlohi(node,num_nodes,start+1,end+1);
      read_array(node,TWODIM,start,end,var.phi_surface.info[0].name,
                 var.phi_surface.info[0].index,var.phi_surface.value,EPIC_FLOAT_ARRAY,nc_id);
    }
    if (IAMNODE == NODE0) {
      /* Close file. */
      nc_close(nc_id);
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);

    return;
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                                       *
   * The required PHI_SURFACE does not exist, so we now construct it.      *
   *                                                                       *
   * We prefer to construct the topography from spherical harmonic data,   *
   * rather than gridded data, because the former provides a good          *
   * way to interpolate smoothly onto our user-defined model grid.         *
   *                                                                       *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Spherical harmonic data are truncated to maximum degree max_l for radius data.
   * Read in max_l.
   */
  if (IAMNODE == NODE0) {
    if (strncmp(planet->name,"Venus",5) == 0) {
      /*-------* 
       * Venus *
       *-------*/

      /*
       * Venus topography data:
       *   Rappaport, N.J., A.S. Konopliv, A.B. Kucinskas, P.G. Ford, 1999, 
       *     An improved 360 degree and order model of Venus topography, Icarus 139, 19-31.
       */
      strcpy(data_file_name,EPIC_PATH"/data/Venus/TOPOCOEF-gtdr3.2");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) data from %s \n",data_file_name);

      fscanf(infile_r,"%d",&max_l);
    }
    else if (strncmp(planet->name,"Earth",5) == 0) {
      /*-------* 
       * Earth *
       *-------*/

      /*
       * Earth topography data:
       *   Wieczorek MA, 2007, Gravity and topography of the terrestrial planets, 
       *     Treatise on Geophysics 10, 165-205.
       */
      strcpy(data_file_name,EPIC_PATH"/data/Earth/srtmp720.msl");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) data from %s \n",data_file_name);

      /* Skip header */
      for (i = 0; i < 7; i++) {
        fgets(header,N_STR,infile_r);
      }
      /* Input max_l */
      fscanf(infile_r,"%d",&max_l);
      /* Skip column heading */
      for (i = 0; i < 2; i++) {
        fgets(header,N_STR,infile_r);
      }
    }
    else if (strncmp(planet->name,"Mars",4) == 0) {
      /*------* 
       * Mars *
       *------*/

      /*
       * Radius data from the Mars Orbiter Laser Altimeter (MOLA) on the Mars Global Surveyor (MGS) spacecraft.
       */  
      strcpy(data_file_name,EPIC_PATH"/data/Mars/gtm090aa_sha.txt");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) radius data from %s \n",data_file_name);

      fscanf(infile_r," %*f, %*f, %*f, %d, %*d, %*d, %*f, %*f",&max_l);

      /* Convert reference radius from km to m. */
      r0 *= 1.e+3;  
    }
    else if (strncmp(planet->name,"Titan",5) == 0) {
      /*-------* 
       * Titan *
       *-------*/

      /*
       * Titan topography data:
       *   Corlies P, Hayes AG, Birch SPD, Lorenz R, Stiles BW, Kirk R, Poggiali V, Zebker H, Iess L, 2017,
       *     Titan's topography and shape at the end of the Cassini mission, Geophys. Res. Letts., 44, 11,754--11,761,
       *     doi:10.1002/2017GL075518. 
       */
      strcpy(data_file_name,EPIC_PATH"/data/Titan/topo_harmonic_coeffs.txt");
      infile_r = fopen(data_file_name,"r");
      if (!infile_r) {
        sprintf(Message,"could not open %s",data_file_name);
        epic_error(dbmsname,Message);
      }
      fprintf(stderr,"Reading spherical harmonic (l,m) radius data from %s \n",data_file_name);

      /* Skip header */
      for (i = 0; i < 9; i++) {
        fgets(header,N_STR,infile_r);
      }
      /* Input max_l */
      fscanf(infile_r,"%d",&max_l);
      /* Skip column heading */
      for (i = 0; i < 3; i++) {
        fgets(header,N_STR,infile_r);
      }
    }
    else {
      sprintf(Message,"failed to open spherical-harmonic file(s) for %s",planet->name);
      epic_error(dbmsname,Message);
    }
  }
  
  /*
   * Allocate memory for (l,m) spherical-harmonic coefficients, p_lm.
   */
  p_lm = (EPIC_FLOAT ***)calloc(JHI-JLO+1,sizeof(EPIC_FLOAT **));
  if (!p_lm) {
    sprintf(Message,"calloc error allocating p_lm");
    epic_error(dbmsname,Message);
  }
  for (J = JLO; J <= JHI; J++) {
    p_lm[J-JLO] = (EPIC_FLOAT **)calloc(max_l+1,sizeof(EPIC_FLOAT *));
    if (!p_lm[J-JLO]) {
      sprintf(Message,"calloc error allocating p_lm[%d-%d]",J,JLO);
      epic_error(dbmsname,Message);
    }
    for (l = 0; l <= max_l; l++) {
      p_lm[J-JLO][l] = fvector(0,l,dbmsname);
    }
  }

  cos_mlon = (EPIC_FLOAT **)calloc(grid.ni,sizeof(EPIC_FLOAT *));
  if (!cos_mlon) {
    sprintf(Message,"calloc error allocating cos_mlon");
    epic_error(dbmsname,Message);
  }
  sin_mlon = (EPIC_FLOAT **)calloc(grid.ni,sizeof(EPIC_FLOAT *));
  if (!sin_mlon) {
    sprintf(Message,"calloc error allocating sin_mlon");
    epic_error(dbmsname,Message);
  }
  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    cos_mlon[ilocal] = fvector(0,max_l,dbmsname);
    sin_mlon[ilocal] = fvector(0,max_l,dbmsname);
  }

  if (max_l >= 0) {
    c_lm_r = (EPIC_FLOAT **)calloc(max_l+1,sizeof(EPIC_FLOAT *));
    if (!c_lm_r) {
      sprintf(Message,"calloc error allocating c_lm_r");
      epic_error(dbmsname,Message);
    }

    s_lm_r = (EPIC_FLOAT **)calloc(max_l+1,sizeof(EPIC_FLOAT *));
    if (!s_lm_r) {
      sprintf(Message,"calloc error allocating s_lm_r");
      epic_error(dbmsname,Message);
    }

    for (l = 0; l <= max_l; l++) {
      c_lm_r[l] = fvector(0,l,dbmsname);
      s_lm_r[l] = fvector(0,l,dbmsname);
    }
  }

  /*
   * Input spherical-harmonic coefficients.
   */
  if (IAMNODE == NODE0) {
    if (strncmp(planet->name,"Venus",5) == 0) {
      /*-------* 
       * Venus *
       *-------*/

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile_r," %lf",&r0);
#else
      fscanf(infile_r," %f",&r0);
#endif

      /* Convert reference radius from km to m. */
      r0 *= 1.e+3;
      for (l = 0; l <= max_l; l++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %lf",c_lm_r[l]);
#else
          fscanf(infile_r," %f",c_lm_r[l]);
#endif

        for (m = 1; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %lf",c_lm_r[l]+m);
          fscanf(infile_r," %lf",s_lm_r[l]+m);
#else
          fscanf(infile_r," %f",c_lm_r[l]+m);
          fscanf(infile_r," %f",s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);
    }
    else if (strncmp(planet->name,"Earth",5) == 0) {
      /*-------* 
       * Earth *
       *-------*/

      /* Data gives elevation rather than radius; signal this with r0 = 0. */
      r0 = 0.;
      for (l = 0; l <= max_l; l++) {
        for (m = 0; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r,"%*d %*d %lf %lf",c_lm_r[l]+m,s_lm_r[l]+m);
#else
          fscanf(infile_r,"%*d %*d %f %f",c_lm_r[l]+m,s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);
    }
    else if (strncmp(planet->name,"Mars",4) == 0) {
      /*------* 
       * Mars *
       *------*/

      for (l = 0; l <= max_l; l++) {
        for (m = 0; m <= l; m++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %*d, %*d, %lf, %lf, %*f, %*f",c_lm_r[l]+m,s_lm_r[l]+m);
#else
          fscanf(infile_r," %*d, %*d, %f, %f, %*f, %*f",c_lm_r[l]+m,s_lm_r[l]+m);
#endif

        }
      }
      fclose(infile_r);
    }
    else if (strncmp(planet->name,"Titan",5) == 0) {
      /*-------* 
       * Titan *
       *-------*/

      EPIC_FLOAT
        tmp1,tmp2,
        sign;

      /*
       * See the supporting material for Corlies et al (2017), p.4, for a discussion
       * of the spherical harmonic coefficients.  The convention for m < 0 is C_lm = S_l|m|.
       */
      for (l = 0; l <= max_l; l++) {
        m    = 0;
        sign = 1.;

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile_r," %*d %*d %lf %*lf",c_lm_r[l]+m);
#else
        fscanf(infile_r," %*d %*d %f %*f",  c_lm_r[l]+m);
#endif

        for (m = 1; m <= l; m++) {
          /*
           * The Corlies et al.(2017) coefficients need the (-1)^m, Condon-Shortley phase factor.
           */
          sign = -sign;


#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(infile_r," %*d %*d %lf %*lf %*d %*d %lf %*lf",&tmp1,&tmp2);
#else
          fscanf(infile_r," %*d %*d %f %*f %*d %*d %f %*f",    &tmp1,&tmp2);
#endif

          c_lm_r[l][m] = sign*tmp1;
          s_lm_r[l][m] = sign*tmp2;
        }
      }
    }
    else {
      sprintf(Message,"failed to read spherical-harmonic topography-coefficient file for %s",planet->name);
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Calculate PHI_SURFACE(J,I).
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"Calculating topography: \n");
  }

  /*
   * Precalculate cos, sin functions.
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Calculating and storing cos and sin factors.\n");
  }
  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    lonr   = grid.lon[2*I+1]*DEG;
    for (m = 0; m <= max_l; m++) {
      cos_mlon[ilocal][m] = cos((EPIC_FLOAT)m*lonr);
      sin_mlon[ilocal][m] = sin((EPIC_FLOAT)m*lonr);
    }
  }

  /*
   * Precalculate associated Legendre functions (ALFs).
   *
   * For the normalized ALFs, use a recursive formula by 
   *   Colombo C, 1981, Numerical methods for harmonic analysis on the sphere. Rep 310, 
   *     Dept. of Geodetic Science and Surveying, The Ohio State University, Columbus
   * which is cited and reproduced as eqn.(11) in Holmes and Featherstone (2002), 
   * and also described on the MITGCM website: http://mitgcm.org/~mlosch/geoidcookbook/node11.html.
   *
   */
  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Calculating and storing normalized Legendre factors:   0%%");
  }

  rerp = planet->re/planet->rp;

  for (J = JLO; J <= JHI; J++) {
    jlocal  = J-JLO;
    latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
    sin_lat = sin(latr);
    /*
     * Start with sectorial (l == m) values. These serve as seeds for the non-sectorial (l != m) recursion.
     * We divide out the cos_lat factor to avoid underflow for high degree and order, following 
     *   Holmes SA, Featherstone WE, 2002,  A unified approach to the Clenshaw summation and the recursive
     *     computation of very high degree and order normalised associated Legendre functions, 
     *     J. Geodesy 76, 279-299   
     */
    p_lm[jlocal][0][0] = 1.;
    p_lm[jlocal][1][1] = sqrt(3.);
    for (m = 2; m <= max_l; m++) {
      /*
       * The cos_lat factor is deliberately left off, and brought in below, following Holmes and Featherstone (2002).
       */
      p_lm[jlocal][m][m] = sqrt(1.+1./(2.*(double)m))*p_lm[jlocal][m-1][m-1];
    }
    /*
     * Calculate non-sectorial values.
     */
    for (l = 1; l <= max_l; l++) {
      two_l = 2.*(double)l;
      m     = l-1;
      p_lm[jlocal][l][m] = sqrt(two_l+1.)*p_lm[jlocal][l-1][m]*sin_lat;
      for (m = 0; m < l-1; m++) {
        lplusm  = (double)(l+m);
        lminusm = (double)(l-m);
        p_lm[jlocal][l][m] = p_lm[jlocal][l-1][m]*sqrt((two_l-1.)*(two_l+1.)/(lminusm*lplusm))*sin_lat
                            -p_lm[jlocal][l-2][m]*sqrt((two_l+1.)*(lplusm-1.)*(lminusm-1.)/(lminusm*lplusm*(two_l-3.)));
      }
    }
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
    }

    if (strncmp(planet->name,"Titan",5) == 0) {
      EPIC_FLOAT
        coef,norm;

      /*
       * The Titan community uses unnormalized Legendre functions, so unnormalize p_lm.
       * See Holmes and Featherstone (2002), Eq.(8).
       *
       * Algorithm by M.E. Bradley, 5 Sept 2018.
       */
      for (l = 1; l <= max_l; l++) {
        m = 0;
        norm                = sqrt(2.*l+1.);
        p_lm[jlocal][l][m] /= norm;

        coef = 1.;
        for (m = 1; m <= l; m++) {
          coef               /= (l-m+1)*(l+m);
          norm                = sqrt((4.*l+2.)*coef);
          p_lm[jlocal][l][m] /= norm;
        }
      }
    }
  }
  if (IAMNODE == NODE0) {
    fprintf(stderr,"\n");
  }

  if (IAMNODE == NODE0) {
    fprintf(stderr,"  Summing to get PHI_SURFACE(J,I):   0%%");
  }
  if (strncmp(planet->name,"Venus",5) == 0) {
    /*-------* 
     * Venus *
     *-------*/

    /*
     * Coefficients yield radius divided by average radius, which we convert to elevation, z, and then to phi.
     */
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        sum = 0.;
        for (m = max_l; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l; l++) {
            sum += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                      +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }
        PHI_SURFACE(J,I) = gravity[J-JLO]*r0*(sum-1.);
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }
  else if (strncmp(planet->name,"Earth",5) == 0) {
    /*-------* 
     * Earth *
     *-------*/

    /*
     * Coefficients yield elevation, z, which we convert to phi.
     */
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        sum = 0.;
        for (m = max_l; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l; l++) {
            sum += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                      +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }
        /*
         * NOTE: Earth topography data set currently includes ocean-floor topography,
         *       which will need to be set to sea level for atmospheric applications.
         */
        PHI_SURFACE(J,I) = gravity[J-JLO]*sum;
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }
  else if (strncmp(planet->name,"Mars",4) == 0) {
    /*------* 
     * Mars *
     *------*/

    /*
     * Using g0*(r-r0) rather than the actual geopotential, to match the OpenMARS model.  Additionally.
     * the Mars modeling community tends to assume Mars is a sphere, so we remove the rotational flattening.
     *
     * NOTE: We tried to use the spherical harmonic model of the geopotential itself, but because
     *       the OpenMARS model did not do this, we were not able to find a self-consistent
     *       way of using it.
     */
    EPIC_FLOAT
      g0         = 3.71,
      reduce_J2  = 0.06,
      lowest_lat = 1.e+20,
      avg;
    int
      J_lowest_lat;

     /* 
      * Remove rotational flattening by reducing J2 = C_20 by 94%,
      * following Wieczorek & Zuber (2004, JGR 109, doi:10.1029/2003JE002153).
      */
    c_lm_r[2][0] *= reduce_J2;


    J_lowest_lat = (JHI-JLO)/2;
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);

      if (fabs(latr) < lowest_lat) {
        lowest_lat   = fabs(latr);
        J_lowest_lat = J;
      }

      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        r = 0.;
        for (m = max_l; m >= 0; m--) {
          r *= cos_lat;
          for (l = m; l <= max_l; l++) {
            r += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                    +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }

        PHI_SURFACE(J,I) = g0*r;
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);

    /*
     * Subtract equatorial average.
     *
     * NOTE: NOT MPI ready.
     */
    avg = 0.;
    J   = J_lowest_lat;
    for (I = ILO; I <= IHI; I++) {
      avg += PHI_SURFACE(J,I);
    }
    avg /= grid.ni;

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        PHI_SURFACE(J,I) -= avg;
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);
    
    /*
     * Restore J2.
     */
    c_lm_r[2][0] /= reduce_J2;

    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }
  else if (strncmp(planet->name,"Titan",5) == 0) {
    /*-------* 
     * Titan *
     *-------*/

    /*
     * Titan reference ellipsoid model from Eq. (12) of Mitri et al. (2014, Icarus 236, 169--177).
     * The values for C_20 = - J_2 and C_22 are from Table 2, SOL1a of Iess et al. (2012, Science 337, 457--459).
     * 
     * NOTE: Corlies et al. and Mitri et al. use unnormalized Legengre polynomials, hence p_lm
     *       has been unnormalized above.
     *
     * NOTE: The second instance of "q_r" in Mitri et al. (2014), Eq. (12) should be "q_t".
     */
    EPIC_FLOAT
      q_r  =   3.9528e-5,
      q_t  =  -1.1858e-4,
      C_20 = -33.5990e-6,
      C_22 =  10.1210e-6,
      reference_ellipsoid;

    /*
     * Coefficients yield radius, which we convert to elevation above the geoid, z, and then to phi.
     */
    for (J = JLO; J <= JHI; J++) {
      jlocal = J-JLO;
      latr    = lat_graphic_to_centric(grid.lat[2*J+1],rerp)*DEG;
      cos_lat = cos(latr);
      for (I = ILO; I <= IHI; I++) {
        ilocal = I-ILO;
        /*
         * The cos_lat factor is folded back in here, following Holmes and Featherstone (2002).
         */
        sum = 0.;
        for (m = max_l; m >= 0; m--) {
          sum *= cos_lat;
          for (l = m; l <= max_l; l++) {
            sum += p_lm[jlocal][l][m]*(c_lm_r[l][m]*cos_mlon[ilocal][m]
                                      +s_lm_r[l][m]*sin_mlon[ilocal][m]);
          }
        }
        reference_ellipsoid = c_lm_r[0][0]*( 1.+(C_20-5.*q_r/6.)*p_lm[jlocal][2][0]
                                               +(C_22+   q_t/4.)*p_lm[jlocal][2][2]*cos(2.*latr) );
        PHI_SURFACE(J,I)    = gravity[J-JLO]*(sum-reference_ellipsoid);
      }
      if (IAMNODE == NODE0) {
        fprintf(stderr,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)(jlocal+1)/(JHI-JLO+1)));
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(var.phi_surface.value,TWODIM);
    if (IAMNODE == NODE0) {
      fprintf(stderr,"\n\n");
    }
  }

  /*
   * Free allocated memory.
   */
  if (max_l >= 0) {
    for (l = max_l; l >= 0; l--) {
      free_fvector(s_lm_r[l],0,l,dbmsname);
      free_fvector(c_lm_r[l],0,l,dbmsname);
    }
    free(s_lm_r);
    free(c_lm_r);
  }

  for (I = ILO; I <= IHI; I++) {
    ilocal = I-ILO;
    free_fvector(sin_mlon[ilocal],0,max_l,dbmsname);
    free_fvector(cos_mlon[ilocal],0,max_l,dbmsname);
  }
  free(sin_mlon);
  free(cos_mlon);

  for (J = JLO; J <= JHI; J++) {
    jlocal = J-JLO;
    for (l = 0; l <= max_l; l++) {
      free_fvector(   p_lm[jlocal][l],0,max_l,dbmsname);
    }
    free(   p_lm[jlocal]);
  }
  free(   p_lm);

  /*
   * Apply high-latitude, low-pass zonal filter to prevent numerical instabilities.
   *
   * NOTE: This must be done right after the topography is assembled, because so much
   *       depends on the mountains. To do this later literally causes an earthquake
   *       in the model.
   *
   *       One practical implication is that zonal_filter() needs set_fmn(),
   *       which needs set_re_rp(), but the latter involves a hydrostatic integration that uses
   *       a reference density profile that not set at this point. This is gas-giant functionality
   *       interfering with terrestrial-planet functionality.  We work around this by
   *       temporarilty assigning grid.re[K] = planet->re, grid.rp[K] = planet->rp here.
   */
  for (K = KLOPAD; K <= KHIPAD; K++) {
    grid.re[K] = planet->re;
    grid.rp[K] = planet->rp;
  }
  set_fmn(planet);
  /*
   * Call with negative index to cause function re-initialization on next call.
   */
  zonal_filter(-PHI_SURFACE_INDEX,var.phi_surface.value);

  /*
   * Write PHI_SURFACE(J,I) data to a phi_surface.nc file.
   */

  if (IAMNODE == NODE0) {
    /*
     * Enter define mode:
     */
    sprintf(phi_surface_nc,"./%s_%02d_phi_surface.nc",planet->name,num_files);
    nc_err = nc_create(phi_surface_nc,NC_CLOBBER,&nc_id);
    if (nc_err != NC_NOERR) {
      fprintf(stderr,"Cannot write %s \n",phi_surface_nc);
      exit(1);
    }
    WRITEC(planet->name,planet_name,32);
    WRITEI(&grid.nj,grid_nj,1);
    WRITEI(&grid.ni,grid_ni,1);
    WRITEC(grid.geometry,grid_geometry,GEOM_STR);
    WRITEF(&grid.globe_lonbot,grid_globe_lonbot,1);
    WRITEF(&grid.globe_lontop,grid_globe_lontop,1);
    WRITEF(&grid.globe_latbot,grid_globe_latbot,1);
    WRITEF(&grid.globe_lattop,grid_globe_lattop,1);

    /*
     * PHI_SURFACE(J,I) resides on the p-grid.
     */
    jdim = grid.nj-grid.jlo+1;
    idim = grid.ni;

    /* 
     * lon (I direction): 
     */
    nc_def_dim(nc_id,"lon_p",idim,&var.phi_surface.info[0].dimid[NETCDF_I_INDEX]);
    nc_def_var(nc_id,"lon_p",nc_float_type,ONEDIM,
               &var.phi_surface.info[0].dimid[NETCDF_I_INDEX],
               &var.phi_surface.info[0].coorid[NETCDF_I_INDEX]);
    nc_put_att_text(nc_id,var.phi_surface.info[0].coorid[NETCDF_I_INDEX],"units",
                    strlen("degrees_east")+1,"degrees_east");
    /* 
     * lat (J direction):
     */
    nc_def_dim(nc_id,"lat_p",jdim,&var.phi_surface.info[0].dimid[NETCDF_J_INDEX]);
    nc_def_var(nc_id,"lat_p",nc_float_type,ONEDIM,
               &var.phi_surface.info[0].dimid[NETCDF_J_INDEX],
               &var.phi_surface.info[0].coorid[NETCDF_J_INDEX]);
    nc_put_att_text(nc_id,var.phi_surface.info[0].coorid[NETCDF_J_INDEX],"units",
                    strlen("degrees_north")+1,"degrees_north");
    nc_put_att_text(nc_id,var.phi_surface.info[0].coorid[NETCDF_J_INDEX],"mapping",
                    strlen("planetographic")+1,"planetographic");

    /* phi_surface */
    nc_def_var(nc_id,"phi_surface",nc_float_type,TWODIM,
               &var.phi_surface.info[0].dimid[NETCDF_J_INDEX],
               &var.phi_surface.info[0].id);
    nc_put_att_text(nc_id,var.phi_surface.info[0].id,"units",
                    strlen("J/kg")+1,"J/kg");

    /*
     * Leave define mode:
     */
    nc_enddef(nc_id);

    /*
     * Assign values to coordinates.
     * Longitude:
     */
    for (I = 1; I <= grid.ni; I++) {
      index[0] = I-1;
#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.phi_surface.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.phi_surface.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
#endif
    }
    /*
     * Latitude:
     */
    for (J = grid.jlo; J <= grid.nj; J++) {
      index[0] = J-grid.jlo;
#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.phi_surface.info[0].coorid[NETCDF_J_INDEX],
                         index,&(grid.lat[2*J+1]));
#else
      nc_put_var1_float(nc_id,var.phi_surface.info[0].coorid[NETCDF_J_INDEX],
                        index,&(grid.lat[2*J+1]));
#endif
    }
  }

  /*
   * Setup write_array().
   */
  num_nodes = setup_write_array();
  /*
   * Load start and end vectors.
   */
  start[0] = 1;
  end[  0] = grid.ni;

  for (node = 0; node < num_nodes; node++) {
    get_jlohi(node,num_nodes,start+1,end+1);
    write_array(node,TWODIM,start,end,0,var.phi_surface.info[0].name,
                var.phi_surface.info[0].index,var.phi_surface.value,EPIC_FLOAT_ARRAY,nc_id);
  }

  if (IAMNODE == NODE0) {
    /* Close file. */
    nc_close(nc_id);
  }

  return;
}

/*====================== end of init_phi_surface() ===========================*/

/*====================== is_phi_surface_file() ===============================*/

/*
 * Returns TRUE if entry->d_name contains the string "phi_surface.nc",
 * otherwise returns FALSE.
 */

int is_phi_surface_file(const struct dirent *entry) 
{
  if (strstr(entry->d_name,"phi_surface.nc")) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

/*====================== end of is_phi_surface_file() ========================*/

/*====================== init_with_u() =======================================*/

/* 
 * Initialize the prognostic variables given meridional-plane
 * zonal wind profile and T(p) sounding.
 *
 * The input parameters floor_tp and ceiling_tp indicate the needed index range of
 * the sounding profile (t_vs_p). 
 *
 * NOTE: Assumes bottom P is already assigned.
 * NOTE: Not MPI ready.
 */
#include <epic_pv_schemes.h>

#define MONT_FINE(iitp,J)  mont_fine[iitp][J-JLO]
#define EXNER_FINE(iitp,J) exner_fine[iitp][J-JLO]
#define P_FINE(iitp,J)     p_fine[iitp][J-JLO]
#define PHI_FINE(iitp,J)   phi_fine[iitp][J-JLO]
#define FP_FINE(iitp,J)    fp_fine[iitp][J-JLO]
#define T_FINE(iitp,J)     t_fine[iitp][J-JLO]
#define THETA_FINE(iitp,J) theta_fine[iitp][J-JLO]
#define RHO_FINE(iitp,J)   rho_fine[iitp][J-JLO]
#define U_FINE(iitp,J)     u_fine[iitp][J-JLO]

#undef  ACOMPACT
#define ACOMPACT(i,j) acompact[(m1+m2+1)*(i)+(j)]
#undef  AORIG
#define AORIG(i,j) aorig[(m1+m2+1)*(i)+(j)]
#undef  AL
#define AL(i,j) al[m1*(i)+(j)]

void init_with_u(planetspec        *planet,
                 int                floor_tp,
                 int                ceiling_tp,
                 init_defaultspec  *def,
                 EPIC_FLOAT       **Buff2D)
{
  int
    K,J,I,
    kk,kaykay,k_isen,jj,iq,
    iitp,i,ii,
    n,m1,m2,
    count,
   *index;
  EPIC_FLOAT
    pressure,temperature,density,
    fpara,mu,theta,
    fgibb,fpe,uoup,
    theta_ortho,theta_para,
    pbot,sigma,sum,
    x,dx,dy,dp,dmont,dphi,
    lnp,lnp_d,
    dz10,dz20,dz30,dz21,dz31,dz32,
    df0,df1,df2,df3,
    neglogp,h,tmp,tmp2;
  EPIC_FLOAT
    *mudat,
    *exnerdat,
    *rhodat,
    *thetadat,
    *fpdat,
    *udy,
    *pvhudy,
    *bern,
    *p,
    *hh,
    **mont_fine,
    **phi_fine,
    **theta_fine,
    **p_fine,
    **exner_fine,
    **rho_fine,
    **u_fine;
  EPIC_FLOAT
    *aorig,*acompact,*al,*b,*borig,d;
  float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_with_u";

  /*
   * Screen for bad inputs.
   */
  if (strcmp(planet->type,"gas-giant") != 0) {
    sprintf(Message,"not implemented for planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }
  if (var.ntp < 1) {
    sprintf(Message,"var.ntp=%d < 1\n",var.ntp);
    epic_error(dbmsname,Message);
  }

  /*
   * Allocate memory:
   */
   if (grid.coord_type == COORD_ISENTROPIC) {
     thetadat = fvector(0,var.ntp-1,dbmsname);
     exnerdat = fvector(0,var.ntp-1,dbmsname);
   }
   else {
     rhodat = fvector(0,var.ntp-1,dbmsname);
     mudat  = fvector(0,var.ntp-1,dbmsname);
   }

   fpdat        = fvector(0,var.ntp-1,dbmsname);
   udy          = fvector(0,grid.nj+1,dbmsname);
   pvhudy       = fvector(0,grid.nj+1,dbmsname);
   bern         = fvector(0,grid.nj+1,dbmsname);
   buff_triplet = ftriplet(0,IMAX(var.ntp,grid.nk+1),dbmsname);

   /*
    * Allocate 2D arrays for mont_fine or phi_fine, etc.
    */
   if (grid.coord_type == COORD_ISENTROPIC) {
     mont_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!mont_fine) {
       sprintf(Message,"unable to allocate mont_fine");
       epic_error(dbmsname,Message);
     }

     p_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!p_fine) {
       sprintf(Message,"unable to allocate p_fine");
       epic_error(dbmsname,Message);
     }

     exner_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!exner_fine) {
       sprintf(Message,"unable to allocate exner_fine");
       epic_error(dbmsname,Message);
     }

     for (iitp = 0; iitp < var.ntp; iitp++) {
       mont_fine[iitp]  = fvector(0,JHI-JLO,dbmsname);
       p_fine[iitp]     = fvector(0,JHI-JLO,dbmsname);
       exner_fine[iitp] = fvector(0,JHI-JLO,dbmsname);
     }
   }
   else {
     phi_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!phi_fine) {
       sprintf(Message,"unable to allocate phi_fine");
       epic_error(dbmsname,Message);
     }

     theta_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!theta_fine) {
       sprintf(Message,"unable to allocate theta_fine");
       epic_error(dbmsname,Message);
     }

     rho_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
     if (!rho_fine) {
       sprintf(Message,"unable to allocate rho_fine");
       epic_error(dbmsname,Message);
     }

     for (iitp = 0; iitp < var.ntp; iitp++) {
       phi_fine[iitp]   = fvector(0,JHI-JLO,dbmsname);
       theta_fine[iitp] = fvector(0,JHI-JLO,dbmsname);
       rho_fine[iitp]   = fvector(0,JHI-JLO,dbmsname);
     }
   }

   u_fine = (EPIC_FLOAT **)calloc(var.ntp,sizeof(EPIC_FLOAT *));
   if (!u_fine) {
     sprintf(Message,"unable to allocate u_fine");
     epic_error(dbmsname,Message);
   }
   for (iitp = 0; iitp < var.ntp; iitp++) {
     u_fine[iitp] = fvector(0,JHI-JLO,dbmsname);
   }

  /*
   * Set fpdat to equilibrium value.
   */
  for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
    fpdat[iitp] = return_fpe(var.tdat[iitp]);
  }

  if (grid.coord_type == COORD_ISENTROPIC) {
    /*
     * Set thetadat, exnerdat.
     */
    for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
      thetadat[iitp] = return_theta(planet,fpdat[iitp],var.pdat[iitp],var.tdat[iitp],&theta_ortho,&theta_para);
      exnerdat[iitp] = planet->cp*var.tdat[iitp]/thetadat[iitp];
    }
  }
  else {
    /*
     * Set rhodat, mudat.  The function mu_p() contains an estimate of the effects of 
     * optional condensables.
     */
    for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
      mudat[iitp]  = mu_p(var.pdat[iitp]);
      rhodat[iitp] = return_density(planet,fpdat[iitp],var.pdat[iitp],var.tdat[iitp],mudat[iitp],PASSING_T);
    }
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                           *
   * 1. Specify T(p) at one latitude and determine MONT or PHI *
   *    via hydrostatic balance.                               *
   *                                                           *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  if (grid.coord_type == COORD_ISENTROPIC) {
    /*
     * Determine fine (high vertical resolution) MONT at J = grid.jtp.
     * Vertically integrate the hydrostatic balance equation for this column.
     * Arbitrarily assign MONT_FINE(floor_tp,J) = 0 for now.
     *
     * NOTE: There is no map factor in the integration to get MONT_FINE.
     */

    J                 = grid.jtp;
    iitp              = floor_tp;
    MONT_FINE(iitp,J) = 0.;

    /*
     * Short centered step up from bottom.
     */
    iitp  = floor_tp+1;
    dmont = .5*(exnerdat[iitp]+exnerdat[iitp-1])*(thetadat[iitp]-thetadat[iitp-1]);
    MONT_FINE(iitp,J) = MONT_FINE(iitp-1,J)+dmont;

    /* 
     * Use the quadrature scheme of Leslie and Purser (1992).
     */
    for (iitp = floor_tp+2; iitp <= ceiling_tp-1; iitp++) {
      /* 
       * This accurate integration scheme is based on the one described and illustrated in Fig. 8(d) of 
       *   Leslie LM,  Purser RJ, 1992,  A comparative study of the performance of various
       *   vertical discretization schemes, Meteorol. Atmos. Phys. 50, 61-73.
       */
      dz10 = thetadat[iitp-1]-thetadat[iitp-2];
      dz21 = thetadat[iitp  ]-thetadat[iitp-1];
      dz32 = thetadat[iitp+1]-thetadat[iitp  ];
      dz20 = dz21+dz10;
      dz30 = dz32+dz20;
      dz31 = dz32+dz21;

      df0  = exnerdat[iitp-2];
      df1  = exnerdat[iitp-1];
      df2  = exnerdat[iitp  ];
      df3  = exnerdat[iitp+1];

      dmont = -(df0*(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)
               +df1*(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)
               +df2*(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)
               +df3*(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)               )/12.;

      if (dmont <= 0.) {
        /*
         * We have found that the Leslie & Purser (1992) integration scheme can be too stiff
         * when one of the dz values is much smaller than the others. In this case, we revert to
         * a less accurate but monotonic scheme (the short-centered one used at the bottom and top).
         */
        dmont = .5*(exnerdat[iitp]+exnerdat[iitp-1])*dz21;
      }

      if (dmont < 0.) {
        sprintf(Message,"dmont=%g <= 0., var.pdat[%d]=%g, var.pdat[%d]=%g hPa; thetadat[%d]=%g, thetadat[%d]=%g",
                         dmont,iitp,var.pdat[iitp]/100.,iitp-1,var.pdat[iitp-1]/100.,
                               iitp,thetadat[iitp],     iitp-1,thetadat[iitp-1]);
        epic_error(dbmsname,Message);
      }

      MONT_FINE(iitp,J) = MONT_FINE(iitp-1,J)+dmont;
    }

    /*
     * Short centered step up to top.
     */
    iitp              = ceiling_tp;
    dmont             = .5*(exnerdat[iitp]+exnerdat[iitp-1])*(thetadat[iitp]-thetadat[iitp-1]);
    MONT_FINE(iitp,J) = MONT_FINE(iitp-1,J)+dmont;
  }
  else {
    /*
     * Determine fine (high vertical resolution) PHI at J = grid.jtp.
     * Vertically integrate the hydrostatic balance equation for this column.
     * Arbitrarily assign PHI_FINE(floor_tp,J) = 0 for now.
     *
     * NOTE: There is no map factor in the integration to get PHI_FINE.
     */

    J                = grid.jtp;
    iitp             = floor_tp;
    PHI_FINE(iitp,J) = 0.;

    /*
     * Short centered step up from bottom.
     */
    iitp = floor_tp+1;
    dphi = sqrt(var.pdat[iitp-1]/rhodat[iitp-1]*var.pdat[iitp]/rhodat[iitp])
          *log(var.pdat[iitp-1]/var.pdat[iitp]);
    PHI_FINE(iitp,J) = PHI_FINE(iitp-1,J)+dphi;

    /* 
     * Use the quadrature scheme of Leslie and Purser (1992).
     */
    for (iitp = floor_tp+2; iitp <= ceiling_tp-1; iitp++) {
      /* 
       * This accurate integration scheme is based on the one described and illustrated in Fig. 8(d) of 
       *   Leslie LM,  Purser RJ, 1992,  A comparative study of the performance of various
       *   vertical discretization schemes, Meteorol. Atmos. Phys. 50, 61-73.
       */
      dz10 = log(var.pdat[iitp-2]/var.pdat[iitp-1]);
      dz21 = log(var.pdat[iitp-1]/var.pdat[iitp  ]);
      dz32 = log(var.pdat[iitp  ]/var.pdat[iitp+1]);
      dz20 = dz21+dz10;
      dz30 = dz32+dz20;
      dz31 = dz32+dz21;

      df0  = var.pdat[iitp-2]/rhodat[iitp-2];
      df1  = var.pdat[iitp-1]/rhodat[iitp-1];
      df2  = var.pdat[iitp  ]/rhodat[iitp  ];
      df3  = var.pdat[iitp+1]/rhodat[iitp+1];

      dphi = -(df0*(dz21*dz21*dz21*(dz31+dz32))/(dz10*dz20*dz30)
              +df1*(dz21*(dz21*(dz21+2.*(dz10-dz31))-6.*dz10*dz31))/(dz10*dz31)
              +df2*(dz21*(dz21*(dz21+2.*(dz32-dz20))-6.*dz20*dz32))/(dz20*dz32)
              +df3*(dz21*dz21*dz21*(dz10+dz20))/(dz30*dz31*dz32)               )/12.;

      if (dphi <= 0.) {
        /*
         * We have found that the Leslie & Purser (1992) integration scheme can be too stiff
         * when one of the dz values is much smaller than the others. In this case, we revert to
         * a less accurate but monotonic scheme (the short-centered one used at the bottom and top).
         */
        dphi = sqrt(var.pdat[iitp-1]/rhodat[iitp-1]*var.pdat[iitp]/rhodat[iitp])*dz21;
      }

      if (dphi < 0.) {
        sprintf(Message,"dphi=%g <= 0., var.pdat[%d]=%g, var.pdat[%d]=%g hPa",
                         dphi,iitp-1,var.pdat[iitp-1]/100.,iitp,var.pdat[iitp]/100.);
        epic_error(dbmsname,Message);
      }

      PHI_FINE(iitp,J) = PHI_FINE(iitp-1,J)+dphi;
    }

    /*
     * Short centered step up to top.
     */
    iitp             = ceiling_tp;
    dphi             = sqrt(var.pdat[iitp-1]/rhodat[iitp-1]*var.pdat[iitp]/rhodat[iitp])
                      *log(var.pdat[iitp-1]/var.pdat[iitp]);
    PHI_FINE(iitp,J) = PHI_FINE(iitp-1,J)+dphi;
  }

  /* * * * * * * * * * * * * * * * * * * *
   *                                     *
   * 2. Specify target U_FINE.           *
   *                                     *
   * * * * * * * * * * * * * * * * * * * */

  if (def->u_scale == 0.) {
    for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
      for (J = JLO; J <= JHI; J++) {
        U_FINE(iitp,J) = 0.;
      }
    }
  }
  else {
    for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
      for (J = JLO; J <= JHI; J++) {
        U_FINE(iitp,J) = def->u_scale*planet->u(var.pdat[iitp],grid.lat[2*J+1]);
      }
    }
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *                                                                   *
   * 3. Integrate the gradient-balance equation horizontally to get    *
   *    MONT_FINE or PHI_FINE.                                         *
   *                                                                   *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /*
   * Temporarily use the KLO plane of the U, V, and PV2 arrays, in order to
   * employ the same algorithms for (zeta+f)*u*dy as in the DVDT code.
   */
  K = KLO;
  for (iitp = ceiling_tp; iitp >= floor_tp; iitp--) {
    if (grid.coord_type == COORD_ISENTROPIC) {
      /*
       * Determine kk for which grid.theta_ref[kk] is closest to thetadat[iitp].
       * Use the corresponding map factors.
       *
       * NOTE: kk is not related to K in this iitp loop.
       */
      tmp2 = FLOAT_MAX;
      for (kaykay = 1; kaykay <= 2*KHI+1; kaykay++) {
        tmp = fabs(thetadat[iitp]-grid.theta_ref[kaykay]);
        if (tmp < tmp2) {
          kk   = kaykay;
          tmp2 = tmp;
        } 
      }
    }
    else {
      /*
       * Determine kk for which grid.p_ref[kk] is closest to var.pdat[iitp].
       * Use the corresponding map factors.
       *
       * NOTE: kk is not related to K in this iitp loop.
       */
      tmp2 = FLOAT_MAX;
      for (kaykay = 1; kaykay <= 2*KHI+1; kaykay++) {
        tmp = fabs(var.pdat[iitp]-grid.p_ref[kaykay]);
        if (tmp < tmp2) {
          kk   = kaykay;
          tmp2 = tmp;
        } 
      }
    }

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        U(grid.it_uv,K,J,I) = U_FINE(iitp,J);
        V(grid.it_uv,K,J,I) = 0.;
      }
    }
    bc_lateral(var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,TWODIM);
    bc_lateral(var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,TWODIM);

    /*
     * Calculate absolute vorticity.
     * Store in PV so that epic_pv_schemes.h macros work.
     * Use the same scheme as for dvdt in timestep().
     *
     * NOTE: kk is not related to K in this iitp loop.
     */
    vorticity(ON_SIGMATHETA,ABSOLUTE,kk,
              var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
              var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
              NULL,
              var.pv2.value+(K-Kshift)*Nelem2d);
    /*
     * Form udy. 
     */
    I = ILO;
    for (J = JLO; J <= JHI; J++) {
      dy     = 1./grid.n[kk][2*J+1];
      udy[J] = U(grid.it_uv,K,J,I)*dy;
    }
    /*
     * Form (zeta+f)*u*dy = pvhudy.
     */
    for (J = JFIRST; J <= JHI; J++) {
      pvhudy[J] = (GA_V*udy[J  ]+DE_V*udy[J  ]
                  +AL_V*udy[J-1]+BE_V*udy[J-1])*PV_COEF;
    }

    if (grid.coord_type == COORD_ISENTROPIC) {
      /*
       * Calculate bern by integrating -pvhudy.
       * Calculate MONT_FINE = bern-kin.
       */
      bern[grid.jtp] = MONT_FINE(iitp,grid.jtp)+get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                               var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,grid.jtp,I);
      for (J = grid.jtp+1; J <= JHI; J++) {
        bern[J]           = bern[J-1]-pvhudy[J];
        MONT_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                   var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,J,I);
      }
      for (J = grid.jtp-1; J >= JLO; J--) {
        bern[J]           = bern[J+1]+pvhudy[J+1];
        MONT_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                   var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,J,I);
      }
    }
    else {
      /*
       * Calculate bern by integrating -pvhudy.
       * Calculate PHI_FINE = bern-kin.
       */
      bern[grid.jtp] = PHI_FINE(iitp,grid.jtp)+get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                              var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,grid.jtp,I);
      for (J = grid.jtp+1; J <= JHI; J++) {
        bern[J]          = bern[J-1]-pvhudy[J];
        PHI_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                  var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,J,I);
      }
      for (J = grid.jtp-1; J >= JLO; J--) {
        bern[J]          = bern[J+1]+pvhudy[J+1];
        PHI_FINE(iitp,J) = bern[J]-get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                                                  var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,kk,J,I);
      }
    }
  }

  if (grid.coord_type == COORD_ISENTROPIC) {
    /*
     * Screen for non-monotonic MONT_FINE.
     */
    for (J = JLO; J <= JHI; J++) {
      for (iitp = ceiling_tp-1; iitp >= floor_tp; iitp--) {
        if (MONT_FINE(iitp,J) >= MONT_FINE(iitp+1,J)) {
          sprintf(Message,"non-monotonic: pdat=%ehPa; MONT_FINE(%d,%d)=%g >= MONT_FINE(%d,%d)=%g",
                           var.pdat[iitp]/100.,iitp,J,MONT_FINE(iitp,J),iitp+1,J,MONT_FINE(iitp+1,J));
          epic_error(dbmsname,Message);
        }
      }
    }
  }
  else {
    /*
     * Screen for non-monotonic PHI_FINE.
     */
    for (J = JLO; J <= JHI; J++) {
      for (iitp = ceiling_tp-1; iitp >= floor_tp; iitp--) {
        if (PHI_FINE(iitp,J) >= PHI_FINE(iitp+1,J)) {
          sprintf(Message,"non-monotonic: pdat=%ehPa; PHI_FINE(%d,%d)=%g >= PHI_FINE(%d,%d)=%g",
                           var.pdat[iitp]/100.,iitp,J,PHI_FINE(iitp,J),iitp+1,J,PHI_FINE(iitp+1,J));
          epic_error(dbmsname,Message);
        }
      }
    }
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
   *                                                                         *
   * 4. Calculate EXNER_FINE(lat,p) from MONT_FINE via hydrostatic balance.  *
   *    Get P_FINE(lat,p) = p(theta,exner).                                  *
   *                                                                         *
   *                                -OR-                                     *
   *                                                                         *
   *    Calculate RHO_FINE(lat,p) from PHI_FINE via hydrostatic balance.     *
   *    Get THETA_FINE(lat,p) = theta(p,rho).                                *
   *                                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /* 
   * NOTE: For convenience, assuming initial fpara on each pressure level is given by fpdat rather than fpe(tdat(fpara,pdat)),
   *       which would require iteration without much gain in accuracy. Likewise, assuming mu = mudat.
   *
   * NOTE: We tried the 4-banded inversion of the Leslie and Purser quadrature here, but had trouble with it.
   */

  if (grid.coord_type == COORD_ISENTROPIC) {
    for (J = JLO; J <= JHI; J++) {
      for (iitp = floor_tp+1; iitp < ceiling_tp; iitp++) {
        EXNER_FINE(iitp,J) = (MONT_FINE(iitp+1,J)-MONT_FINE(iitp-1,J))/(thetadat[iitp+1]-thetadat[iitp-1]);
      }

      iitp               = ceiling_tp;
      EXNER_FINE(iitp,J) = exnerdat[iitp];

      iitp               = floor_tp;
      dmont              = MONT_FINE(iitp+1,J)-MONT_FINE(iitp,J);
      EXNER_FINE(iitp,J) = 2.*dmont/(thetadat[iitp+1]-thetadat[iitp])-EXNER_FINE(iitp+1,J);
    }

    /*
     * Calculate P_FINE.
     */
    for (J = JLO; J <= JHI; J++) {
      for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
        temperature    = EXNER_FINE(iitp,J)*thetadat[iitp]/planet->cp;
        P_FINE(iitp,J) = return_press(planet,fpdat[iitp],temperature,thetadat[iitp]);
      }
    }
  }
  else {
    for (J = JLO; J <= JHI; J++) {
      for (iitp = floor_tp+1; iitp < ceiling_tp; iitp++) {
        tmp              = (PHI_FINE(iitp+1,J)-PHI_FINE(iitp-1,J))/(var.pdat[iitp]*log(var.pdat[iitp-1]/var.pdat[iitp+1]));
        RHO_FINE(iitp,J) = 1./tmp;
      }

      iitp             = ceiling_tp;
      RHO_FINE(iitp,J) = rhodat[iitp];

      iitp             = floor_tp;
      dphi             = PHI_FINE(iitp+1,J)-PHI_FINE(iitp,J);
      RHO_FINE(iitp,J) = var.pdat[iitp]/(2.*dphi/log(var.pdat[iitp]/var.pdat[iitp+1])-var.pdat[iitp+1]/RHO_FINE(iitp+1,J));
    }

    /*
     * Calculate THETA_FINE.
     */
    for (J = JLO; J <= JHI; J++) {
      for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
        temperature        = alt_return_temp(planet,fpdat[iitp],var.pdat[iitp],mudat[iitp],RHO_FINE(iitp,J));
        THETA_FINE(iitp,J) = return_theta(planet,fpdat[iitp],var.pdat[iitp],temperature,&theta_ortho,&theta_para);
      }
    }
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   *                                                         *
   * 5. The fine-resolution variables are now ready to be    *
   *    interpolated onto the model's vertical coordinate.   *
   *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      I = ILO;
      for (J = JLO; J <= JHI; J++) {
        if (def->spacing_type == SPACING_P) {
          count = 0;
          for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
            buff_triplet[count].x = thetadat[iitp];
            buff_triplet[count].y = P_FINE(iitp,J);
            count++;
          }

          spline_pchip(count,buff_triplet);

          ii = -2;
          for (K = KLO-1; K <= KHI; K++) {
            x            = grid.sigmatheta[2*K+1];
            THETA(K,J,I) = x;
            ii           = hunt_place_in_table(count,buff_triplet,x,&dx,ii);
            P3(K,J,I)    = splint_pchip(x,buff_triplet+ii,dx);
          }
        }
        else {
          count = 0;
          for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
            buff_triplet[count].x = thetadat[iitp];
            buff_triplet[count].y = log(P_FINE(iitp,J));
            count++;
          }

          spline_pchip(count,buff_triplet);

          ii = -2;
          for (K = KLO-1; K <= KHI; K++) {
            x            = grid.sigmatheta[2*K+1];
            THETA(K,J,I) = x;
            ii           = hunt_place_in_table(count,buff_triplet,x,&dx,ii);
            P3(K,J,I)    = exp(splint_pchip(x,buff_triplet+ii,dx));
          }
        }
        /*
         * K = KHI+1 is the abyssal layer, used for gas giants.
         * Assume an adiabatic temperature lapse rate (constant THETA).
         */
        K            = KHI+1;
        THETA(K,J,I) = THETA(K-1,J,I);
      }

      /*
       * Extend in I direction.
       */
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO+1; I <= IHI; I++) {
            P3(K,J,I)    = P3(K,J,ILO);
            THETA(K,J,I) = THETA(K,J,ILO);
          }
        }
      }
    break;
    case COORD_ISOBARIC:
      I     = ILO;
      for (J = JLO; J <= JHI; J++) {
        if (def->spacing_type == SPACING_P) {
          count = 0;
          for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
            buff_triplet[count].x = -var.pdat[iitp];
            buff_triplet[count].y = THETA_FINE(iitp,J);
            count++;
          }

          spline_pchip(count,buff_triplet);

          ii = -2;
          for (K = KLO-1; K <= KHI; K++) {
            x            = -grid.p_ref[2*K+1];
            P3(K,J,I)    = -x;
            ii           = hunt_place_in_table(count,buff_triplet,x,&dx,ii);
            THETA(K,J,I) = splint_pchip(x,buff_triplet+ii,dx);
          }
        }
        else {
          count = 0;
          for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
            buff_triplet[count].x = -log(var.pdat[iitp]);
            buff_triplet[count].y = THETA_FINE(iitp,J);
            count++;
          }

          spline_pchip(count,buff_triplet);

          ii = -2;
          for (K = KLO-1; K <= KHI; K++) {
            x            = -log(grid.p_ref[2*K+1]);
            P3(K,J,I)    = grid.p_ref[2*K+1];
            ii           = hunt_place_in_table(count,buff_triplet,x,&dx,ii);
            THETA(K,J,I) = splint_pchip(x,buff_triplet+ii,dx);
          }
        }
        /*
         * K = KHI+1 is the abyssal layer, used for gas giants.
         * Assume an adiabatic temperature lapse rate (constant THETA).
         */
        K            = KHI+1;
        THETA(K,J,I) = THETA(K-1,J,I);
      }

      /*
       * Extend in I direction.
       */
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO+1; I <= IHI; I++) {
            P3(K,J,I)    = P3(K,J,ILO);
            THETA(K,J,I) = THETA(K,J,ILO);
          }
        }
      }
    break;
    case COORD_HYBRID:
      k_isen = grid.k_sigma-2;
      fprintf(stderr,"Interpolating variables onto hybrid-coordinate surfaces.\n");
      for (J = JLO; J <= JHI; J++) {
        /*
         * Loop over I to handle terrain-following coordinate.
         *
         * Find P3 using sigmatheta as the independent variable.
         */
        for (I = ILO; I <= IHI; I++) {
          pbot = P3(KHI,J,I);
          /*
           * Start with p=pbot table entry.
           */
          buff_triplet[0].x = grid.zeta0;
          if (def->spacing_type == SPACING_P) {
            buff_triplet[0].y = pbot;
          }
          else {
            buff_triplet[0].y = log(pbot);
          }
          count = 1;
          for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
            if (fcmp(var.pdat[iitp],pbot) < 0) {
              buff_triplet[count].x = return_sigmatheta(THETA_FINE(iitp,J),var.pdat[iitp],pbot);
              /*
               * Check that x is monotonic.
               */
              if (buff_triplet[count].x <= buff_triplet[count-1].x) {
                sprintf(Message,"J=%2d, iitp=%2d, sgth=%g %g not monotonically increasing (maybe try a different p_sigma)",
                                 J,iitp,buff_triplet[count-1].x,buff_triplet[count].x);
                epic_error(dbmsname,Message);
              }

              if (def->spacing_type == SPACING_P) {
                buff_triplet[count].y = var.pdat[iitp];
              }
              else {
                buff_triplet[count].y = log(var.pdat[iitp]);
              }
              count++;
            }
          }

          spline_pchip(count,buff_triplet);

          /*
           * NOTE: A cubic spline may be too stiff, reverting to linear interpolation.
           */
          ii = -2;
          for (K = KLO-1; K < KHI; K++) {
            x  = grid.sigmatheta[2*K+1];
            ii = hunt_place_in_table(count,buff_triplet,x,&dx,ii);

            if (def->spacing_type == SPACING_P) {
              P3(K,J,I) = linint(x,buff_triplet+ii,dx);
            }
            else {
              P3(K,J,I) = exp(linint(x,buff_triplet+ii,dx));
            }
          }

          /*
           * Set THETA to its diagnostic value.
           */
          for (K = KLO; K <= k_isen; K++) {
            sigma        = get_sigma(pbot,P3(K,J,I));
            THETA(K,J,I) = (grid.sigmatheta[2*K+1]-f_sigma(sigma))/g_sigma(sigma);
          }

          /*
           * Interpolate to get THETA in the lower region, using the appropriate function of P3
           * as the independent variable.
           */
          if (def->spacing_type == SPACING_P) {
            for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
              buff_triplet[iitp-floor_tp].x = -var.pdat[iitp];
              buff_triplet[iitp-floor_tp].y = THETA_FINE(iitp,J);
            }

            spline_pchip(ceiling_tp-floor_tp+1,buff_triplet);

            ii = -2;
            for (K = k_isen+1; K <= KHI; K++) {
              x            = -P3(K,J,I);
              ii           = hunt_place_in_table(ceiling_tp-floor_tp+1,buff_triplet,x,&dx,ii);
              THETA(K,J,I) = splint_pchip(x,buff_triplet+ii,dx);
            }
          }
          else {
            for (iitp = floor_tp; iitp <= ceiling_tp; iitp++) {
              /*
               * Need independent variable to be increasing for spline_pchip().
               */
              buff_triplet[iitp-floor_tp].x = -log(var.pdat[iitp]);
              buff_triplet[iitp-floor_tp].y = THETA_FINE(iitp,J);
            }

            spline_pchip(ceiling_tp-floor_tp+1,buff_triplet);

            ii = -2;
            for (K = k_isen+1; K <= KHI; K++) {
              x            = -log(P3(K,J,I));
              ii           = hunt_place_in_table(ceiling_tp-floor_tp+1,buff_triplet,x,&dx,ii);
              THETA(K,J,I) = splint_pchip(x,buff_triplet+ii,dx);
            }
          }

          /*
           * K = KHI+1 is the abyssal layer, used for gas giants.
           * Assume an adiabatic temperature lapse rate (constant THETA).
           */
          K            = KHI+1;
          THETA(K,J,I) = THETA(K-1,J,I);

        } /* loop over I */
      } /* loop over J */
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  bc_lateral(var.p3.value,   THREEDIM);
  bc_lateral(var.theta.value,THREEDIM);

  /* * * * * * * * * *
   *                 *
   * 5. Set H, U     *
   *                 *
   * * * * * * * * * */

  /*
   * Set H.
   */
  p  = fvector(0,2*grid.nk+1,dbmsname);
  hh = fvector(0,2*grid.nk+1,dbmsname);
  for (J = JLOPAD; J <= JHIPAD; J++) {
    jj = 2*J+1;
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        p[kk] = get_p(planet,P2_INDEX,kk,J,I);
      }
      calc_h(jj,p,hh);

      for (K = KLO; K <= KHI; K++) {
        H(K,J,I) = hh[2*K];
      }
      K = 0;
      H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
      K = KHI+1;
      H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
    }
  }
  free_fvector(p, 0,2*grid.nk+1,dbmsname);
  free_fvector(hh,0,2*grid.nk+1,dbmsname);

  /*
   * Set P2, HDRY, THETA2, etc.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  /*
   * Set U.
   */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      for (K = KLO; K <= KHI; K++) {
        U(grid.it_uv,K,J,I) = def->u_scale*planet->u(P2(K,J,I),grid.lat[2*J+1]);
      }
      K = 0;
      U(grid.it_uv,K,J,I) = def->u_scale*planet->u(P3(K,  J,I),grid.lat[2*J+1]);
      K = KHI+1;
      U(grid.it_uv,K,J,I) = def->u_scale*planet->u(P3(K-1,J,I),grid.lat[2*J+1]);
    }
  }
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

  /*
   * Initial V is zero.
   */

  /*
   * Free allocated memory:
   */
  free_ftriplet(buff_triplet,0,IMAX(var.ntp,grid.nk+1),dbmsname);

  if (grid.coord_type == COORD_ISENTROPIC) {
    for (iitp = 0; iitp < var.ntp; iitp++) {
      free_fvector(exner_fine[iitp],0,JHI-JLO,dbmsname);
      free_fvector(p_fine[iitp],    0,JHI-JLO,dbmsname);
      free_fvector(mont_fine[iitp], 0,JHI-JLO,dbmsname);
      free_fvector(u_fine[iitp],    0,JHI-JLO,dbmsname);
    }
    free(exner_fine);
    free(p_fine);
    free(mont_fine);

    free_fvector(thetadat,0,var.ntp-1,dbmsname);
    free_fvector(exnerdat,0,var.ntp-1,dbmsname);
  }
  else {
    for (iitp = 0; iitp < var.ntp; iitp++) {
      free_fvector(rho_fine[iitp],  0,JHI-JLO,dbmsname);
      free_fvector(theta_fine[iitp],0,JHI-JLO,dbmsname);
      free_fvector(phi_fine[iitp],  0,JHI-JLO,dbmsname);
      free_fvector(u_fine[iitp],    0,JHI-JLO,dbmsname);
    }

    free(rho_fine);
    free(theta_fine);
    free(phi_fine);

    free_fvector(rhodat,0,var.ntp-1,dbmsname);
    free_fvector(mudat, 0,var.ntp-1,dbmsname);
  }
  free(u_fine);
  free_fvector(bern,    0,grid.nj+1,dbmsname);
  free_fvector(pvhudy,  0,grid.nj+1,dbmsname);
  free_fvector(udy,     0,grid.nj+1,dbmsname);
  free_fvector(fpdat,   0,var.ntp-1,dbmsname);

  return;
}

#undef MONT_FINE
#undef EXNER_FINE
#undef P_FINE
#undef PHI_FINE
#undef FP_FINE
#undef T_FINE
#undef THETA_FINE
#undef RHO_FINE
#undef U_FINE

/*====================== end of init_with_u() ================================*/

/*====================== init_with_ref() =====================================*/
    
/* 
 * For the hybrid vertical coordinate, initialize p in the upper region
 * using grid.p_ref[kk], and THETA in the lower region using grid.theta_ref[kk].
 * Use the diagnostic values for P and THETA in the complementary regions.
 *
 * Assumes the surface pressure has already been set.
 */

void init_with_ref(planetspec  *planet,
                   EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kk,ki,jj,iq;
  EPIC_FLOAT
    pbot,sigma,gsg,
    s,s0,s1,
    p_diag,theta_diag,
    sum,lnp,lnp_d,
   *p,*h;
  float_triplet
    theta_table[grid.nk+1];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_with_ref";

  if (grid.coord_type == COORD_HYBRID) {
    s1 = grid.sigmatheta[2*grid.k_sigma-1];
    s0 = grid.sigmatheta[2*KHI+1];
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        K         = 0;
        P3(K,J,I) = grid.ptop;
        pbot      = P3(grid.nk,J,I);
        /*
         * Use P3 = p_ref in hybrid portion of model.
         */
        for (K = KLO; K < grid.k_sigma; K++) {
          P3(K,J,I) = grid.p_ref[2*K+1];
        }
        /*
         *  Transition to terrain-following in lower portion of model.
         */
        for (K = grid.k_sigma; K < KHI; K++) {
          s         = (grid.sigmatheta[2*K+1]-s0)/(s1-s0);
          P3(K,J,I) = pbot*exp(s*log(P3(grid.k_sigma-1,J,I)/pbot));
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        K         = 0;
        P3(K,J,I) = grid.ptop;
        for (K = KLO; K < KHI; K++) {
          P3(K,J,I) = grid.p_ref[2*K+1];
        }
        K = KHI+1;
        P3(K,J,I) = grid.p_ref[2*K+1];
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * Set H.
   */
  p = fvector(0,2*grid.nk+1,dbmsname);
  h = fvector(0,2*grid.nk+1,dbmsname);
  for (J = JLOPAD; J <= JHIPAD; J++) {
    jj = 2*J+1;
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (kk = 1; kk <= 2*KHI+1; kk++) {
        p[kk] = get_p(planet,P2_INDEX,kk,J,I);
      }
      calc_h(jj,p,h);

      for (K = KLO; K <= KHI; K++) {
        H(K,J,I) = h[2*K];
      }
      K = 0;
      H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
      K = KHI+1;
      H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
    }
  }
  free_fvector(p,0,2*grid.nk+1,dbmsname);
  free_fvector(h,0,2*grid.nk+1,dbmsname);

  if (grid.coord_type == COORD_HYBRID) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        /*
         * Set THETA by interpolating theta_ref with respect to log(p).
         */
        for (K = KLOPAD; K <= KHIPAD; K++) {
          theta_table[K].x = log(grid.p_ref[2*K+1]);
          theta_table[K].y = grid.theta_ref[2*K+1];
        }
        spline_pchip(grid.nk+1,theta_table);
        ki = -2;
        for (K = KLOPAD; K <= KHIPAD; K++) {
          if (P3(K,J,I) <= grid.p_ref[2*grid.nk+1]) {
            lnp          = log(P3(K,J,I));
            ki           = hunt_place_in_table(grid.nk+1,theta_table,lnp,&lnp_d,ki);
            THETA(K,J,I) = splint_pchip(lnp,theta_table+ki,lnp_d);
          }
          else {
            /*
             * Below the end of the reference values, assume the lapse rate is adiabatic.
             */
            THETA(K,J,I) = THETA(K-1,J,I);
          }
        }
        /*
         * K = KHI+1 is the abyssal layer, used for gas giants.
         * Assume an adiabatic temperature lapse rate (constant THETA).
         */
        K            = KHI+1;
        THETA(K,J,I) = THETA(K-1,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  else {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          THETA(K,J,I) = grid.theta_ref[2*K+1];
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * Set P2, THETA2, etc.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  /*
   * Start with U = 0, V = 0, which requires no action.
   */

  return;
}

/*====================== end of init_with_ref() ==============================*/

/*====================== init_fpara_as_fpe() =================================*/

/*
 * Use a root finder to solve for fp using the implicit equation:
 *   fp = fpe(T(p,theta,fp))
 * Global variables used to communicate with fpe_minus_fpe():
 */

void init_fpara_as_fpe(planetspec *planet)
{
  int
    K,J,I;
  EPIC_FLOAT
    fpe,fptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_fpara_as_fpe";

  /* Check that fpara is turned on. */
  if (!var.fpara.on) {
    fprintf(stderr,"Warning: init_fpara_as_fpe(): var.fpara.on is off \n");
    return;
  }

  FPEMFPE_planet = planet;

  fptol = pow(machine_epsilon(),2./3.);

  for (K = 0; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /*
         * Global variables used to communicate with fpe_minus_fpe().
         */
        FPEMFPE_p     = P3(   K,J,I);
        FPEMFPE_theta = THETA(K,J,I);
        find_root(0.25,1.00,fptol,&fpe,fpe_minus_fpe);

        FPARA(K,J,I) = fpe;
      }
    }
  }

  return;
}

/*====================== end of init_fpara_as_fpe() ==========================*/

/*====================== init_species() ======================================*/

/*
 * Top function for initializing optional species (condensables) in the atmosphere.
 *
 * Species are carried as Q(is,ip,K,J,I), where Q is the mass mixing ratio,
 * equal to mass_density_i/mass_density_dry. The indices "is" and "ip" are the species index
 * (e.g. CH_4_INDEX) and phase index (e.g. VAPOR), respectively.
 *
 * Synchronized with Q is X(is,ip,K,J,I), the mole fraction, equal to
 * number_density/number_density_total.  Notice the different denominators, 'dry' vs. 'total',
 * between Q and X.
 *
 * NOTE: Many textbooks use "q" for "specific humidity", which is density_i/density_total,
 *       and not identical to the mass mixing ratio.
 *
 * NOTE: Need to calculate required diagnostic variables like T here,
 *       since this function is called before they are set.
 */

void init_species(planetspec        *planet,
                  init_defaultspec  *def,
                  int                prompt_mode,
                  EPIC_FLOAT       **Buff2D)
{
  int
    K,J,I,
    is,ip,iq;
  EPIC_FLOAT
    factor,
    mu_dry,
    pressure,theta,fpara,
    fpe,fptol;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_species";

  /*
   * Make sure pressure-related diagnostic variables are set.
   */
  set_p2_etc(planet,DONT_UPDATE_THETA,Buff2D);

  /*
   * Calculate T3.
   *
   * Global variable used to communicate with fpe_minus_fpe():
   */
  FPEMFPE_planet = planet;

  fptol = pow(machine_epsilon(),2./3.);

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (K = 0; K <= KHI; K++) {
        pressure = P3(K,J,I);
        theta    = THETA(K,J,I);
        if (var.fpara.on) {
          fpara = FPARA(K,J,I);
        }
        else {
          /*
           * Use a root finder to solve for fp using the implicit equation:
           *   fp = fpe(T(p,theta,fp))
           * Global variables used to communicate with fpe_minus_fpe():
           */
          FPEMFPE_p     = pressure;
          FPEMFPE_theta = theta;
          find_root(0.25,1.00,fptol,&fpe,fpe_minus_fpe);
          fpara = fpe;
        }
        T3(K,J,I) = return_temp(planet,fpara,pressure,theta);
      }
    }
  }

  /*
   * Calculate mole fractions, X (aka volume mixing ratios, vmr).
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      switch (is) {
        case NH_3_INDEX:
          switch (planet->index) {
            case JUPITER_INDEX:
            case SATURN_INDEX:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        case CH_4_INDEX:
          switch (planet->index) {
            case JUPITER_INDEX:
            case SATURN_INDEX:
            case URANUS_INDEX:
              init_vmr_via_obs(planet,is,&def->mole_fraction[is],prompt_mode);
            break;
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        case C_2H_2_INDEX:
          switch (planet->index) {
            case JUPITER_INDEX:
            case SATURN_INDEX:
              init_vmr_via_obs(planet,is,NULL,prompt_mode);
            break;
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        case C_2H_4_INDEX:
          switch (planet->index) {
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        case C_2H_6_INDEX:
          switch (planet->index) {
            case JUPITER_INDEX:
            case SATURN_INDEX:
              init_vmr_via_obs(planet,is,NULL,prompt_mode);
            break;
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        case PH_3_INDEX:
          switch (planet->index) {
            default:
              init_vmr_via_deep_value(planet,is,&def->mole_fraction[is],
                                      &def->mole_fraction_over_solar[is],&def->rh_max[is],prompt_mode);
            break;
          }
        break;
        default:
          sprintf(Message,"%s not yet implemented",var.species[is].info[0].name);
          epic_error(dbmsname,Message);
        break;
      }
      /*
       * Set abyssal layer to equal KHI layer.
       */
      K = KHI+1;
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (var.species[is].phase[ip].on) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              X(is,ip,K,J,I) = X(is,ip,K-1,J,I);
            }
          }
          /* No need to call bc_lateral() here. */
        }
      }
    }
  }

  /*
   * Synchronize mass mixing ratios, Q, to X.
   */
  for (K = KLOPAD; K <= KHIPAD; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        mu_dry = R_GAS/planet->rgas;
        factor  = 1.;
        for (iq = 0; iq < grid.nq; iq++) {
          factor -= X(grid.is[iq],grid.ip[iq],K,J,I);
        }
        factor = 1./(mu_dry*factor);
        for (iq = 0; iq < grid.nq; iq++) {
          Q(grid.is[iq],grid.ip[iq],K,J,I) = X(grid.is[iq],grid.ip[iq],K,J,I)*var.species[grid.is[iq]].molar_mass*factor;
        }
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Update pressure-related diagnostic variables.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  return;
}

/*====================== end of init_species() ===============================*/

/*====================== init_vmr_via_deep_value() ===========================*/

/*
 * Assign volume mixing ratio (vmr, aka mole fraction) starting with a deep value,
 * which is assumed to be bounded by the saturation curve and a cold trap effect.
 */

void init_vmr_via_deep_value(planetspec *planet,
                             int         is,
                             EPIC_FLOAT *mole_fraction,
                             EPIC_FLOAT *mole_fraction_over_solar,
                             EPIC_FLOAT *rh_max,
                             int         prompt_mode)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    solar,
    sat_vapor_p,
    x_sat,
    mole_frac0;
  char
    min_element[4],
    prompt[64];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_vmr_via_deep_value";

  /*
   * Establish deep-atmosphere mole fractions and store in *mole_fraction.
   */
  if (prompt_mode == USE_PROMPTS) {
    solar = solar_fraction(var.species[is].info[0].name,MOLAR,min_element);
    fprintf(stdout,"Solar mole fraction of %s is %e \n",var.species[is].info[0].name,solar);
    sprintf(prompt,"Input %s mole fraction relative to solar [1.=solar]:\n",var.species[is].info[0].name);
    *mole_fraction_over_solar = input_float(prompt,*mole_fraction_over_solar);
    *mole_fraction            = *mole_fraction_over_solar*solar;
    *rh_max                   = input_float("Input maximum initial relative humidity [1.=sat]\n",*rh_max);
  }

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      /* 
       * Start at lowest layer with mole fractions set to specified value. 
       */
      mole_frac0 = *mole_fraction;

      for (K = KHI; K >= 0; K--) {
        /* 
         * Calculate saturation mole fraction.
         * Assume the mole fraction is bounded by the saturation curve,
         * and store it in the mass variable, Q.
         */
        sat_vapor_p       = var.species[is].sat_vapor_p(T3(K,J,I));
        x_sat             = sat_vapor_p/P3(K,J,I);
        X(is,VAPOR,K,J,I) = MIN(mole_frac0,x_sat*(*rh_max));

        /*
         * Assume that as we travel upwards, mole fraction stays 
         * reduced once it has been limited by saturation (cold trap effect).
         */
        mole_frac0 = X(is,VAPOR,K,J,I);
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  return;
}

/*====================== end of init_vmr_via_deep_value() ====================*/

/*====================== init_vmr_via_obs() ==================================*/

/*
 * Assign volume mixing ratio (vmr, aka mole fraction) using observations.
 * The intent is to continually update this function to use state-of-the-art
 * observations for each planet and each species. 
 */

void init_vmr_via_obs(planetspec *planet,
                      int         is,
                      EPIC_FLOAT *mole_fraction,
                      int         prompt_mode)
{
  int
    K,J,I,
    k;
  EPIC_FLOAT
    p,logp,dp,lat,temp,rh_trad,
    x_tmp,mole_frac0;
  char
    prompt[64],
    header[N_STR];
  static int
    np,
    initialized = FALSE;
  static float_triplet
    *methane_table;
  FILE
    *infile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_vmr_via_obs";

  if (!initialized) {
    if (var.species[CH_4_INDEX].on) {
      switch(planet->index) {
        case JUPITER_INDEX:
          infile = fopen(EPIC_PATH"/data/Jupiter/CH_4.Jupiter.ModelC","r");

          /* Skip 9-line header */
          for(k = 0; k < 9; k++) {
            fgets(header,N_STR,infile);
          }

          /* Input number of pressure levels. */
          fscanf(infile,"%d",&np);
          fgets(header,N_STR,infile);
          fgets(header,N_STR,infile);

          /* Allocate memory */
          methane_table = ftriplet(0,np-1,dbmsname);

          /* Input methane mole fraction data. */
          for (k = np-1; k >= 0; k--) {
            fscanf(infile,"%lf %lf",&p,&x_tmp);
            /* Convert pressure from [hPa] to [Pa], store as log(p) */
            methane_table[k].x = log(p*100.);
  
            /* Store mole fraction as log(x) */
            methane_table[k].y = log(x_tmp);
          }

          /* Set up spline. */
          spline_pchip(np,methane_table);

          fclose(infile);
        break;
        case SATURN_INDEX:
          infile = fopen(EPIC_PATH"/data/Saturn/CH_4.Saturn.ModelC","r");

          /* Skip 9-line header */
          for(K = 0; K < 9; K++) {
            fgets(header,N_STR,infile);
          }

          /* Input number of pressure levels. */
          fscanf(infile,"%d",&np);
          fgets(header,N_STR,infile);
          fgets(header,N_STR,infile);

          /* Allocate memory */
          methane_table = ftriplet(0,np-1,dbmsname);

          /* Input methane mole fraction data. */
          for (k = np-1; k >= 0; k--) {
            fscanf(infile,"%lf %lf",&p,&x_tmp);
            /* Convert pressure from [hPa] to [Pa], store as log(p) */
            methane_table[k].x = log(p*100.);

            /* Store mole fraction as log(x) */
            methane_table[k].y = log(x_tmp);
          }

          /* Set up spline. */
          spline_pchip(np,methane_table);

          fclose(infile);
        break;
        case URANUS_INDEX:
          /* 
           * CH_4 for Uranus is handled as a function, below.
           */
        break;
        default:
          sprintf(Message,"CH_4 not yet implemented for %s",planet->name);
          epic_error(dbmsname,Message);
        break;
      }
    }

    initialized = TRUE;
  }

  switch(planet->index) {
    case JUPITER_INDEX:
      switch(is) {
        case CH_4_INDEX:
          /*
           * The spline table is set up in the initialization above.
           */
          for (K = KLO-1; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              k = -2;
              for (I = ILO; I <= IHI; I++) {
                logp              = log(P3(K,J,I));
                k                 = hunt_place_in_table(np,methane_table,logp,&dp,k);
                X(is,VAPOR,K,J,I) = exp(splint_pchip(logp,methane_table+k,dp));
              }
            }
          }
          /* Need to call bc_lateral() here. */
          bc_lateral(var.species[is].phase[VAPOR].x,THREEDIM);
        break;
        case C_2H_2_INDEX:
        case C_2H_6_INDEX:
          for (K = KLO-1; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              for (I = ILO; I <= IHI; I++) {
                X(is,VAPOR,K,J,I) = CIRS_data(planet,is,P3(K,J,I),grid.lat[2*J+1]);
              }
            }
          }
          /* Need to call bc_lateral() here. */
          bc_lateral(var.species[is].phase[VAPOR].x,THREEDIM);
        break;
        default:
          sprintf(Message,"species %d not yet implemented",is);
          epic_error(dbmsname,Message);
        break;
      }
    break;
    case SATURN_INDEX:
      switch(is) {
        case CH_4_INDEX:
          /*
           * The spline table is set up in the initialization above.
           */
          for (K = KLO-1; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              k = -2;
              for (I = ILO; I <= IHI; I++) {
                logp              = log(P3(K,J,I));
                k                 = hunt_place_in_table(np,methane_table,logp,&dp,k);
                X(is,VAPOR,K,J,I) = exp(splint_pchip(logp,methane_table+k,dp));
              }
            }
          }
          /* Need to call bc_lateral() here. */
          bc_lateral(var.species[is].phase[VAPOR].x,THREEDIM);
        break;
        case C_2H_2_INDEX:
        case C_2H_6_INDEX:
          for (K = KLO-1; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              for (I = ILO; I <= IHI; I++) {
                X(is,VAPOR,K,J,I) = CIRS_data(planet,is,P3(K,J,I),grid.lat[2*J+1]);
              }
            }
          }
          /* Need to call bc_lateral() here. */

          bc_lateral(var.species[is].phase[VAPOR].x,THREEDIM);
        break;
        default:
          sprintf(Message,"species %d not yet implemented",is);
          epic_error(dbmsname,Message);
        break;
      }
    break;
    case URANUS_INDEX:
      switch(is) {
        case CH_4_INDEX:
          /*
           * Methane profile derived from Karkoschka & Tomasko, 2009, "The haze
           * and methane distributions on Uranus from HST-STIS spectroscopy".
           *
           * Relative humidity is 48% below 1 bar, and variable via the function:
           *   RH = 48% * [1 - (1-P)^2]
           * above 1 bar where P is pressure in bar.  This results in a relative
           * humidity of ~10% at the tropopause. (1 bar = 1.e5 Pa)
           */

          if (prompt_mode == USE_PROMPTS) {
            /*
             * Establish deep-atmosphere mole fractions and store in *mole_fraction.
             */
            sprintf(prompt,"Input %s deep mole fraction:\n",var.species[is].info[0].name);
            *mole_fraction = input_float(prompt,*mole_fraction);
          }
          /*
           * NOTE: We are using grid.aux_fa to record input value of deep mole fraction.
           */
          grid.aux_fa = *mole_fraction;
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              /* 
               * Start at lowest layer with mole fractions set to specified value.
               */
              mole_frac0  = *mole_fraction;
              for (K = KHI; K >= 0; K--) {
                /* rh_trad = (partial pressure)/(saturation partial pressure) */
                rh_trad = 0.48;
                if (P3(K,J,I) < 1.e5) {
                  rh_trad *= (1.-SQR(1.-P3(K,J,I)/1.e5));
                }
                x_tmp             = rh_trad*var.species[CH_4_INDEX].sat_vapor_p(T3(K,J,I))/P3(K,J,I);
                X(is,VAPOR,K,J,I) = MIN(x_tmp,mole_frac0);
                /*
                 * Assume that as we travel upwards, mole fraction stays 
                 * reduced once it has been limited by saturation (cold trap effect).
                 */
                mole_frac0 = X(is,VAPOR,K,J,I);
              }
            }
          }
          /* No need to call bc_lateral() here */
        break;
        default:
          sprintf(Message,"species %d not yet implemented",is);
          epic_error(dbmsname,Message);
        break;
      }
    break;
    default:
      sprintf(Message,"case planet->index=%d not yet implemented",planet->index);
      epic_error(dbmsname,Message);
    break;
  }
}

/*====================== end of init_vmr_via_obs() ============================*/

/*====================== CIRS_data() ==========================================*/

/*
 * Returns the volume mixing ratio (vmr) or temperature [K] 
 * from meridional-plane Cassini CIRS data, for the given
 * pressure [Pa] and latitude [deg planetographic].
 * Uses bilinear interpolation.
 */

#define DATA_T(k,j)      data_t[     j+(k)*nlat_t     ]
#define DATA_C_2H_2(k,j) data_C_2H_2[j+(k)*nlat_C_2H_2]
#define DATA_C_2H_6(k,j) data_C_2H_6[j+(k)*nlat_C_2H_6]

EPIC_FLOAT CIRS_data(planetspec *planet,
                     int         index,
                     EPIC_FLOAT  p,
                     EPIC_FLOAT  lat)
{
  int
    j,k;
  EPIC_FLOAT
    logp,frac_logp,frac_lat,
    ans;
  static int
    np_t, nlat_t,
    np_C_2H_2,nlat_C_2H_2,
    np_C_2H_6,nlat_C_2H_6,
    initialized = FALSE;
  static EPIC_FLOAT
    *logp_t,     *lat_t,     *data_t,
    *logp_C_2H_2,*lat_C_2H_2,*data_C_2H_2,
    *logp_C_2H_6,*lat_C_2H_6,*data_C_2H_6;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="CIRS_data()";

  if (!initialized) {
    switch(planet->index) {
      case JUPITER_INDEX:
        /* 
         * Determine size of data, allocate memory, and read in data.
         */
        if (var.t2.on || var.t3.on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/T.Jupiter.CIRS",SIZE_DATA,
                                &np_t,&nlat_t,NULL,NULL,NULL);
          logp_t = fvector(0,np_t   -1,dbmsname);
          lat_t  = fvector(0,nlat_t -1,dbmsname);
          data_t = fvector(0,np_t*nlat_t  -1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/T.Jupiter.CIRS",POST_SIZE_DATA,
                                &np_t,&nlat_t,logp_t,lat_t,data_t);
        }

        if (var.species[C_2H_2_INDEX].on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/C_2H_2.Jupiter.CIRS",SIZE_DATA,
                                &np_C_2H_2,&nlat_C_2H_2,NULL,NULL,NULL);
          logp_C_2H_2 = fvector(0,np_C_2H_2  -1,dbmsname);
          lat_C_2H_2  = fvector(0,nlat_C_2H_2-1,dbmsname);
          data_C_2H_2 = fvector(0,np_C_2H_2*nlat_C_2H_2-1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/C_2H_2.Jupiter.CIRS",POST_SIZE_DATA,
                                &np_C_2H_2,&nlat_C_2H_2,logp_C_2H_2,lat_C_2H_2,data_C_2H_2);
        }

        if (var.species[C_2H_6_INDEX].on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/C_2H_6.Jupiter.CIRS",SIZE_DATA,
                                &np_C_2H_6,&nlat_C_2H_6,NULL,NULL,NULL);
          logp_C_2H_6 = fvector(0,np_C_2H_6  -1,dbmsname);
          lat_C_2H_6  = fvector(0,nlat_C_2H_6-1,dbmsname);
          data_C_2H_6 = fvector(0,np_C_2H_6*nlat_C_2H_6-1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/C_2H_6.Jupiter.CIRS",POST_SIZE_DATA,
                                &np_C_2H_6,&nlat_C_2H_6,logp_C_2H_6,lat_C_2H_6,data_C_2H_6);
        }
      break;
      case SATURN_INDEX:
        /* 
         * Determine size of data, allocate memory, and read in data.
         */
        if (var.t2.on || var.t3.on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Saturn/T.Saturn.CIRS",SIZE_DATA,
                                &np_t,&nlat_t,NULL,NULL,NULL);
          logp_t = fvector(0,np_t   -1,dbmsname);
          lat_t  = fvector(0,nlat_t -1,dbmsname);
          data_t = fvector(0,np_t*nlat_t  -1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/aturn/T.Saturn.CIRS",POST_SIZE_DATA,
                                &np_t,&nlat_t,logp_t,lat_t,data_t);
        }

        if (var.species[C_2H_2_INDEX].on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Saturn/C_2H_2.Saturn.CIRS",SIZE_DATA,
                                &np_C_2H_2,&nlat_C_2H_2,NULL,NULL,NULL);
          logp_C_2H_2 = fvector(0,np_C_2H_2  -1,dbmsname);
          lat_C_2H_2  = fvector(0,nlat_C_2H_2-1,dbmsname);
          data_C_2H_2 = fvector(0,np_C_2H_2*nlat_C_2H_2-1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/Saturn/C_2H_2.Saturn.CIRS",POST_SIZE_DATA,
                                &np_C_2H_2,&nlat_C_2H_2,logp_C_2H_2,lat_C_2H_2,data_C_2H_2);
        }

        if (var.species[C_2H_6_INDEX].on) {
          read_meridional_plane(planet,EPIC_PATH"/data/Saturn/C_2H_6.Saturn.CIRS",SIZE_DATA,
                                &np_C_2H_6,&nlat_C_2H_6,NULL,NULL,NULL);
          logp_C_2H_6 = fvector(0,np_C_2H_6  -1,dbmsname);
          lat_C_2H_6  = fvector(0,nlat_C_2H_6-1,dbmsname);
          data_C_2H_6 = fvector(0,np_C_2H_6*nlat_C_2H_6-1,dbmsname);
          read_meridional_plane(planet,EPIC_PATH"/data/Saturn/C_2H_6.Saturn.CIRS",POST_SIZE_DATA,
                                &np_C_2H_6,&nlat_C_2H_6,logp_C_2H_6,lat_C_2H_6,data_C_2H_6);
        }
      break;
      default:
        sprintf(Message,"not yet implemented for %s",planet->name);
        epic_error(dbmsname,Message);
      break;
    }

    initialized = TRUE;
  }

  logp = log(p);

  /*
   * Use bilinear interpolation
   *   T(lat,logp)
   *   logX(lat,logp)
   */

  switch(index) {
    case T2_INDEX:
    case T3_INDEX:
      /*
       * Temperature
       */
      if (lat > lat_t[nlat_t-1]) {
        j        = nlat_t-2;
        frac_lat = (lat-lat_t[j])/(lat_t[j+1]-lat_t[j]);
      }
      else {
        for (j = 1; j < nlat_t; j++) {
          if (lat <= lat_t[j]) {
            j--;
            frac_lat = (lat-lat_t[j])/(lat_t[j+1]-lat_t[j]);
            break;
          }
        }  
      }

      if (logp > logp_t[np_t-1]) {
        k         = np_t-2;
        frac_logp = (logp-logp_t[k])/(logp_t[k+1]-logp_t[k]);
      }
      else {
        for (k = 1; k < np_t; k++) {
          if (logp <= logp_t[k]) {
            k--;
            frac_logp = (logp-logp_t[k])/(logp_t[k+1]-logp_t[k]);
            break;
          }
        }
      }

      ans = DATA_T(k,  j  )*(1.-frac_logp)*(1.-frac_lat)
           +DATA_T(k+1,j  )*(   frac_logp)*(1.-frac_lat)
           +DATA_T(k,  j+1)*(1.-frac_logp)*(   frac_lat)
           +DATA_T(k+1,j+1)*(   frac_logp)*(   frac_lat);
    break;
    case C_2H_2_INDEX:
      /*
       * C_2H_2
       */
      if (lat > lat_C_2H_2[nlat_C_2H_2-1]) {
        j        = nlat_C_2H_2-2;
        frac_lat = (lat-lat_C_2H_2[j])/(lat_C_2H_2[j+1]-lat_C_2H_2[j]);
      }
      else {
        for (j = 1; j < nlat_C_2H_2; j++) {
          if (lat <= lat_C_2H_2[j]) {
            j--;
            frac_lat = (lat-lat_C_2H_2[j])/(lat_C_2H_2[j+1]-lat_C_2H_2[j]);
            break;
          }
        }  
      }

      if (logp > logp_C_2H_2[np_C_2H_2-1]) {
        k         = np_C_2H_2-2;
        frac_logp = (logp-logp_C_2H_2[k])/(logp_C_2H_2[k+1]-logp_C_2H_2[k]);
      }
      else {
        for (k = 1; k < np_C_2H_2; k++) {
          if (logp <= logp_C_2H_2[k]) {
            k--;
            frac_logp = (logp-logp_C_2H_2[k])/(logp_C_2H_2[k+1]-logp_C_2H_2[k]);
            break;
          }
        }
      }
 
      ans = log(DATA_C_2H_2(k,  j  ))*(1.-frac_logp)*(1.-frac_lat)
           +log(DATA_C_2H_2(k+1,j  ))*(   frac_logp)*(1.-frac_lat)
           +log(DATA_C_2H_2(k,  j+1))*(1.-frac_logp)*(   frac_lat)
           +log(DATA_C_2H_2(k+1,j+1))*(   frac_logp)*(   frac_lat);
      ans = exp(ans);
    break;
    case C_2H_6_INDEX:
      /*
       * C_2H_6
       */
      if (lat > lat_C_2H_6[nlat_C_2H_6-1]) {
        j        = nlat_C_2H_6-2;
        frac_lat = (lat-lat_C_2H_6[j])/(lat_C_2H_6[j+1]-lat_C_2H_6[j]);
      }
      else {
        for (j = 1; j < nlat_C_2H_6; j++) {
          if (lat <= lat_C_2H_6[j]) {
            j--;
            frac_lat = (lat-lat_C_2H_6[j])/(lat_C_2H_6[j+1]-lat_C_2H_6[j]);
            break;
          }
        }  
      }

      if (logp > logp_C_2H_6[np_C_2H_6-1]) {
        k         = np_C_2H_6-2;
        frac_logp = (logp-logp_C_2H_6[k])/(logp_C_2H_6[k+1]-logp_C_2H_6[k]);
      }
      else {
        for (k = 1; k < np_C_2H_6; k++) {
          if (logp <= logp_C_2H_6[k]) {
            k--;
            frac_logp = (logp-logp_C_2H_6[k])/(logp_C_2H_6[k+1]-logp_C_2H_6[k]);
            break;
          }
        } 
      }

      ans = log(DATA_C_2H_6(k,  j  ))*(1.-frac_logp)*(1.-frac_lat)
           +log(DATA_C_2H_6(k+1,j  ))*(   frac_logp)*(1.-frac_lat)
           +log(DATA_C_2H_6(k,  j+1))*(1.-frac_logp)*(   frac_lat)
           +log(DATA_C_2H_6(k+1,j+1))*(   frac_logp)*(   frac_lat);
      ans = exp(ans);
    break;
    default:
      sprintf(Message,"index %d not yet implemented",index);
      epic_error(dbmsname,Message);
    break;
  }

  return ans;
}

/*====================== end of CIRS_data() ==================================*/

/*====================== setup_mu_p() ========================================*/

/*
 * Global declarations for mu_p().
 */
static float_triplet
  *mu_p_table;
static int
   n_mu_p;

/*
 * Assumes a preliminary call to init_with_ref() and init_species() has been made. 
 * Used to establish an approximate function relating molar mass to pressure
 * to improve init_with_u()'s gradient balance for models with condensables.
 */
void setup_mu_p(planetspec *planet)
{
  int
    K,kk,ii;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_mu_p";

  n_mu_p = 2*grid.nk+1;

  /*
   * Allocate memory.
   */
  mu_p_table = ftriplet(0,n_mu_p-1,dbmsname);

  K  = 0;
  kk = 2*K+1;
  ii = 0;
  mu_p_table[ii].x = P3(K,JLO,ILO);
  mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);
  for (K = 1; K <= KHI; K++) {
    kk = 2*K;
    ii++;
    mu_p_table[ii].x = P2(K,JLO,ILO);
    mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);

    kk++;
    ii++;
    mu_p_table[ii].x = P3(K,JLO,ILO);
    mu_p_table[ii].y = avg_molar_mass(planet,kk,JLO,ILO);
  }
  spline_pchip(n_mu_p,mu_p_table);

  return;
}

/*====================== end of setup_mu_p() =================================*/

/*====================== mu_p() ==============================================*/

/*
 * Relate molar mass to pressure to help improve init_with_u()'s 
 * gradient balance for models with condensables.
 *
 * Requires initial call to setup_mu_p().
 */

EPIC_FLOAT mu_p(EPIC_FLOAT pressure)
{
  static int
    i = -2;
  EPIC_FLOAT
    mu,p,dp;

  p  = LIMIT_RANGE(mu_p_table[0].x,pressure,mu_p_table[n_mu_p-1].x);
  i  = hunt_place_in_table(n_mu_p,mu_p_table,p,&dp,i);
  mu = splint_pchip(p,mu_p_table+i,dp);

  return mu;
}

/*====================== end of mu_p() =======================================*/

/*====================== t_yp() ==============================================*/

/*
 * Interpolates T(y,p) data table using a bicubic spline.
 */

#define T_DATA(k,j)    t_data[   j+(k)*nlat]
#define T_DATA_Y2(k,j) t_data_y2[j+(k)*nlat]

EPIC_FLOAT t_yp(EPIC_FLOAT        p, 
                EPIC_FLOAT        latr, 
                int               mode,
                init_defaultspec *def)
{
  char
    string[N_STR];
  int
    k,j;
  static int
    initialized=FALSE,
    nlat,
    npress;
  EPIC_FLOAT 
    z,temp,lat_d,
    log_p,p_d,
    a4 = 2.;
  static EPIC_FLOAT
    *lat,
    *press,
    *t_data,
    *t_data_y2,
    *t_k,
    *t_k_y2;
  static float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="t_yp";

  if (!initialized) {
    char
      header[80];
    FILE
      *infile;

    sprintf(string,EPIC_PATH"/data/%s/iris/temperature.dat",planet->name);
    infile = fopen(string,"r");
    if (!infile) {
      fprintf(stderr,"Cannot find input file %s \n",string);
      exit(1);
    }
    fscanf(infile,"%d %d",&nlat,&npress);
    fscanf(infile,"%s",header);
    /* Allocate memory: */
    lat          = fvector( 0,nlat-1,             dbmsname);
    press        = fvector( 0,npress-1,           dbmsname);
    t_k          = fvector( 0,npress-1,           dbmsname);
    t_k_y2       = fvector( 0,npress-1,           dbmsname);
    t_data       = fvector( 0,nlat*npress-1,      dbmsname);
    t_data_y2    = fvector( 0,nlat*npress-1,      dbmsname);
    buff_triplet = ftriplet(0,IMAX(nlat,npress)-1,dbmsname);

    for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",press+k);
#else
      fscanf(infile,"%f",press+k);
#endif

      /* convert press to mks */
      press[k] *= 100.;
      /* spline on log p */
      press[k] = log(press[k]);
    }
    fscanf(infile,"%s",header);
    for (j = 0; j < nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",lat+j);
#else
      fscanf(infile,"%f",lat+j);
#endif

      /* convert lat to radians */
      lat[j] *= DEG;
      for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile,"%lf",&(T_DATA(k,j)));
#else
        fscanf(infile,"%f",&(T_DATA(k,j)));
#endif

      }
    }
    fclose(infile);
    /* Compute spline information for each row: */

    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }

    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer.
       */
      for (j = 0; j < nlat; j++) {
        buff_triplet[j].y = T_DATA(k,j);
      }
      spline_pchip(nlat,buff_triplet);
      /* 
       * Unpack y2. 
       */
      for (j = 0; j < nlat; j++) {
        T_DATA_Y2(k,j) = buff_triplet[j].z;
      }
    }
    if (IAMNODE == NODE0) {
      /* Echo range limits: */
      fprintf(stderr,"Range limits: lat[%2.1f,%2.1f]  p[%2.1f,%2.1f] \n",
                     lat[0]/DEG,lat[nlat-1]/DEG,
                     exp(press[npress-1])/100.,exp(press[0])/100.);
    }

    if (def) {
      /* Set default ranges */
      /* latbot, lattop in degrees */
      def->globe_latbot = lat[0     ]/DEG;
      def->globe_lattop = lat[nlat-1]/DEG;
      /* ptop, pbot in mbar */
      def->ptop         = exp(press[0       ])/100.;
      def->pbot         = exp(press[npress-1])/100.;
    }

    initialized = TRUE;
  }
  /* End initialization. */

  if (p == 0.) return 0.;

  if (mode == 1) {
    /* Synthetic data: */
    z    = -log(p/planet->p0);
    temp = a4*z*cos(4.*latr);
  }
  else if (mode == 0) {
    /* Data */
    if (latr < lat[0]) {
      sprintf(Message,"lat = %f below data range = %f",latr/DEG,lat[0]/DEG);
      epic_error(dbmsname,Message);
    }
    else if (latr > lat[nlat-1]) {
      sprintf(Message,"lat = %f above data range = %f",latr/DEG,lat[nlat-1]/DEG);
      epic_error(dbmsname,Message);
    }
    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    j = find_place_in_table(nlat,buff_triplet,latr,&lat_d);
    /* Make column at latr */
    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer segment.
       */
      buff_triplet[j  ].y = T_DATA(   k,j  );
      buff_triplet[j+1].y = T_DATA(   k,j+1);
      buff_triplet[j  ].z = T_DATA_Y2(k,j  );
      buff_triplet[j+1].z = T_DATA_Y2(k,j+1);

      t_k[k] = splint_pchip(latr,buff_triplet+j,lat_d);
    }
    /*
     * Pack float_triplet buffer.
     */
    for (k = 0; k < npress; k++) {
      buff_triplet[k].x = press[k];
      buff_triplet[k].y = t_k[k];
    }
    spline_pchip(npress,buff_triplet);

    log_p = log(p);
    if (log_p >= press[0] && log_p <= press[npress-1]) {
      k    = find_place_in_table(npress,buff_triplet,log_p,&p_d);
      temp = splint_pchip(log_p,buff_triplet+k,p_d);
    }
    else if (log_p > press[0]) {
      sprintf(Message,"p(%5.1f) = %e out of bounds (increase floor_tp)",
                      latr/DEG,p/100.);
      epic_error(dbmsname,Message);
    }
    else {
      sprintf(Message,"p(%5.1f) = %e out of bounds (decrease ceiling_tp)",
                      latr/DEG,p/100.);
      epic_error(dbmsname,Message);
    }
  }
  else {
    sprintf(Message,"unrecognized mode");
    epic_error(dbmsname,Message);
  }

  return temp;
}

/*====================== end of t_yp() ========================================*/

/*====================== fp_yp() ==============================================*/

/*
 * Interpolates fp data table using a bicubic spline.
 * Uses data-boundary value when pressure out of range.
 */

#define FP_DATA(k,j)    fp_data[   j+(k)*nlat]
#define FP_DATA_Y2(k,j) fp_data_y2[j+(k)*nlat]

EPIC_FLOAT fp_yp(EPIC_FLOAT p, 
                 EPIC_FLOAT latr, 
                 int        mode)
{
  char
    string[N_STR];
  int
    k,j;
  static int
    initialized=FALSE,
    nlat,
    npress;
  EPIC_FLOAT 
    z,fp,lat_d,
    log_p,p_d,
    a1 =  0.5,
    a2 =  0.1,
    a3 = -0.05;
  static EPIC_FLOAT
    *lat,
    *press,
    *fp_data,
    *fp_data_y2,
    *fp_k,
    *fp_k_y2;
  static float_triplet
    *buff_triplet;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="fp_yp";

  if (!initialized) {
    char
      header[80];
    FILE
      *infile;

    sprintf(string,EPIC_PATH"/data/%s/iris/fpara.dat",planet->name);
    infile = fopen(string,"r");
    fscanf(infile,"%d %d",&nlat,&npress);
    fscanf(infile,"%s",header);
    /* 
     * Allocate memory 
     */
    lat          = fvector( 0,nlat-1,             dbmsname);
    press        = fvector( 0,npress-1,           dbmsname);
    fp_k         = fvector( 0,npress-1,           dbmsname);
    fp_k_y2      = fvector( 0,npress-1,           dbmsname);
    fp_data      = fvector( 0,nlat*npress-1,      dbmsname);
    fp_data_y2   = fvector( 0,nlat*npress-1,      dbmsname);
    buff_triplet = ftriplet(0,IMAX(nlat,npress)-1,dbmsname);

    for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",press+k);
#else
      fscanf(infile,"%f",press+k);
#endif

      /* convert press to mks */
      press[k] *= 100.;
      /* spline on log p */
      press[k] = log(press[k]);
    }
    fscanf(infile,"%s",header);
    for (j = 0; j < nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
      fscanf(infile,"%lf",lat+j);
#else
      fscanf(infile,"%f",lat+j);
#endif

      /* convert lat to radians */
      lat[j] *= DEG;
      for (k = 0; k < npress; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile,"%lf",&(FP_DATA(k,j)));
#else
        fscanf(infile,"%f",&(FP_DATA(k,j)));
#endif

      }
    }
    fclose(infile);
    /* Compute spline information for each row: */
    /*
     * Pack float_triplet buffer. 
     */
    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    for (k = 0; k < npress; k++) {
      for (j = 0; j < nlat; j++) {
        buff_triplet[j].y = FP_DATA(k,j);
      }
      spline_pchip(nlat,buff_triplet);
      for (j = 0; j < nlat; j++) {
        FP_DATA_Y2(k,j) = buff_triplet[j].z;
      }
    }
    initialized = TRUE;
  }
  /* End initialization. */

  if (mode == 1) {
    /* Synthetic data */
    z  = -log(p/planet->p0);
    fp = a1+a2*z+a3*(1.-cos(4.*latr));
  }
  else if (mode == 0) {
    /* Data */
    if (latr < lat[0]) {
      sprintf(Message,"lat = %f below data range = %f",latr/DEG,lat[0]/DEG);
      epic_error(dbmsname,Message);
    }
    else if (latr > lat[nlat-1]) {
      sprintf(Message,"lat = %f above data range = %f",latr/DEG,lat[nlat-1]/DEG);
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < nlat; j++) {
      buff_triplet[j].x = lat[j];
    }
    j = find_place_in_table(nlat,buff_triplet,latr,&lat_d);
    /* Make column at latr */
    for (k = 0; k < npress; k++) {
      /*
       * Pack float_triplet buffer segment.
       */
      buff_triplet[j  ].y = FP_DATA(   k,j  );
      buff_triplet[j+1].y = FP_DATA(   k,j+1);
      buff_triplet[j  ].z = FP_DATA_Y2(k,j  );
      buff_triplet[j+1].z = FP_DATA_Y2(k,j+1);

      fp_k[k] = splint_pchip(latr,buff_triplet+j,lat_d);
    }
    /*
     * Pack float_triplet buffer.
     */
    for (k = 0; k < npress; k++) {
      buff_triplet[k].x = press[k];
      buff_triplet[k].y = fp_k[k];
    }
    spline_pchip(npress,buff_triplet);

    log_p = log(p);
    if (log_p < press[0]) {
      log_p = press[0];
    }
    else if (log_p > press[npress-1]) {
      log_p = press[npress-1];
    }
    k  = find_place_in_table(npress,buff_triplet,log_p,&p_d);
    fp = splint_pchip(log_p,buff_triplet+k,p_d);
  }
  else {
    sprintf(Message,"unrecognized mode");
    epic_error(dbmsname,Message);
  }

  return MIN(fp,1.);

}

/*====================== end of fp_yp() ======================================*/

/*====================== p_phi() =============================================*/
/*
 * Returns pressure at the bottom of the model given phi,
 * based on p vs z data file for specified planet.
 * For gas giants or the pure-isentropic vertical coordiante,
 * just return grid.pbot;
 */

EPIC_FLOAT p_phi(planetspec *planet,
                 int         J,
                 EPIC_FLOAT  phi)
{
  int
    i,
    ntp;
  static int
    ndat,
    initialized=FALSE;
  char
    header[N_STR],
    data_file[FILE_STR];
  EPIC_FLOAT
    p,
    temperature,
    tmp,
    phi_d,
    rttop,
    gravity,ge,gp,
    spin_factor,sinlat2,coslat2;
  static EPIC_FLOAT
    pbot0,
    *pdat,
    *tdat,
    *zdat,
    *neglogpdat;
  static float_triplet
    *buff_triplet;
  FILE
    *infile;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="p_phi";

  if (!initialized) {
    if (strcmp(planet->name,"Earth") == 0) {
      /*-------* 
       * Earth *
       *-------*/

      pbot0 = 1013.25*100.;
    }
    else if (strcmp(planet->name,"Held_Suarez") == 0) {
      /*-------------* 
       * Held-Suarez *
       *-------------*/

      pbot0 = 1000.*100.;
    }
    else if (strncmp(planet->name,"Venus",5) == 0) {
      /*-------* 
       * Venus *
       *-------*/

      /*
       * This case includes all systems with a name that starts with "Venus".
       */
      sprintf(data_file,EPIC_PATH"/data/Venus/VIRA.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&temperature,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&temperature,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert mbar to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }
    else if (strncmp(planet->name,"Mars",4) == 0) {
      /*------* 
       * Mars *
       *------*/

      /*
       * This case includes all systems with a name that starts with "Mars".
       */
      sprintf(data_file,EPIC_PATH"/data/Mars/LMD.Hellas.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&temperature,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&temperature,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert hPa to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }
    else if (strcmp(planet->name,"Titan") == 0) {
      /*-------* 
       * Titan *
       *-------*/

      /*
       * Set the reference surface pressure to be 1.5 bar.
       */
      pbot0 = 1500.*100;
      sprintf(data_file,EPIC_PATH"/data/Titan/Titan_Huygens_sounding.dat");
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 6; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);

      for (i = 0; i < 3; i++) {
        fgets(header,N_STR,infile);
      }
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = ndat-1; i >= 0; i--) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %*lf %*lf",zdat+i,pdat+i);
#else
        fscanf(infile," %f %f %*f %*f",zdat+i,pdat+i);
#endif

        /* Interpolate on -log(p). */
        neglogpdat[i] = -log(pdat[i]);
      }
      fclose(infile);
    }
    else {
      /*---------* 
       * Default *
       *---------*/

      /*
       * Generic sounding.dat.
       */
      sprintf(data_file,EPIC_PATH"/data/%s/sounding.dat",planet->name);
      infile = fopen(data_file,"r");
      if (!infile) {
        sprintf(Message,"unable to read %s",data_file);
        epic_error(dbmsname,Message);
      }
      /* Skip header. */
      for (i = 0; i < 5; i++) {
        fgets(header,N_STR,infile);
      }
      /* Input number of data points. */
      fscanf(infile,"%d",&ndat);
      /*
       * Allocate memory.
       */
      pdat         = fvector( 0,ndat-1,dbmsname);
      neglogpdat   = fvector( 0,ndat-1,dbmsname);
      zdat         = fvector( 0,ndat-1,dbmsname);
      buff_triplet = ftriplet(0,ndat-1,dbmsname);
      /*
       * Input data.
       */
      for (i = 0; i < ndat; i++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
        fscanf(infile," %lf %lf %lf %lf",zdat+i,&tmp,pdat+i,&tmp);
#else
        fscanf(infile," %f %f %f %f",zdat+i,&tmp,pdat+i,&tmp);
#endif

        /* Convert km to m. */
        zdat[i] *= 1000.;
        /* Convert hPa to Pa. */
        pdat[i] *= 100.;
        /* Interpolate on log density. */
        neglogpdat[i] = -log(pdat[i]);
        /* Discard the rest of the line. */
        fgets(header,N_STR,infile);
      }
      fclose(infile);
    }

    initialized = TRUE;
  }
  /* End of initialization. */

  /*
   * For gas-giant case, just return pbot.
   */
  if (strcmp(planet->type,"gas-giant") == 0) {
    return grid.pbot;
  }

  /*
   * Calculate surface gravity.
   *
   * NOTE: grid.g[kk][jj] is not available.
   */
  spin_factor = (planet->omega_sidereal*planet->re)*(planet->omega_sidereal*planet->re)*(planet->re/planet->GM);
  ge          = (1.+1.5*planet->J2-spin_factor)*planet->GM/(planet->re*planet->re);
  gp          = (1.-3.*(planet->re/planet->rp)*(planet->re/planet->rp)*planet->J2)*planet->GM/(planet->rp*planet->rp);
  sinlat2     = sin(grid.lat[2*J+1]*DEG);
  sinlat2    *= sinlat2;
  coslat2     = cos(grid.lat[2*J+1]*DEG);
  coslat2    *= coslat2;
  gravity     = (planet->re*ge*coslat2+planet->rp*gp*sinlat2)/sqrt(planet->re*planet->re*coslat2+planet->rp*planet->rp*sinlat2);

  if (strcmp(planet->name,"Earth") == 0) {
    /*
     * Currently, PHI_SURFACE(J,I) = 0. everywhere.
     */
    if (phi == 0.) {
      p = pbot0;
    }
    else {
      sprintf(Message,"phi=%f; case phi != 0. not yet implemented for Earth\n",phi);
      epic_error(dbmsname,Message);
    }
  }
  else if (strcmp(planet->name,"Held_Suarez") == 0) {
    /*
     * PHI_SURFACE(J,I) = 0. everywhere.
     */
    if (phi == 0.) {
      p = pbot0;
    }
    else {
      sprintf(Message,"phi=%f; case phi != 0. not implemented for Held_Suarez\n",phi);
      epic_error(dbmsname,Message);
    }
  }
  else {
    if (phi < gravity*zdat[0]) {
      sprintf(Message,"phi=%g < g*zdat[0]=%g",phi,gravity*zdat[0]);
      epic_error(dbmsname,Message);
    }
    else if (phi >= gravity*zdat[ndat-1]) {
      /* Assume a constant scale height above data table. */
      rttop = grid.g[2*grid.nk+1][2*J+1]*(zdat[ndat-1]-zdat[ndat-2])/(neglogpdat[ndat-1]-neglogpdat[ndat-2]);
      p     = pdat[ndat-1]*exp(-(phi-grid.g[2*grid.nk+1][2*J+1]*zdat[ndat-1])/rttop);
    }
    else {
      /*
       * Use smooth, monotonic interpolation.
       */
      for (i = 0; i < ndat; i++) {
        buff_triplet[i].x = gravity*zdat[i];
        buff_triplet[i].y = neglogpdat[i];
      }
      spline_pchip(ndat,buff_triplet);
      i   = find_place_in_table(ndat,buff_triplet,phi,&phi_d);
      tmp = splint_pchip(phi,buff_triplet+i,phi_d);
      p   = exp(-tmp);
    }
  }

  /*
   * Screen for NaN.
   */
  if (!isfinite(p) || !isfinite(phi)) {
    sprintf(Message,"phi=%g, p=%g",phi,p);
    epic_error(dbmsname,Message);
  }

  return p;
}

/*====================== end of p_phi() =====================================*/

/*======================= openmars_var_read() ===============================*/

/*
 * Read in OpenMARS data.
 * portion = SIZE_DATA:       read in grid dimensions
 * portion = POST_SIZE_DATA:  read in constant parameters
 * portion = VAR_DATA:        read in variables
 */

#define FBUFF2D(j,i)           fbuffer[i+openmars_grid->ni*(j)]
#define FBUFF3D(k,j,i)         fbuffer[i+openmars_grid->ni*(j+openmars_grid->nj*(k))]
#define EPIC_K_OPENMARS_j(K,j) epic_k_openmars_j[j+openmars_grid->nj*(K-KLOPAD)]

void openmars_var_read(planetspec        *planet,
                       openmars_gridspec *openmars_grid,
                       char              *openmars_infile,
                       int                portion,
                       int                input_itime)
{
  int
    K,J,I,
    kk,k,j,i,
    itime,
    is,
    nelem2d,nelem3d;
  char
    expected_project[] = "UPWARDS OpenMARS reanalysis",
    start_season[32],
    end_season[32];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  int
    nc_err,nc_id,
    nc_dimid,nc_varid;
  double
    g0 = 3.71;   /* gravity standard for Mars, m/s^2 */
  int
    *ibuffer;
  float
    *fbuffer;
  double
    *dbuffer;
  size_t
    dimlen;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="openmars_var_read";

  nc_err = lookup_netcdf(openmars_infile,&nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
    epic_error(dbmsname,Message);
  }

  /*
   * Verify data file has expected content by checking its global attribute "Project".
   */
  READC(openmars_grid->project,Project,64);
  if (strcmp(openmars_grid->project,expected_project) != 0) {
    sprintf(Message,"Expecting data file's Project to be:\n%s, but it is:\n%s\n",
                    expected_project,openmars_grid->project);
    epic_error(dbmsname,Message);
  }

  /*
   * Store netcdf file history for later use.
   */
  READC(openmars_grid->infile_history,TEXT,128);

  if (portion == SIZE_DATA) {
    fprintf(stdout,"Reading SIZE_DATA from %s \n",openmars_infile);
  }
  else if (portion == POST_SIZE_DATA) {
    fprintf(stdout,"Reading POST_SIZE_DATA from %s \n",openmars_infile);
  }
  else if (portion == VAR_DATA) {
    ;
  }
  else {
    sprintf(Message,"unrecognized portion = %d",portion);
    epic_error(dbmsname,Message);
  }
  fflush(stderr);

  /*
   * Read in size of model and set parameters needed by openmars_make_arrays().
   */
  if (portion == SIZE_DATA) {
    /* 
     * Number of longitude points
     */
    nc_err = nc_inq_dimid(nc_id,"lon",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }
    openmars_grid->ni = (int)dimlen;

    /* 
     * Number of latitude points
     */
    nc_err = nc_inq_dimid(nc_id,"lat",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }
    openmars_grid->nj = (int)dimlen;

    /* 
     * Number of sigma levels.
     */
    nc_err = nc_inq_dimid(nc_id,"lev",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }
    openmars_grid->nk = (int)dimlen;

    /* 
     * Number of times.
     */
    nc_err = nc_inq_dimid(nc_id,"time",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_infile);
      epic_error(dbmsname,Message);
    }
    openmars_grid->ntime = (int)dimlen;
  }
  else if (portion == POST_SIZE_DATA) {
    /*
     * NOTE: POST_SIZE_DATA is interpreted here to mean constant parameters,
     *       including the constant 1D dimension arrays, lon, lat, etc.
     */

    /*-------------------------*
     * Read 1D lon array [deg] *
     *-------------------------*/

    fbuffer = calloc(openmars_grid->ni,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"lon",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (i = 0; i < openmars_grid->ni; i++) {
       openmars_grid->lon[i] = (EPIC_FLOAT)fbuffer[i];
    }

    free(fbuffer);

    /*
     * If anything involving a zonal FFT is to be used,
     * then grid.ni needs to be a power of 2. In that case, use:
     *   grid.ni = 1 << NINT(log(openmars_grid->ni)/log(2.));
     */
    grid.ni = openmars_grid->ni;

    /*-------------------------*
     * Read 1D lat array [deg] *
     *-------------------------*/

    fbuffer = calloc(openmars_grid->nj,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"lat",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: OpenMARS latitudes are stored in reverse index order from EPIC convention,
     *       so the ordering is changed to the EPIC convention here.
     */
    for (j = 0; j < openmars_grid->nj; j++) {
       openmars_grid->lat[j] = (EPIC_FLOAT)fbuffer[openmars_grid->nj-1-j];
    }

    free(fbuffer);

    grid.nj = NINT(180./(openmars_grid->lat[1]-openmars_grid->lat[0])-1.-sqrt(2.));

    /*------------------------------------------*
     * Read 1D sigma array, sigma = p/p_surface *
     *------------------------------------------*/

    fbuffer = calloc(openmars_grid->nk,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"lev",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: OpenMARS sigma are stored in reverse index order from EPIC convention,
     *       so the ordering is changed to the EPIC convention upon input.
     */
    for (k = 0; k < openmars_grid->nk; k++) {
      openmars_grid->sigma[k] = (EPIC_FLOAT)fbuffer[openmars_grid->nk-1-k];
    }

    free(fbuffer);

    /*
     * OpenMARS and EMARS have 35 and 28 levels, with different
     * placements, and use sigma and hybrid sigma-pressure coordinates,
     * respectively. Despite these differences, the vertical coverage is similar
     * enough for us to map them to the same potential temperature, theta, surfaces.
     *
     * NOTE: We hardwire grid.nk to yield rounded theta levels. This needs to be adjusted for each new application.
     */
    grid.nk = 29;

    /*---------------------------*
     * Read 1D time array [sols] *
     *---------------------------*/

    fbuffer = calloc(openmars_grid->ntime,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"time",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    for (itime = 0; itime < openmars_grid->ntime; itime++) {
       openmars_grid->time[itime] = (EPIC_FLOAT)fbuffer[itime];
    }
    free(fbuffer);

    /*-----------------------------------------*
     * Read 1D Ls array, Solar Longitude [deg] *
     *-----------------------------------------*/

    fbuffer = calloc(openmars_grid->ntime,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"Ls",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Montabone et al (2014, Geosci. Data J. 1, p. 134) suggest
     *       subtracting 0.12deg from Ls to correct for a near-constant bias,
     *       which we implement here (Ls matches precisely in the MACDA and OpenMARS data sets).
     */
    for (itime = 0; itime < openmars_grid->ntime; itime++) {
      openmars_grid->Ls[itime] = (EPIC_FLOAT)fbuffer[itime]-0.12;
      openmars_grid->Ls[itime] = (openmars_grid->Ls[itime] < 0.) ? openmars_grid->Ls[itime]+360. : openmars_grid->Ls[itime];
    }

    season_string(openmars_grid->Ls[0                     ],start_season);
    season_string(openmars_grid->Ls[openmars_grid->ntime-1],  end_season);
    if (strcmp(start_season,end_season) == 0) {
      fprintf(stdout,"  Input data span L_s = %6.2f to %6.2f (%s)\n",
                     openmars_grid->Ls[0],openmars_grid->Ls[openmars_grid->ntime-1],start_season);
    }
    else {
      fprintf(stdout,"   Input data span L_s = %6.2f to %6.2f (%s to %s)\n",
                     openmars_grid->Ls[0],openmars_grid->Ls[openmars_grid->ntime-1],start_season,end_season);
    }

    free(fbuffer);

    /*----------------------------------------------------*
     * Read 1D MY array, Martian Year (Clancy et al 2000) *
     *----------------------------------------------------*/

    ibuffer = calloc(openmars_grid->ntime,sizeof(int));
    if (!ibuffer) {
      sprintf(Message,"calloc error allocating ibuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"MY",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_int(nc_id,nc_varid,ibuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (itime = 0; itime < openmars_grid->ntime; itime++) {
      openmars_grid->MY[itime] = ibuffer[itime];
    }
    free(ibuffer);
  }
  else if (portion == VAR_DATA) {
    int
      count,
      excluded_count;
    EPIC_FLOAT
      p,tmp,x,x_d,
     *epic_k_openmars_j,
     *openmars_phi_surface,
     *openmars_phi;
    size_t
      nc_start2d[3],nc_start3d[4],
      nc_count2d[3],nc_count3d[4];
    float_triplet
      *k_triplet,
      *j_triplet,
      *i_triplet,
      *epic_j_triplet;
    FILE
      *excluded_file;

    /*
     * NOTE: VAR_DATA is interpreted here to mean one timeframe, as specified
     *       by the input argument input_itime, of the variable fields, like ps,
     *       temp, etc.
     */

    nelem2d       = openmars_grid->nj*openmars_grid->ni;
    nc_start2d[0] = input_itime;
    nc_start2d[1] = 0;
    nc_start2d[2] = 0;
    nc_count2d[0] = 1;
    nc_count2d[1] = openmars_grid->nj;
    nc_count2d[2] = openmars_grid->ni;

    /*--------------------------------------*
     * Read ps array, surface pressure [Pa] *
     *--------------------------------------*/

    fbuffer = calloc(nelem2d,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"ps",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start2d,nc_count2d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < openmars_grid->nj; j++) {
      for (i = 0; i < openmars_grid->ni; i++) {
        OPENMARS_PS(j,i) = (EPIC_FLOAT)FBUFF2D(openmars_grid->nj-1-j,i);
      }
    }

    /*-------------------------------------------*
     * Read tsurf array, surface temperature [K] *
     *-------------------------------------------*/

    nc_err = nc_inq_varid(nc_id,"tsurf",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start2d,nc_count2d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < openmars_grid->nj; j++) {
      for (i = 0; i < openmars_grid->ni; i++) {
        OPENMARS_TSURF(j,i) = (EPIC_FLOAT)FBUFF2D(openmars_grid->nj-1-j,i);
      }
    }

    free(fbuffer);

    /*-------------------------------------------*
     * Read temp array, 3D temperature field [K] *
     *-------------------------------------------*/

    nelem3d       = nelem2d*openmars_grid->nk;
    nc_start3d[0] = input_itime;
    nc_start3d[1] = 0;
    nc_start3d[2] = 0;
    nc_start3d[3] = 0;
    nc_count3d[0] = 1;
    nc_count3d[1] = openmars_grid->nk;
    nc_count3d[2] = openmars_grid->nj;
    nc_count3d[3] = openmars_grid->ni;

    fbuffer = calloc(nelem3d,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"temp",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < openmars_grid->nk; k++) {
      for (j = 0; j < openmars_grid->nj; j++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          OPENMARS_TEMP(k,j,i) = (EPIC_FLOAT)FBUFF3D(openmars_grid->nk-1-k,openmars_grid->nj-1-j,i);
        }
      }
    }

    /*-----------------------------------------*
     * Read u array, 3D zonal wind field [m/s] *
     *-----------------------------------------*/

    nc_err = nc_inq_varid(nc_id,"u",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < openmars_grid->nk; k++) {
      for (j = 0; j < openmars_grid->nj; j++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          OPENMARS_U(k,j,i) = (EPIC_FLOAT)FBUFF3D(openmars_grid->nk-1-k,openmars_grid->nj-1-j,i);
        }
      }
    }

    /*----------------------------------------------*
     * Read v array, 3D meridional wind field [m/s] *
     *----------------------------------------------*/

    nc_err = nc_inq_varid(nc_id,"v",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < openmars_grid->nk; k++) {
      for (j = 0; j < openmars_grid->nj; j++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          OPENMARS_V(k,j,i) = (EPIC_FLOAT)FBUFF3D(openmars_grid->nk-1-k,openmars_grid->nj-1-j,i);
        }
      }
    }

    free(fbuffer);

    /*---------------------------------------------------------*
     * Compute theta array, 3D potential temperature field [K] *
     *                                                         *
     * NOTE: For Mars, planet->p0 = 610 Pa                     *
     *---------------------------------------------------------*/

    for (k = 0; k < openmars_grid->nk; k++) {
      for (j = 0; j < openmars_grid->nj; j++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          p                     = openmars_grid->sigma[k]*OPENMARS_PS(j,i);
          OPENMARS_THETA(k,j,i) = OPENMARS_TEMP(k,j,i)*pow(planet->p0/p,planet->kappa);
        }
      }
    }

    /*------------------------------------------------------*
     * Interpolate OpenMARS data onto EPIC staggered grids. *
     *------------------------------------------------------*/

    /* Allocate memory */
    k_triplet            = ftriplet(0,openmars_grid->nk-1,               dbmsname);
    j_triplet            = ftriplet(0,openmars_grid->nj+1,               dbmsname);
    i_triplet            = ftriplet(0,openmars_grid->ni,                 dbmsname);
    epic_j_triplet       = ftriplet(0,grid.nj,                           dbmsname);
    epic_k_openmars_j    = fvector( 0,(KHI-KLOPAD+1)*openmars_grid->nj-1,dbmsname);
    openmars_phi_surface = fvector( 0,openmars_grid->nj-1,               dbmsname);
    openmars_phi         = fvector( 0,openmars_grid->nk-1,               dbmsname);

    /*--------------------------------------------*
     * Zonal wind, U                              *
     * Longitude can be done by direct insertion. * 
     *--------------------------------------------*/

    excluded_file = fopen("excluded_thetas.dat","w");
    if (!excluded_file) {
      sprintf(Message,"file %s not found","excluded_thetas.dat");
      epic_error(dbmsname,Message);
    }
    fprintf(excluded_file," Points in column excluded because of lack of monotonicity of theta.\n");
    fprintf(excluded_file,"  Lon    Lat  Excluded thetas [K]\n");

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping OpenMARS latitudes.
       */
      for (j = 0; j < openmars_grid->nj; j++) {
        k              = openmars_grid->nk-1;
        k_triplet[0].x = OPENMARS_THETA(k,j,i);
        k_triplet[0].y = OPENMARS_U(k,j,i);
        count          = 1;
        excluded_count = 0;
        for (k = openmars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (OPENMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = OPENMARS_THETA(k,j,i);
            k_triplet[count].y = OPENMARS_U(k,j,i);
            count++;
          }
          else {
            /*
             * Keep track of number of excluded points caused by non-monotonicity of theta.
             */
            excluded_count++;
            if (excluded_count == 1) {
              fprintf(excluded_file," %6.1f %5.1f  %5.1f",openmars_grid->lon[i],openmars_grid->lat[j],OPENMARS_THETA(k,j,i));
            }
            else {
              fprintf(excluded_file," %5.1f",OPENMARS_THETA(k,j,i));
            }
          }
        }

        if (excluded_count > 0) {
          fprintf(excluded_file,"\n");
        }

        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLO; K <= KHI; K++) {
          /*
           * U resides in the layer.
           */
          x = grid.sigmatheta[2*K];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_OPENMARS_j(K,j) = OPENMARS_U(openmars_grid->nk-1,j,i);
          }
          else {
            k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_OPENMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       * We add U = 0 pole points.
       */
      for (K = KLO; K <= KHI; K++) {
        j_triplet[0].x = -90.;
        j_triplet[0].y =   0.;
        for (j = 0; j < openmars_grid->nj; j++) {
          j_triplet[j+1].x = openmars_grid->lat[j];
          j_triplet[j+1].y = EPIC_K_OPENMARS_j(K,j);
        }
        j_triplet[j+1].x = 90.;
        j_triplet[j+1].y =  0.;
        spline_pchip(openmars_grid->nj+2,j_triplet);
        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x                   = grid.lat[2*J+1];
          j                   = hunt_place_in_table(openmars_grid->nj+2,j_triplet,x,&x_d,j);
          U(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

    fclose(excluded_file);

    /*-------------------------------------------------------------------*
     * Meridonal wind, V                                                 *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping OpenMARS latitudes.
       */
      for (j = 0; j < openmars_grid->nj; j++) {
        k              = openmars_grid->nk-1;
        k_triplet[0].x = OPENMARS_THETA(k,j,i);
        k_triplet[0].y = OPENMARS_V(k,j,i);
        count          = 1;
        for (k = openmars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (OPENMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = OPENMARS_THETA(k,j,i);
            k_triplet[count].y = OPENMARS_V(k,j,i);
            count++;
          }
        }
        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLO; K <= KHI; K++) {
          x = grid.sigmatheta[2*K];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_OPENMARS_j(K,j) = OPENMARS_V(openmars_grid->nk-1,j,i);
          }
          else {
            k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_OPENMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       * We add V = 0 pole points.
       */
      for (K = KLO; K <= KHI; K++) {
        j_triplet[0].x = -90.;
        j_triplet[0].y =   0.;
        for (j = 0; j < openmars_grid->nj; j++) {
          j_triplet[j+1].x = openmars_grid->lat[j];
          j_triplet[j+1].y = EPIC_K_OPENMARS_j(K,j);
        }
        j_triplet[j+1].x = 90.;
        j_triplet[j+1].y =  0.;
        spline_pchip(openmars_grid->nj+2,j_triplet);
        j = -2;
        for (J = JFIRST; J <= JHI; J++) {
          x                   = grid.lat[2*J];
          j                   = hunt_place_in_table(openmars_grid->nj+2,j_triplet,x,&x_d,j);
          V(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JFIRST; J <= JHI; J++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = openmars_grid->lon[i];
          i_triplet[i].y = V(grid.it_uv,K,J,I);
        }
        i_triplet[i].x = openmars_grid->lon[i-1]+openmars_grid->lon[1]-openmars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(openmars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x                   = grid.lon[2*I+1];
          i                   = hunt_place_in_table(openmars_grid->ni+1,i_triplet,x,&x_d,i);
          V(grid.it_uv,K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

    /*-------------------------------------------------------------------*
     * Pressure on the layer interfaces, P3.                             *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping OpenMARS latitudes.
       */
      for (j = 0; j < openmars_grid->nj; j++) {
        k              = openmars_grid->nk-1;
        k_triplet[0].x = OPENMARS_THETA(k,j,i);
        k_triplet[0].y = openmars_grid->sigma[k]*OPENMARS_PS(j,i);
        count          = 1;
        for (k = openmars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (OPENMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = OPENMARS_THETA(k,j,i);
            k_triplet[count].y = openmars_grid->sigma[k]*OPENMARS_PS(j,i);
            count++;
          }
        }
        spline_pchip(count,k_triplet);

        k = -2;
        for (K = KLOPAD; K <= KHI; K++) {
          x = grid.sigmatheta[2*K+1];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_OPENMARS_j(K,j) = OPENMARS_PS(j,i);
          }
          else {
            k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_OPENMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      for (K = KLOPAD; K <= KHI; K++) {
        for (j = 0; j < openmars_grid->nj; j++) {
          j_triplet[j].x = openmars_grid->lat[j];
          j_triplet[j].y = EPIC_K_OPENMARS_j(K,j);
        }
        spline_pchip(openmars_grid->nj,j_triplet);
        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x         = grid.lat[2*J+1];
          j         = hunt_place_in_table(openmars_grid->nj,j_triplet,x,&x_d,j);
          P3(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLOPAD; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = openmars_grid->lon[i];
          i_triplet[i].y = P3(K,J,I);
        }
        i_triplet[i].x = openmars_grid->lon[i-1]+openmars_grid->lon[1]-openmars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(openmars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x         = grid.lon[2*I+1];
          i         = hunt_place_in_table(openmars_grid->ni+1,i_triplet,x,&x_d,i);
          P3(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */ 
    bc_lateral(var.p3.value,THREEDIM);

    /*---------------------------------------------------------------------*
     * Potential temperature, THETA.                                       *
     * Using isentropic coordinates, so THETA takes the coordinate value,  *
     * except when the layer is underground, in which case the surface     *
     * value is used.                                                      *
     *---------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      for (j = 0; j < openmars_grid->nj; j++) {
        for (K = KLOPAD; K <= KHI; K++) {
          /*
           * NOTE: Using the bottom position of EPIC layers, 2*K+1 rather than 2*K, 
           *       to make it easier to specify the output THETA values as whole numbers
           *       (we are not setting up variables for an EPIC model run).
           */
          EPIC_K_OPENMARS_j(K,j) = MAX(grid.sigmatheta[2*K+1],OPENMARS_THETA(openmars_grid->nk-1,j,i));
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      for (K = KLOPAD; K <= KHI; K++) {
        for (j = 0; j < openmars_grid->nj; j++) {
          j_triplet[j].x = openmars_grid->lat[j];
          j_triplet[j].y = EPIC_K_OPENMARS_j(K,j);
        }
        spline_pchip(openmars_grid->nj,j_triplet);

        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x            = grid.lat[2*J+1];
          j            = hunt_place_in_table(openmars_grid->nj,j_triplet,x,&x_d,j);
          THETA(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLOPAD; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (i = 0; i < openmars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = openmars_grid->lon[i];
          i_triplet[i].y = THETA(K,J,I);
        }
        i_triplet[i].x = openmars_grid->lon[i-1]+openmars_grid->lon[1]-openmars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(openmars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x            = grid.lon[2*I+1];
          i            = hunt_place_in_table(openmars_grid->ni+1,i_triplet,x,&x_d,i);
          THETA(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */ 
    bc_lateral(var.theta.value,THREEDIM);


    /*--------------------------------------------------------------------*
     * Geopotential at the bottom of the model, PHI3(KHI,J,I). This is    *
     * the geopotential on the bottom isentropic surface, grid.thetabot,  *
     * not on the planet's surface, PHI_SURFACE(J,I).                     *
     * Longitude needs to be staggered, but start with direct insertion.  *
     *--------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * Interpolate EPIC's PHI_SURFACE(J,I) for Mars, which is calculated
       * on the planet's surface, onto the OpenMARS latitude grid.
       *
       * NOTE: The southernmost and northermost point assigned to openmars_phi_surface[j]
       *       are each an extrapolation beyond the range of the EPIC grid.  This issue
       *       could be alleviated by importing more Mars topography data to cover the
       *       poles.  However, the output itself is confined to the EPIC grid, so
       *       the issue is minor.
       */
      for (J = JLO; J <= JHI; J++) {
        epic_j_triplet[J].x = grid.lat[2*J+1];
        epic_j_triplet[J].y = PHI_SURFACE(J,I);
      }
      spline_pchip(JHI-JLO+1,epic_j_triplet);
      J = -2;
      for (j = 0; j < openmars_grid->nj; j++) {
        x                       = openmars_grid->lat[j];
        J                       = hunt_place_in_table(JHI-JLO+1,epic_j_triplet,x,&x_d,J);
        openmars_phi_surface[j] = splint_pchip(x,epic_j_triplet+J,x_d);
      }

      for (j = 0; j < openmars_grid->nj; j++) {
        /*
         * Calculate the OpenMARS geopotential vertically by integrating
         * the hydrostatic balance equation up from the surface.
         */
        openmars_phi[openmars_grid->nk-1] = openmars_phi_surface[j];
        for (k = openmars_grid->nk-2; k >= 0; k--) {
          /*
           *  Hydrostatic integration algorithm (K increases towards greater pressure): 
           *    phi[k] = phi[k+1]+(p[k+1]-p[k])*.5*(1/rho[k]+1/rho[k+1]),
           *  with the ideal-gas law, 1/rho = rgas*T/p.
           */
          openmars_phi[k] = openmars_phi[k+1]+(openmars_grid->sigma[k+1]-openmars_grid->sigma[k])
                                       *.5*planet->rgas*(OPENMARS_TEMP(k,  j,i)/openmars_grid->sigma[k  ]
                                                        +OPENMARS_TEMP(k+1,j,i)/openmars_grid->sigma[k+1]);
        }

        /*
         * Next, interpolate onto EPIC theta surfaces, keeping OpenMARS latitudes.
         */
        k              = openmars_grid->nk-1;
        k_triplet[0].x = OPENMARS_THETA(k,j,i);
        k_triplet[0].y = openmars_phi[k];
        count          = 1;
        for (k = openmars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (OPENMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = OPENMARS_THETA(k,j,i);
            k_triplet[count].y = openmars_phi[k];
            count++;
          }
        }
        spline_pchip(count,k_triplet);

        /* Only need the geopotential at grid.thetabot. */
        K = KHI;
        x = grid.sigmatheta[2*K+1];

        if (x < k_triplet[0].x) {
          /* 
           * Layer is below the ground, so place it on the ground.
           */
          EPIC_K_OPENMARS_j(K,j) = openmars_phi[openmars_grid->nk];
        }
        else {
          k                      = -2; 
          k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
          EPIC_K_OPENMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
        }
      } /* j loop */

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      K = KHI;
      for (j = 0; j < openmars_grid->nj; j++) {
        j_triplet[j].x = openmars_grid->lat[j];
        j_triplet[j].y = EPIC_K_OPENMARS_j(K,j);
      }
      spline_pchip(openmars_grid->nj,j_triplet);
      j = -2;
      for (J = JLO; J <= JHI; J++) {
        x           = grid.lat[2*J+1];
        j           = hunt_place_in_table(openmars_grid->nj,j_triplet,x,&x_d,j);
        PHI3(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
      }
    }  /* I loop */

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    K = KHI;
    for (J = JLO; J <= JHI; J++) {
      for (i = 0; i < openmars_grid->ni; i++) {
        I              = i+1;
        i_triplet[i].x = openmars_grid->lon[i];
        i_triplet[i].y = PHI3(K,J,I);
      }
      i_triplet[i].x = openmars_grid->lon[i-1]+openmars_grid->lon[1]-openmars_grid->lon[0];
      i_triplet[i].y = i_triplet[0].y;
      periodic_spline_pchip(openmars_grid->ni+1,i_triplet);
      i = -2;
      for (I = ILO; I <= IHI; I++) {
        x           = grid.lon[2*I+1];
        i           = hunt_place_in_table(openmars_grid->ni+1,i_triplet,x,&x_d,i);
        PHI3(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.phi3.value,THREEDIM);

    /* Free allocated memory */
    free_ftriplet(k_triplet,          0,openmars_grid->nk-1,               dbmsname);
    free_ftriplet(j_triplet,          0,openmars_grid->nj-1,               dbmsname);
    free_ftriplet(i_triplet,          0,openmars_grid->ni,                 dbmsname);
    free_ftriplet(epic_j_triplet,     0,grid.nj,                           dbmsname);
    free_fvector(epic_k_openmars_j,   0,(KHI-KLOPAD+1)*openmars_grid->nj-1,dbmsname);
    free_fvector(openmars_phi_surface,0,openmars_grid->nj-1,               dbmsname);
    free_fvector(openmars_phi,        0,openmars_grid->nk-1,               dbmsname);

  } /* end portion == VAR_DATA */

  return;
}

#undef FBUFF2D
#undef FBUFF3D
#undef EPIC_K_OPENMARS_j

/*======================= end of openmars_var_read() ========================*/

/*======================= openmars_make_arrays() ============================*/

void openmars_make_arrays(planetspec        *planet,
                          openmars_gridspec *openmars_grid)
{
  int
    nelem2d,nelem3d;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="openmars_make_arrays";

  nelem2d = openmars_grid->nj*openmars_grid->ni;
  nelem3d = nelem2d*openmars_grid->nk;

  openmars_grid->lon       = fvector(0,openmars_grid->ni-1,   dbmsname);
  openmars_grid->lat       = fvector(0,openmars_grid->nj-1,   dbmsname);
  openmars_grid->sigma     = fvector(0,openmars_grid->nk-1,   dbmsname);
  openmars_grid->time      = fvector(0,openmars_grid->ntime-1,dbmsname);
  openmars_grid->Ls        = fvector(0,openmars_grid->ntime-1,dbmsname);
  openmars_grid->MY        = ivector(0,openmars_grid->ntime-1,dbmsname);
  openmars_grid->ps        = fvector(0,nelem2d-1,             dbmsname);
  openmars_grid->tsurf     = fvector(0,nelem2d-1,             dbmsname);
  openmars_grid->temp      = fvector(0,nelem3d-1,             dbmsname);
  openmars_grid->u         = fvector(0,nelem3d-1,             dbmsname);
  openmars_grid->v         = fvector(0,nelem3d-1,             dbmsname);
  openmars_grid->theta     = fvector(0,nelem3d-1,             dbmsname);

  return;
}

/*======================= end of openmars_make_arrays() =====================*/

/*======================= openmars_free_arrays() ============================*/

void openmars_free_arrays(planetspec        *planet,
                          openmars_gridspec *openmars_grid)
{
  int
    nelem2d,nelem3d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="openmars_free_arrays";

  nelem2d = openmars_grid->ni*openmars_grid->nj;
  nelem3d = nelem2d*openmars_grid->nk;

  free_fvector(openmars_grid->lon,        0,openmars_grid->ni-1,   dbmsname);
  free_fvector(openmars_grid->lat,        0,openmars_grid->nj-1,   dbmsname);
  free_fvector(openmars_grid->sigma,      0,openmars_grid->nk-1,   dbmsname);
  free_fvector(openmars_grid->time,       0,openmars_grid->ntime-1,dbmsname);
  free_fvector(openmars_grid->Ls,         0,openmars_grid->ntime-1,dbmsname);
  free_ivector(openmars_grid->MY,         0,openmars_grid->ntime-1,dbmsname);
  free_fvector(openmars_grid->ps,         0,nelem2d-1,             dbmsname);
  free_fvector(openmars_grid->tsurf,      0,nelem2d-1,             dbmsname);
  free_fvector(openmars_grid->temp,       0,nelem3d-1,             dbmsname);
  free_fvector(openmars_grid->u,          0,nelem3d-1,             dbmsname);
  free_fvector(openmars_grid->v,          0,nelem3d-1,             dbmsname);
  free_fvector(openmars_grid->theta,      0,nelem3d-1,             dbmsname);

  return;
}

/*======================= end of openmars_free_arrays() =====================*/

/*======================= openmars_epic_nc() ================================*/

/*
 * Create openmars_epic.nc by running initial.
 *
 * The input parameters are "hard wired" and need to be adjusted for each new application.
 */

void openmars_epic_nc(planetspec *planet)
{
  pid_t
    pid;
  int
    commpipe[2],
    ierr,status,
    saved_stdout;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="openmars_epic_nc";

  /*
   * The code for executing initial from change is adapted from examples at following websites:
   *    https://www.gidforums.com/t-3369.html
   *    https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/forkexec.html
   *    http://www.unix.com/programming/173811-c-execl-pipes.html
   *    http://stackoverflow.com/questions/11042218/c-restore-stdout-to-terminal
   *
   * Set up a communications pipeline between the parent and child processes, 
   * which here are change and initial, respectively.
   */
  ierr = pipe(commpipe);
  if (ierr) {
    perror("pipe");
    sprintf(Message,"error calling pipe()");
    epic_error(dbmsname,Message);
  }

  /*
   * Run initial by forking a child process.
   */
  switch (pid = fork()) {
    case -1:
      sprintf(Message,"error calling fork()");
      epic_error(dbmsname,Message);
    break;
    case 0:
      /*
       * The child process.
       * Connect its stdin to the pipe, then specify the process to be initial.
       */
      dup2(commpipe[0],STDIN_FILENO);  /* set stdin                           */
      close(commpipe[1]);              /* close unused end of pipe            */
      /*
       * Specify the child to be initial.
       */
      ierr = execl(EPIC_PATH"/bin/initial",EPIC_PATH"/bin/initial",NULL);
      if (ierr) {
        perror("execl");
        sprintf(Message,"error calling execl()");
        epic_error(dbmsname,Message);
      }
    break;
    default:
      /*
       * The parent process, which is change.
       */
      saved_stdout = dup(STDOUT_FILENO);
      dup2(commpipe[1],STDOUT_FILENO);       /* set stdout                         */
      close(commpipe[0]);                    /* close unused end of pipe           */
      setvbuf(stdout,(char *)NULL,_IOLBF,0); /* set line-buffered output on stdout */
      /*
       * Send inputs to initial.
       *
       * NOTE: These are not 'smart' and hence will need to be updated if the
       *       inputs to initial are modified.
       */
      printf("Mars\n");
      printf("%d\n",COORD_ISENTROPIC);
      printf("2\n");     /*prompt for starting date values */
      /*
       * Converter between Earth and Mars calendar time:
       * http://www-mars.lmd.jussieu.fr/mars/time/mars_date_to_earth_date.html
       *
       * OpenMARS data are referred to Mars Year 24, L_s = 0, which according to the 
       * above website is Julian Date 2451009.27883, give or take.  This corresponds
       * to late in the day UTC on July 14, 1998; the value below is a best guess
       * at this instant in time---an official precise value is elusive
       * to track down.
       */ 
      printf("1998\n");  /* year             */
      printf("7\n");     /* month            */
      printf("14\n");    /* day              */
      printf("18\n");    /* hour             */
      printf("41\n");    /* minute           */
      printf("31\n");    /* second           */

      printf("1\n");                /* initial-wind scaling factor            */
      printf("0\n");                /* radiation scheme                       */
      printf("0\n");                /* turbulence scheme off                  */
      printf("-1\n");               /* sponge off                             */
      printf("%d\n",SPACING_THETA);  /* layer spacing even in potential temp. */

      /*
       * Adjust theta at the model top and bottom for each new application.
       */
      printf("410.0\n");            /* theta at model top [K]                 */
      printf("120.0\n");            /* theta at model bottom [K]              */
      printf("%d\n",grid.nk);       /* number of vertical layers              */
      printf("globe\n");            /* geometry                               */
      printf("-90\n");              /* latbot                                 */
      printf("90\n");               /* lattop                                 */
      printf("-180\n");             /* lonbot                                 */
      printf("180\n");              /* lontop                                 */
      printf("none\n");             /* optional prognostic variables          */
      printf("%d\n",grid.nj);       /* number of latitude gridpoints          */
      printf("%d\n",grid.ni);       /* number of longitude gridpoints         */
      printf("\n");                 /* timestep, use default                  */
      /*
       * NOTE: numerical damping must be turned off if grid.ni is not a
       *       power of 2, which it is not in the OpenMARS dataset.
       */
      printf("0.\n");               /* divergence damping                     */
      printf("0\n");                /* hyperviscosity                         */
      printf("none\n");             /* extract variables set in initial       */

      wait(&status);                /* wait for initial to end                */

      /*
       * Restore stdout for change.
       */
      dup2(saved_stdout,STDOUT_FILENO);
      close(saved_stdout);
      close(commpipe[1]);
    break;
  }
     
  system("mv ./epic.nc openmars_epic.nc");
}

/*======================= end of openmars_epic_nc() ==========================*/

/*======================= openmars_conversion() ==============================*/

/*
 * Input OpenMARS data into EPIC, compute diagnostic fields on isentropic surfaces,
 * and output results.
 */
void openmars_conversion(planetspec        *planet,
                         openmars_gridspec *openmars_grid,
                         char              *openmars_infile,
                         char              *openmars_outfile_qb,
                         char              *openmars_outfile_uvpt,
                         EPIC_FLOAT       **Buff2D)
{
  int
    jj,kk,i,
    K,J,I,
    openmars_itime,
    nc_err,
    nc_id_qb,
    nc_id_uvpt,
    dimid[4],
    coorid[4];
  double
    g,
    c2factor,
   *buffer;
  register double
    tmp;
  time_t
    now;
  struct tm
   *today;
  size_t
    nc_index[1],
    nc_start[4], 
    nc_count[4];
  char
    date[16];
  id_information
    L_s_info,
    MY_info,
    pv_info,bsf_info,NH2_info,kin_info,
    u_info,ug_info,v_info,p_info,t_info,phi_info,ma_info;
  EPIC_FLOAT
   *p,
   *h;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="openmars_conversion";

  /*
   * Generate current-date string.
   */
  time(&now);
  today = localtime(&now);
  strftime(date,12,"%Y-%m-%d",today);

  /*
   * Create OpenMARS output .nc files.
   */
  sprintf(openmars_outfile_qb,  "./ooQB");
  strcat(openmars_outfile_qb,strstr(openmars_infile,"_my"));
  nc_err = nc_create(openmars_outfile_qb,NC_CLOBBER,&nc_id_qb);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_outfile_qb);
    epic_error(dbmsname,Message);
  }

  sprintf(openmars_outfile_uvpt,  "./ooUVPT");
  strcat(openmars_outfile_uvpt,strstr(openmars_infile,"_my"));
  nc_err = nc_create(openmars_outfile_uvpt,NC_CLOBBER,&nc_id_uvpt);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),openmars_outfile_uvpt);
    epic_error(dbmsname,Message);
  }

  /*
   * Write global attributes to OpenMARS outfiles.
   */

  /* QB */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model calculations of Q vs. B on Theta Surfaces; input data are OpenMARS U, V, P and T");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s; \n"
                  "%s: Q, B and NH on isentropic surfaces uses the EPIC GCM v%4.2f",
                  openmars_grid->infile_history,date,grid.epic_version);
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"OpenMARS: PI Stephen R. Lewis (s.r.lewis@open.ac.uk); "
                   "EPIC: PI Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains potential vorticity, Q, Bernoulli Streamfunction, B, "
                  "and the product of buoyancy frequency times pressure scale height, NH, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Mars.");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /* UVPT */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model re-gridding of U, V, P and T on Theta Surfaces; input data are OpenMARS U, V, P and T on Sigma Levels");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s; \n"
                  "%s: U, V, P and T on isentropic surfaces uses the EPIC GCM v%4.2f",
                  openmars_grid->infile_history,date,grid.epic_version);
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"OpenMARS: PI Stephen R. Lewis (s.r.lewis@open.ac.uk); "
                   "EPIC: PI Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains zonal and meridional winds, u and v, pressure, p, "
                  "temperature, T, and geopotential, phi, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Mars.");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /*
   * Define coordinates and variables.
   */

  /*
   * lon (I direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_qb,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_uvpt,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /*
   * lat (J direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_qb,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_uvpt,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /*
   * theta (K direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_qb,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_uvpt,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /*
   * time
   */

  /* QB */
  nc_def_dim(nc_id_qb,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_qb,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"long_name",strlen("Time")+1,"Time");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"units",strlen("days since 0000-00-0 00:00:00")+1,"days since 0000-00-0 00:00:00");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_uvpt,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"long_name",strlen("Time")+1,"Time");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"units",strlen("days since 0000-00-0 00:00:00")+1,"days since 0000-00-0 00:00:00");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");

  /*
   * L_s
   */
  L_s_info.name = (char *)malloc(strlen("L_s")+1);
  sprintf(L_s_info.name,"%s","L_s");
  L_s_info.standard_name = (char *)malloc(strlen("solar_longitude")+1);
  sprintf(L_s_info.standard_name,"%s","solar_longitude");
  L_s_info.long_name = (char *)malloc(strlen("Solar longitude")+1);
  sprintf(L_s_info.long_name,"%s","Solar longitude");
  L_s_info.units = (char *)malloc(strlen("degrees")+1);
  sprintf(L_s_info.units,"%s","degrees");
  L_s_info.dimid[NETCDF_T_INDEX] = dimid[NETCDF_T_INDEX];

  /* QB */
  nc_def_var(nc_id_qb,L_s_info.name,NC_DOUBLE,1,&L_s_info.dimid[NETCDF_T_INDEX],&L_s_info.id);
  nc_put_att_text(nc_id_qb,L_s_info.id,"standard_name",strlen(L_s_info.standard_name)+1,L_s_info.standard_name);
  nc_put_att_text(nc_id_qb,L_s_info.id,"long_name",strlen(L_s_info.long_name)+1,L_s_info.long_name);
  nc_put_att_text(nc_id_qb,L_s_info.id,"units",strlen(L_s_info.units)+1,L_s_info.units);

  /* UVPT */
  nc_def_var(nc_id_uvpt,L_s_info.name,NC_DOUBLE,1,&L_s_info.dimid[NETCDF_T_INDEX],&L_s_info.id);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"standard_name",strlen(L_s_info.standard_name)+1,L_s_info.standard_name);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"long_name",strlen(L_s_info.long_name)+1,L_s_info.long_name);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"units",strlen(L_s_info.units)+1,L_s_info.units);

  /*
   * MY
   */
  MY_info.name = (char *)malloc(strlen("MY")+1);
  sprintf(MY_info.name,"%s","MY");
  MY_info.long_name = (char *)malloc(strlen("Martian year")+1);
  sprintf(MY_info.long_name,"%s","Martian year");
  MY_info.units = (char *)malloc(strlen("1")+1);
  sprintf(MY_info.units,"%s","1");
  MY_info.dimid[NETCDF_T_INDEX] = dimid[NETCDF_T_INDEX];

  /* QB */
  nc_def_var(nc_id_qb,MY_info.name,NC_SHORT,1,&MY_info.dimid[NETCDF_T_INDEX],&MY_info.id);
  nc_put_att_text(nc_id_qb,MY_info.id,"long_name",strlen(MY_info.long_name)+1,MY_info.long_name);
  nc_put_att_text(nc_id_qb,MY_info.id,"units",strlen(MY_info.units)+1,MY_info.units);

  /* UVPT */
  nc_def_var(nc_id_uvpt,MY_info.name,NC_SHORT,1,&MY_info.dimid[NETCDF_T_INDEX],&MY_info.id);
  nc_put_att_text(nc_id_uvpt,MY_info.id,"long_name",strlen(MY_info.long_name)+1,MY_info.long_name);
  nc_put_att_text(nc_id_uvpt,MY_info.id,"units",strlen(MY_info.units)+1,MY_info.units);

  /*
   * QB specific fields
   */

  /*
   * Potential vorticity (PV)
   *
   * NOTE: For isentropic coordinates, and for the isentropic region of hybrid coordinates,
   *       this PV is Ertel's isentropic potential vorticity. For isobaric coordinates, it is
   *       proportional to absolute vorticity on pressure levels, because h is constant.
   */
  pv_info.name = (char *)malloc(strlen("pv")+1);
  sprintf(pv_info.name,"%s","pv");
  pv_info.standard_name = (char *)malloc(strlen("potential_vorticity")+1);
  sprintf(pv_info.standard_name,"%s","potential_vorticity");
  pv_info.long_name = (char *)malloc(strlen("potential vorticity")+1);
  sprintf(pv_info.long_name,"%s","potential vorticity");
  pv_info.units = (char *)malloc(strlen("K m2 kg-1 s-1")+1);
  sprintf(pv_info.units,"%s","K m2 kg-1 s-1");
  for (i = 0; i < 4; i++) {
    pv_info.dimid[ i] = dimid[ i];
    pv_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,pv_info.name,NC_DOUBLE,4,pv_info.dimid,&pv_info.id);
  nc_put_att_text(nc_id_qb,pv_info.id,"standard_name",strlen(pv_info.standard_name)+1,pv_info.standard_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"long_name",strlen(pv_info.long_name)+1,pv_info.long_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"units",strlen(pv_info.units)+1,pv_info.units);

  /*
   * Bernoulli streamfunction
   */
  bsf_info.name = (char *)malloc(strlen("bsf")+1);
  sprintf(bsf_info.name,"%s","bsf");
  bsf_info.standard_name = (char *)malloc(strlen("specific_dry_energy_of_air")+1);
  sprintf(bsf_info.standard_name,"%s","specific_dry_energy_of_air");
  bsf_info.long_name = (char *)malloc(strlen("Bernoulli streamfunction")+1);
  sprintf(bsf_info.long_name,"%s","Bernoulli streamfunction");
  bsf_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(bsf_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    bsf_info.dimid[ i] = dimid[ i];
    bsf_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,bsf_info.name,NC_DOUBLE,4,bsf_info.dimid,&bsf_info.id);
  nc_put_att_text(nc_id_qb,bsf_info.id,"standard_name",strlen(bsf_info.standard_name)+1,bsf_info.standard_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"long_name",strlen(bsf_info.long_name)+1,bsf_info.long_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"units",strlen(bsf_info.units)+1,bsf_info.units);

  /*
   * NH2
   */
  NH2_info.name = (char *)malloc(strlen("NHsquared")+1);
  sprintf(NH2_info.name,"%s","NHsquared");
  NH2_info.long_name = (char *)malloc(strlen("(NH)^2")+1);
  sprintf(NH2_info.long_name,"%s","(NH)^2");
  NH2_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(NH2_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    NH2_info.dimid[ i] = dimid[ i];
    NH2_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,NH2_info.name,NC_DOUBLE,4,NH2_info.dimid,&NH2_info.id);
  nc_put_att_text(nc_id_qb,NH2_info.id,"long_name",strlen(NH2_info.long_name)+1,NH2_info.long_name);
  nc_put_att_text(nc_id_qb,NH2_info.id,"units",strlen(NH2_info.units)+1,NH2_info.units);
  nc_put_att_text(nc_id_qb,NH2_info.id,"comment",strlen(
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H")+1,
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H");

  /*
   * Kinetic energy
   */
  kin_info.name = (char *)malloc(strlen("kin")+1);
  sprintf(kin_info.name,"%s","kin");
  kin_info.standard_name = (char *)malloc(strlen("specific_kinetic_energy_of_air")+1);
  sprintf(kin_info.standard_name,"%s","specific_kinetic_energy_of_air");
  kin_info.long_name = (char *)malloc(strlen("Kinetic energy")+1);
  sprintf(kin_info.long_name,"%s","Kinetic energy");
  kin_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(kin_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    kin_info.dimid[ i] = dimid[ i];
    kin_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,kin_info.name,NC_DOUBLE,4,kin_info.dimid,&kin_info.id);
  nc_put_att_text(nc_id_qb,kin_info.id,"long_name",strlen(kin_info.long_name)+1,kin_info.long_name);
  nc_put_att_text(nc_id_qb,kin_info.id,"units",strlen(kin_info.units)+1,kin_info.units);
  nc_put_att_text(nc_id_qb,kin_info.id,"comment",strlen(
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)")+1,
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)");

  /*
   * UVPT specific fields
   */

  /*
   * Zonal wind
   */
  u_info.name = (char *)malloc(strlen("u")+1);
  sprintf(u_info.name,"%s","u");
  u_info.standard_name = (char *)malloc(strlen("eastward_wind")+1);
  sprintf(u_info.standard_name,"%s","eastward_wind");
  u_info.long_name = (char *)malloc(strlen("Zonal wind")+1);
  sprintf(u_info.long_name,"%s","Zonal wind");
  u_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(u_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    u_info.dimid[ i] = dimid[ i];
    u_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,u_info.name,NC_DOUBLE,4,u_info.dimid,&u_info.id);
  nc_put_att_text(nc_id_uvpt,u_info.id,"standard_name",strlen(u_info.standard_name)+1,u_info.standard_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"long_name",strlen(u_info.long_name)+1,u_info.long_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"units",strlen(u_info.units)+1,u_info.units);

//   /*
//    * Geostrophic zonal wind (assuming variable-f geostrophy)
//    */
//   ug_info.name = (char *)malloc(strlen("ug")+1);
//   sprintf(ug_info.name,"%s","ug");
//   ug_info.standard_name = (char *)malloc(strlen("geostrophic_eastward_wind")+1);
//   sprintf(ug_info.standard_name,"%s","geostrophic_eastward_wind");
//   ug_info.long_name = (char *)malloc(strlen("Geostrophic zonal wind")+1);
//   sprintf(ug_info.long_name,"%s","Geostrophic zonal wind");
//   ug_info.units = (char *)malloc(strlen("m s-1")+1);
//   sprintf(ug_info.units,"%s","m s-1");
//   for (i = 0; i < 4; i++) {
//     ug_info.dimid[ i] = dimid[ i];
//     ug_info.coorid[i] = coorid[i];
//   }
// 
//   nc_def_var(nc_id_uvpt,ug_info.name,NC_DOUBLE,4,ug_info.dimid,&ug_info.id);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"standard_name",strlen(ug_info.standard_name)+1,ug_info.standard_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"long_name",strlen(ug_info.long_name)+1,ug_info.long_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"units",strlen(ug_info.units)+1,ug_info.units);

  /*
   * Meridional wind
   */
  v_info.name = (char *)malloc(strlen("v")+1);
  sprintf(v_info.name,"%s","v");
  v_info.standard_name = (char *)malloc(strlen("northward_wind")+1);
  sprintf(v_info.standard_name,"%s","northward_wind");
  v_info.long_name = (char *)malloc(strlen("Meridional wind")+1);
  sprintf(v_info.long_name,"%s","Meridional wind");
  v_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(v_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    v_info.dimid[ i] = dimid[ i];
    v_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,v_info.name,NC_DOUBLE,4,v_info.dimid,&v_info.id);
  nc_put_att_text(nc_id_uvpt,v_info.id,"standard_name",strlen(v_info.standard_name)+1,v_info.standard_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"long_name",strlen(v_info.long_name)+1,v_info.long_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"units",strlen(v_info.units)+1,v_info.units);

  /*
   * Pressure
   */
  p_info.name = (char *)malloc(strlen("p")+1);
  sprintf(p_info.name,"%s","p");
  p_info.standard_name = (char *)malloc(strlen("air_pressure")+1);
  sprintf(p_info.standard_name,"%s","air_pressure");
  p_info.long_name = (char *)malloc(strlen("Pressure")+1);
  sprintf(p_info.long_name,"%s","Pressure");
  p_info.units = (char *)malloc(strlen("Pa")+1);
  sprintf(p_info.units,"%s","Pa");
  for (i = 0; i < 4; i++) {
    p_info.dimid[ i] = dimid[ i];
    p_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,p_info.name,NC_DOUBLE,4,p_info.dimid,&p_info.id);
  nc_put_att_text(nc_id_uvpt,p_info.id,"standard_name",strlen(p_info.standard_name)+1,p_info.standard_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"long_name",strlen(p_info.long_name)+1,p_info.long_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"units",strlen(p_info.units)+1,p_info.units);

  /*
   * Temperature
   */
  t_info.name = (char *)malloc(strlen("T")+1);
  sprintf(t_info.name,"%s","T");
  t_info.standard_name = (char *)malloc(strlen("air_temperature")+1);
  sprintf(t_info.standard_name,"%s","air_temperature");
  t_info.long_name = (char *)malloc(strlen("Temperature")+1);
  sprintf(t_info.long_name,"%s","Temperature");
  t_info.units = (char *)malloc(strlen("K")+1);
  sprintf(t_info.units,"%s","K");
  for (i = 0; i < 4; i++) {
    t_info.dimid[ i] = dimid[ i];
    t_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,t_info.name,NC_DOUBLE,4,t_info.dimid,&t_info.id);
  nc_put_att_text(nc_id_uvpt,t_info.id,"standard_name",strlen(t_info.standard_name)+1,t_info.standard_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"long_name",strlen(t_info.long_name)+1,t_info.long_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"units",strlen(t_info.units)+1,t_info.units);

  /*
   * Geopotential
   */
  phi_info.name = (char *)malloc(strlen("phi")+1);
  sprintf(phi_info.name,"%s","phi");
  phi_info.standard_name = (char *)malloc(strlen("geopotential")+1);
  sprintf(phi_info.standard_name,"%s","geopotential");
  phi_info.long_name = (char *)malloc(strlen("Geopotential")+1);
  sprintf(phi_info.long_name,"%s","Geopotential");
  phi_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(phi_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    phi_info.dimid[ i] = dimid[ i];
    phi_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,phi_info.name,NC_DOUBLE,4,phi_info.dimid,&phi_info.id);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"standard_name",strlen(phi_info.standard_name)+1,phi_info.standard_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"long_name",strlen(phi_info.long_name)+1,phi_info.long_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"units",strlen(phi_info.units)+1,phi_info.units);

  /*
   * Mach number
   */
  ma_info.name = (char *)malloc(strlen("Ma")+1);
  sprintf(ma_info.name,"%s","Ma");
  ma_info.long_name = (char *)malloc(strlen("Mach number")+1);
  sprintf(ma_info.long_name,"%s","Mach number");
  ma_info.units = (char *)malloc(strlen("1")+1);
  sprintf(ma_info.units,"%s","1");
  for (i = 0; i < 4; i++) {
    ma_info.dimid[ i] = dimid[ i];
    ma_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,ma_info.name,NC_DOUBLE,4,ma_info.dimid,&ma_info.id);
  nc_put_att_text(nc_id_uvpt,ma_info.id,"long_name",strlen(ma_info.long_name)+1,ma_info.long_name);
  nc_put_att_text(nc_id_uvpt,ma_info.id,"units",strlen(ma_info.units)+1,ma_info.units);

  /*---------------------------*
   * Leave netcdf define mode: *
   *---------------------------*/
  nc_enddef(nc_id_qb);
  nc_enddef(nc_id_uvpt);

  /*
   * Assign values to X, Y and Z coordinates.
   */
  /* lon */
  for (I = ILO; I <= IHI; I++) {
    nc_index[0] = I-ILO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
  }

  /* lat */
  for (J = JLO; J <= JHI; J++) {
    nc_index[0] = J-JLO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
  }

  /* theta */
  for (K = KHI; K >= KLO; K--) {
    nc_index[0] = KHI-K;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
  }

  fprintf(stdout,"Computing variables on theta surfaces for OpenMARS QB and UVPT files...  0%%");


  /* Allocate memory */
  p      = fvector(0,2*grid.nk+1,dbmsname);
  h      = fvector(0,2*grid.nk+1,dbmsname);
  buffer = (double *)calloc(grid.ni,sizeof(double));

  for (openmars_itime = 0; openmars_itime < openmars_grid->ntime; openmars_itime++) {
    /* Show progress */
    fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(double)(openmars_itime+1)/openmars_grid->ntime));
    fflush(stdout);

    /*
     * Input OpenMARS data for timeframe openmars_itime, and interpolate into EPIC variables.
     */
    openmars_var_read(planet,openmars_grid,openmars_infile,VAR_DATA,openmars_itime);

    /*
     * Write time.
     */

    /* QB */
    nc_index[0] = openmars_itime;
    nc_put_var1_double(nc_id_qb,coorid[NETCDF_T_INDEX],nc_index,         &openmars_grid->time[openmars_itime]);
    nc_put_var1_double(nc_id_qb,L_s_info.id,           nc_index,         &openmars_grid->Ls[  openmars_itime]);
    nc_put_var1_short( nc_id_qb,MY_info.id,            nc_index,(short *)&openmars_grid->MY[  openmars_itime]);

    /* UVPT */
    nc_index[0] = openmars_itime;
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_T_INDEX],nc_index,         &openmars_grid->time[openmars_itime]);
    nc_put_var1_double(nc_id_uvpt,L_s_info.id,           nc_index,         &openmars_grid->Ls[  openmars_itime]);
    nc_put_var1_short( nc_id_uvpt,MY_info.id,            nc_index,(short *)&openmars_grid->MY[  openmars_itime]);

    /*
     * Compute Q and B on theta surfaces by running the normal EPIC model
     * diagnostic-variable calculations.
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      for (I = ILO; I <= IHI; I++) {
        for (kk = 1; kk <= 2*KHI+1; kk++) {
          p[kk] = get_p(planet,P2_INDEX,kk,J,I);
        }
        calc_h(jj,p,h);

        for (K = KLO; K <= KHI; K++) {
          H(K,J,I) = h[2*K];
        }
        K = 0;
        H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
        K = KHI+1;
        H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
      }
    }
    bc_lateral(var.h.value,THREEDIM);

    /*
     * PHI3(KHI,J,I) is set above by openmars_var_read(VAR_DATA).
     */
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,PASSING_PHI3NK);
    store_diag(planet);

    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;
    nc_count[NETCDF_I_INDEX] = grid.ni;

    nc_start[NETCDF_T_INDEX] = openmars_itime;
    nc_start[NETCDF_I_INDEX] = 0;

    for (K = KHI; K >= KLO; K--) {
      kk                       = 2*K;
      nc_start[NETCDF_K_INDEX] = KHI-K;
      for (J = JLO; J <= JHI; J++) {
        jj = 2*J;

        g  = GRAVITY2(K,J);

        nc_start[NETCDF_J_INDEX] = J-JLO;

        /*--------------------*
         * QB specific fields *
         *--------------------*/

        /*
         * Write potential vorticity (pv, aka Q).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .25*(PV2(K,J,I)+PV2(K,J,I+1)+PV2(K,J+1,I)+PV2(K,J+1,I+1));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,pv_info.id,nc_start,nc_count,buffer);

        /*
         * Write Bernoulli streamfunction (bsf).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = MONT2(K,J,I)
                         +get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,bsf_info.id,nc_start,nc_count,buffer);

        /*
         * Write (NH)^2, the square of the product of the buoyancy (Brunt-Vaisala) frequency, N, 
         * and the pressure scale height, H = RT/g.
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_brunt2(planet,kk,J,I)
                         *SQR(planet->rgas*T2(K,J,I)/g);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,NH2_info.id,nc_start,nc_count,buffer);

        /*
         * Write kinetic energy per mass (kin), .5*(u*u+v*v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,kin_info.id,nc_start,nc_count,buffer);

        /*----------------------*
         * UVPT specific fields *
         *----------------------*/

        /*
         * Write zonal wind (u).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(U(grid.it_uv,K,J,I)+U(grid.it_uv,K,J,I+1));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,u_info.id,nc_start,nc_count,buffer);

//         /*
//          * Calculate and write geostrophic zonal wind (ug).
//          */
// 
//         if (J == JLO) {
//           /*
//            * Southern edge: take forward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         else if (J == JHI) {
//           /*
//            * Northern edge: take backward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J,I)-MONT2(K,J-1,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         else {
//           /*
//            * Interior point: take long step to center back onto h-grid.
//            */
//           tmp = -.5*grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J-1,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         nc_put_vara_double(nc_id_uvpt,ug_info.id,nc_start,nc_count,buffer);

        /*
         * Write meridional wind (v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(V(grid.it_uv,K,J,I)+V(grid.it_uv,K,J+1,I));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,v_info.id,nc_start,nc_count,buffer);

        /*
         * Write pressure (p).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = P2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,p_info.id,nc_start,nc_count,buffer);

        /*
         * Write temperature (T).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = T2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,t_info.id,nc_start,nc_count,buffer);

        /*
         * Write geopotential (phi).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = PHI2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,phi_info.id,nc_start,nc_count,buffer);

        /*
         * Write Mach number, Ma.
         */
        for (I = ILO; I <= IHI; I++) {
          c2factor = return_cp(planet,0.,0.,T2(K,J,I));
          c2factor = planet->rgas*c2factor/(c2factor-planet->rgas);
          buffer[I-ILO] = sqrt(2.*get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                                 var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I)
                               /(T2(K,J,I)*c2factor));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,ma_info.id,nc_start,nc_count,buffer);
      } /* J loop */
    } /* K loop */
  } /* openmars_itime loop */

  nc_close(nc_id_qb);
  nc_close(nc_id_uvpt);

  fprintf(stdout,"\b\b\b\b100%%\n");
  fflush(stdout);
  fprintf(stdout,"Output written to %s\n",openmars_outfile_qb);

  /*
   * The conversion is complete, so clean up (write default values,
   * free allocated memory, etc.) and exit.
   */
  free_fvector(p,0,2*grid.nk+1,dbmsname);
  free_fvector(h,0,2*grid.nk+1,dbmsname);
  free(buffer);
}

/*======================= end of openmars_conversion() =======================*/

/*======================= emars_var_read() ===================================*/

/*
 * Read in EMARS data.
 * portion = SIZE_DATA:       read in grid dimensions
 * portion = POST_SIZE_DATA:  read in constant parameters
 * portion = VAR_DATA:        read in variables
 */

#define FBUFF2D(j,i)        fbuffer[i+emars_grid->ni*(j)]
#define DBUFF2D(j,i)        dbuffer[i+emars_grid->ni*(j)]
#define FBUFF3D(k,j,i)      fbuffer[i+emars_grid->ni*(j+emars_grid->nj*(k))]
#define EPIC_K_EMARS_j(K,j) epic_k_emars_j[j+emars_grid->nj*(K-KLOPAD)]

void emars_var_read(planetspec     *planet,
                    emars_gridspec *emars_grid,
                    char           *emars_infile,
                    int             portion,
                    int             input_itime)
{
  int
    K,J,I,
    kk,k,j,i,ii,
    itime,
    is,
    nelem2d,nelem3d;
  char
    start_season[32],
    end_season[32];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0,
    i_seam;
  int
    nc_err,nc_id,
    nc_dimid,nc_varid;
  double
    g0 = 3.71,   /* gravity standard for Mars, m/s^2 */
    lowest_lon;
  int
    *ibuffer;
  float
    *fbuffer;
  double
    *dbuffer;
  size_t
    dimlen;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="emars_var_read";

  nc_err = lookup_netcdf(emars_infile,&nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
    epic_error(dbmsname,Message);
  }

  if (portion == SIZE_DATA) {
    fprintf(stdout,"Reading SIZE_DATA from %s \n",emars_infile);
  }
  else if (portion == POST_SIZE_DATA) {
    fprintf(stdout,"Reading POST_SIZE_DATA from %s \n",emars_infile);
  }
  else if (portion == VAR_DATA) {
    ;
  }
  else {
    sprintf(Message,"unrecognized portion = %d",portion);
    epic_error(dbmsname,Message);
  }
  fflush(stderr);

  /*
   * Read in size of model and set parameters needed by emars_make_arrays().
   */
  if (portion == SIZE_DATA) {
    /* 
     * Number of longitude points
     */
    nc_err = nc_inq_dimid(nc_id,"lon",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }
    emars_grid->ni = (int)dimlen;

    /* 
     * Number of latitude points
     */
    nc_err = nc_inq_dimid(nc_id,"lat",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }
    emars_grid->nj = (int)dimlen;

    /* 
     * Number of levels.
     *
     * NOTE: u, v, t are on "pfull" and the hybrid sigma-p variables, ak and bk, are on "phalf",
     *       which has one more position than pfull. We are using the pfull number of levels.
     */
    nc_err = nc_inq_dimid(nc_id,"pfull",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }
    emars_grid->nk = (int)dimlen;

    /* 
     * Number of times.
     *
     * NOTE: EMARS records every Mars hour, whereas OpenMARS records every 2 hours.
     *       We drop back to every 2 hours for EMARS.
     */
    nc_err = nc_inq_dimid(nc_id,"time",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_infile);
      epic_error(dbmsname,Message);
    }
    emars_grid->ntime = ((int)dimlen+1)/2;
  }
  else if (portion == POST_SIZE_DATA) {
    /*
     * NOTE: POST_SIZE_DATA is interpreted here to mean constant parameters,
     *       including the constant 1D dimension arrays, lon, lat, etc.
     *       and the 2D surface geopotential.
     */

    /*-------------------------*
     * Read 1D lon array [deg] *
     *-------------------------*/

    dbuffer = calloc(emars_grid->ni,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"lon",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    lowest_lon = 360.;
    for (i = 0; i < emars_grid->ni; i++) {
      /* 
       * Map EMARS [180,360] east longitude to [-180,0] east longitude, to match OpenMARS and EPIC.
       */
      emars_grid->lon[i] = (dbuffer[i] < 180.) ? dbuffer[i] : dbuffer[i]-360.;
      /*
       * Find longitude-seam index value, i_seam.
       */
      if (emars_grid->lon[i] < lowest_lon) {
        i_seam     = i;
        lowest_lon = emars_grid->lon[i];
      }
    }

    /*
     * Roll lon array to put the lon = -180 seam at the beginning, to match OpenMARS and EPIC.
     *
     * NOTE: This must be done consistently for each EMARS array that has a longitude dimension.
     */
    ii = 0;
    for (i = i_seam; i < emars_grid->ni; i++) {
      dbuffer[ii++] = emars_grid->lon[i];
    }
    for (i = 0; i < i_seam; i++) {
      dbuffer[ii++] = emars_grid->lon[i];
    }
    for (i = 0; i < emars_grid->ni; i++) {
      emars_grid->lon[i] = dbuffer[i];
    }

    free(dbuffer);

    /*
     * If anything involving a zonal FFT is to be used,
     * then grid.ni needs to be a power of 2, for example:
     *   grid.ni = 1 << NINT(log(emars_grid->ni)/log(2.));
     *
     * NOTE: EMARS has ni = 60 and OpenMARS has ni = 72. We are keeping
     *       these native values, rather than interpolating one onto the other.
     */
    grid.ni = emars_grid->ni;

    /*-------------------------*
     * Read 1D lat array [deg] *
     *-------------------------*/

    dbuffer = calloc(emars_grid->nj,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"lat",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Both OpenMARS and EMARS have 36 values for latitude.  However, the EMARS spacing
     *       is 5.14286 deg, except at the polar edges where it is 3.85714 deg, whereas the OpenMARS spacing
     *       is 5.0 deg everywhere. We set EPIC's grid.nj the same for both.
     *
     * NOTE: Unlike OpenMARS, the EMARS latitudes are stored in the same order as the EPIC convention.
     */
    for (j = 0; j < emars_grid->nj; j++) {
       emars_grid->lat[j] = (EPIC_FLOAT)dbuffer[j];
    }

    free(dbuffer);

    grid.nj = NINT(180./(NINT(emars_grid->lat[2]-emars_grid->lat[1]))-1.-sqrt(2.));

    /*------------------------------------------*
     * Read 1D ak and bk arrays,                *
     *  interface p := phalf = ak + bk*psurface *
     *------------------------------------------*/

    /*
     * NOTE: ak and bk are on the phalf levels, which have one more value than the pfull levels.
     */
    fbuffer = calloc(emars_grid->nk+1,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"ak",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Unlike OpenMARS sigma (lev), the EMARS ak and bk arrays are stored
     *       in the same top-down order (the planet's surface is last) as the EPIC convention.
     */
    for (k = 0; k < emars_grid->nk+1; k++) {
      emars_grid->ak[k] = (EPIC_FLOAT)fbuffer[k];
    }

    nc_err = nc_inq_varid(nc_id,"bk",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_float(nc_id,nc_varid,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < emars_grid->nk+1; k++) {
      emars_grid->bk[k] = (EPIC_FLOAT)fbuffer[k];
    }

    free(fbuffer);

    /*
     * OpenMARS and EMARS have 35 and 28 levels, with different
     * placements, and use sigma and hybrid sigma-pressure coordinates,
     * respectively. Despite these differences, the vertical coverage is similar
     * enough for us to map them to the same potential temperature, theta, surfaces.
     *
     * NOTE: We hardwire grid.nk to yield rounded theta levels. This needs to be adjusted for each new application.
     */
    grid.nk = 29;

    /*---------------------------*
     * Read 1D time array [sols] *
     *---------------------------*/

    /*
     * EMARS includes macda_sol (graciously), which when combined with mars_hour/24.
     * yields the time variable used in OpenMARS, as constructed here.
     *
     * NOTE: EMARS saves every hour, so we are sampling every two hours
     *       to match OpenMARS.
     */

    dbuffer = calloc(emars_grid->ntime*2,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"macda_sol",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    for (itime = 0; itime < emars_grid->ntime; itime++) {
      emars_grid->time[itime] = (EPIC_FLOAT)dbuffer[itime*2];
    }

    nc_err = nc_inq_varid(nc_id,"mars_hour",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
    for (itime = 0; itime < emars_grid->ntime; itime++) {
       emars_grid->time[itime] += (EPIC_FLOAT)(dbuffer[itime*2]/24.);
    }

    free(dbuffer);

    /*-----------------------------------------*
     * Read 1D Ls array, Solar Longitude [deg] *
     *-----------------------------------------*/

    dbuffer = calloc(emars_grid->ntime*2,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"Ls",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Montabone et al (2014, Geosci. Data J. 1, p. 134) suggest subtracting 0.12deg from the MACDA Ls
     *       to correct for a near-constant bias (midnight vs. actual start of spring), which we implement
     *       for OpenMARS (Ls matches precisely in the MACDA and OpenMARS data sets). We compared Ls in
     *       EMARS versus OpenMARS, by spot checking at time = 4351.083 (sol since start of MY 24) and the
     *       EMARS Ls is nearly the same as the corrected value for OpenMARS, 161.660 vs. 161.795-0.12 = 161.675,
     *       respectively. Hence, we do not do a correction on the EMARS Ls values.
     */
    for (itime = 0; itime < emars_grid->ntime; itime++) {
      emars_grid->Ls[itime] = (EPIC_FLOAT)dbuffer[itime*2];
    }

    season_string(emars_grid->Ls[0                  ],start_season);
    season_string(emars_grid->Ls[emars_grid->ntime-1],  end_season);
    if (strcmp(start_season,end_season) == 0) {
      fprintf(stdout,"  Input data span L_s = %6.2f to %6.2f (%s)\n",
                     emars_grid->Ls[0],emars_grid->Ls[emars_grid->ntime-1],start_season);
    }
    else {
      fprintf(stdout,"   Input data span L_s = %6.2f to %6.2f (%s to %s)\n",
                     emars_grid->Ls[0],emars_grid->Ls[emars_grid->ntime-1],start_season,end_season);
    }

    free(dbuffer);

    /*----------------------------------------------------*
     * Read 1D MY array, Martian Year (Clancy et al 2000) *
     *----------------------------------------------------*/

    dbuffer = calloc(emars_grid->ntime*2,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"MY",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (itime = 0; itime < emars_grid->ntime; itime++) {
      emars_grid->MY[itime] = (int)dbuffer[itime*2];
    }
    free(dbuffer);

    /*----------------------------------------------------*
     * Read 2D Surface_geopotential array                 *
     *                                                    *
     * NOTE: This is obtained from emars_v1.0_anal_*.nc   *
     *----------------------------------------------------*/
    {
      char
        **anal_gattname=NULL,
        **anal_varname =NULL;
      int
        anal_ngatts    =0,
        anal_num_progs =0,
        anal_nc_id;
      size_t
        anal_nc_start2d[2],
        anal_nc_count2d[2];
      char
        anal_emars_infile[]="emars_v1.0_anal_mean_MY24_Ls090-120.nc";

      nc_err = lookup_netcdf(anal_emars_infile,&anal_nc_id,&anal_ngatts,&anal_gattname,&anal_num_progs,&anal_varname);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),anal_emars_infile);
        epic_error(dbmsname,Message);
      }

      nelem2d            = emars_grid->nj*emars_grid->ni;
      anal_nc_start2d[0] = 0;
      anal_nc_start2d[1] = 0;
      anal_nc_count2d[0] = emars_grid->nj;
      anal_nc_count2d[1] = emars_grid->ni;

      dbuffer = calloc(nelem2d,sizeof(double));
      if (!dbuffer) {
        sprintf(Message,"calloc error allocating dbuffer");
        epic_error(dbmsname,Message);
      }

      nc_err = nc_inq_varid(anal_nc_id,"Surface_geopotential",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      nc_err = nc_get_vara_double(anal_nc_id,nc_varid,anal_nc_start2d,anal_nc_count2d,dbuffer);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      for (j = 0; j < emars_grid->nj; j++) {
        /*
         * Roll longitude -180 seam to beginning of array.
         */
        ii = 0;
        for (i = i_seam; i < emars_grid->ni; i++) {
          EMARS_PHI_SURFACE(j,ii++) = (EPIC_FLOAT)DBUFF2D(j,i);
        }
        for (i = 0; i < i_seam; i++) {
          EMARS_PHI_SURFACE(j,ii++) = (EPIC_FLOAT)DBUFF2D(j,i);
        }
      }

      free(dbuffer);
      nc_close(anal_nc_id);
    }
  }
  else if (portion == VAR_DATA) {
    int
      count,
      excluded_count;
    EPIC_FLOAT
      p1,p2,p3,
      tmp,x,x_d,
     *epic_k_emars_j,
     *emars_phi;
    size_t
      nc_start2d[3],nc_start3d[4],
      nc_count2d[3],nc_count3d[4];
    float_triplet
      *k_triplet,
      *j_triplet,
      *i_triplet,
      *epic_j_triplet;
    FILE
      *excluded_file;

    /*
     * NOTE: VAR_DATA is interpreted here to mean one timeframe, as specified
     *       by the input argument input_itime, of the variable fields like ps, etc.
     *       We are sampling every other timeframe in EMARS.
     */

    nelem2d       = emars_grid->nj*emars_grid->ni;
    nc_start2d[0] = input_itime*2;                   /* sampling every other timeframe */
    nc_start2d[1] = 0;
    nc_start2d[2] = 0;
    nc_count2d[0] = 1;
    nc_count2d[1] = emars_grid->nj;
    nc_count2d[2] = emars_grid->ni;

    /*--------------------------------------*
     * Read ps array, surface pressure [Pa] *
     *--------------------------------------*/

    /*
     * NOTE: The EMARS and EPIC 2D and 3D array layouts are the same, but the data types
     *       are typically different, i.e. input floats vs. internal doubles (EPIC_FLOAT is usually
     *       specified to be double). Thus, we use point-by-point assignment instead of memcpy().
     */

    fbuffer = calloc(nelem2d,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"ps",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start2d,nc_count2d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < emars_grid->nj; j++) {
      /*
       * Roll longitude -180 seam to beginning of array.
       */
      ii = 0;
      for (i = i_seam; i < emars_grid->ni; i++) {
        EMARS_PS(j,ii++) = (EPIC_FLOAT)FBUFF2D(j,i);
      }
      for (i = 0; i < i_seam; i++) {
        EMARS_PS(j,ii++) = (EPIC_FLOAT)FBUFF2D(j,i);
      }
    }

    /*-------------------------------------------*
     * Read tsurf array, surface temperature [K] *
     *-------------------------------------------*/

    /*
     * EMARS:ts = OpenMARS:tsurf; we call it tsurf in EPIC.
     */
    nc_err = nc_inq_varid(nc_id,"ts",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start2d,nc_count2d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (j = 0; j < emars_grid->nj; j++) {
      /*
       * Roll longitude -180 seam to beginning of array.
       */
      ii = 0;
      for (i = i_seam; i < emars_grid->ni; i++) {
        EMARS_TSURF(j,ii++) = (EPIC_FLOAT)FBUFF2D(j,i);
      }
      for (i = 0; i < i_seam; i++) {
        EMARS_TSURF(j,ii++) = (EPIC_FLOAT)FBUFF2D(j,i);
      }
    }

    free(fbuffer);

    /*-------------------------------------------*
     * Read temp array, 3D temperature field [K] *
     *-------------------------------------------*/

    /*
     * EMARS:t = OpenMARS:temp; we call it temp in EPIC.
     */
    nelem3d       = nelem2d*emars_grid->nk;
    nc_start3d[0] = input_itime*2;
    nc_start3d[1] = 0;
    nc_start3d[2] = 0;
    nc_start3d[3] = 0;
    nc_count3d[0] = 1;
    nc_count3d[1] = emars_grid->nk;
    nc_count3d[2] = emars_grid->nj;
    nc_count3d[3] = emars_grid->ni;

    fbuffer = calloc(nelem3d,sizeof(float));
    if (!fbuffer) {
      sprintf(Message,"calloc error allocating fbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"t",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < emars_grid->nk; k++) {
      for (j = 0; j < emars_grid->nj; j++) {
        /*
         * Roll longitude -180 seam to beginning of array.
         */
        ii = 0;
        for (i = i_seam; i < emars_grid->ni; i++) {
          EMARS_TEMP(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
        for (i = 0; i < i_seam; i++) {
          EMARS_TEMP(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
      }
    }

    /*-----------------------------------------*
     * Read u array, 3D zonal wind field [m/s] *
     *-----------------------------------------*/

    nc_err = nc_inq_varid(nc_id,"u",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < emars_grid->nk; k++) {
      for (j = 0; j < emars_grid->nj; j++) {
        /*
         * Roll longitude -180 seam to beginning of array.
         */
        ii = 0;
        for (i = i_seam; i < emars_grid->ni; i++) {
          EMARS_U(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
        for (i = 0; i < i_seam; i++) {
          EMARS_U(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
      }
    }

    /*----------------------------------------------*
     * Read v array, 3D meridional wind field [m/s] *
     *----------------------------------------------*/

    nc_err = nc_inq_varid(nc_id,"v",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start3d,nc_count3d,fbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < emars_grid->nk; k++) {
      for (j = 0; j < emars_grid->nj; j++) {
        /*
         * Roll longitude -180 seam to beginning of array.
         */
        ii = 0;
        for (i = i_seam; i < emars_grid->ni; i++) {
          EMARS_V(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
        for (i = 0; i < i_seam; i++) {
          EMARS_V(k,j,ii++) = (EPIC_FLOAT)FBUFF3D(k,j,i);
        }
      }
    }

    free(fbuffer);

    /*---------------------------------------------------------*
     * Compute pressure and potential temperature arrays,      *
     *                                                         *
     * NOTE: For Mars, planet->p0 = 610 Pa                     *
     *---------------------------------------------------------*/


    for (k = 0; k < emars_grid->nk; k++) {
      for (j = 0; j < emars_grid->nj; j++) {
        for (i = 0; i < emars_grid->ni; i++) {
          /*
           * Construct layer pressure. See Greybush et al. (2019), p. 142.
           */
          p1 = emars_grid->ak[k  ]+emars_grid->bk[k  ]*EMARS_PS(j,i);
          p3 = emars_grid->ak[k+1]+emars_grid->bk[k+1]*EMARS_PS(j,i);
          EMARS_P(k,j,i)     = (p3-p1)/log(p3/p1);
          EMARS_THETA(k,j,i) = EMARS_TEMP(k,j,i)*pow(planet->p0/EMARS_P(k,j,i),planet->kappa);
        }
      }
    }


    /*---------------------------------------------------*
     * Interpolate EMARS data onto EPIC staggered grids. *
     *---------------------------------------------------*/

    /* Allocate memory */
    k_triplet         = ftriplet(0,emars_grid->nk-1,               dbmsname);
    j_triplet         = ftriplet(0,emars_grid->nj+1,               dbmsname);
    i_triplet         = ftriplet(0,emars_grid->ni,                 dbmsname);
    epic_j_triplet    = ftriplet(0,grid.nj,                        dbmsname);
    epic_k_emars_j    = fvector( 0,(KHI-KLOPAD+1)*emars_grid->nj-1,dbmsname);
    emars_phi         = fvector( 0,emars_grid->nk,                 dbmsname);  /* on interfaces, so gets an extra position */

    /*--------------------------------------------*
     * Zonal wind, U                              *
     * Longitude can be done by direct insertion. * 
     *--------------------------------------------*/

    excluded_file = fopen("excluded_thetas.dat","w");
    if (!excluded_file) {
      sprintf(Message,"file %s not found","excluded_thetas.dat");
      epic_error(dbmsname,Message);
    }
    fprintf(excluded_file," Points in column excluded because of lack of monotonicity of theta.\n");
    fprintf(excluded_file,"  Lon    Lat  Excluded thetas [K]\n");

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping EMARS latitudes.
       */
      for (j = 0; j < emars_grid->nj; j++) {
        k              = emars_grid->nk-1;
        k_triplet[0].x = EMARS_THETA(k,j,i);
        k_triplet[0].y = EMARS_U(k,j,i);
        count          = 1;
        excluded_count = 0;
        for (k = emars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (EMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = EMARS_THETA(k,j,i);
            k_triplet[count].y = EMARS_U(k,j,i);
            count++;
          }
          else {
            /*
             * Keep track of number of excluded points caused by non-monotonicity of theta.
             */
            excluded_count++;
            if (excluded_count == 1) {
              fprintf(excluded_file," %6.1f %5.1f  %5.1f",emars_grid->lon[i],emars_grid->lat[j],EMARS_THETA(k,j,i));
            }
            else {
              fprintf(excluded_file," %5.1f",EMARS_THETA(k,j,i));
            }
          }
        }

        if (excluded_count > 0) {
          fprintf(excluded_file,"\n");
        }

        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLO; K <= KHI; K++) {
          /*
           * U resides in the layer.
           */
          x = grid.sigmatheta[2*K];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_EMARS_j(K,j) = EMARS_U(emars_grid->nk-1,j,i);
          }
          else {
            k                   = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_EMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       * We add U = 0 pole points.
       */
      for (K = KLO; K <= KHI; K++) {
        j_triplet[0].x = -90.;
        j_triplet[0].y =   0.;
        for (j = 0; j < emars_grid->nj; j++) {
          j_triplet[j+1].x = emars_grid->lat[j];
          j_triplet[j+1].y = EPIC_K_EMARS_j(K,j);
        }
        j_triplet[j+1].x = 90.;
        j_triplet[j+1].y =  0.;

        spline_pchip(emars_grid->nj+2,j_triplet);
        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x                   = grid.lat[2*J+1];
          j                   = hunt_place_in_table(emars_grid->nj+2,j_triplet,x,&x_d,j);
          U(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

    fclose(excluded_file);

    /*-------------------------------------------------------------------*
     * Meridonal wind, V                                                 *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping EMARS latitudes.
       */
      for (j = 0; j < emars_grid->nj; j++) {
        k              = emars_grid->nk-1;
        k_triplet[0].x = EMARS_THETA(k,j,i);
        k_triplet[0].y = EMARS_V(k,j,i);
        count          = 1;
        for (k = emars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (EMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = EMARS_THETA(k,j,i);
            k_triplet[count].y = EMARS_V(k,j,i);
            count++;
          }
        }
        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLO; K <= KHI; K++) {
          x = grid.sigmatheta[2*K];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_EMARS_j(K,j) = EMARS_V(emars_grid->nk-1,j,i);
          }
          else {
            k                   = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_EMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       * We add V = 0 pole points.
       */
      for (K = KLO; K <= KHI; K++) {
        j_triplet[0].x = -90.;
        j_triplet[0].y =   0.;
        for (j = 0; j < emars_grid->nj; j++) {
          j_triplet[j+1].x = emars_grid->lat[j];
          j_triplet[j+1].y = EPIC_K_EMARS_j(K,j);
        }
        j_triplet[j+1].x = 90.;
        j_triplet[j+1].y =  0.;
        spline_pchip(emars_grid->nj+2,j_triplet);
        j = -2;
        for (J = JFIRST; J <= JHI; J++) {
          x                   = grid.lat[2*J];
          j                   = hunt_place_in_table(emars_grid->nj+2,j_triplet,x,&x_d,j);
          V(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JFIRST; J <= JHI; J++) {
        for (i = 0; i < emars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = emars_grid->lon[i];
          i_triplet[i].y = V(grid.it_uv,K,J,I);
        }
        i_triplet[i].x = emars_grid->lon[i-1]+emars_grid->lon[1]-emars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(emars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x                   = grid.lon[2*I+1];
          i                   = hunt_place_in_table(emars_grid->ni+1,i_triplet,x,&x_d,i);
          V(grid.it_uv,K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

    /*-------------------------------------------------------------------*
     * Pressure on the layer interfaces, P3.                             *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate EMARS in-layer values onto EPIC THETA surfaces, keeping EMARS latitudes.
       */
      for (j = 0; j < emars_grid->nj; j++) {
        k              = emars_grid->nk-1;
        k_triplet[0].x = EMARS_THETA(k,j,i);
        k_triplet[0].y = EMARS_P(    k,j,i);
        count          = 1;
        for (k = emars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (EMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = EMARS_THETA(k,j,i);
            k_triplet[count].y = EMARS_P(    k,j,i);
            count++;
          }
        }
        spline_pchip(count,k_triplet);

        k = -2;
        for (K = KLOPAD; K <= KHI; K++) {
          x = grid.sigmatheta[2*K+1];

          if (x < k_triplet[0].x) {
            /* 
             * Layer is below the ground, so place it on the ground.
             */
            EPIC_K_EMARS_j(K,j) = EMARS_PS(j,i);
          }
          else {
            k                   = hunt_place_in_table(count,k_triplet,x,&x_d,k);
            EPIC_K_EMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
          }
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      for (K = KLOPAD; K <= KHI; K++) {
        for (j = 0; j < emars_grid->nj; j++) {
          j_triplet[j].x = emars_grid->lat[j];
          j_triplet[j].y = EPIC_K_EMARS_j(K,j);
        }
        spline_pchip(emars_grid->nj,j_triplet);

        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x         = grid.lat[2*J+1];
          j         = hunt_place_in_table(emars_grid->nj,j_triplet,x,&x_d,j);
          P3(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLOPAD; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (i = 0; i < emars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = emars_grid->lon[i];
          i_triplet[i].y = P3(K,J,I);
        }
        i_triplet[i].x = emars_grid->lon[i-1]+emars_grid->lon[1]-emars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(emars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x         = grid.lon[2*I+1];
          i         = hunt_place_in_table(emars_grid->ni+1,i_triplet,x,&x_d,i);
          P3(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */ 
    bc_lateral(var.p3.value,THREEDIM);

    /*---------------------------------------------------------------------*
     * Potential temperature, THETA.                                       *
     * Using isentropic coordinates, so THETA takes the coordinate value,  *
     * except when the layer is underground, in which case the surface     *
     * value is used.                                                      *
     *---------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      for (j = 0; j < emars_grid->nj; j++) {
        for (K = KLOPAD; K <= KHI; K++) {
          /*
           * NOTE: Using the bottom position of EPIC layers, 2*K+1 rather than 2*K, 
           *       to make it easier to specify the output THETA values as whole numbers
           *       (we are not setting up variables for an EPIC model run).
           */
          EPIC_K_EMARS_j(K,j) = MAX(grid.sigmatheta[2*K+1],EMARS_THETA(emars_grid->nk-1,j,i));
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      for (K = KLOPAD; K <= KHI; K++) {
        for (j = 0; j < emars_grid->nj; j++) {
          j_triplet[j].x = emars_grid->lat[j];
          j_triplet[j].y = EPIC_K_EMARS_j(K,j);
        }
        spline_pchip(emars_grid->nj,j_triplet);

        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x            = grid.lat[2*J+1];
          j            = hunt_place_in_table(emars_grid->nj,j_triplet,x,&x_d,j);
          THETA(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLOPAD; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (i = 0; i < emars_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = emars_grid->lon[i];
          i_triplet[i].y = THETA(K,J,I);
        }
        i_triplet[i].x = emars_grid->lon[i-1]+emars_grid->lon[1]-emars_grid->lon[0];
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(emars_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x            = grid.lon[2*I+1];
          i            = hunt_place_in_table(emars_grid->ni+1,i_triplet,x,&x_d,i);
          THETA(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */ 
    bc_lateral(var.theta.value,THREEDIM);

    /*--------------------------------------------------------------------*
     * Geopotential at the bottom of the model, PHI3(KHI,J,I). This is    *
     * the geopotential on the bottom isentropic surface, grid.thetabot,  *
     * not on the planet's surface, PHI_SURFACE(J,I).                     *
     * Longitude needs to be staggered, but start with direct insertion.  *
     *--------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * Unlike for OpenMARS, for EMARS the surface geopotential is provided,
       * so we do not need to use EPIC's version.
       */
      for (j = 0; j < emars_grid->nj; j++) {
        /*
         * Calculate the EMARS geopotential vertically by integrating
         * the hydrostatic balance equation up from the surface.
         */
        emars_phi[emars_grid->nk] = EMARS_PHI_SURFACE(j,i);
        for (k = emars_grid->nk-1; k >= 0; k--) {
          /*
           *  Hydrostatic integration algorithm (k increases towards greater pressure): 
           *    phi1 = phi3+(p3-p1)/rho2; with the ideal-gas law, 1/rho2 = rgas*temp2/p2.
           *
           *    1 ---phi[0]-------------
           *    2    p[0], temp[0]       top layer, k = 0
           *    3 ---phi[1]-------------
           *    ...
           *    1 ---phi[nk-1]----------
           *    2    p[nk-1], temp[nk-1] bottom layer, k = nk-1
           *    3 ---phi[nk]------------ surface
           */
          p1 = emars_grid->ak[k  ]+emars_grid->bk[k  ]*EMARS_PS(j,i);
          p3 = emars_grid->ak[k+1]+emars_grid->bk[k+1]*EMARS_PS(j,i);

          emars_phi[k] = emars_phi[k+1]+(p3-p1)*planet->rgas*EMARS_TEMP(k,j,i)/EMARS_P(k,j,i);
        }

        /*
         * Interpolate onto EPIC theta surfaces, keeping EMARS latitudes.
         */
        k              = emars_grid->nk-1;
        k_triplet[0].x = EMARS_THETA(k,j,i);
        k_triplet[0].y = .5*(emars_phi[k]+emars_phi[k+1]);
        count          = 1;
        for (k = emars_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (EMARS_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = EMARS_THETA(k,j,i);
            k_triplet[count].y = .5*(emars_phi[k]+emars_phi[k+1]);
            count++;
          }
        }
        spline_pchip(count,k_triplet);

        /* 
         * Only need the geopotential at grid.thetabot. 
         */
        K = KHI;
        x = grid.sigmatheta[2*K+1];

        if (x < k_triplet[0].x) {
          /* 
           * Layer is below the ground, so place it on the ground.
           */
          EPIC_K_EMARS_j(K,j) = emars_phi[emars_grid->nk];
        }
        else {
          k                   = -2; 
          k                   = hunt_place_in_table(count,k_triplet,x,&x_d,k);
          EPIC_K_EMARS_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
        }
      } /* j loop */

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      K = KHI;
      for (j = 0; j < emars_grid->nj; j++) {
        j_triplet[j].x = emars_grid->lat[j];
        j_triplet[j].y = EPIC_K_EMARS_j(K,j);
      }
      spline_pchip(emars_grid->nj,j_triplet);
      j = -2;
      for (J = JLO; J <= JHI; J++) {
        x           = grid.lat[2*J+1];
        j           = hunt_place_in_table(emars_grid->nj,j_triplet,x,&x_d,j);
        PHI3(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
      }
    }  /* I loop */

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    K = KHI;
    for (J = JLO; J <= JHI; J++) {
      for (i = 0; i < emars_grid->ni; i++) {
        I              = i+1;
        i_triplet[i].x = emars_grid->lon[i];
        i_triplet[i].y = PHI3(K,J,I);
      }
      i_triplet[i].x = emars_grid->lon[i-1]+emars_grid->lon[1]-emars_grid->lon[0];
      i_triplet[i].y = i_triplet[0].y;
      periodic_spline_pchip(emars_grid->ni+1,i_triplet);
      i = -2;
      for (I = ILO; I <= IHI; I++) {
        x           = grid.lon[2*I+1];
        i           = hunt_place_in_table(emars_grid->ni+1,i_triplet,x,&x_d,i);
        PHI3(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.phi3.value,THREEDIM);

    /* Free allocated memory */
    free_ftriplet(k_triplet,     0,emars_grid->nk-1,               dbmsname);
    free_ftriplet(j_triplet,     0,emars_grid->nj-1,               dbmsname);
    free_ftriplet(i_triplet,     0,emars_grid->ni,                 dbmsname);
    free_ftriplet(epic_j_triplet,0,grid.nj,                        dbmsname);
    free_fvector(epic_k_emars_j, 0,(KHI-KLOPAD+1)*emars_grid->nj-1,dbmsname);
    free_fvector(emars_phi,      0,emars_grid->nk,                 dbmsname);

  } /* end portion == VAR_DATA */

  return;
}

#undef FBUFF2D
#undef DBUFF2D
#undef FBUFF3D
#undef EPIC_K_EMARS_j

/*======================= end of emars_var_read() ===========================*/

/*======================= emars_make_arrays() ===============================*/

void emars_make_arrays(planetspec     *planet,
                       emars_gridspec *emars_grid)
{
  int
    nelem2d,nelem3d;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="emars_make_arrays";

  nelem2d = emars_grid->nj*emars_grid->ni;
  nelem3d = nelem2d*emars_grid->nk;

  emars_grid->lon         = fvector(0,emars_grid->ni-1,   dbmsname);
  emars_grid->lat         = fvector(0,emars_grid->nj-1,   dbmsname);
  emars_grid->ak          = fvector(0,emars_grid->nk,     dbmsname);
  emars_grid->bk          = fvector(0,emars_grid->nk,     dbmsname);
  emars_grid->time        = fvector(0,emars_grid->ntime-1,dbmsname);
  emars_grid->Ls          = fvector(0,emars_grid->ntime-1,dbmsname);
  emars_grid->MY          = ivector(0,emars_grid->ntime-1,dbmsname);
  emars_grid->ps          = fvector(0,nelem2d-1,          dbmsname);
  emars_grid->tsurf       = fvector(0,nelem2d-1,          dbmsname);
  emars_grid->phi_surface = fvector(0,nelem2d-1,          dbmsname);
  emars_grid->temp        = fvector(0,nelem3d-1,          dbmsname);
  emars_grid->p           = fvector(0,nelem3d-1,          dbmsname);
  emars_grid->u           = fvector(0,nelem3d-1,          dbmsname);
  emars_grid->v           = fvector(0,nelem3d-1,          dbmsname);
  emars_grid->theta       = fvector(0,nelem3d-1,          dbmsname);

  return;
}

/*======================= end of emars_make_arrays() ========================*/

/*======================= emars_free_arrays() ===============================*/

void emars_free_arrays(planetspec     *planet,
                       emars_gridspec *emars_grid)
{
  int
    nelem2d,nelem3d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="emars_free_arrays";

  nelem2d = emars_grid->ni*emars_grid->nj;
  nelem3d = nelem2d*emars_grid->nk;

  free_fvector(emars_grid->lon,        0,emars_grid->ni-1,   dbmsname);
  free_fvector(emars_grid->lat,        0,emars_grid->nj-1,   dbmsname);
  free_fvector(emars_grid->ak,         0,emars_grid->nk,     dbmsname);
  free_fvector(emars_grid->bk,         0,emars_grid->nk,     dbmsname);
  free_fvector(emars_grid->time,       0,emars_grid->ntime-1,dbmsname);
  free_fvector(emars_grid->Ls,         0,emars_grid->ntime-1,dbmsname);
  free_ivector(emars_grid->MY,         0,emars_grid->ntime-1,dbmsname);
  free_fvector(emars_grid->ps,         0,nelem2d-1,          dbmsname);
  free_fvector(emars_grid->tsurf,      0,nelem2d-1,          dbmsname);
  free_fvector(emars_grid->phi_surface,0,nelem2d-1,          dbmsname);
  free_fvector(emars_grid->temp,       0,nelem3d-1,          dbmsname);
  free_fvector(emars_grid->u,          0,nelem3d-1,          dbmsname);
  free_fvector(emars_grid->v,          0,nelem3d-1,          dbmsname);
  free_fvector(emars_grid->theta,      0,nelem3d-1,          dbmsname);

  return;
}

/*======================= end of emars_free_arrays() ========================*/

/*======================= emars_epic_nc() ===================================*/

/*
 * Create emars_epic.nc by running initial.
 *
 * The input parameters are "hard wired" and need to be adjusted for each new application.
 */

void emars_epic_nc(planetspec *planet)
{
  pid_t
    pid;
  int
    commpipe[2],
    ierr,status,
    saved_stdout;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="emars_epic_nc";

  /*
   * The code for executing initial from change is adapted from examples at following websites:
   *    https://www.gidforums.com/t-3369.html
   *    https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/forkexec.html
   *    http://www.unix.com/programming/173811-c-execl-pipes.html
   *    http://stackoverflow.com/questions/11042218/c-restore-stdout-to-terminal
   *
   * Set up a communications pipeline between the parent and child processes, 
   * which here are change and initial, respectively.
   */
  ierr = pipe(commpipe);
  if (ierr) {
    perror("pipe");
    sprintf(Message,"error calling pipe()");
    epic_error(dbmsname,Message);
  }

  /*
   * Run initial by forking a child process.
   */
  switch (pid = fork()) {
    case -1:
      sprintf(Message,"error calling fork()");
      epic_error(dbmsname,Message);
    break;
    case 0:
      /*
       * The child process.
       * Connect its stdin to the pipe, then specify the process to be initial.
       */
      dup2(commpipe[0],STDIN_FILENO);  /* set stdin                           */
      close(commpipe[1]);              /* close unused end of pipe            */
      /*
       * Specify the child to be initial.
       */
      ierr = execl(EPIC_PATH"/bin/initial",EPIC_PATH"/bin/initial",NULL);
      if (ierr) {
        perror("execl");
        sprintf(Message,"error calling execl()");
        epic_error(dbmsname,Message);
      }
    break;
    default:
      /*
       * The parent process, which is change.
       */
      saved_stdout = dup(STDOUT_FILENO);
      dup2(commpipe[1],STDOUT_FILENO);       /* set stdout                         */
      close(commpipe[0]);                    /* close unused end of pipe           */
      setvbuf(stdout,(char *)NULL,_IOLBF,0); /* set line-buffered output on stdout */
      /*
       * Send inputs to initial.
       *
       * NOTE: These are not 'smart' and hence will need to be updated if the
       *       inputs to initial are modified.
       */
      printf("Mars\n");
      printf("%d\n",COORD_ISENTROPIC);
      printf("2\n");     /*prompt for starting date values */
      /*
       * Converter between Earth and Mars calendar time:
       * http://www-mars.lmd.jussieu.fr/mars/time/mars_date_to_earth_date.html
       *
       * EMARS data are referenced to the start of Mars Year 24 via macda_sol, which according to the 
       * above website is Julian Date 2451009.27883, give or take.  This corresponds
       * to late in the day UTC on July 14, 1998; the value below is a best guess
       * at this instant in time---an official precise value is elusive
       * to track down.
       */ 
      printf("1998\n");  /* year             */
      printf("7\n");     /* month            */
      printf("14\n");    /* day              */
      printf("18\n");    /* hour             */
      printf("41\n");    /* minute           */
      printf("31\n");    /* second           */

      printf("1\n");                /* initial-wind scaling factor            */
      printf("0\n");                /* radiation scheme                       */
      printf("0\n");                /* turbulence scheme off                  */
      printf("-1\n");               /* sponge off                             */
      printf("%d\n",SPACING_THETA);  /* layer spacing even in potential temp. */

      /*
       * Adjust theta at the model top and bottom for each new application.
       */
      printf("410.0\n");            /* theta at model top [K]                 */
      printf("120.0\n");            /* theta at model bottom [K]              */
      printf("%d\n",grid.nk);       /* number of vertical layers              */
      printf("globe\n");            /* geometry                               */
      printf("-90\n");              /* latbot                                 */
      printf("90\n");               /* lattop                                 */
      printf("-180\n");             /* lonbot                                 */
      printf("180\n");              /* lontop                                 */
      printf("none\n");             /* optional prognostic variables          */
      printf("%d\n",grid.nj);       /* number of latitude gridpoints          */
      printf("%d\n",grid.ni);       /* number of longitude gridpoints         */
      printf("\n");                 /* timestep, use default                  */
      /*
       * NOTE: numerical damping must be turned off if grid.ni is not a
       *       power of 2, which it is not in the EMARS dataset.
       */
      printf("0.\n");               /* divergence damping                     */
      printf("0\n");                /* hyperviscosity                         */
      printf("none\n");             /* extract variables set in initial       */

      wait(&status);                /* wait for initial to end                */

      /*
       * Restore stdout for change.
       */
      dup2(saved_stdout,STDOUT_FILENO);
      close(saved_stdout);
      close(commpipe[1]);
    break;
  }
     
  system("mv ./epic.nc emars_epic.nc");
}

/*======================= end of emars_epic_nc() =============================*/

/*======================= emars_conversion() =================================*/

/*
 * Input EMARS data into EPIC, compute diagnostic fields on isentropic surfaces,
 * and output results.
 */
void emars_conversion(planetspec     *planet,
                      emars_gridspec *emars_grid,
                      char           *emars_infile,
                      char           *emars_outfile_qb,
                      char           *emars_outfile_uvpt,
                      EPIC_FLOAT    **Buff2D)
{
  int
    jj,kk,i,
    K,J,I,
    emars_itime,
    nc_err,
    nc_id_qb,
    nc_id_uvpt,
    dimid[4],
    coorid[4];
  double
    g,
    c2factor,
   *buffer;
  register double
    tmp;
  time_t
    now;
  struct tm
   *today;
  size_t
    nc_index[1],
    nc_start[4], 
    nc_count[4];
  char
    date[16];
  id_information
    L_s_info,
    MY_info,
    pv_info,bsf_info,NH2_info,kin_info,
    u_info,v_info,ug_info,p_info,t_info,phi_info,ma_info;
  EPIC_FLOAT
   *p,
   *h;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="emars_conversion";

  /*
   * Generate current-date string.
   */
  time(&now);
  today = localtime(&now);
  strftime(date,12,"%Y-%m-%d",today);

  /*
   * Create EMARS output .nc files.
   */
  sprintf(emars_outfile_qb,  "./eeQB");
  strcat(emars_outfile_qb,strstr(emars_infile,"_back"));
  nc_err = nc_create(emars_outfile_qb,NC_CLOBBER,&nc_id_qb);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_outfile_qb);
    epic_error(dbmsname,Message);
  }

  sprintf(emars_outfile_uvpt,  "./eeUVPT");
  strcat(emars_outfile_uvpt,strstr(emars_infile,"_back"));
  nc_err = nc_create(emars_outfile_uvpt,NC_CLOBBER,&nc_id_uvpt);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),emars_outfile_uvpt);
    epic_error(dbmsname,Message);
  }

  /*
   * Write global attributes to EMARS outfiles.
   */

  /* QB */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model calculations of Q vs. B on Theta Surfaces; input data are EMARS U, V, P and T");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s: Q, B and NH on isentropic surfaces uses the EPIC GCM v%4.2f",
                  date,grid.epic_version);
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"EMARS: PI Steven J. Greybush (sjg213@psu.edu); "
                   "EPIC: PI Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains potential vorticity, Q, Bernoulli Streamfunction, B, "
                  "and the square of buoyancy frequency times pressure scale height, (NH)^2, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Mars.");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /* UVPT */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model re-gridding of U, V, P and T on Theta Surfaces; input data are EMARS U, V, P and T on hybrid sigma-P Levels");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s: U, V, P and T on isentropic surfaces uses the EPIC GCM v%4.2f",
                  date,grid.epic_version);
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"EMARS: PI Steven J. Greybush (sjg213@psu.edu); "
                   "EPIC: PI Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains zonal and meridional winds, u and v, pressure, p, "
                  "temperature, T, and geopotential, phi, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Mars.");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /*
   * Define coordinates and variables.
   */

  /*
   * lon (I direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_qb,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_uvpt,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /*
   * lat (J direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_qb,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_uvpt,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /*
   * theta (K direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_qb,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_uvpt,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /*
   * time
   */

  /* QB */
  nc_def_dim(nc_id_qb,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_qb,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"long_name",strlen("Time")+1,"Time");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"units",strlen("days since 0000-00-0 00:00:00")+1,"days since 0000-00-0 00:00:00");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_uvpt,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"long_name",strlen("Time")+1,"Time");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"units",strlen("days since 0000-00-0 00:00:00")+1,"days since 0000-00-0 00:00:00");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");

  /*
   * L_s
   */
  L_s_info.name = (char *)malloc(strlen("L_s")+1);
  sprintf(L_s_info.name,"%s","L_s");
  L_s_info.standard_name = (char *)malloc(strlen("solar_longitude")+1);
  sprintf(L_s_info.standard_name,"%s","solar_longitude");
  L_s_info.long_name = (char *)malloc(strlen("Solar longitude")+1);
  sprintf(L_s_info.long_name,"%s","Solar longitude");
  L_s_info.units = (char *)malloc(strlen("degrees")+1);
  sprintf(L_s_info.units,"%s","degrees");
  L_s_info.dimid[NETCDF_T_INDEX] = dimid[NETCDF_T_INDEX];

  /* QB */
  nc_def_var(nc_id_qb,L_s_info.name,NC_DOUBLE,1,&L_s_info.dimid[NETCDF_T_INDEX],&L_s_info.id);
  nc_put_att_text(nc_id_qb,L_s_info.id,"standard_name",strlen(L_s_info.standard_name)+1,L_s_info.standard_name);
  nc_put_att_text(nc_id_qb,L_s_info.id,"long_name",strlen(L_s_info.long_name)+1,L_s_info.long_name);
  nc_put_att_text(nc_id_qb,L_s_info.id,"units",strlen(L_s_info.units)+1,L_s_info.units);

  /* UVPT */
  nc_def_var(nc_id_uvpt,L_s_info.name,NC_DOUBLE,1,&L_s_info.dimid[NETCDF_T_INDEX],&L_s_info.id);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"standard_name",strlen(L_s_info.standard_name)+1,L_s_info.standard_name);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"long_name",strlen(L_s_info.long_name)+1,L_s_info.long_name);
  nc_put_att_text(nc_id_uvpt,L_s_info.id,"units",strlen(L_s_info.units)+1,L_s_info.units);

  /*
   * MY
   */
  MY_info.name = (char *)malloc(strlen("MY")+1);
  sprintf(MY_info.name,"%s","MY");
  MY_info.long_name = (char *)malloc(strlen("Martian year")+1);
  sprintf(MY_info.long_name,"%s","Martian year");
  MY_info.units = (char *)malloc(strlen("1")+1);
  sprintf(MY_info.units,"%s","1");
  MY_info.dimid[NETCDF_T_INDEX] = dimid[NETCDF_T_INDEX];

  /* QB */
  nc_def_var(nc_id_qb,MY_info.name,NC_SHORT,1,&MY_info.dimid[NETCDF_T_INDEX],&MY_info.id);
  nc_put_att_text(nc_id_qb,MY_info.id,"long_name",strlen(MY_info.long_name)+1,MY_info.long_name);
  nc_put_att_text(nc_id_qb,MY_info.id,"units",strlen(MY_info.units)+1,MY_info.units);

  /* UVPT */
  nc_def_var(nc_id_uvpt,MY_info.name,NC_SHORT,1,&MY_info.dimid[NETCDF_T_INDEX],&MY_info.id);
  nc_put_att_text(nc_id_uvpt,MY_info.id,"long_name",strlen(MY_info.long_name)+1,MY_info.long_name);
  nc_put_att_text(nc_id_uvpt,MY_info.id,"units",strlen(MY_info.units)+1,MY_info.units);

  /*
   * QB specific fields
   */

  /*
   * Potential vorticity (PV)
   *
   * NOTE: For isentropic coordinates, and for the isentropic region of hybrid coordinates,
   *       this PV is Ertel's isentropic potential vorticity. For isobaric coordinates, it is
   *       proportional to absolute vorticity on pressure levels, because in that case h is constant.
   */
  pv_info.name = (char *)malloc(strlen("pv")+1);
  sprintf(pv_info.name,"%s","pv");
  pv_info.standard_name = (char *)malloc(strlen("potential_vorticity")+1);
  sprintf(pv_info.standard_name,"%s","potential_vorticity");
  pv_info.long_name = (char *)malloc(strlen("potential vorticity")+1);
  sprintf(pv_info.long_name,"%s","potential vorticity");
  pv_info.units = (char *)malloc(strlen("K m2 kg-1 s-1")+1);
  sprintf(pv_info.units,"%s","K m2 kg-1 s-1");
  for (i = 0; i < 4; i++) {
    pv_info.dimid[ i] = dimid[ i];
    pv_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,pv_info.name,NC_DOUBLE,4,pv_info.dimid,&pv_info.id);
  nc_put_att_text(nc_id_qb,pv_info.id,"standard_name",strlen(pv_info.standard_name)+1,pv_info.standard_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"long_name",strlen(pv_info.long_name)+1,pv_info.long_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"units",strlen(pv_info.units)+1,pv_info.units);

  /*
   * Bernoulli streamfunction
   */
  bsf_info.name = (char *)malloc(strlen("bsf")+1);
  sprintf(bsf_info.name,"%s","bsf");
  bsf_info.standard_name = (char *)malloc(strlen("specific_dry_energy_of_air")+1);
  sprintf(bsf_info.standard_name,"%s","specific_dry_energy_of_air");
  bsf_info.long_name = (char *)malloc(strlen("Bernoulli streamfunction")+1);
  sprintf(bsf_info.long_name,"%s","Bernoulli streamfunction");
  bsf_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(bsf_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    bsf_info.dimid[ i] = dimid[ i];
    bsf_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,bsf_info.name,NC_DOUBLE,4,bsf_info.dimid,&bsf_info.id);
  nc_put_att_text(nc_id_qb,bsf_info.id,"standard_name",strlen(bsf_info.standard_name)+1,bsf_info.standard_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"long_name",strlen(bsf_info.long_name)+1,bsf_info.long_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"units",strlen(bsf_info.units)+1,bsf_info.units);

  /*
   * NH2
   */
  NH2_info.name = (char *)malloc(strlen("NHsquared")+1);
  sprintf(NH2_info.name,"%s","NHsquared");
  NH2_info.long_name = (char *)malloc(strlen("(NH)^2")+1);
  sprintf(NH2_info.long_name,"%s","(NH)^2");
  NH2_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(NH2_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    NH2_info.dimid[ i] = dimid[ i];
    NH2_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,NH2_info.name,NC_DOUBLE,4,NH2_info.dimid,&NH2_info.id);
  nc_put_att_text(nc_id_qb,NH2_info.id,"long_name",strlen(NH2_info.long_name)+1,NH2_info.long_name);
  nc_put_att_text(nc_id_qb,NH2_info.id,"units",strlen(NH2_info.units)+1,NH2_info.units);
  nc_put_att_text(nc_id_qb,NH2_info.id,"comment",strlen(
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H")+1,
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H");

  /*
   * Kinetic energy
   */
  kin_info.name = (char *)malloc(strlen("kin")+1);
  sprintf(kin_info.name,"%s","kin");
  kin_info.standard_name = (char *)malloc(strlen("specific_kinetic_energy_of_air")+1);
  sprintf(kin_info.standard_name,"%s","specific_kinetic_energy_of_air");
  kin_info.long_name = (char *)malloc(strlen("Kinetic energy")+1);
  sprintf(kin_info.long_name,"%s","Kinetic energy");
  kin_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(kin_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    kin_info.dimid[ i] = dimid[ i];
    kin_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,kin_info.name,NC_DOUBLE,4,kin_info.dimid,&kin_info.id);
  nc_put_att_text(nc_id_qb,kin_info.id,"long_name",strlen(kin_info.long_name)+1,kin_info.long_name);
  nc_put_att_text(nc_id_qb,kin_info.id,"units",strlen(kin_info.units)+1,kin_info.units);
  nc_put_att_text(nc_id_qb,kin_info.id,"comment",strlen(
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)")+1,
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)");

  /*
   * UVPT specific fields
   */

  /*
   * Zonal wind
   */
  u_info.name = (char *)malloc(strlen("u")+1);
  sprintf(u_info.name,"%s","u");
  u_info.standard_name = (char *)malloc(strlen("eastward_wind")+1);
  sprintf(u_info.standard_name,"%s","eastward_wind");
  u_info.long_name = (char *)malloc(strlen("Zonal wind")+1);
  sprintf(u_info.long_name,"%s","Zonal wind");
  u_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(u_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    u_info.dimid[ i] = dimid[ i];
    u_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,u_info.name,NC_DOUBLE,4,u_info.dimid,&u_info.id);
  nc_put_att_text(nc_id_uvpt,u_info.id,"standard_name",strlen(u_info.standard_name)+1,u_info.standard_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"long_name",strlen(u_info.long_name)+1,u_info.long_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"units",strlen(u_info.units)+1,u_info.units);

//   /*
//    * Geostrophic zonal wind (assuming variable-f geostrophy)
//    */
//   ug_info.name = (char *)malloc(strlen("ug")+1);
//   sprintf(ug_info.name,"%s","ug");
//   ug_info.standard_name = (char *)malloc(strlen("geostrophic_eastward_wind")+1);
//   sprintf(ug_info.standard_name,"%s","geostrophic_eastward_wind");
//   ug_info.long_name = (char *)malloc(strlen("Geostrophic zonal wind")+1);
//   sprintf(ug_info.long_name,"%s","Geostrophic zonal wind");
//   ug_info.units = (char *)malloc(strlen("m s-1")+1);
//   sprintf(ug_info.units,"%s","m s-1");
//   for (i = 0; i < 4; i++) {
//     ug_info.dimid[ i] = dimid[ i];
//     ug_info.coorid[i] = coorid[i];
//   }
// 
//   nc_def_var(nc_id_uvpt,ug_info.name,NC_DOUBLE,4,ug_info.dimid,&ug_info.id);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"standard_name",strlen(ug_info.standard_name)+1,ug_info.standard_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"long_name",strlen(ug_info.long_name)+1,ug_info.long_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"units",strlen(ug_info.units)+1,ug_info.units);

  /*
   * Meridional wind
   */
  v_info.name = (char *)malloc(strlen("v")+1);
  sprintf(v_info.name,"%s","v");
  v_info.standard_name = (char *)malloc(strlen("northward_wind")+1);
  sprintf(v_info.standard_name,"%s","northward_wind");
  v_info.long_name = (char *)malloc(strlen("Meridional wind")+1);
  sprintf(v_info.long_name,"%s","Meridional wind");
  v_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(v_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    v_info.dimid[ i] = dimid[ i];
    v_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,v_info.name,NC_DOUBLE,4,v_info.dimid,&v_info.id);
  nc_put_att_text(nc_id_uvpt,v_info.id,"standard_name",strlen(v_info.standard_name)+1,v_info.standard_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"long_name",strlen(v_info.long_name)+1,v_info.long_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"units",strlen(v_info.units)+1,v_info.units);

  /*
   * Pressure
   */
  p_info.name = (char *)malloc(strlen("p")+1);
  sprintf(p_info.name,"%s","p");
  p_info.standard_name = (char *)malloc(strlen("air_pressure")+1);
  sprintf(p_info.standard_name,"%s","air_pressure");
  p_info.long_name = (char *)malloc(strlen("Pressure")+1);
  sprintf(p_info.long_name,"%s","Pressure");
  p_info.units = (char *)malloc(strlen("Pa")+1);
  sprintf(p_info.units,"%s","Pa");
  for (i = 0; i < 4; i++) {
    p_info.dimid[ i] = dimid[ i];
    p_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,p_info.name,NC_DOUBLE,4,p_info.dimid,&p_info.id);
  nc_put_att_text(nc_id_uvpt,p_info.id,"standard_name",strlen(p_info.standard_name)+1,p_info.standard_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"long_name",strlen(p_info.long_name)+1,p_info.long_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"units",strlen(p_info.units)+1,p_info.units);

  /*
   * Temperature
   */
  t_info.name = (char *)malloc(strlen("T")+1);
  sprintf(t_info.name,"%s","T");
  t_info.standard_name = (char *)malloc(strlen("air_temperature")+1);
  sprintf(t_info.standard_name,"%s","air_temperature");
  t_info.long_name = (char *)malloc(strlen("Temperature")+1);
  sprintf(t_info.long_name,"%s","Temperature");
  t_info.units = (char *)malloc(strlen("K")+1);
  sprintf(t_info.units,"%s","K");
  for (i = 0; i < 4; i++) {
    t_info.dimid[ i] = dimid[ i];
    t_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,t_info.name,NC_DOUBLE,4,t_info.dimid,&t_info.id);
  nc_put_att_text(nc_id_uvpt,t_info.id,"standard_name",strlen(t_info.standard_name)+1,t_info.standard_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"long_name",strlen(t_info.long_name)+1,t_info.long_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"units",strlen(t_info.units)+1,t_info.units);

  /*
   * Geopotential
   */
  phi_info.name = (char *)malloc(strlen("phi")+1);
  sprintf(phi_info.name,"%s","phi");
  phi_info.standard_name = (char *)malloc(strlen("geopotential")+1);
  sprintf(phi_info.standard_name,"%s","geopotential");
  phi_info.long_name = (char *)malloc(strlen("Geopotential")+1);
  sprintf(phi_info.long_name,"%s","Geopotential");
  phi_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(phi_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    phi_info.dimid[ i] = dimid[ i];
    phi_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,phi_info.name,NC_DOUBLE,4,phi_info.dimid,&phi_info.id);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"standard_name",strlen(phi_info.standard_name)+1,phi_info.standard_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"long_name",strlen(phi_info.long_name)+1,phi_info.long_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"units",strlen(phi_info.units)+1,phi_info.units);

  /*
   * Mach number
   */
  ma_info.name = (char *)malloc(strlen("Ma")+1);
  sprintf(ma_info.name,"%s","Ma");
  ma_info.long_name = (char *)malloc(strlen("Mach number")+1);
  sprintf(ma_info.long_name,"%s","Mach number");
  ma_info.units = (char *)malloc(strlen("1")+1);
  sprintf(ma_info.units,"%s","1");
  for (i = 0; i < 4; i++) {
    ma_info.dimid[ i] = dimid[ i];
    ma_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,ma_info.name,NC_DOUBLE,4,ma_info.dimid,&ma_info.id);
  nc_put_att_text(nc_id_uvpt,ma_info.id,"long_name",strlen(ma_info.long_name)+1,ma_info.long_name);
  nc_put_att_text(nc_id_uvpt,ma_info.id,"units",strlen(ma_info.units)+1,ma_info.units);

  /*---------------------------*
   * Leave netcdf define mode: *
   *---------------------------*/
  nc_enddef(nc_id_qb);
  nc_enddef(nc_id_uvpt);

  /*
   * Assign values to X, Y and Z coordinates.
   */
  /* lon */
  for (I = ILO; I <= IHI; I++) {
    nc_index[0] = I-ILO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
  }

  /* lat */
  for (J = JLO; J <= JHI; J++) {
    nc_index[0] = J-JLO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
  }

  /* theta */
  for (K = KHI; K >= KLO; K--) {
    nc_index[0] = KHI-K;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
  }

  fprintf(stdout,"Computing variables on theta surfaces for EMARS QB and UVPT files...  0%%");

  /* Allocate memory */
  p      = fvector(0,2*grid.nk+1,dbmsname);
  h      = fvector(0,2*grid.nk+1,dbmsname);
  buffer = (double *)calloc(grid.ni,sizeof(double));

  for (emars_itime = 0; emars_itime < emars_grid->ntime; emars_itime++) {
    /* Show progress */
    fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(double)(emars_itime+1)/emars_grid->ntime));
    fflush(stdout);

    /*
     * Input EMARS data for timeframe emars_itime, and interpolate into EPIC variables.
     */
    emars_var_read(planet,emars_grid,emars_infile,VAR_DATA,emars_itime);

    /*
     * Write time.
     */

    /* QB */
    nc_index[0] = emars_itime;
    nc_put_var1_double(nc_id_qb,coorid[NETCDF_T_INDEX],nc_index,         &emars_grid->time[emars_itime]);
    nc_put_var1_double(nc_id_qb,L_s_info.id,           nc_index,         &emars_grid->Ls[  emars_itime]);
    nc_put_var1_short( nc_id_qb,MY_info.id,            nc_index,(short *)&emars_grid->MY[  emars_itime]);

    /* UVPT */
    nc_index[0] = emars_itime;
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_T_INDEX],nc_index,         &emars_grid->time[emars_itime]);
    nc_put_var1_double(nc_id_uvpt,L_s_info.id,           nc_index,         &emars_grid->Ls[  emars_itime]);
    nc_put_var1_short( nc_id_uvpt,MY_info.id,            nc_index,(short *)&emars_grid->MY[  emars_itime]);

    /*
     * Compute Q and B on theta surfaces by running the normal EPIC model
     * diagnostic-variable calculations.
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      for (I = ILO; I <= IHI; I++) {
        for (kk = 1; kk <= 2*KHI+1; kk++) {
          p[kk] = get_p(planet,P2_INDEX,kk,J,I);
        }
        calc_h(jj,p,h);

        for (K = KLO; K <= KHI; K++) {
          H(K,J,I) = h[2*K];
        }
        K = 0;
        H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
        K = KHI+1;
        H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
      }
    }
    bc_lateral(var.h.value,THREEDIM);

    /*
     * PHI3(KHI,J,I) is set above by emars_var_read(VAR_DATA).
     */
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,PASSING_PHI3NK);
    store_diag(planet);

    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;
    nc_count[NETCDF_I_INDEX] = grid.ni;

    nc_start[NETCDF_T_INDEX] = emars_itime;
    nc_start[NETCDF_I_INDEX] = 0;

    for (K = KHI; K >= KLO; K--) {
      kk                       = 2*K;
      nc_start[NETCDF_K_INDEX] = KHI-K;
      for (J = JLO; J <= JHI; J++) {
        jj = 2*J;

        g  = GRAVITY2(K,J);

        nc_start[NETCDF_J_INDEX] = J-JLO;

        /*--------------------*
         * QB specific fields *
         *--------------------*/

        /*
         * Write potential vorticity (pv, aka Q).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .25*(PV2(K,J,I)+PV2(K,J,I+1)+PV2(K,J+1,I)+PV2(K,J+1,I+1));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,pv_info.id,nc_start,nc_count,buffer);

        /*
         * Write Bernoulli streamfunction (bsf).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = MONT2(K,J,I)
                         +get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,bsf_info.id,nc_start,nc_count,buffer);

        /*
         * Write (NH)^2, the square of the product of the buoyancy (Brunt-Vaisala) frequency, N, 
         * and the pressure scale height, H = RT/g.
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_brunt2(planet,kk,J,I)
                         *SQR(planet->rgas*T2(K,J,I)/g);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,NH2_info.id,nc_start,nc_count,buffer);

        /*
         * Write kinetic energy per mass (kin), .5*(u*u+v*v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_qb,kin_info.id,nc_start,nc_count,buffer);

        /*----------------------*
         * UVPT specific fields *
         *----------------------*/

        /*
         * Write zonal wind (u).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(U(grid.it_uv,K,J,I)+U(grid.it_uv,K,J,I+1));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,u_info.id,nc_start,nc_count,buffer);

//         /*
//          * Calculate and write geostrophic zonal wind (ug).
//          */
//         if (J == JLO) {
//           /*
//            * Southern edge: take forward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         else if (J == JHI) {
//           /*
//            * Northern edge: take backward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J,I)-MONT2(K,J-1,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         else {
//           /*
//            * Interior point: take long step to center back onto h-grid.
//            */
//           tmp = -.5*grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J-1,I));
//             if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
//               /* Flag underground value with NAN */
//               buffer[I-ILO] = NAN;
//             }
//           }
//         }
//         nc_put_vara_double(nc_id_uvpt,ug_info.id,nc_start,nc_count,buffer);

        /*
         * Write meridional wind (v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(V(grid.it_uv,K,J,I)+V(grid.it_uv,K,J+1,I));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,v_info.id,nc_start,nc_count,buffer);

        /*
         * Write pressure (p).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = P2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,p_info.id,nc_start,nc_count,buffer);

        /*
         * Write temperature (T).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = T2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,t_info.id,nc_start,nc_count,buffer);

        /*
         * Write geopotential (phi).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = PHI2(K,J,I);
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,phi_info.id,nc_start,nc_count,buffer);

        /*
         * Write Mach number, Ma.
         */
        for (I = ILO; I <= IHI; I++) {
          c2factor = return_cp(planet,0.,0.,T2(K,J,I));
          c2factor = planet->rgas*c2factor/(c2factor-planet->rgas);
          buffer[I-ILO] = sqrt(2.*get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                                 var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I)
                               /(T2(K,J,I)*c2factor));
          if (grid.sigmatheta[2*K+1] < THETA(K,J,I)) {
            /* Flag underground value with NAN */
            buffer[I-ILO] = NAN;
          }
        }
        nc_put_vara_double(nc_id_uvpt,ma_info.id,nc_start,nc_count,buffer);
      } /* J loop */
    } /* K loop */
  } /* emars_itime loop */

  nc_close(nc_id_qb);
  nc_close(nc_id_uvpt);

  fprintf(stdout,"\b\b\b\b100%%\n");
  fflush(stdout);
  fprintf(stdout,"Output written to %s\n",emars_outfile_qb);

  /*
   * The conversion is complete, so clean up (write default values,
   * free allocated memory, etc.) and exit.
   */
  free_fvector(p,0,2*grid.nk+1,dbmsname);
  free_fvector(h,0,2*grid.nk+1,dbmsname);
  free(buffer);
}

/*======================= end of emars_conversion() ==========================*/

/*======================= weizmann_var_read() ================================*/

/*
 * Read in Weizmann Institute (Jupiter) data.
 * portion = SIZE_DATA:       read in grid dimensions
 * portion = POST_SIZE_DATA:  read in constant parameters
 * portion = VAR_DATA:        read in variables
 */

#define DBUFF3D(k,j,i)         dbuffer[k+weizmann_grid->nk*(j+weizmann_grid->nj*(i))]
#define EPIC_K_WEIZMANN_j(K,j) epic_k_weizmann_j[j+weizmann_grid->nj*(K-KLOPAD)]

void weizmann_var_read(planetspec        *planet,
                       weizmann_gridspec *weizmann_grid,
                       char              *weizmann_infile,
                       int                portion,
                       int                itime)
{
  int
    K,J,I,
    kk,k,j,i,
    is,
    nelem2d,nelem3d;
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  int
    nc_err,nc_id,
    nc_dimid,nc_varid;
  double
    *dbuffer;
  size_t
    dimlen;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="weizmann_var_read";

  nc_err = lookup_netcdf(weizmann_infile,&nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_infile);
    epic_error(dbmsname,Message);
  }

  /*
   * Store netcdf file history for later use.
   */
  /*
  READC(weizmann_grid->infile_history,history,128);
  */
  sprintf(weizmann_grid->infile_history,"[No input-file history]");

  if (portion == SIZE_DATA) {
    fprintf(stdout,"Reading SIZE_DATA from %s \n",weizmann_infile);
  }
  else if (portion == POST_SIZE_DATA) {
    fprintf(stdout,"Reading POST_SIZE_DATA from %s \n",weizmann_infile);
  }
  else if (portion == VAR_DATA) {
    ;
  }
  else {
    sprintf(Message,"unrecognized portion = %d",portion);
    epic_error(dbmsname,Message);
  }
  fflush(stderr);

  /*
   * Read in size of model and set parameters needed by weizmann_make_arrays().
   */
  if (portion == SIZE_DATA) {
    /* 
     * Number of longitude points
     */
    weizmann_grid->ni = 1;

    /* 
     * Number of latitude points
     */
    nc_err = nc_inq_dimid(nc_id,"Lat",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_infile);
      epic_error(dbmsname,Message);
    }
    weizmann_grid->nj = (int)dimlen;

    /* 
     * Number of pressure levels.
     */
    nc_err = nc_inq_dimid(nc_id,"P",&nc_dimid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_infile);
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_dimlen(nc_id,nc_dimid,&dimlen);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_infile);
      epic_error(dbmsname,Message);
    }
    weizmann_grid->nk = (int)dimlen;

    /* 
     * Number of times.
     */
    weizmann_grid->ntime = 1;
  }
  else if (portion == POST_SIZE_DATA) {
    /*
     * NOTE: POST_SIZE_DATA is interpreted here to mean constant parameters,
     *       including the constant 1D dimension arrays, Lat, P, etc.
     */

    /*-------------------------*
     * Read 1D lon array [deg] *
     *-------------------------*/

    i = 0;
    weizmann_grid->lon[i] = (EPIC_FLOAT)0.;

    grid.ni = weizmann_grid->ni;

    /*-------------------------*
     * Read 1D Lat array [deg] *
     *-------------------------*/

    dbuffer = calloc(weizmann_grid->nj,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"Lat",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Weizmann latitudes are stored in the same index order as the EPIC convention.
     */
    for (j = 0; j < weizmann_grid->nj; j++) {
       weizmann_grid->lat[j] = (EPIC_FLOAT)dbuffer[j];
    }

    free(dbuffer);

    grid.nj = NINT(180./(weizmann_grid->lat[1]-weizmann_grid->lat[0])-1.-sqrt(2.));

    /*---------------------------*
     * Read 1D pressure array, P *
     *---------------------------*/

    dbuffer = calloc(weizmann_grid->nk,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }

    nc_err = nc_inq_varid(nc_id,"P",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_var_double(nc_id,nc_varid,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    /*
     * NOTE: Weizmann P are stored in the same index order as the EPIC convention.
     * The units are converted from [bar] to [Pa].
     */
    for (k = 0; k < weizmann_grid->nk; k++) {
      weizmann_grid->p[k] = (EPIC_FLOAT)dbuffer[k]*1.e+5;
    }

    free(dbuffer);

    /*
     * We hardwire grid.nk.
     */
    grid.nk = 8;

    /*---------------------*
     * Read 1D time array  *
     *---------------------*/

    weizmann_grid->time[0] = 0.;
  }
  else if (portion == VAR_DATA) {
    int
      count,
      excluded_count;
    EPIC_FLOAT
      p,tmp,x,x_d,
      phi0,mont0,
      fgibb,fpe,uoup,
     *u1d,
     *mont1d,
     *epic_k_weizmann_j;
    size_t
      nc_start3d[2],
      nc_count3d[2];
    float_triplet
      *k_triplet,
      *j_triplet,
      *i_triplet,
      *epic_j_triplet;
    FILE
      *excluded_file;

    /*
     * "3d" is used for the meridional-plane data
     */
    nelem3d       = weizmann_grid->nj*weizmann_grid->nk;
    nc_start3d[0] = 0;
    nc_start3d[1] = 0;
    nc_count3d[0] = weizmann_grid->nj;
    nc_count3d[1] = weizmann_grid->nk;

    /*------------------*
     * Read T array [K] *
     *------------------*/

    dbuffer = calloc(nelem3d,sizeof(double));
    if (!dbuffer) {
      sprintf(Message,"calloc error allocating dbuffer");
      epic_error(dbmsname,Message);
    }
    nc_err = nc_inq_varid(nc_id,"T",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_double(nc_id,nc_varid,nc_start3d,nc_count3d,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < weizmann_grid->nk; k++) {
      for (j = 0; j < weizmann_grid->nj; j++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          WEIZMANN_TEMP(k,j,i) = (EPIC_FLOAT)DBUFF3D(k,j,i);
        }
      }
    }

    /*--------------------*
     * Read U array [m/s] *
     *--------------------*/

    nc_err = nc_inq_varid(nc_id,"U",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_double(nc_id,nc_varid,nc_start3d,nc_count3d,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < weizmann_grid->nk; k++) {
      for (j = 0; j < weizmann_grid->nj; j++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          WEIZMANN_U(k,j,i) = (EPIC_FLOAT)DBUFF3D(k,j,i);
        }
      }
    }

    /*--------------------*
     * Read V array [m/s] *
     *--------------------*/

    for (k = 0; k < weizmann_grid->nk; k++) {
      for (j = 0; j < weizmann_grid->nj; j++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          WEIZMANN_V(k,j,i) = 0.;
        }
      }
    }

    /*-------------------------*
     * Read Rho array [kg/m^3] *
     *-------------------------*/

    nc_err = nc_inq_varid(nc_id,"Rho",&nc_varid);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    nc_err = nc_get_vara_double(nc_id,nc_varid,nc_start3d,nc_count3d,dbuffer);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"%s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }

    for (k = 0; k < weizmann_grid->nk; k++) {
      for (j = 0; j < weizmann_grid->nj; j++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          WEIZMANN_RHO(k,j,i) = (EPIC_FLOAT)DBUFF3D(k,j,i);
        }
      }
    }

    /*---------------------------------------------------------*
     * Compute theta array, 3D potential temperature field [K] *
     *                                                         *
     *---------------------------------------------------------*/

    for (k = 0; k < weizmann_grid->nk; k++) {
      for (j = 0; j < weizmann_grid->nj; j++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          WEIZMANN_THETA(k,j,i) = WEIZMANN_TEMP(k,j,i)*pow(planet->p0/weizmann_grid->p[k],planet->kappa);
        }
      }
    }

    /*------------------------------------------------------*
     * Interpolate Weizmann data onto EPIC staggered grids. *
     *------------------------------------------------------*/

    /* Allocate memory */
    k_triplet         = ftriplet(0,weizmann_grid->nk-1,                  dbmsname);
    j_triplet         = ftriplet(0,weizmann_grid->nj-1,                  dbmsname);
    i_triplet         = ftriplet(0,weizmann_grid->ni,                    dbmsname);
    epic_j_triplet    = ftriplet(0,grid.nj,                              dbmsname);
    epic_k_weizmann_j = fvector( 0,(KHIPAD-KLOPAD+1)*weizmann_grid->nj-1,dbmsname);

    /*---------------------------------------------------------------------*
     * Potential temperature, THETA.                                       *
     * Using isentropic coordinates, so THETA takes the coordinate values. *
     *---------------------------------------------------------------------*/

    for (K = KLOPAD; K <= KHIPAD; K++) {
      tmp = grid.sigmatheta[2*K+1];
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          THETA(K,J,I) = tmp;
        }
      }
    }
    /* No need to apply bc_lateral() here */

    /*--------------------------------------------*
     * Zonal wind, U                              *
     * Longitude can be done by direct insertion. * 
     *--------------------------------------------*/

    excluded_file = fopen("excluded_thetas.dat","w");
    if (!excluded_file) {
      sprintf(Message,"file %s not found","excluded_thetas.dat");
      epic_error(dbmsname,Message);
    }
    fprintf(excluded_file," Points in column excluded because of lack of monotonicity of theta.\n");
    fprintf(excluded_file,"  Lon    Lat  Excluded thetas [K]\n");

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping Weizmann latitudes.
       */
      for (j = 0; j < weizmann_grid->nj; j++) {
        k              = weizmann_grid->nk-1;
        k_triplet[0].x = WEIZMANN_THETA(k,j,i);
        k_triplet[0].y = WEIZMANN_U(k,j,i);
        count          = 1;
        excluded_count = 0;
        for (k = weizmann_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (WEIZMANN_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = WEIZMANN_THETA(k,j,i);
            k_triplet[count].y = WEIZMANN_U(k,j,i);
            count++;
          }
          else {
            /*
             * Keep track of number of excluded points caused by non-monotonicity of theta.
             */
            excluded_count++;
            if (excluded_count == 1) {
              fprintf(excluded_file," %6.1f %5.1f  %5.1f",weizmann_grid->lon[i],weizmann_grid->lat[j],WEIZMANN_THETA(k,j,i));
            }
            else {
              fprintf(excluded_file," %5.1f",WEIZMANN_THETA(k,j,i));
            }
          }
        }

        if (excluded_count > 0) {
          fprintf(excluded_file,"\n");
        }

        spline_pchip(count,k_triplet);

        k = -2;
        for (K = KLO; K <= KHIPAD; K++) {
          x                      = grid.sigmatheta[2*K];
          k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
          EPIC_K_WEIZMANN_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       *
       * NOTE: Pole values are included in the Weizmann latitude range.
       */
      for (K = KLO; K <= KHIPAD; K++) {
        for (j = 0; j < weizmann_grid->nj; j++) {
          j_triplet[j].x = weizmann_grid->lat[j];
          j_triplet[j].y = EPIC_K_WEIZMANN_j(K,j);
        }
        spline_pchip(weizmann_grid->nj,j_triplet);
        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x                   = grid.lat[2*J+1];
          j                   = hunt_place_in_table(weizmann_grid->nj,j_triplet,x,&x_d,j);
          U(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);

    fclose(excluded_file);

    /*-------------------------------------------------------------------*
     * Meridonal wind, V                                                 *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping Weizmann latitudes.
       */
      for (j = 0; j < weizmann_grid->nj; j++) {
        k              = weizmann_grid->nk-1;
        k_triplet[0].x = WEIZMANN_THETA(k,j,i);
        k_triplet[0].y = WEIZMANN_V(k,j,i);
        count          = 1;
        for (k = weizmann_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (WEIZMANN_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = WEIZMANN_THETA(k,j,i);
            k_triplet[count].y = WEIZMANN_V(k,j,i);
            count++;
          }
        }
        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLO; K <= KHIPAD; K++) {
          x                      = grid.sigmatheta[2*K];
          k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
          EPIC_K_WEIZMANN_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       * We include V = 0 pole points.
       */
      for (K = KLO; K <= KHIPAD; K++) {
        for (j = 0; j < weizmann_grid->nj; j++) {
          j_triplet[j].x = weizmann_grid->lat[j];
          j_triplet[j].y = EPIC_K_WEIZMANN_j(K,j);
        }
        spline_pchip(weizmann_grid->nj,j_triplet);
        j = -2;
        for (J = JFIRST; J <= JHI; J++) {
          x                   = grid.lat[2*J];
          j                   = hunt_place_in_table(weizmann_grid->nj,j_triplet,x,&x_d,j);
          V(grid.it_uv,K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLO; K <= KHIPAD; K++) {
      for (J = JFIRST; J <= JHI; J++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = weizmann_grid->lon[i];
          i_triplet[i].y = V(grid.it_uv,K,J,I);
        }
        if (weizmann_grid->ni > 1) {
          i_triplet[i].x = weizmann_grid->lon[i-1]+weizmann_grid->lon[1]-weizmann_grid->lon[0];
        }
        else {
          /* There is no longitude dimension, so just put in a unit increment. */
          i_triplet[i].x = weizmann_grid->lon[i-1]+1.;
        }
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(weizmann_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x                   = grid.lon[2*I+1];
          i                   = hunt_place_in_table(weizmann_grid->ni+1,i_triplet,x,&x_d,i);
          V(grid.it_uv,K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

    /*-------------------------------------------------------------------*
     * Pressure on the layer interfaces, P3.                             *
     * Longitude needs to be staggered, but start with direct insertion. *
     *-------------------------------------------------------------------*/

    for (I = ILO; I <= IHI; I++) {
      i = I-1;
      /*
       * First, interpolate onto EPIC theta surfaces, keeping Weizmann latitudes.
       */
      for (j = 0; j < weizmann_grid->nj; j++) {
        k              = weizmann_grid->nk-1;
        k_triplet[0].x = WEIZMANN_THETA(k,j,i);
        k_triplet[0].y = weizmann_grid->p[k];
        count          = 1;
        for (k = weizmann_grid->nk-2; k >= 0; k--) {
          /*
           * Following Du, Dowling and Bradley (2015, JAMES), we exclude from
           * the source interpolation table points where theta is not larger
           * than the value on the level below.
           */
          if (WEIZMANN_THETA(k,j,i) > k_triplet[count-1].x) {
            k_triplet[count].x = WEIZMANN_THETA(k,j,i);
            k_triplet[count].y = weizmann_grid->p[k];
            count++;
          }
        }
        spline_pchip(count,k_triplet);
        k = -2;
        for (K = KLOPAD; K <= KHI; K++) {
          x                      = grid.sigmatheta[2*K+1];
          k                      = hunt_place_in_table(count,k_triplet,x,&x_d,k);
          EPIC_K_WEIZMANN_j(K,j) = splint_pchip(x,k_triplet+k,x_d);
        }
      }

      /*
       * Next, interpolate onto EPIC latitudes.
       */
      for (K = KLOPAD; K <= KHI; K++) {
        for (j = 0; j < weizmann_grid->nj; j++) {
          j_triplet[j].x = weizmann_grid->lat[j];
          j_triplet[j].y = EPIC_K_WEIZMANN_j(K,j);
        }
        spline_pchip(weizmann_grid->nj,j_triplet);
        j = -2;
        for (J = JLO; J <= JHI; J++) {
          x         = grid.lat[2*J+1];
          j         = hunt_place_in_table(weizmann_grid->nj,j_triplet,x,&x_d,j);
          P3(K,J,I) = splint_pchip(x,j_triplet+j,x_d);
        }
      }
    }

    /*
     * Finally, interpolate onto staggered longitude grid.
     */
    for (K = KLOPAD; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (i = 0; i < weizmann_grid->ni; i++) {
          I              = i+1;
          i_triplet[i].x = weizmann_grid->lon[i];
          i_triplet[i].y = P3(K,J,I);
        }
        if (weizmann_grid->ni > 1) {
          i_triplet[i].x = weizmann_grid->lon[i-1]+weizmann_grid->lon[1]-weizmann_grid->lon[0];
        }
        else {
          /* There is no longitude dimension, so just put in a unit increment. */
          i_triplet[i].x = weizmann_grid->lon[i-1]+1.;
        }
        i_triplet[i].y = i_triplet[0].y;
        periodic_spline_pchip(weizmann_grid->ni+1,i_triplet);
        i = -2;
        for (I = ILO; I <= IHI; I++) {
          x         = grid.lon[2*I+1];
          i         = hunt_place_in_table(weizmann_grid->ni+1,i_triplet,x,&x_d,i);
          P3(K,J,I) = splint_pchip(x,i_triplet+i,x_d);
        }
      }
    }
    /* Need to call bc_lateral() here */ 
    bc_lateral(var.p3.value,THREEDIM);

    /*--------------------------------------------------------------------*
     * Geopotential at the bottom of the model, PHI3(KHI,J,I). This is    *
     * the geopotential on the bottom isentropic surface, grid.thetabot.  *
     * Integrate the gradient-balance equation assuming zonal symmetry.   *
     *--------------------------------------------------------------------*/

    /* Allocate memory */
    u1d    = fvector(0,JADIM-1,dbmsname);
    mont1d = fvector(0,JADIM-1,dbmsname);

    for (J = JLOPAD; J <= JHIPAD; J++) {
      /*
       * Calculate zonal averages.
       */
      U1D(J) = 0.;
      for (I = ILO; I <= IHI; I++) {
        U1D(J) += U(grid.it_uv,KHI+1,J,I);
      }
      U1D(J) /= grid.ni;
    }

    K      = KHI;
    phi0   = 0.;
    mont0 = phi0+planet->cp*grid.t_ref[2*K+1];

    mont_from_u(planet,2*KHI+1,u1d,mont1d,grid.jtp,mont0);

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        T3(K,J,I)   = THETA(K,J,I)*pow(P3(K,J,I)/planet->p0,planet->kappa);
        PHI3(K,J,I) = MONT1D(J)-planet->cp*T3(K,J,I);
      }
    }
    /* Need to call bc_lateral() here */
    bc_lateral(var.phi3.value+(K-Kshift)*Nelem2d,TWODIM);

    /* Free allocated memory */
    free_ftriplet(k_triplet,       0,weizmann_grid->nk-1,            dbmsname);
    free_ftriplet(j_triplet,       0,weizmann_grid->nj-1,            dbmsname);
    free_ftriplet(i_triplet,       0,weizmann_grid->ni,              dbmsname);
    free_ftriplet(epic_j_triplet,  0,grid.nj,                        dbmsname);
    free_fvector(epic_k_weizmann_j,0,(KHIPAD-KLOPAD+1)*weizmann_grid->nj-1,dbmsname);
    free_fvector(u1d,              0,JADIM-1,dbmsname);
    free_fvector(mont1d,           0,JADIM-1,dbmsname);
  } /* end portion == VAR_DATA */

  return;
}

#undef DBUFF3D
#undef EPIC_K_WEIZMANN_j

/*======================= end of weizmann_var_read() =========================*/

/*======================= weizmann_make_arrays() =============================*/

void weizmann_make_arrays(planetspec        *planet,
                          weizmann_gridspec *weizmann_grid)
{
  int
    nelem2d,nelem3d;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="weizmann_make_arrays";

  nelem2d = weizmann_grid->nj*weizmann_grid->ni;
  nelem3d = nelem2d*weizmann_grid->nk;

  weizmann_grid->lon       = fvector(0,weizmann_grid->ni-1,   dbmsname);
  weizmann_grid->lat       = fvector(0,weizmann_grid->nj-1,   dbmsname);
  weizmann_grid->p         = fvector(0,weizmann_grid->nk-1,   dbmsname);
  weizmann_grid->time      = fvector(0,weizmann_grid->ntime-1,dbmsname);
  weizmann_grid->temp      = fvector(0,nelem3d-1,             dbmsname);
  weizmann_grid->u         = fvector(0,nelem3d-1,             dbmsname);
  weizmann_grid->v         = fvector(0,nelem3d-1,             dbmsname);
  weizmann_grid->theta     = fvector(0,nelem3d-1,             dbmsname);
  weizmann_grid->rho       = fvector(0,nelem3d-1,             dbmsname);

  return;
}

/*======================= end of weizmann_make_arrays() =====================*/

/*======================= weizmann_free_arrays() ============================*/

void weizmann_free_arrays(planetspec        *planet,
                          weizmann_gridspec *weizmann_grid)
{
  int
    nelem2d,nelem3d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="weizmann_free_arrays";

  nelem2d = weizmann_grid->ni*weizmann_grid->nj;
  nelem3d = nelem2d*weizmann_grid->nk;

  free_fvector(weizmann_grid->lon,  0,weizmann_grid->ni-1,   dbmsname);
  free_fvector(weizmann_grid->lat,  0,weizmann_grid->nj-1,   dbmsname);
  free_fvector(weizmann_grid->p,    0,weizmann_grid->nk-1,   dbmsname);
  free_fvector(weizmann_grid->time, 0,weizmann_grid->ntime-1,dbmsname);
  free_fvector(weizmann_grid->temp, 0,nelem3d-1,             dbmsname);
  free_fvector(weizmann_grid->u,    0,nelem3d-1,             dbmsname);
  free_fvector(weizmann_grid->v,    0,nelem3d-1,             dbmsname);
  free_fvector(weizmann_grid->theta,0,nelem3d-1,             dbmsname);
  free_fvector(weizmann_grid->rho,  0,nelem3d-1,             dbmsname);

  return;
}

/*======================= end of weizmann_free_arrays() ======================*/

/*======================= weizmann_epic_nc() =================================*/

/*
 * Bootstrap weizmann_epic.nc by running initial.
 */

void weizmann_epic_nc(planetspec *planet)
{
  pid_t
    pid;
  int
    commpipe[2],
    ierr,status,
    saved_stdout;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="weizmann_epic_nc";

  /*
   * The code for executing initial from change is adapted from examples at following websites:
   *    https://www.gidforums.com/t-3369.html
   *    https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/forkexec.html
   *    http://www.unix.com/programming/173811-c-execl-pipes.html
   *    http://stackoverflow.com/questions/11042218/c-restore-stdout-to-terminal
   *
   * Set up a communications pipeline between the parent and child processes, 
   * which here are change and initial, respectively.
   */
  ierr = pipe(commpipe);
  if (ierr) {
    perror("pipe");
    sprintf(Message,"error calling pipe()");
    epic_error(dbmsname,Message);
  }

  /*
   * Run initial by forking a child process.
   */
  switch (pid = fork()) {
    case -1:
      sprintf(Message,"error calling fork()");
      epic_error(dbmsname,Message);
    break;
    case 0:
      /*
       * The child process.
       * Connect its stdin to the pipe, then specify the process to be initial.
       */
      dup2(commpipe[0],STDIN_FILENO);  /* set stdin                           */
      close(commpipe[1]);              /* close unused end of pipe            */
      /*
       * Specify the child to be initial.
       */
      ierr = execl(EPIC_PATH"/bin/initial",EPIC_PATH"/bin/initial",NULL);
      if (ierr) {
        perror("execl");
        sprintf(Message,"error calling execl()");
        epic_error(dbmsname,Message);
      }
    break;
    default:
      /*
       * The parent process, which is change.
       */
      saved_stdout = dup(STDOUT_FILENO);
      dup2(commpipe[1],STDOUT_FILENO);       /* set stdout                         */
      close(commpipe[0]);                    /* close unused end of pipe           */
      setvbuf(stdout,(char *)NULL,_IOLBF,0); /* set line-buffered output on stdout */
      /*
       * Send inputs to initial.
       *
       * NOTE: These are not 'smart' and hence will need to be updated if the
       *       inputs to initial are modified.
       */
      printf("Jupiter\n");
      printf("%d\n",COORD_ISENTROPIC);
      printf("0\n");                     /*default starting date and time          */
      printf("1\n");                     /* initial-wind scaling factor            */
      printf("1.0\n");                   /* Vertical wind profile factor           */
      printf("0\n");                     /* radiation scheme off                   */
      printf("0\n");                     /* turbulence scheme off                  */
      printf("-1\n");                    /* sponge off                             */
      printf("%d\n",SPACING_LOGP);       /* layer spacing even in log p            */
      printf("100.\n");                  /* pressure at model top [hPa]            */
      printf("1600.\n");                 /* pressure at the model bottom [hPa]     */
      printf("%d\n",grid.nk);            /* number of vertical layers              */
      printf("globe\n");                 /* geometry                               */
      printf("-90\n");                   /* latbot                                 */
      printf("90\n");                    /* lattop                                 */
      printf("-180\n");                  /* lonbot                                 */
      printf("180\n");                   /* lontop                                 */
      printf("none\n");                  /* optional prognostic variables          */
      printf("%d\n",grid.nj);            /* number of latitude gridpoints          */
      printf("%d\n",grid.ni);            /* number of longitude gridpoints         */
      printf("\n");                      /* timestep, use default                  */
      printf("0.");                      /* lat. to apply sounding, use 0deg       */
      printf("0.\n");                    /* divergence damping                     */
      printf("0\n");                     /* hyperviscosity                         */
      printf("none\n");                  /* extract variables set in initial       */
      printf("\n");                      /* This extra return prevents hanging.    */

      wait(&status);                     /* wait for initial to end                */

      /*
       * Restore stdout for change.
       */
      dup2(saved_stdout,STDOUT_FILENO);
      close(saved_stdout);
      close(commpipe[1]);
    break;
  }
     
  system("mv ./epic.nc weizmann_epic.nc");
}

/*======================= end of weizmann_epic_nc() ==========================*/

/*======================= weizmann_conversion() ==============================*/

/*
 * Input Weizmann Institute data into EPIC, compute diagnostic fields on isentropic surfaces,
 * and output results.
 */

void weizmann_conversion(planetspec        *planet,
                         weizmann_gridspec *weizmann_grid,
                         char              *weizmann_infile,
                         char              *weizmann_outfile_qb,
                         char              *weizmann_outfile_uvpt,
                         EPIC_FLOAT       **Buff2D)
{
  int
    kk,jj,i,
    K,J,I,
    weizmann_itime,
    nc_err,
    nc_id_qb,
    nc_id_uvpt,
    dimid[4],
    coorid[4];
  double
    g,
    c2factor,
   *buffer;
  register double
    tmp;
  EPIC_FLOAT
   *p,
   *h;
  time_t
    now;
  struct tm
   *today;
  size_t
    nc_index[1],
    nc_start[4], 
    nc_count[4],
    message_len;
  char
    date[16];
  id_information
    pv_info,bsf_info,NH2_info,kin_info,
    u_info,ug_info,v_info,p_info,t_info,phi_info;

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="weizmann_conversion";

  /*
   * Generate current-date string.
   */
  time(&now);
  today = localtime(&now);
  strftime(date,12,"%Y-%m-%d",today);

  /*
   * Create weizmann output .nc files.
   */
  sprintf(weizmann_outfile_qb,"./QB_weizmann.nc");
  nc_err = nc_create(weizmann_outfile_qb,NC_CLOBBER,&nc_id_qb);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_outfile_qb);
    epic_error(dbmsname,Message);
  }

  sprintf(weizmann_outfile_uvpt,  "./UVPT_weizmann.nc");
  nc_err = nc_create(weizmann_outfile_uvpt,NC_CLOBBER,&nc_id_uvpt);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),weizmann_outfile_uvpt);
    epic_error(dbmsname,Message);
  }

  /*
   * Write global attributes to weizmann outfiles.
   */

  /* QB */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model calculations of Q vs. B on Theta Surfaces; input data are Weizmann Institute U and T and Pressure Levels");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s; \n"
                  "%s: Q, B and NH on isentropic surfaces uses the EPIC GCM v%4.2f",
                  weizmann_grid->infile_history,date,grid.epic_version);
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"Weizmann Institute: Eli Galanti (eli.galanti@weizmann@ac.il), Yohai Kaspi (yohai.kaspi@weizmann.ac.il), "
                  "EPIC: Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains potential vorticity, Q, Bernoulli Streamfunction, B, "
                  "and the product of buoyancy frequency times pressure scale height, NH, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Jupiter.");
  nc_put_att_text(nc_id_qb,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /* UVPT */
  sprintf(Message,"CF-1.4");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
  sprintf(Message,"EPIC Model re-gridding of U, V, P, and T on Theta Surfaces; input data are Weizmann Institute U and T on Pressure Levels");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"title",strlen(Message)+1,Message);
  sprintf(Message,"%s; \n"
                  "%s: U, V, P and T on isentropic surfaces uses the EPIC GCM v%4.2f",
                  weizmann_grid->infile_history,date,grid.epic_version);
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"history",strlen(Message)+1,Message);
  sprintf(Message,"Weizmann Institute:  Eli Galanti (eli.galanti@weizmann@ac.il), Yohai Kaspi (yohai.kaspi@weizmann.ac.il), "
                  "EPIC: Timothy E. Dowling (dowling@louiville.edu)");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"contact",strlen(Message)+1,Message);
  sprintf(Message,"This file contains zonal and meridional winds, u and v, pressure, p, "
                  "temperature, T, and geopotential, phi, on potential-temperature surfaces, "
                  "theta, for the atmosphere of Jupiter.");
  nc_put_att_text(nc_id_uvpt,NC_GLOBAL,"comment",strlen(Message)+1,Message);

  /*
   * Define coordinates and variables.
   */

  /*
   * lon (I direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_qb,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lon",grid.ni,&dimid[NETCDF_I_INDEX]);
  nc_def_var(nc_id_uvpt,"lon",NC_DOUBLE,1,&dimid[NETCDF_I_INDEX],&coorid[NETCDF_I_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"standard_name",strlen("longitude")+1,"longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"long_name",strlen("Longitude")+1,"Longitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"units",strlen("degrees_east")+1,"degrees_east");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_I_INDEX],"axis",strlen("X")+1,"X");

  /*
   * lat (J direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_qb,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"lat",grid.nj+1,&dimid[NETCDF_J_INDEX]);
  nc_def_var(nc_id_uvpt,"lat",NC_DOUBLE,1,&dimid[NETCDF_J_INDEX],&coorid[NETCDF_J_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"standard_name",strlen("latitude")+1,"latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"long_name",strlen("Latitude")+1,"Latitude");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"units",strlen("degrees_north")+1,"degrees_north");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_J_INDEX],"axis",strlen("Y")+1,"Y");

  /*
   * theta (K direction)
   */

  /* QB */
  nc_def_dim(nc_id_qb,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_qb,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"theta",grid.nk,&dimid[NETCDF_K_INDEX]);
  nc_def_var(nc_id_uvpt,"theta",NC_DOUBLE,1,&dimid[NETCDF_K_INDEX],&coorid[NETCDF_K_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"standard_name",strlen("air_potential_temperature")+1,"air_potential_temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"long_name",strlen("Potential temperature")+1,"Potential temperature");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"units",strlen("K")+1,"K");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"axis",strlen("Z")+1,"Z");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_K_INDEX],"positive",strlen("up")+1,"up");

  /*
   * time
   */

  /* QB */
  nc_def_dim(nc_id_qb,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_qb,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"long_name",strlen("elapsed time")+1,"elapsed time");
  message_len = strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S 0",gmtime(&var.start_time));
  nc_err = nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"units",message_len+1,Message);
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");
  nc_put_att_text(nc_id_qb,coorid[NETCDF_T_INDEX],"calendar",strlen("julian")+1,"julian");

  /* UVPT */
  nc_def_dim(nc_id_uvpt,"time",NC_UNLIMITED,&dimid[NETCDF_T_INDEX]);
  nc_def_var(nc_id_uvpt,"time",NC_DOUBLE,1,&dimid[NETCDF_T_INDEX],&coorid[NETCDF_T_INDEX]);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"standard_name",strlen("time")+1,"time");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"long_name",strlen("Time")+1,"Time");
  message_len = strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S 0",gmtime(&var.start_time));
  nc_err = nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"units",message_len+1,Message);
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"axis",strlen("T")+1,"T");
  nc_put_att_text(nc_id_uvpt,coorid[NETCDF_T_INDEX],"calendar",strlen("julian")+1,"julian");

  /*
   * QB specific fields
   */

  /*
   * Potential vorticity (PV)
   *
   * NOTE: For isentropic coordinates, and for the isentropic region of hybrid coordinates,
   *       this PV is Ertel's isentropic potential vorticity. For isobaric coordinates, it is
   *       proportional to absolute vorticity on pressure levels, because h is constant.
   */
  pv_info.name = (char *)malloc(strlen("pv")+1);
  sprintf(pv_info.name,"%s","pv");
  pv_info.standard_name = (char *)malloc(strlen("potential_vorticity")+1);
  sprintf(pv_info.standard_name,"%s","potential_vorticity");
  pv_info.long_name = (char *)malloc(strlen("potential vorticity")+1);
  sprintf(pv_info.long_name,"%s","potential vorticity");
  pv_info.units = (char *)malloc(strlen("K m2 kg-1 s-1")+1);
  sprintf(pv_info.units,"%s","K m2 kg-1 s-1");
  for (i = 0; i < 4; i++) {
    pv_info.dimid[ i] = dimid[ i];
    pv_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,pv_info.name,NC_DOUBLE,4,pv_info.dimid,&pv_info.id);
  nc_put_att_text(nc_id_qb,pv_info.id,"standard_name",strlen(pv_info.standard_name)+1,pv_info.standard_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"long_name",strlen(pv_info.long_name)+1,pv_info.long_name);
  nc_put_att_text(nc_id_qb,pv_info.id,"units",strlen(pv_info.units)+1,pv_info.units);

  /*
   * Bernoulli streamfunction
   */
  bsf_info.name = (char *)malloc(strlen("bsf")+1);
  sprintf(bsf_info.name,"%s","bsf");
  bsf_info.standard_name = (char *)malloc(strlen("specific_dry_energy_of_air")+1);
  sprintf(bsf_info.standard_name,"%s","specific_dry_energy_of_air");
  bsf_info.long_name = (char *)malloc(strlen("Bernoulli streamfunction")+1);
  sprintf(bsf_info.long_name,"%s","Bernoulli streamfunction");
  bsf_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(bsf_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    bsf_info.dimid[ i] = dimid[ i];
    bsf_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,bsf_info.name,NC_DOUBLE,4,bsf_info.dimid,&bsf_info.id);
  nc_put_att_text(nc_id_qb,bsf_info.id,"standard_name",strlen(bsf_info.standard_name)+1,bsf_info.standard_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"long_name",strlen(bsf_info.long_name)+1,bsf_info.long_name);
  nc_put_att_text(nc_id_qb,bsf_info.id,"units",strlen(bsf_info.units)+1,bsf_info.units);

  /*
   * NH2
   */
  NH2_info.name = (char *)malloc(strlen("NHsquared")+1);
  sprintf(NH2_info.name,"%s","NHsquared");
  NH2_info.long_name = (char *)malloc(strlen("(NH)^2")+1);
  sprintf(NH2_info.long_name,"%s","(NH)^2");
  NH2_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(NH2_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    NH2_info.dimid[ i] = dimid[ i];
    NH2_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,NH2_info.name,NC_DOUBLE,4,NH2_info.dimid,&NH2_info.id);
  nc_put_att_text(nc_id_qb,NH2_info.id,"long_name",strlen(NH2_info.long_name)+1,NH2_info.long_name);
  nc_put_att_text(nc_id_qb,NH2_info.id,"units",strlen(NH2_info.units)+1,NH2_info.units);
  nc_put_att_text(nc_id_qb,NH2_info.id,"comment",strlen(
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H")+1,
    "square of product of Brunt-Vailsala frequency, N, and pressure scale height, H");

  /*
   * Kinetic energy
   */
  kin_info.name = (char *)malloc(strlen("kin")+1);
  sprintf(kin_info.name,"%s","kin");
  kin_info.standard_name = (char *)malloc(strlen("specific_kinetic_energy_of_air")+1);
  sprintf(kin_info.standard_name,"%s","specific_kinetic_energy_of_air");
  kin_info.long_name = (char *)malloc(strlen("Kinetic energy")+1);
  sprintf(kin_info.long_name,"%s","Kinetic energy");
  kin_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(kin_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    kin_info.dimid[ i] = dimid[ i];
    kin_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_qb,kin_info.name,NC_DOUBLE,4,kin_info.dimid,&kin_info.id);
  nc_put_att_text(nc_id_qb,kin_info.id,"long_name",strlen(kin_info.long_name)+1,kin_info.long_name);
  nc_put_att_text(nc_id_qb,kin_info.id,"units",strlen(kin_info.units)+1,kin_info.units);
  nc_put_att_text(nc_id_qb,kin_info.id,"comment",strlen(
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)")+1,
    "The horizontal kinetic energy per mass, .5*(u*u+v*v)");

  /*
   * UVPT specific fields
   */

  /*
   * Zonal wind
   */
  u_info.name = (char *)malloc(strlen("u")+1);
  sprintf(u_info.name,"%s","u");
  u_info.standard_name = (char *)malloc(strlen("eastward_wind")+1);
  sprintf(u_info.standard_name,"%s","eastward_wind");
  u_info.long_name = (char *)malloc(strlen("Zonal wind")+1);
  sprintf(u_info.long_name,"%s","Zonal wind");
  u_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(u_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    u_info.dimid[ i] = dimid[ i];
    u_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,u_info.name,NC_DOUBLE,4,u_info.dimid,&u_info.id);
  nc_put_att_text(nc_id_uvpt,u_info.id,"standard_name",strlen(u_info.standard_name)+1,u_info.standard_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"long_name",strlen(u_info.long_name)+1,u_info.long_name);
  nc_put_att_text(nc_id_uvpt,u_info.id,"units",strlen(u_info.units)+1,u_info.units);

//   /*
//    * Geostrophic zonal wind (assuming variable-f geostrophy)
//    */
//   ug_info.name = (char *)malloc(strlen("ug")+1);
//   sprintf(ug_info.name,"%s","ug");
//   ug_info.standard_name = (char *)malloc(strlen("geostrophic_eastward_wind")+1);
//   sprintf(ug_info.standard_name,"%s","geostrophic_eastward_wind");
//   ug_info.long_name = (char *)malloc(strlen("Geostrophic zonal wind")+1);
//   sprintf(ug_info.long_name,"%s","Geostrophic zonal wind");
//   ug_info.units = (char *)malloc(strlen("m s-1")+1);
//   sprintf(ug_info.units,"%s","m s-1");
//   for (i = 0; i < 4; i++) {
//     ug_info.dimid[ i] = dimid[ i];
//     ug_info.coorid[i] = coorid[i];
//   }
// 
//   nc_def_var(nc_id_uvpt,ug_info.name,NC_DOUBLE,4,ug_info.dimid,&ug_info.id);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"standard_name",strlen(ug_info.standard_name)+1,ug_info.standard_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"long_name",strlen(ug_info.long_name)+1,ug_info.long_name);
//   nc_put_att_text(nc_id_uvpt,ug_info.id,"units",strlen(ug_info.units)+1,ug_info.units);

  /*
   * Meridional wind
   */
  v_info.name = (char *)malloc(strlen("v")+1);
  sprintf(v_info.name,"%s","v");
  v_info.standard_name = (char *)malloc(strlen("northward_wind")+1);
  sprintf(v_info.standard_name,"%s","northward_wind");
  v_info.long_name = (char *)malloc(strlen("Meridional wind")+1);
  sprintf(v_info.long_name,"%s","Meridional wind");
  v_info.units = (char *)malloc(strlen("m s-1")+1);
  sprintf(v_info.units,"%s","m s-1");
  for (i = 0; i < 4; i++) {
    v_info.dimid[ i] = dimid[ i];
    v_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,v_info.name,NC_DOUBLE,4,v_info.dimid,&v_info.id);
  nc_put_att_text(nc_id_uvpt,v_info.id,"standard_name",strlen(v_info.standard_name)+1,v_info.standard_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"long_name",strlen(v_info.long_name)+1,v_info.long_name);
  nc_put_att_text(nc_id_uvpt,v_info.id,"units",strlen(v_info.units)+1,v_info.units);

  /*
   * Pressure
   */
  p_info.name = (char *)malloc(strlen("p")+1);
  sprintf(p_info.name,"%s","p");
  p_info.standard_name = (char *)malloc(strlen("air_pressure")+1);
  sprintf(p_info.standard_name,"%s","air_pressure");
  p_info.long_name = (char *)malloc(strlen("Pressure")+1);
  sprintf(p_info.long_name,"%s","Pressure");
  p_info.units = (char *)malloc(strlen("Pa")+1);
  sprintf(p_info.units,"%s","Pa");
  for (i = 0; i < 4; i++) {
    p_info.dimid[ i] = dimid[ i];
    p_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,p_info.name,NC_DOUBLE,4,p_info.dimid,&p_info.id);
  nc_put_att_text(nc_id_uvpt,p_info.id,"standard_name",strlen(p_info.standard_name)+1,p_info.standard_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"long_name",strlen(p_info.long_name)+1,p_info.long_name);
  nc_put_att_text(nc_id_uvpt,p_info.id,"units",strlen(p_info.units)+1,p_info.units);

  /*
   * Temperature
   */
  t_info.name = (char *)malloc(strlen("T")+1);
  sprintf(t_info.name,"%s","T");
  t_info.standard_name = (char *)malloc(strlen("air_temperature")+1);
  sprintf(t_info.standard_name,"%s","air_temperature");
  t_info.long_name = (char *)malloc(strlen("Temperature")+1);
  sprintf(t_info.long_name,"%s","Temperature");
  t_info.units = (char *)malloc(strlen("K")+1);
  sprintf(t_info.units,"%s","K");
  for (i = 0; i < 4; i++) {
    t_info.dimid[ i] = dimid[ i];
    t_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,t_info.name,NC_DOUBLE,4,t_info.dimid,&t_info.id);
  nc_put_att_text(nc_id_uvpt,t_info.id,"standard_name",strlen(t_info.standard_name)+1,t_info.standard_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"long_name",strlen(t_info.long_name)+1,t_info.long_name);
  nc_put_att_text(nc_id_uvpt,t_info.id,"units",strlen(t_info.units)+1,t_info.units);

  /*
   * Geopotential
   */
  phi_info.name = (char *)malloc(strlen("phi")+1);
  sprintf(phi_info.name,"%s","phi");
  phi_info.standard_name = (char *)malloc(strlen("geopotential")+1);
  sprintf(phi_info.standard_name,"%s","geopotential");
  phi_info.long_name = (char *)malloc(strlen("Geopotential")+1);
  sprintf(phi_info.long_name,"%s","Geopotential");
  phi_info.units = (char *)malloc(strlen("m2 s-2")+1);
  sprintf(phi_info.units,"%s","m2 s-2");
  for (i = 0; i < 4; i++) {
    phi_info.dimid[ i] = dimid[ i];
    phi_info.coorid[i] = coorid[i];
  }

  nc_def_var(nc_id_uvpt,phi_info.name,NC_DOUBLE,4,phi_info.dimid,&phi_info.id);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"standard_name",strlen(phi_info.standard_name)+1,phi_info.standard_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"long_name",strlen(phi_info.long_name)+1,phi_info.long_name);
  nc_put_att_text(nc_id_uvpt,phi_info.id,"units",strlen(phi_info.units)+1,phi_info.units);

  /*---------------------------*
   * Leave netcdf define mode: *
   *---------------------------*/
  nc_enddef(nc_id_qb);
  nc_enddef(nc_id_uvpt);

  /*
   * Assign values to X, Y and Z coordinates.
   */
  /* lon */
  for (I = ILO; I <= IHI; I++) {
    nc_index[0] = I-ILO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_I_INDEX],nc_index,&(grid.lon[2*I+1]));
  }

  /* lat */
  for (J = JLO; J <= JHI; J++) {
    nc_index[0] = J-JLO;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_J_INDEX],nc_index,&(grid.lat[2*J+1]));
  }

  /* theta */
  for (K = KHI; K >= KLO; K--) {
    nc_index[0] = KHI-K;
    nc_put_var1_double(nc_id_qb,  coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_K_INDEX],nc_index,&(grid.sigmatheta[2*K]));
  }

  fprintf(stdout,"Computing Q, B and (NH)^2 on theta surfaces...  0%%");

  /* Allocate memory */
  p      = fvector(0,2*grid.nk+1,dbmsname);
  h      = fvector(0,2*grid.nk+1,dbmsname);
  buffer = (double *)calloc(grid.ni,sizeof(double));

  for (weizmann_itime = 0; weizmann_itime < weizmann_grid->ntime; weizmann_itime++) {
    /* Show progress */
    fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(double)(weizmann_itime+1)/weizmann_grid->ntime));
    fflush(stdout);

    /*
     * Input Weizmann data for timeframe weizmann_itime, and interpolate into EPIC variables.
     */
    weizmann_var_read(planet,weizmann_grid,weizmann_infile,VAR_DATA,weizmann_itime);

    /*
     * Write time.
     */

    /* QB */
    nc_index[0] = weizmann_itime;
    nc_put_var1_double(nc_id_qb,coorid[NETCDF_T_INDEX],nc_index,&weizmann_grid->time[weizmann_itime]);

    /* UVPT */
    nc_index[0] = weizmann_itime;
    nc_put_var1_double(nc_id_uvpt,coorid[NETCDF_T_INDEX],nc_index,&weizmann_grid->time[weizmann_itime]);

    /*
     * Compute Q and B on theta surfaces by running the normal EPIC model
     * diagnostic-variable calculations.
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      for (I = ILO; I <= IHI; I++) {
        for (kk = 1; kk <= 2*KHI+1; kk++) {
          p[kk] = get_p(planet,P2_INDEX,kk,J,I);
        }
        calc_h(jj,p,h);

        for (K = KLO; K <= KHI; K++) {
          H(K,J,I) = h[2*K];
        }
        K = 0;
        H(K,J,I) = SQR(H(K+1,J,I))/H(K+2,J,I);
        K = KHI+1;
        H(K,J,I) = SQR(H(K-1,J,I))/H(K-2,J,I);
      }
    }
    bc_lateral(var.h.value,THREEDIM);

    /*
     * PHI3(KHI,J,I) is set above by weizmann_var_read(VAR_DATA).
     */
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,PASSING_PHI3NK);
    store_diag(planet);

    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;
    nc_count[NETCDF_I_INDEX] = grid.ni;

    nc_start[NETCDF_T_INDEX] = weizmann_itime;
    nc_start[NETCDF_I_INDEX] = 0;

    for (K = KHI; K >= KLO; K--) {
      kk                       = 2*K;
      nc_start[NETCDF_K_INDEX] = KHI-K;
      for (J = JLO; J <= JHI; J++) {
        jj = 2*J;

        g  = GRAVITY2(K,J);

        nc_start[NETCDF_J_INDEX] = J-JLO;

        /*--------------------*
         * QB specific fields *
         *--------------------*/

        /*
         * Write potential vorticity (pv, aka Q).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .25*(PV2(K,J,I)+PV2(K,J,I+1)+PV2(K,J+1,I)+PV2(K,J+1,I+1));
        }
        nc_put_vara_double(nc_id_qb,pv_info.id,nc_start,nc_count,buffer);

        /*
         * Write Bernoulli streamfunction (bsf).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = MONT2(K,J,I)
                         +get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
        }
        nc_put_vara_double(nc_id_qb,bsf_info.id,nc_start,nc_count,buffer);

        /*
         * Write (NH)^2, the square of the product of the buoyancy (Brunt-Vaisala) frequency, N, 
         * and the pressure scale height, H = RT/g.
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_brunt2(planet,kk,J,I)
                         *SQR(planet->rgas*T2(K,J,I)/g);
        }
        nc_put_vara_double(nc_id_qb,NH2_info.id,nc_start,nc_count,buffer);

        /*
         * Write kinetic energy per mass (kin), .5*(u*u+v*v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                         var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
        }
        nc_put_vara_double(nc_id_qb,kin_info.id,nc_start,nc_count,buffer);

        /*----------------------*
         * UVPT specific fields *
         *----------------------*/

        /*
         * Write zonal wind (u).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(U(grid.it_uv,K,J,I)+U(grid.it_uv,K,J,I+1));
        }
        nc_put_vara_double(nc_id_uvpt,u_info.id,nc_start,nc_count,buffer);

//         /*
//          * Calculate and write geostrophic zonal wind (ug).
//          */
//         if (J == JLO) {
//           /*
//            * Southern edge: take forward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J,I));
//           }
//         }
//         else if (J == JHI) {
//           /*
//            * Northern edge: take backward step.
//            */
//           tmp = -grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J,I)-MONT2(K,J-1,I));
//           }
//         }
//         else {
//           /*
//            * Interior point: take long step to center back onto h-grid.
//            */
//           tmp = -.5*grid.n[kk][jj+1]/grid.f[jj+1];
//           for (I = ILO; I <= IHI; I++) {
//             buffer[I-ILO] = tmp*(MONT2(K,J+1,I)-MONT2(K,J-1,I));
//           }
//         }
//         nc_put_vara_double(nc_id_uvpt,ug_info.id,nc_start,nc_count,buffer);

        /*
         * Write meridional wind (v).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = .5*(V(grid.it_uv,K,J,I)+V(grid.it_uv,K,J+1,I));
        }
        nc_put_vara_double(nc_id_uvpt,v_info.id,nc_start,nc_count,buffer);

        /*
         * Write pressure (p).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = P2(K,J,I);
        }
        nc_put_vara_double(nc_id_uvpt,p_info.id,nc_start,nc_count,buffer);

        /*
         * Write temperature (T).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = T2(K,J,I);
        }
        nc_put_vara_double(nc_id_uvpt,t_info.id,nc_start,nc_count,buffer);

        /*
         * Write geopotential (phi).
         */
        for (I = ILO; I <= IHI; I++) {
          buffer[I-ILO] = PHI2(K,J,I);
        }
        nc_put_vara_double(nc_id_uvpt,phi_info.id,nc_start,nc_count,buffer);

      } /* J loop */
    } /* K loop */
  } /* weizmann_itime loop */

  /* Free allocated memory */
  free_fvector(p,0,2*grid.nk+1,dbmsname);
  free_fvector(h,0,2*grid.nk+1,dbmsname);
  free(buffer);

  nc_close(nc_id_qb);
  nc_close(nc_id_uvpt);

  fprintf(stdout,"\b\b\b\b100%%\n");
  fflush(stdout);
  fprintf(stdout,"Output written to %s\n",weizmann_outfile_qb);

}

/*======================= end of weizmann_conversion() =======================*/

/* * * * * * * * * * * * end of epic_funcs_init.c * * * * * * * * * * * * * * */
