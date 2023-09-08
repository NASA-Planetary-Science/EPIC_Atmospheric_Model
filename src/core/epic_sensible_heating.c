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

/* * * * * * * * epic_sensible_heating.c * * * * * * * * * * * * * * 
 *                                                                 *
 *  This file includes functions used to calculate contributions   *
 *  to the heating rate, HEAT3(K,J,I), from processes involving    *
 *  temperature differences (sensible heating).                    *
 *                                                                 *
 *  The heating-rate units are [W/kg].                             *
 *                                                                 *
 *       radiative_heating()                                       *
 *       heating_from_file()                                       *
 *       newtonian_cooling()                                       *
 *         t_rad()                                                 *
 *         temp_eq()                                               *
 *       rt_heating()                                              *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>
/*======================= radiative_heating() ===============================*/

/*
 * Wrapper function for radiative-heating function specified 
 * by grid.radiation_scheme.
 *
 * Values of input parameter action are:
 *   EPIC_ALLOC     (allocate memory before applying radiative heating)
 *   EPIC_APPLY     (apply radiative heating)
 *   EPIC_FREE      (free allocated memory when done)
 */
void radiative_heating(planetspec *planet,
                       int         action)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="radiative_heating";

  if (strcmp(grid.radiation_scheme,"off") == 0) {
    return;
  }
  else if (strcmp(grid.radiation_scheme,"Correlated k") == 0) {
    rt_heating(planet,action);
    return;
  }
  else if (strcmp(grid.radiation_scheme,"Newtonian") == 0) {
    newtonian_cooling(planet,action);
    return;
  }
  else if (strcmp(grid.radiation_scheme,"Heating from file") == 0) {
    heating_from_file(planet,action);
    return;
  }
  else {
    sprintf(Message,"grid.radiation_scheme=%s not recognized",grid.radiation_scheme);
    epic_error(dbmsname,Message);
  }
}

/*====================== end radiative_heating() ============================*/

/*====================== heating_from_file() ================================*/

/*
 * Apply heating/cooling based on input from a file.  This heating is 
 * one way---there is no feedback from the model's dynamics and chemistry.
 *
 * NOTE: For seasonal variation, we could use a truncated Fourier series with
 *       complex coefficients. For an example, see EPIC Version 3.85, the
 *       funtion temp_eq().
 */

#define FILE_HEATING(k,j) file_heating[j+(k)*nlat]

void heating_from_file(planetspec *planet,
                       int         action)
{
  int
    K,J,I,
    jj,j,k;
  EPIC_FLOAT
    frac_lat,frac_logp,
    logp3;
  static int
    np,nlat;
  static EPIC_FLOAT
    *logp,
    *lat,
    *file_heating;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="heating_from_file";

  switch(action) {
    case EPIC_ALLOC:
      /*
       * Determine size, allocate memory, and read in data.
       */
      if (strcmp(planet->name,"Jupiter") == 0) {
        /* Determine size of data file. */
        read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/stratospheric_netheating.Jupiter.Allen_etal",SIZE_DATA,
                              &np,&nlat,NULL,NULL,NULL);

        /* Screen for bad input values. */
        if (np < 1) {
          sprintf(Message,"np=%d < 1",np);
          epic_error(dbmsname,Message);
        }

        /* Allocate memory */
        logp         = fvector(0,     np-1,dbmsname);
        lat          = fvector(0,   nlat-1,dbmsname);
        file_heating = fvector(0,np*nlat-1,dbmsname);

        /* Read in data */
        read_meridional_plane(planet,EPIC_PATH"/data/Jupiter/stratospheric_netheating.Jupiter.Allen_etal",POST_SIZE_DATA,
                              &np,&nlat,logp,lat,file_heating);
      }
      else {
        sprintf(Message,"action = %d not yet implemented for %s",action,planet->name);
        epic_error(dbmsname,Message);
      }
    break;

    case EPIC_FREE:
      /*
       * Free allocated memory.
       */
      if (strcmp(planet->name,"Jupiter") == 0) {
        free_fvector(logp,        0,np-1,     dbmsname);
        free_fvector(lat,         0,nlat-1,   dbmsname);
        free_fvector(file_heating,0,np*nlat-1,dbmsname);
      }
      else {
        sprintf(Message,"action = %d not yet implemented for %s",action,planet->name);
        epic_error(dbmsname,Message);
      }
    break;

    case EPIC_APPLY:
      /*
       * Apply heating
       */
      if (strcmp(planet->name,"Jupiter") == 0) {
        /*
         * Use bilinear interpolation.
         */
        for (J = JLO; J <= JHI; J++) {
          jj = 2*J+1;

          if (grid.lat[jj] > lat[nlat-1]) {
            j        = nlat-2;
            frac_lat = (grid.lat[jj]-lat[j])/(lat[j+1]-lat[j]);
          }
          else {
            for (j = 1; j < nlat; j++) {
              if (grid.lat[jj] <= lat[j]) {
                j--;
                frac_lat = (grid.lat[jj]-lat[j])/(lat[j+1]-lat[j]);
                break;
              }
            }  
          }
          for (I = ILO; I <= IHI; I++) {
            for (K = KLO; K <= KHI; K++) {
              /*
               * Add no contribution to HEAT if outside pressure range of heating data.
               */
              logp3 = log(P3(K,J,I));
              if (logp3 < logp[0]) {
                continue;
              }
              else if (logp3 > logp[np-1]) {
                continue;
              }
              else {
                for (k = 1; k < np; k++) {
                  if (logp3 <= logp[k]) {
                    k--;
                    frac_logp = (logp3-logp[k])/(logp[k+1]-logp[k]);
                    break;
                  }
                }
              }  
              HEAT3(K,J,I) += FILE_HEATING(k,  j  )*(1.-frac_logp)*(1.-frac_lat)
                             +FILE_HEATING(k+1,j  )*(   frac_logp)*(1.-frac_lat)
                             +FILE_HEATING(k,  j+1)*(1.-frac_logp)*(   frac_lat)
                             +FILE_HEATING(k+1,j+1)*(   frac_logp)*(   frac_lat);
            }
          }
        }
        /* HEAT does not need bc_lateral() here. */
      }
      else {
        sprintf(Message,"need to implement for %s",planet->name);
        epic_error(dbmsname,Message);
      }
    break;

    default:
      sprintf(Message,"action = %d not yet defined",action);
      epic_error(dbmsname,Message);
    break;
  }

  return;
}

/*====================== end of heating_from_file() =========================*/

/*====================== newtonian_cooling() ================================*/

/*
 * Force temperature profile to the given radiative-equilibrium profile
 * at an appropriate rate.  Keep the layer-average heating zero if requested.
 *
 *
 * NOTE: Lateral boundary conditions, bc_lateral(), are not applied to HEAT array here.
 */

void newtonian_cooling(planetspec *planet,
                       int         action)
{
  register int
    K,J,I,
    ki;
  register EPIC_FLOAT
    lat,
    pressure,
    temperature,
    fpara,
    cp_over_time,
    da,heat_avg;
  static EPIC_FLOAT
   *buffji;
  EPIC_FLOAT
    t_eq,
    area,heat_area,
    tmp;
  char
    header[N_STR],
    infile[FILE_STR];
  FILE
    *t_cool_vs_p; 

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="newtonian_cooling";

  switch(action) {
    case EPIC_ALLOC:
      /* 
       * Allocate memory.
       */
      buffji = fvector(0,Nelem2d-1,dbmsname);

      /*
       * Set up t_cool_table.
       * Input t_cool_vs_p, looking first in the local directory.
       *
       * grid.newt_cool_adjust == TRUE: sets layer avg. of newt. cooling to 0
       */

      /* Look in local directory first. */
      sprintf(infile,"./t_cool_vs_p.%s",planet->name);
      t_cool_vs_p = fopen(infile,"r");
      if (!t_cool_vs_p) {
        sprintf(infile,EPIC_PATH"/data/%s/t_cool_vs_p.%s",planet->name,planet->name);
        t_cool_vs_p = fopen(infile,"r");
      }
      if (t_cool_vs_p) {
        FILE
          *relax_times;

        for (ki = 0; ki < 6; ki++) {
          fgets(header,N_STR,t_cool_vs_p);
        }
        /* input number of data points */
        fscanf(t_cool_vs_p, "%d", &var.n_t_cool); 

        /* 
         * Allocate memory for var.t_cool_table.
         */
        var.t_cool_table = ftriplet(0,var.n_t_cool-1,dbmsname);

        relax_times = fopen("./relax_times.dat","w");
        fprintf(relax_times," %s\n",planet->name);
        fprintf(relax_times,"   p[hPa]        trad[s] \n");

        /* stored in order of decreasing p (increasing -log p) */
        for (ki = var.n_t_cool-1; ki >= 0;  ki--) {  

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(t_cool_vs_p,"%lf %*f %lf",&var.t_cool_table[ki].x,&var.t_cool_table[ki].y);
#else
          fscanf(t_cool_vs_p,"%f %*f %f",&var.t_cool_table[ki].x,&var.t_cool_table[ki].y);
#endif

          /* convert from hPa to Pa */
          var.t_cool_table[ki].x *= 100.;

          /* spline on -log p  */
          var.t_cool_table[ki].x = -log(var.t_cool_table[ki].x);  

          /* Output table of relaxation times vs pressure: */
          fprintf(relax_times," %e  %e \n",
                  exp(-var.t_cool_table[ki].x)/100.,var.t_cool_table[ki].y);
        }
        fclose(t_cool_vs_p);
        fclose(relax_times);
      }
      else {
        /*
         * Unable to read t_cool_vs_p data file.
         */
        var.n_t_cool = 0;
        if (strcmp(planet->name,"Held_Suarez") == 0) {
          /*
           * The Held-Suarez reference case does not use var.t_cool_table.
           */
          ;
        }
        else {
          sprintf(Message,"cannot find %s",infile);
          epic_error(dbmsname,Message);
        }
      }
    break;

    case EPIC_FREE:
      /*
       * Free allocated memory.
       */
      free_fvector(buffji,0,Nelem2d-1,dbmsname);
      free_ftriplet(var.t_cool_table,0,var.n_t_cool-1,dbmsname);
    break;

    case EPIC_APPLY:
      for (K = KLO; K <= KHI; K++) {
        for (J = JLO; J <= JHI; J++) {
          lat = grid.lat[2*J+1];
          for (I = ILO; I <= IHI; I++) {
            pressure     = P3(K,J,I);
            temperature  = T3(K,J,I);
            if (var.fpara.on) {
              fpara = FPARA(K,J,I);
            }
            else {
              fpara = return_fpe(temperature);
            }
            cp_over_time = return_cp(planet,fpara,pressure,temperature)/t_rad(planet,K,J,I);
            /*
             * Determine T_eq:
             */
            switch(planet->index) {
              case HELD_SUAREZ_INDEX:
                t_eq = temp_eq(planet,lat,pressure,0.);
              break;
              case VENUS_LLR05_INDEX:
                t_eq = temp_eq(planet,lat,pressure,0.);
              break;
              case TITAN_INDEX:
                t_eq = temp_eq(planet,lat,pressure,0.);
              break;
              default:
                /*
                 * Use t_vs_p data for radiative equilibrium profile.
                 */
                get_sounding(planet,pressure,"temperature",&t_eq);
              break;
            }
            BUFFJI(J,I) = -cp_over_time*(temperature-t_eq);
          }
        }
        /* No need to apply bc_lateral() here. */

        /*
         * Make layer average in sponge zero, and do this for the rest
         * of the model if requested.
         */
        if (K <= grid.k_sponge || grid.newt_cool_adjust) {
          /*
           * Ensure that layer average is zero.
           */
          area      = 0.;
          heat_area = 0.;
          for (J = JLO; J <= JHI; J++) {
            da = 1./grid.mn[2*K+1][2*J+1];
            for (I = ILO; I <= IHI; I++) {
              area      += da;
              heat_area += BUFFJI(J,I)*da;
            }
          }

#if defined(EPIC_MPI)
          tmp = area;
          MPI_Allreduce(&tmp,&area,1,float_type,MPI_SUM,para.comm);
          tmp = heat_area;
          MPI_Allreduce(&tmp,&heat_area,1,float_type,MPI_SUM,para.comm);
#endif

          heat_avg = heat_area/area;

          for (J = JLO; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              BUFFJI(J,I) -= heat_avg;
            }
          }
        }

        /*
         * Apply result to HEAT3 array.
         */
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            HEAT3(K,J,I) += BUFFJI(J,I); 
          }
        }
      }
    break;

    default:
      sprintf(Message,"action = %d not yet defined",action);
      epic_error(dbmsname,Message);
    break;
  }

  return;
}

/*====================== end of newtonian_cooling() =========================*/

/*======================= t_rad() ===========================================*/

/*
 * Return radiative cooling time [s], given position K,J,I.
 * The value is for the bottom interface of layer K.
 */

EPIC_FLOAT t_rad(planetspec *planet,
                 int         K,
                 int         J,
                 int         I)
{
  register int
    ki;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    ka,ks;
  EPIC_FLOAT
    pressure,
    t_cool,
    p_t_cool_d,   
    neglogp;

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="t_rad";

  /*
   * Check that Newtonian cooling is turned on.
   */
  if (strcmp(grid.radiation_scheme,"Newtonian") != 0) {
    sprintf(Message,"grid.radiation_scheme=%s",grid.radiation_scheme);
    epic_error(dbmsname,Message);
  }

  if (!initialized) {
    if (strcmp(planet->name,"Held_Suarez") == 0) {
      /*
       * Set parameters ka, ks:
       */
      ka = 1./(40.*60.*60.*24.);
      ks = 1./( 4.*60.*60.*24.);
    }
    else {
      spline_pchip(var.n_t_cool,var.t_cool_table);
    }
    initialized = TRUE;
  }
  /* End of initialization. */

  /*
   * Use interface pressure.
   */
  pressure = P3(K,J,I);

  if (strcmp(planet->name,"Held_Suarez") == 0) {
    EPIC_FLOAT
      kt,amp,
      pbot,cos_lat;

    kt   = ka;
    pbot = P3(grid.nk,J,I);
    amp  = (pressure/pbot-0.7)/(1.-0.7);
    if (amp > 0.) {
      cos_lat = cos(grid.lat[2*J+1]*DEG);
      kt     += (ks-ka)*pow(cos_lat,4.)*amp;
    }
    t_cool = 1./kt;
  }
  else {
    /*
     *  Interpolate to get Newtonian cooling time:
     */
    neglogp = -log(pressure);
    if (neglogp >= var.t_cool_table[var.n_t_cool-1].x) {
      /* past range of t_cool data */
      t_cool = var.t_cool_table[var.n_t_cool-1].y;
    }
    else if (neglogp <= var.t_cool_table[0].x) {
      /* past range of t_cool data */
      t_cool = var.t_cool_table[0].y;
    }
    else {
      ki     = find_place_in_table(var.n_t_cool,var.t_cool_table,neglogp,&p_t_cool_d);
      t_cool = splint_pchip(neglogp,var.t_cool_table+ki,p_t_cool_d);
    }
  }

  return t_cool;
}

/*======================= end of t_rad() ====================================*/

/*======================= temp_eq() =========================================*/

/*
 * Inputs: pressure [Pa], latitude [deg], time [sec].
 *
 * time = 0.0 corresponds to northern spring equinox.
 */
EPIC_FLOAT temp_eq(planetspec *planet,
                   EPIC_FLOAT  latitude,
                   EPIC_FLOAT  pressure,
                   EPIC_FLOAT  time)
{
  int
    ki;
  static int
    initialized=FALSE;
  register EPIC_FLOAT
    t_tp,dt_tp;
  EPIC_FLOAT
    ans,
    neglogp,p_d,
    sin_y_j,sin_y;
  static EPIC_FLOAT
    dt_y,dth_z,p0;
  static float_triplet
    *t_table,
    *dt_table;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="temp_eq";

  /*
   * Check that Newtonian cooling is turned on.
   */
  if (strcmp(grid.radiation_scheme,"Newtonian") != 0) {
    sprintf(Message,"grid.radiation_scheme=%s",grid.radiation_scheme);
    epic_error(dbmsname,Message);
  }

  if (!initialized) {
    if (strcmp(planet->name,"Held_Suarez") == 0) {
      /*
       * Set parameters for the held_suarez test case.
       */
      dt_y  = 60.;
      dth_z = 10.;
      p0    = 1000.*100.;   /* Held-Suarez case */
    }
    else {
      if (var.ntp <= 0) {
        sprintf(Message,"var.ntp=%d",var.ntp);
        epic_error(dbmsname,Message);
      }

      /* Allocate memory */
      t_table  = ftriplet(0,var.ntp-1,dbmsname);
      dt_table = ftriplet(0,var.ntp-1,dbmsname);

      /* Assign table values */
      for (ki = 0; ki < var.ntp; ki++) { 
        t_table[ ki].x = dt_table[ ki].x = -log(var.pdat[ki]);
        t_table[ ki].y = var.tdat[ ki];
        dt_table[ki].y = var.dtdat[ki];
      }
      spline_pchip(var.ntp, t_table);
      spline_pchip(var.ntp,dt_table); 
    }

    initialized = TRUE;
  }
  /* end of initialization */

  if (strcmp(planet->name,"Held_Suarez") == 0) {
    EPIC_FLOAT
      cos2_lat,
      sin2_lat,
      p_p0;

    cos2_lat  = cos(latitude*DEG);
    cos2_lat *= cos2_lat;
    sin2_lat  = 1.-cos2_lat;
    p_p0      = pressure/p0;
    ans       = (315.-dt_y*sin2_lat-dth_z*log(p_p0)*cos2_lat)*pow(p_p0,planet->kappa);
    ans       = MAX(200.,ans);
  }
  else if (strcmp(planet->name,"Venus_LLR05") == 0) {
    /*
     *  Interpolate on data table:
     */
    neglogp = -log(pressure);
    if (neglogp >= t_table[var.ntp-1].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[ var.ntp-1].y;
      dt_tp = dt_table[var.ntp-1].y;
    }
    else if (neglogp <= t_table[0].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[0].y;
      dt_tp = dt_table[0].y;
    }
    else {
      ki   = find_place_in_table(var.ntp,t_table,neglogp,&p_d);
      t_tp  = splint_pchip(neglogp, t_table+ki,p_d);
      dt_tp = splint_pchip(neglogp,dt_table+ki,p_d);
    }
    /*
     * Venus forcing of Lee, Lewis, and Read (2005, Adv. Space Res. 36, 2142-2145).
     */
    ans = t_tp+dt_tp*(cos(latitude*DEG)-0.642);
  }
  else if (strcmp(planet->name,"Titan") == 0) {
    /*
     *  Interpolate on data table:
     */
    neglogp = -log(pressure);
    if (neglogp >= t_table[var.ntp-1].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[ var.ntp-1].y;
      dt_tp = dt_table[var.ntp-1].y;
    }
    else if (neglogp <= t_table[0].x) {
      /* past range of t_vs_p data */
      t_tp  = t_table[0].y;
      dt_tp = dt_table[0].y;
    }
    else {
      ki   = find_place_in_table(var.ntp,t_table,neglogp,&p_d);
      t_tp  = splint_pchip(neglogp, t_table+ki,p_d);
      dt_tp = splint_pchip(neglogp,dt_table+ki,p_d);
    }
    /*
     * Titan forcing as in Flasar et al. (1981, Nature 292, 693-698).
     */
    ans = t_tp+dt_tp*cos(latitude*DEG);
  }
  else {
    sprintf(Message,"not yet implemented for %s",planet->name);
    epic_error(dbmsname,Message);
  }
   
  return ans;
}

/*======================= end of temp_eq() ==================================*/

/*====================== rt_heating() =======================================*/

/*
 * This function connects the EPIC model to the full radiative transfer code
 * of Greathouse et al., which is located in the subdirectory epic/src/rt.
 *
 * NOTE: The heating arrays take into account the net effects of both
 *       heating and cooling processes, and have units [W/kg].
 */

#undef  SHORTWAVE_HEATING
#define SHORTWAVE_HEATING(k,j,i) shortwave_heating[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#undef  LONGWAVE_HEATING
#define LONGWAVE_HEATING(k,j,i)  longwave_heating[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/*
 * Set N_ZONAL_SAMPLE, the number of evenly strided samples in the zonal direction,
 * for economical zonal averaging when grid.zonal_average_rt is turned on.
 * 
 * NOTE: Set to grid.ni for no effect.
 */
#undef  N_ZONAL_SAMPLE
#define N_ZONAL_SAMPLE 32

void rt_heating(planetspec *planet,
                int         action)
{
  int
    K,J,I,
    i_stride,i_count,
    update_shortwave,
    update_longwave;
  EPIC_FLOAT
    avg;
  static int
    warned = FALSE;
  static EPIC_FLOAT
   *shortwave_heating,
   *longwave_heating;
  extern void rt_shortwave(planetspec *planet,
                           EPIC_FLOAT *heating,
                           int         i_stride,
                           int         action);
  extern void rt_longwave(planetspec *planet,
                          EPIC_FLOAT *heating,
                          int         i_stride,
                          int         action);
#if defined(EPIC_MPI)
  EPIC_FLOAT
    mpi_tmp;
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="rt_heating";

  switch(action) {
    case EPIC_ALLOC:
      /*
       * Allocate memory
       */

      /* Start printing progress to stdout */
      if (IAMNODE == NODE0) {
        fprintf(stdout,"Setting up radiative transfer:"); fflush(stdout);
      }

      rt_longwave(planet,NULL,0,action);
      longwave_heating = fvector(0,Nelem3d-1,dbmsname);

      rt_shortwave(planet,NULL,0,action);
      shortwave_heating = fvector(0,Nelem3d-1,dbmsname);

      /* Finish printing progress to stdout */
      if (IAMNODE == NODE0) {
        fprintf(stdout,"; RT ready\n"); fflush(stdout);
      }
    break;

    case EPIC_FREE:
      /*
       * Free allocated memory
       */
      rt_longwave( planet,NULL,0,action);
      free_fvector(longwave_heating, 0,Nelem3d-1,dbmsname);

      rt_shortwave(planet,NULL,0,action);
      free_fvector(shortwave_heating,0,Nelem3d-1,dbmsname);

    break;

    case EPIC_APPLY:
      /*
       * Determine whether to update the shortwave and longwave
       * heating arrays, or use the values from the previous timestep.
       */
      update_shortwave = TRUE;
      update_longwave  = TRUE;

      /*
       * Set i_stride, and count number of points along I direction to be sampled in zonal averaging.
       */
      i_stride = IMAX(1,grid.ni/N_ZONAL_SAMPLE);

      i_count = 0;
      for (I = grid.ilo; I <= grid.ni; I += i_stride) {
        i_count++;
      }

      /*
       * Apply radiative-transfer heating to timestep.
       */
      if (update_shortwave) {
        if (grid.zonal_average_rt) {
          /*
           * Sample N_ZONAL_SAMPLE longitudes spanning the domain, evenly strided.
           *
           * NOTE: We are averaging over longitude of a snapshot, rather than over time.
           */

          if (grid.ni < N_ZONAL_SAMPLE && warned == FALSE && IAMNODE == NODE0) {
            sprintf(Message,"grid.ni=%d < %d with grid.zonal_average_rt == TRUE; recommend grid.ni >= %d",grid.ni,N_ZONAL_SAMPLE,N_ZONAL_SAMPLE);
            epic_warning(dbmsname,Message);

            warned = TRUE;
          }

          rt_shortwave(planet,shortwave_heating,i_stride,action);

          /*
           * Average over i_stride points.
           */
          for (K = 0; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              /*
               * Start at grid.ilo so that i_stride is corret.
               */

              avg = 0.;
              for (I = grid.ilo; I <= IHI; I += i_stride) {
                if (I < ILO) {
                  continue;
                }
                avg += SHORTWAVE_HEATING(K,J,I);
              }

#if defined(EPIC_MPI)
              mpi_tmp = avg;
              MPI_Allreduce(&mpi_tmp,&avg,1,float_type,MPI_SUM,para.comm_JLO);
#endif

              avg /= i_count;

              for (I = ILO; I <= IHI; I++) {
                SHORTWAVE_HEATING(K,J,I) = avg;
              }
            }
          }
          /* No need to apply bc_lateral() here. */
        }
        else {
          /*
           * Do every longitude point.
           */
          i_stride = 1;
          rt_shortwave(planet,shortwave_heating,i_stride,action);
        }
      }

      if (update_longwave) {
        if (grid.zonal_average_rt) {
          /*
           * Sample N_ZONAL_SAMPLE longitudes spanning the domain, evenly strided.
           *
           * NOTE: We are averaging over longitude of a snapshot, rather than over time.
           */
          rt_longwave(planet,longwave_heating,i_stride,action);

          /*
           * Average over i_stride points.
           */
          for (K = 0; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {

              avg = 0.;
              for (I = grid.ilo; I <= IHI; I += i_stride) {
                if (I < ILO) {
                  continue;
                }
                avg += LONGWAVE_HEATING(K,J,I);
              }

#if defined(EPIC_MPI)
              mpi_tmp = avg;
              MPI_Allreduce(&mpi_tmp,&avg,1,float_type,MPI_SUM,para.comm_JLO);
#endif

              avg /= i_count;

              for (I = ILO; I <= IHI; I++) {
                LONGWAVE_HEATING(K,J,I) = avg;
              }
            }
          }
          /* No need to apply bc_lateral() here. */
        }
        else {
          /*
           * Do every longitude point.
           */
          i_stride = 1;
          rt_longwave(planet,longwave_heating,i_stride,action);
        }
      }

      for (K = 0; K <= KHI; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            HEAT3(K,J,I) += SHORTWAVE_HEATING(K,J,I)+LONGWAVE_HEATING(K,J,I);
          }
        }
      }
      /* Apply bc_lateral() */
      bc_lateral(var.heat3.value,THREEDIM);
    break;

    default:
      sprintf(Message,"unrecognized action=%d",action);
      epic_error(dbmsname,Message);
    break;
  }
  
  return;
}

/*====================== end of rt_heating() ================================*/

/************************ end of epic_sensible_heating.c *********************/
