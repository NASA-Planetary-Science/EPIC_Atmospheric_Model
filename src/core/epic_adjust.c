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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 *  This file contains the following functions:                    *
 *      restore_mass()                                             *
 *      zonal_filter()                                             *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*=================== restore_mass() ========================================*/

/*
 * We currently just maintain a floor of grid.h_min[K] for H and 0.0
 * for mixing-ratio variables like H_2O_liquid.
 *
 * Returns TRUE if any adjustment was made, otherwise returns FALSE.
 *
 * NOTE: Previously we distributed mass vertically, but we found that
 *       this tends to generate unwanted gravity waves.
 *       
 * NOTE: This function should not call set_p2_etc() itself.
 */

int restore_mass(planetspec *planet,
                 int         species_index,
                 int         phase_index)
{
  register int
    K,J,I;
  int
    changed=FALSE;
  register EPIC_FLOAT
    tol;
  EPIC_FLOAT
    *pt;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="restore_mass";
    
  if (species_index >= FIRST_SPECIES && species_index <= LAST_SPECIES) {
    tol = Q_MIN;
    for (K = KLOPAD; K <= KHIPAD; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pt = &Q(species_index,phase_index,K,J,I);
          /*
           * Introducing the truncation gradually proves to be beneficial.
           */
          if (*pt < 0.) {
            *pt     = tol;
            changed = TRUE;
          }
          else if (*pt < 2.*tol) {
            *pt     = tol+.5*(*pt);
            changed = TRUE;
          }
        }
      }
    }
  }
  else {
    switch(species_index) {
      case H_INDEX:
        for (K = KLOPAD; K <= KHIPAD; K++) {
          tol = grid.h_min[K];
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              pt = &H(K,J,I);
              /*
               * Introducing the truncation gradually proves to be beneficial.
               */
              if (*pt < 0.) {
                *pt     = tol;
                changed = TRUE;
              }
              else if (*pt < 2.*tol) {
                *pt     = tol+.5*(*pt);
                changed = TRUE;
              }
            }
          }
        }
      break;
      case NU_TURB_INDEX:
        tol = 1.e-4*planet->kinvisc;
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              pt = &NU_TURB(K,J,I);
              /*
               * Introducing the truncation gradually proves to be beneficial.
               */
              if (*pt < 0.) {
                *pt     = tol;
                changed = TRUE;
              }
              else if (*pt < 2.*tol) {
                *pt     = tol+.5*(*pt);
                changed = TRUE;
              }
            }
          }
        }
        /* No need to apply bc_lateral() here. */
      break;
      case THETA_INDEX:
        for (K = KLO; K <= KHI; K++) {
          tol = 1.e-4*grid.theta_ref[2*K+1];
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              pt = &THETA(K,J,I);
              /*
               * Introducing the truncation gradually proves to be beneficial.
               */
              if (*pt < 0.) {
                *pt     = tol;
                changed = TRUE;
              }
              else if (*pt < 2.*tol) {
                *pt     = tol+.5*(*pt);
                changed = TRUE;
              }
            }
          }
        }
      break;
      case FPARA_INDEX:
        tol = 0.;
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              pt = &FPARA(K,J,I);
              if (*pt < tol) {
                *pt     = tol;
                changed = TRUE;
              }
            }
          }
        }
      break;
      default:
        sprintf(Message,"unimplemented index=%d",species_index);
        epic_error(dbmsname,Message);
      break;
    }
  }

  return changed;
}

/*=================== end of restore_mass() =================================*/

/*======================= zonal_filter() ====================================*/

/* 
 * Filter to remove the CFL violation at the poles.
 * Input pointer buffji points to the relevant 2D plane.
 *
 * A negative index sets up re-initialization on the next call.
 *
 * NOTE: ni must be an integer power of 2, in order to use realft().
 *
 * NOTE: We used to have the zonal hyperviscosity here, but took
 *       it out because it was prone to Gibbs-effect ringing.
 */

/*
 * In the high-lat filter arrays, the index j is local, spanning JLO to JHI,
 * whereas the index i spans 1 to grid.ni/2+1.
 */
#define HIGH_LAT_H(j,i)   high_lat_h[  (i-1)+(j-JLO)*(grid.ni/2+1)]
#define HIGH_LAT_V(j,i)   high_lat_v[  (i-1)+(j-JLO)*(grid.ni/2+1)]
#define HIGH_LAT_DIV(j,i) high_lat_div[(i-1)+(j-JLO)*(grid.ni/2+1)]

void zonal_filter(int         index,
                  EPIC_FLOAT *buffji)
{
  register int 
    Kay,J,I,
    kk,jj;
  int
    itmp;
  static int 
    initialized = FALSE,
    has_pole    = FALSE,
    warned      = TRUE;   /* Set warned to FALSE here if a one-time warning is desired */
  static EPIC_FLOAT 
    *data,
    *sendbuf,
    *high_lat_h,
    *high_lat_v,
    *high_lat_div;
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
    dbmsname[]="zonal_filter";

  /*
   * No need to filter if ni = 1, or using cartesian geometry.
   */
  if (grid.ni == 1 || 
      (strcmp(grid.geometry,"f-plane")==0 && strcmp(grid.f_plane_map,"cartesian") == 0)) {
    return;
  }

  /*
   * Check whether grid.ni is a power of 2,
   * if not, return (and warn once).
   */
  if (frexp((double)grid.ni,&itmp) != 0.5) {
    /*
     * Warn once if grid.ni is not a power of 2.
     */
    if (warned == FALSE) {
      if (IAMNODE == NODE0) {
        sprintf(Message,"ni=%d is not a power of 2, %s() not applied",grid.ni,dbmsname);
        epic_warning(dbmsname,Message);
      }
      warned = TRUE;
    }
    return;
  }

  if (initialized && index < 0) {
    /*
     * Free allocated memory
     */
    free_fvector(data,   1,grid.ni,dbmsname);
    free_fvector(sendbuf,1,grid.ni,dbmsname);

    free_fvector(high_lat_h,  0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    free_fvector(high_lat_v,  0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    free_fvector(high_lat_div,0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);

    initialized = FALSE;
  }

  if (!initialized) {
    register int
      n,nn;
    EPIC_FLOAT 
      m0,n0,rln,rlt,lat0,
      re,rp,r,
      tmp0,tmp1,fac,al1,al2;
   
    /*
     * Allocate memory
     */
    data       = fvector(1,grid.ni,dbmsname);
    sendbuf    = fvector(1,grid.ni,dbmsname);

    high_lat_h   = fvector(0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    high_lat_v   = fvector(0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    high_lat_div = fvector(0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);

    lat0 = LAT0*DEG;

    /*
     * NOTE: The map factors (m,n) depend on K, but only weakly, so 
     *       here we flatten the filters to 2D by using the settings
     *       for K = grid.nk.
     */
    Kay = grid.nk;
    kk  = 2*Kay;

    if (strcmp(grid.geometry,"globe") == 0) {
      re   = grid.re[Kay];
      rp   = grid.rp[Kay];
      /* Filter is applied poleward of LAT0 */
      rln  = re/sqrt( 1.+ pow(rp/re*tan(lat0),2.) );
      rlt  = rln/( cos(lat0)*( pow(sin(lat0),2.) +
                   pow(re/rp*cos(lat0),2.) ) );
      m0 = 1./(rln*grid.dln*DEG);
      n0 = 1./(rlt*grid.dlt*DEG);

      has_pole = TRUE;
    }
    else if (strcmp(grid.geometry,"f-plane")  == 0 &&
             strcmp(grid.f_plane_map,"polar") == 0) {
      /* nj = 2*(nj+1-1)/2 */
      m0 = grid.m[kk][2*grid.nj/2];
      n0 = grid.n[kk][2*grid.nj/2];
      has_pole = TRUE;
    }
    else {
      sprintf(Message,"geometry not yet implemented");
      epic_error(dbmsname,Message);
    }

    /*
     * Increasing the latitude parameter, al1, and/or decreasing 
     * the wavenumber parameter,al2, makes the filter stronger.
     */

    /*
     * Filter for H-grid and U-grid.
     */
    al1 = 2.0;
    al2 = 1.0;
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      HIGH_LAT_H(J,1) = 1.;    
      if (has_pole) {   
        /* r is roughly 1.5 at a pole, 1 otherwise: */  
        r    = grid.mn[kk][jj]/(grid.m[kk][jj]*grid.n[kk][jj]);
        tmp0 = (grid.n[kk][jj]/n0)/(grid.m[kk][jj]/m0);
        tmp0 = pow(tmp0,al1);
        for (I = 2; I <= grid.ni/2+1; I++) {
          tmp1  = sin((I-1)*(grid.dln)*DEG/2.);
          tmp1  = pow(tmp1,-al2);
          tmp1 *= r*tmp0;
          tmp1  = (tmp1 < 1.) ? tmp1 : 1.;
          HIGH_LAT_H(J,I) = tmp1;
        }
      }
      else {
        for (I = 2; I <= grid.ni/2+1; I++) {
          HIGH_LAT_H(J,I) = 1.;
        }
      }
    }

    /*
     * Filter for V-grid and PV-grid.
     */
    al1 = 2.0;
    al2 = 1.0;
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      HIGH_LAT_V(J,1) = 1.; 
      if (has_pole) {      
        /* r is roughly 1.5 at a pole, 1 otherwise: */  
        r    = grid.mn[kk][jj]/(grid.m[kk][jj]*grid.n[kk][jj]);
        tmp0 = (grid.n[kk][jj]/n0)/(grid.m[kk][jj]/m0);
        tmp0 = pow(tmp0,al1);
        for (I = 2; I <= grid.ni/2+1; I++) {
          tmp1  = sin((I-1)*(grid.dln)*DEG/2.);
          tmp1  = pow(tmp1,-al2);
          tmp1 *= r*tmp0;
          tmp1  = (tmp1 < 1.) ? tmp1 : 1.;
          HIGH_LAT_V(J,I) = tmp1;
        }
      }
      else {
        for (I = 2; I <= grid.ni/2+1; I++) {
          HIGH_LAT_V(J,I) = 1.;
        }
      }
    }

    /*
     * Filter for divergence (H-grid).
     *
     * NOTE: We do this case separately because the numerics associated
     *       with keeping divergence damping stable in the polar regions may be different.
     *       The weaker the zonal filter, the more the damping is applied to the actual
     *       divergence field instead of a filtered version.
     */
    al1 = 2.0;
    al2 = 1.0;
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      HIGH_LAT_DIV(J,1) = 1.;    
      if (has_pole) {   
        /* r is roughly 1.5 at a pole, 1 otherwise: */  
        r    = grid.mn[kk][jj]/(grid.m[kk][jj]*grid.n[kk][jj]);
        tmp0 = (grid.n[kk][jj]/n0)/(grid.m[kk][jj]/m0);
        tmp0 = pow(tmp0,al1);
        for (I = 2; I <= grid.ni/2+1; I++) {
          tmp1  = sin((I-1)*(grid.dln)*DEG/2.);
          tmp1  = pow(tmp1,-al2);
          tmp1 *= r*tmp0;
          tmp1  = (tmp1 < 1.) ? tmp1 : 1.;
          HIGH_LAT_DIV(J,I) = tmp1;
        }
      }
      else {
        for (I = 2; I <= grid.ni/2+1; I++) {
          HIGH_LAT_DIV(J,I) = 1.;
        }
      }
    }

    if (index > 0) {
      initialized = TRUE;
    }
    else {
      /*
       * Negative index sets up re-initialization on next call.
       */
      index       = -index;
      initialized = FALSE;
    }
  } 
  /* End of initialization. */

  switch(index) {
    case V_INDEX:
    case PV2_INDEX:
    case REL_VORT2_INDEX:
    case ABS_VORT2_INDEX:
      /*
       * Variable defined on whole-integer J grid (v-grid or pv-grid).
       */
      for (J = JFIRST; J <= JHI; J++) {
        if (fabs(grid.lat[2*J]) < LAT0) {
          /*
           * Only apply for lat >= |LAT0|.
           */
          continue;
        }

#if defined(EPIC_MPI)
        /* zero sendbuf[] array, which is global in span, running from 1 to grid.ni */
        memset(sendbuf+1,0,grid.ni*sizeof(EPIC_FLOAT));

        /* Assign local values to sendbuf[] */
        for (I = ILO; I <= IHI; I++) {
          sendbuf[I] = BUFFJI(J,I);
        }

        /* 
         * Coadd across communicator para.comm_JLO to assign data[] on each processor,
         * which is global in span, running from 1 to grid.ni.
         *
         * NOTE: This is part of the inefficiency regarding cutting in the I direction with a zonal FFT.
         */
        MPI_Allreduce(sendbuf+1,data+1,grid.ni,float_type,MPI_SUM,para.comm_JLO);
#else
        for (I = 1; I <= grid.ni; I++) {
          data[I] = BUFFJI(J,I);
        }
#endif

        /*
         * NOTE: Doing the zonal fft redundantly on each processor is part of the inefficiency
         *       mentioned above. At least we do not have to broadcast the results.
         */
        realft(data,grid.ni,1);

        data[1] *= HIGH_LAT_V(J,1);
        for (I = 2; I <= (grid.ni)/2; I++) {
          data[2*I-1] *= HIGH_LAT_V(J,I);
          data[2*I  ] *= HIGH_LAT_V(J,I);
        }
        data[2] *= HIGH_LAT_V(J,grid.ni/2+1);

        realft(data,grid.ni,-1);

        for (I = ILO; I <= IHI; I++) {
          BUFFJI(J,I) = data[I]*(2./grid.ni);
        }
      }
    break;

    case DIV_UV2_INDEX:
      for (J = JLO; J <= JHI; J++) {
        if (fabs(grid.lat[2*J+1]) < LAT0) {
          /*
           * Only apply for lat >= |LAT0|.
           */
          continue;
        }

#if defined(EPIC_MPI)
        /* zero sendbuf[] array, which is global in span, running from 1 to grid.ni */
        memset(sendbuf+1,0,grid.ni*sizeof(EPIC_FLOAT));

        /* Assign local values to sendbuf[] */
        for (I = ILO; I <= IHI; I++) {
          sendbuf[I] = BUFFJI(J,I);
        }

        /* 
         * Coadd across communicator para.comm_JLO to assign the global data[] array on each processor.
         *
         * NOTE: This is part of the inefficiency regarding cutting in the I direction with a zonal FFT.
         */
        MPI_Allreduce(sendbuf+1,data+1,grid.ni,float_type,MPI_SUM,para.comm_JLO);
#else
        for (I = 1; I <= grid.ni; I++) {
          data[I] = BUFFJI(J,I);
        }
#endif

        /*
         * NOTE: Doing the zonal fft redundantly on each processor is part of the inefficiency
         *       mentioned above. At least we do not have to broadcast the results.
         */
        realft(data,grid.ni,1);

        data[1] *= HIGH_LAT_DIV(J,1);
        for (I = 2; I <= (grid.ni)/2; I++) {
          data[2*I-1] *= HIGH_LAT_DIV(J,I);
          data[2*I  ] *= HIGH_LAT_DIV(J,I);
        }
        data[2] *= HIGH_LAT_DIV(J,grid.ni/2+1);

        realft(data,grid.ni,-1);

        for (I = ILO; I <= IHI; I++) {
          BUFFJI(J,I) = data[I]*(2./grid.ni);
        }
      }
    break;

    case H_INDEX:
    case U_INDEX:
    default:
      for (J = JLO; J <= JHI; J++) {
        if (fabs(grid.lat[2*J+1]) < LAT0) {
          /*
           * Only apply for lat >= |LAT0|.
           */
          continue;
        }

#if defined(EPIC_MPI)
        /* zero sendbuf[] array, which is global in span, running from 1 to grid.ni */
        memset(sendbuf+1,0,grid.ni*sizeof(EPIC_FLOAT));

        /* Assign local values to sendbuf[] */
        for (I = ILO; I <= IHI; I++) {
          sendbuf[I] = BUFFJI(J,I);
        }

        /* 
         * Coadd across communicator para.comm_JLO to assign the global data[] array on each processor.
         *
         * NOTE: This is part of the inefficiency regarding cutting in the I direction with a zonal FFT.
         */
        MPI_Allreduce(sendbuf+1,data+1,grid.ni,float_type,MPI_SUM,para.comm_JLO);
#else
        for (I = 1; I <= grid.ni; I++) {
          data[I] = BUFFJI(J,I);
        }
#endif

        /*
         * NOTE: Doing the zonal fft redundantly on each processor is part of the inefficiency
         *       mentioned above. At least we do not have to broadcast the results.
         */
        realft(data,grid.ni,1);

        data[1] *= HIGH_LAT_H(J,1);
        for (I = 2; I <= (grid.ni)/2; I++) {
          data[2*I-1] *= HIGH_LAT_H(J,I);
          data[2*I  ] *= HIGH_LAT_H(J,I);
        }
        data[2] *= HIGH_LAT_H(J,grid.ni/2+1);

        realft(data,grid.ni,-1);

        for (I = ILO; I <= IHI; I++) {
          BUFFJI(J,I) = data[I]*(2./grid.ni);
        }
      }
    break;
  }

  /* Need to call bc_lateral() here. */
  bc_lateral(buffji,TWODIM);

  if (initialized == FALSE) {
    /*
     * Free allocated memory
     */
    free_fvector(data,   1,grid.ni,dbmsname);
    free_fvector(sendbuf,1,grid.ni,dbmsname);

    free_fvector(high_lat_h,  0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    free_fvector(high_lat_v,  0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
    free_fvector(high_lat_div,0,(JHI-JLO+1)*(grid.ni/2+1)-1,dbmsname);
  }

  return;
}

#undef HIGH_LAT_H
#undef HIGH_LAT_V
#undef HIGH_LAT_DIV

/*======================= end of zonal_filter() ==========================*/

/* * * * * * * * * * * * end of epic_adjust.c * * * * * * * * * * * * * * */













