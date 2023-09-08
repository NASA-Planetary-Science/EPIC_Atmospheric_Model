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

/* * * * * * * * * * epic_timestep.c * * * * * * * * * * * * * * * * 
 *                                                                 *
 * Integrate the prognostic variables ahead one timestep.          *
 *                                                                 *
 * This file includes the following:                               *
 *                                                                 *
 *           timestep()                                            *
 *           adams_bashforth_step()                                *
 *           leapfrog_step()                                       *
 *           uv_core()                                             *
 *           uv_pgrad()                                            *
 *           uv_sponge()                                           *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>
#include <epic_pv_schemes.h>

/*
 * Local function prototypes.
 */
void uv_pgrad_traditional(planetspec  *planet,
                          int          Kstart,
                          int          Kend,
                          EPIC_FLOAT **Buff2D);

void uv_pgrad_green_gauss(planetspec *planet,
                          int         Kstart,
                          int         Kend);


/*======================= timestep() ========================================*/

/*
 * March the prognostic variables ahead one timestep using the 
 * specified schemes.
 *
 * If action == STEP_PROGS_AND_UPDATE_DIAGS, we assume the diagnostic variables
 * are all up-to-date upon entry.
 *
 * If action == UPDATE_DIAGS_ONLY, we bring the diagnostic variables up-to-date
 * with the prognostic variables, but do not advance the progostic variables.
 * This is used to prime the pump at startup.
 */

void timestep(planetspec  *planet,
              int          action,
              EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    i,itmp,
    is,ip,iq,
    shift;
  register EPIC_FLOAT
    tmp;
  static unsigned long
    nbytes_2d,
    nbytes_3d;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="timestep";

  if (!initialized) {
    /* Allocate memory. */
    nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);
    nbytes_3d = Nelem3d*sizeof(EPIC_FLOAT);

    initialized = TRUE;
  }

  if (grid.itime%100 == 0) {
    /*
     * Occasionally screen model variables for nan.
     * The input NULL implies silent mode.
     */
    check_nan(NULL);
  }

  if (action == STEP_PROGS_AND_SYNC_DIAGS) {
    /*
     *  Clear u,v tendencies for current time.
     */
    memset(var.u.tendency+grid.it_uv_tend*Nelem3d,0,nbytes_3d);
    memset(var.v.tendency+grid.it_uv_tend*Nelem3d,0,nbytes_3d);

    /*
     * Apply sources and sinks to h-grid variables.
     */
    source_sink(planet,Buff2D);

    /*
     * The HEAT array is now complete with the current latent heating.
     * Update the hybrid vertical velocity, W3, and DZDT2.
     */
    calc_w(planet,action);

    /*
     * Apply sources and sinks to turbulence variables.
     */
    if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
        strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
      source_sink_turb(planet,Buff2D);
    }
    else if (strcmp(grid.turbulence_scheme,"off") == 0) {
      ;
    }
    else {
      sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
      epic_error(dbmsname,Message);
    }

    /* * * * * * * * * * * * * * * * * * * * * * * *
     *                                             *
     *  Diagnostic variables are now synchronized. *
     *  Write to extract.nc if appropriate.        *
     *                                             *
     * * * * * * * * * * * * * * * * * * * * * * * */

    if (var.extract_on) {
      /* 
       * Write data to extract file periodically. 
       */
      if ((grid.itime)%(grid.itextract) == 0) {
        /*
         * The counter var.extract_time_index is incremented
         * after each use here.
         */
        if (grid.extract_append[0] != '\0') {
          if (grid.itime > 0) {
            var_write(planet,grid.extract_append,EXTRACT_DATA,var.extract_time_index++,0);
          }
        }
        else {
          var_write(planet,"extract.nc",EXTRACT_DATA,var.extract_time_index++,0);
        }
      }
      /*
       * Return after writing last extract frame.
       */
      if (grid.itime == grid.itrun) return;
    }

    /*
     * Calculate Coriolis and advection tendencies for (u,v).
     */
    uv_core(planet,Buff2D);

    /*
     * Add sponge-layer drag to (u,v) tendencies.
     * For the leapfrog timestep, lagging is used for stability.
     */
    uv_sponge(planet);

    /*
     * Apply subgrid-scale model to scalar variables.
     */
    scalar_vertical_subgrid(planet,Buff2D);
    scalar_horizontal_subgrid(planet,Buff2D);

    /*
     * Advect scalar prognostic variables forward one timestep.
     */
    advection(planet,Buff2D);

    /*
     * Apply high-latitude, low-pass filter to scalar prognostic variables.
     * Ensure positive-definite status as appropriate.
     */
    for (K = KLO; K <= KHI; K++){
      shift = (K-Kshift)*Nelem2d;

      if (var.h.on) {
        zonal_filter(H_INDEX,var.h.value+shift);
        restore_mass(planet,H_INDEX,NO_PHASE);
      }

      if (var.theta.on) {
        zonal_filter(THETA_INDEX,var.theta.value+shift);
        restore_mass(planet,THETA_INDEX,NO_PHASE);
      }

      if (var.fpara.on) {
        zonal_filter(FPARA_INDEX,var.fpara.value+shift);
        restore_mass(planet,FPARA_INDEX,NO_PHASE);
      }

      if (var.nu_turb.on) {
        zonal_filter(NU_TURB_INDEX,var.nu_turb.value+shift);
        restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
      }

      for (iq = 0; iq < grid.nq; iq++) {
        zonal_filter(grid.is[iq],var.species[grid.is[iq]].phase[grid.ip[iq]].q+shift);
        restore_mass(planet,grid.is[iq],grid.ip[iq]);
      }
    }

    /*
     * Update the diagnostic variables needed for uv_pgrad().
     *
     * NOTE: Calling store_pgrad_vars() starts the calculation of DZDT2, which
     *       is finished in calc_w().  This means it would be wrong to call
     *       store_pgrad_vars() again in a timestep after DZDT2 is finished.
     *
     * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
     *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
     *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
     */
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,action,CALC_PHI3NK);
    uv_pgrad(planet,Buff2D);

    /*
     * Add horizontal subgrid-scale model to wind tendencies.
     * The vertical turbulence model is applied implicitly on (u,v) below.
     */
    uv_horizontal_subgrid(planet,Buff2D);

    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      /*
       * March u,v forward one timestep.
       */
      adams_bashforth_step(planet,U_INDEX);
      adams_bashforth_step(planet,V_INDEX);

      /*
       *  Cycle time index backwards:
       */
      itmp      = IT_MINUS2;
      IT_MINUS2 = IT_MINUS1;
      IT_MINUS1 = IT_ZERO;
      IT_ZERO   = itmp;
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      /*
       * March u,v forward one timestep.
       */
      leapfrog_step(planet,U_INDEX);
      leapfrog_step(planet,V_INDEX);

      /* 
       * Cycle time index backwards.
       */
      itmp      = IT_MINUS1;
      IT_MINUS1 = IT_ZERO;
      IT_ZERO   = itmp;
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme = %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }

    timeplane_bookkeeping();

    /*
     * Apply vertical turbulence model to U,V.
     * This is done implicitly to handle thin layers.
     */
    uv_vertical_subgrid(planet,Buff2D);

    /*
     * Apply hyperviscosity to U,V.
     */
    uv_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D);

    /*
     * Apply high-latitude, low-pass filter to prevent numerical instability.
     */
    for(K = KLO;K <= KHI; K++){
      zonal_filter(U_INDEX,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d);
      zonal_filter(V_INDEX,var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d);
    }

    /*
     *  Advance time:
     */
    var.model_time += (time_t)grid.dt;

    /*
     * Update solar longitude, L_s.
     */
    L_s = solar_longitude(planet,var.model_time);

    /*
     * Store commonly used diagnostic variables not calculated elsewhere.
     */
    store_diag(planet);

    /*
     * Start the calculation of HEAT, in W/kg, on the layer interfaces (all except
     * the latent heating).
     */
    calc_heating(planet);
  } 
  else if (action == SYNC_DIAGS_ONLY) {
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,action,CALC_PHI3NK);

    /*
     * Store commonly used diagnostic variables not calculated elsewhere.
     */
    store_diag(planet);

    /*
     * Calculate HEAT, in W/kg, on the layer interfaces.
     *
     * NOTE: This currently does not include latent heating from cloud microphysics.
     */
    calc_heating(planet);

    /*
     * Calculate the hybrid vertical velocity, W, on layer interfaces.
     */
    calc_w(planet,action);
  }
  else {
    sprintf(Message,"unrecognized action=%d",action);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of timestep() =================================*/

/*======================= adams_bashforth_step() ============================*/

/*
 * Use the 3rd-order Adams-Bashforth timestep.
 * This timestep is discussed by D. Durran (1991, MWR 119, 702-720).
 * It is appropriate for dissipative terms as well as for
 * conservative terms, does not suffer from the time-splitting
 * numerical instability of the leap-frog timestep, and is more
 * accurate than the leapfrog timestep.  Its main drawback is that
 * it requires more memory because it uses two previous time derivatives.
 */

void adams_bashforth_step(planetspec *planet,
                          int         index)
{
  register int
    K,J,I,
    jlo;
  wind_variable
    *wind;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="adams_bashforth_step";

  /*
   * The index jlo, as opposed to JLO, arises because U and V are on different grids
   * in the Arakawa C-grid scheme.
   */
  if (index == U_INDEX) {
    wind = &var.u;
    jlo  = JLO;
  }
  else if (index == V_INDEX) {
    wind = &var.v;
    jlo  = JFIRST;
  }
  else {
    sprintf(Message,"not set up to handle index=%d",index);
    epic_error(dbmsname,Message);
  }

  /* 
   * Specify Adams-Bashforth coefficients.
   * The FLOAT_MAX flags are set in epic_initial.c. 
   *
   * NOTE: Use jlo rather than JLO for this test.
   */
  if (DWINDDT(wind,IT_MINUS2,KLO,jlo,ILO) == FLOAT_MAX) {
    if (DWINDDT(wind,IT_MINUS1,KLO,jlo,ILO) == FLOAT_MAX) {
      /* 
       * Use 1st-order Adams-Bashforth, aka forward (Euler) difference,
       * for the initial step.
       */
      grid.ab[0] = 1.;
      grid.ab[1] = 0.;
      grid.ab[2] = 0.;
    }
    else {
      /*
       * Use 2nd-order Adams-Bashforth for the second step.
       */
      grid.ab[0] =  3./2.;
      grid.ab[1] = -1./2.;
      grid.ab[2] =     0.;
    }
  }
  else {
   /*
    * Use 3rd-order Adams-Bashforth for the third step and beyond.
    */
    grid.ab[0] =  23./12.;
    grid.ab[1] = -16./12.;
    grid.ab[2] =   5./12.;
  }

  /*
   * Advance variable.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = jlo; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        WIND(wind,grid.it_uv,K,J,I) += DT*( grid.ab[0]*DWINDDT(wind,IT_ZERO,  K,J,I)
                                           +grid.ab[1]*DWINDDT(wind,IT_MINUS1,K,J,I)
                                           +grid.ab[2]*DWINDDT(wind,IT_MINUS2,K,J,I) );
      }
    }
  }
  /*
   * Tie K = KLO-1 values to K = KLO values.
   * NOTE: Assumes K dimension is not cut.
   */
  K = KLO-1;
  for (J = jlo; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      WIND(wind,grid.it_uv,K,J,I) = WIND(wind,grid.it_uv,K+1,J,I);
    }
  }
  /* Need to apply bc_lateral() here. */
  if (index == U_INDEX) {
    bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  }
  else if (index == V_INDEX) {
    bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);
  }
  else {
    sprintf(Message,"unrecognized index=%d",index);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of adams_bashforth_step() =====================*/

/*======================= leapfrog_step() ===================================*/

/*
 * Apply the leapfrog step. 
 * Use Asselin filtering to remove the odd-even instability.
 *
 * Durran (1991) cites Williamson (1983) for the value gamma = 0.06
 * for the Asselin filter parameter.
 */

#define GAMMA_ASSELIN 0.06

void leapfrog_step(planetspec *planet,
                   int         index)
{
  register int
    K,J,I,
    jlo;
  register EPIC_FLOAT
    twodt,
    old_filtered,
    present,
    new_leap;
  wind_variable
    *wind;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="leapfrog_step";

  /*
   * The index jlo, as opposed to JLO, arises because U and V are on different grids
   * in the Arakawa C-grid scheme.
   */
  if (index == U_INDEX) {
    wind = &var.u;
    jlo  = JLO;
  }
  else if (index == V_INDEX) {
    wind = &var.v;
    jlo  = JFIRST;
  }
  else {
    sprintf(Message,"not set up to handle index=%d",index);
    epic_error(dbmsname,Message);
  }

  /* 
   * Use a forward (Euler) difference for the initial step.
   */

  /*
   * NOTE: Need jlo, not JLO for this test, otherwise the flagged V=FLOAT_MAX will be V=0.
   */
  if (WIND(wind,IT_MINUS1,KLO,jlo,ILO) == FLOAT_MAX) {
    for (K = KLO; K <= KHI; K++) {
      for (J = jlo; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          present                    = WIND(wind,IT_ZERO,K,J,I);
          WIND(wind,IT_MINUS1,K,J,I) = present+DT*DWINDDT(wind,grid.it_uv_tend,K,J,I);
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    if (index == U_INDEX) {
      bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
    }
    else if (index == V_INDEX) {
      bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized index=%d",index);
      epic_error(dbmsname,Message);
    }

    return;
  }

  twodt = 2.*DT;

  /*
   * Take leapfrog step, store in old timeframe. 
   * Asselin-filter the current step.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = jlo; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        old_filtered               = WIND(wind,IT_MINUS1,K,J,I);
        present                    = WIND(wind,IT_ZERO,  K,J,I);
        new_leap                   = old_filtered+twodt*DWINDDT(wind,grid.it_uv_tend,K,J,I);
        WIND(wind,IT_MINUS1,K,J,I) = new_leap;
        WIND(wind,IT_ZERO,  K,J,I) = present*(1.-2.*GAMMA_ASSELIN)+GAMMA_ASSELIN*(old_filtered+new_leap);
      }
    }
  }
  /*
   * Tie K = KLO-1 values to K = KLO values.
   * NOTE: Assumes K dimension is not cut.
   */
  K = KLO-1;
  for (J = jlo; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      WIND(wind,IT_MINUS1,K,J,I) = WIND(wind,IT_MINUS1,K+1,J,I);
      WIND(wind,IT_ZERO,  K,J,I) = WIND(wind,IT_ZERO,  K+1,J,I);
    }
  }
  /* Need to apply bc_lateral() here. */
  if (index == U_INDEX) {
    bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
    bc_lateral(var.u.value+IT_ZERO*Nelem3d,  THREEDIM);
  }
  else if (index == V_INDEX) {
    bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
    bc_lateral(var.v.value+IT_ZERO*Nelem3d,  THREEDIM);
  }
  else {
    sprintf(Message,"unrecognized index=%d",index);
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of leapfrog_step() ============================*/

/*======================= uv_core() =========================================*/

/*
 * Calculate core tendencies for u,v for the current state.
 * These include advection and Coriolis acceleration.
 * The pressure-gradient terms are calculated separately as part
 * of the economical explicit scheme. 
 */

void uv_core(planetspec  *planet,
             EPIC_FLOAT **Buff2D)
{
  register int    
    K,J,I,
    kk,jj;
  unsigned long
    nbytes_2d;
  register EPIC_FLOAT
    al, be,       /* Used in AL_U, BE_U, etc. macros.         */
    ga, de,       /* See Arakawa and Lamb (1981) eqn. (3.34)  */
    ep1,ep2,      /*       "                     "            */
    ph1,ph2,      /*       "                     "            */
    m_2jp1,
    m_2j_inv,
    n_2j,n_2jp1_inv,
    havg,
    d1,d2,d1d2,davg;
  EPIC_FLOAT
    *uh,*vh,*kin;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_core";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

 /*
  * Calculate the horizontal advection and Coriolis terms for U and V.
  *
  * NOTE: In order to get the correct Coriolis terms, use the same type of 
  *       hybrid density to form UH and VH here as is used in the 
  *       potential vorticity.
  */ 
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;

    memset(Buff2D[0],0,nbytes_2d);
    memset(Buff2D[1],0,nbytes_2d);
    memset(Buff2D[2],0,nbytes_2d);
    uh  = Buff2D[0];      
    vh  = Buff2D[1];      
    kin = Buff2D[2];  

    for (J = JLO; J <= JHI; J++) {
      n_2jp1_inv = 1./grid.n[kk][2*J+1];
      for (I = ILO; I <= IHI; I++) {
        havg    = .5*(H(K,J,I)+H(K,J,I-1));
        UH(J,I) = U(grid.it_uv,K,J,I)*havg*n_2jp1_inv;
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(uh,TWODIM);

    for (J = JFIRST; J <= JHI; J++) {
      m_2j_inv  = 1./grid.m[kk][2*J];
      for (I = ILO; I <= IHI; I++) {
        havg    = .5*(H(K,J,I)+H(K,J-1,I));
        VH(J,I) = V(grid.it_uv,K,J,I)*havg*m_2j_inv;
      }
    }
    /* Need bc_lateral() here. */
    bc_lateral(vh,TWODIM);

    /*
     *  Calculate kinetic energy per unit mass, kin:
     */
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        KIN(J,I) = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                  var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,kk,J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(kin,TWODIM);

    /*
     * Add Coriolis, horizontal advection, and kin terms to dvdt.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      n_2j = grid.n[kk][jj];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) += ((-GA_V*UH(J  ,I+1)-DE_V*UH(J  ,I  )
                                         -AL_V*UH(J-1,I  )-BE_V*UH(J-1,I+1)
                                         +PH2*VH( J-1,I  )-PH1*VH( J+1,I  ))*PV_COEF
                                       +KIN(J-1,I)-KIN(J,I))*n_2j;
      }
      if (grid.include_nontrad_accel) {
        for (I = ILO; I <= IHI; I++) {
          DVDT(grid.it_uv_tend,K,J,I) -= .5*(DZDT2(K,J,I)+DZDT2(K,J-1,I))*(V(grid.it_uv,K,J,I)/grid.rlt[kk][jj]);
        }
      }
    }

    /*
     * Add Coriolis, horizontal advection, and kin terms to dudt.
     *
     * The coefficients AL_U, BE_U, etc. are linear combinations
     * of PV(K,J,I) and are defined in $EPIC_PATH/include/epic_pv_schemes.h.
     *
     * grid.f2 is 2*Omega*cos(lat).
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      m_2jp1 = grid.m[kk][jj];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) += ( (AL_U*VH(J+1,I  )+BE_U*VH(J+1,I-1)
                                         +GA_U*VH(J  ,I-1)+DE_U*VH(J,  I  )
                                         +EP2*UH( J,  I-1)-EP1*UH( J,  I+1))*PV_COEF
                                       +KIN(J,I-1)-KIN(J,I) )*m_2jp1;
      }
      if (grid.include_nontrad_accel) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(grid.it_uv_tend,K,J,I) -= .5*(DZDT2(K,J,I)+DZDT2(K,J,I-1))*(grid.f2[jj]+U(grid.it_uv,K,J,I)/grid.rlt[kk][jj]);
        }
      }
    }
  }

  /*
   * Calculate vertical advection terms for U and V.
   */
  uv_vertical_advection(planet);

  return;
}

/*======================= end of uv_core() ======================================*/

/*======================= uv_pgrad() ============================================*/

/*
 * Calculate pressure-gradient tendencies for (u,v).
 *
 * A choice of algorithms is available for the horizontal pressure-gradient
 * force: uv_pgrad_traditional() and uv_pgrad_green_gauss().
 */

void uv_pgrad(planetspec  *planet,
              EPIC_FLOAT **Buff2D)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    K,
    idbms=0;
  static char
    dbmsname[]="uv_pgrad";

  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
    case COORD_ISOBARIC:
      uv_pgrad_traditional(planet,KLO,KHI,Buff2D);
    break;
    case COORD_HYBRID:
      /* --gurn green-gauss--
       * Ideally, we would use the green-gauss (finite volume) algorithm, since
       * it is designed to reduce spurious velocities in steep terrain. However, 
       * in this version of EPIC it is associated with the formation of unstable
       * polar ripples in V and U that are approximately barotropic. We have not
       * yet been able to isolate what is causing this problem, and so have
       * reverted to using the traditional pressure-gradient force algorithm.
       */
      /* uv_pgrad_green_gauss(planet,KLO,KHI); */
      uv_pgrad_traditional(planet,KLO,KHI,Buff2D);
    break;
    default:
      sprintf(Message,"Need to specify a pressure-gradient force algorithm for grid.coord_type=%d",
                       grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }  
}

/*====================== end of uv_pgrad() ========================================*/

/*======================= uv_pgrad_traditional() ==================================*/

/*
 * An implementation of the horizontal pressure-gradient force (PGF) terms
 * via a finite-difference discretization of the PDE form of the PGF, applied
 * to layers Kstart to Kend.
 */

#define FPARA2(j,i) fpara2[i+(j)*Iadim-Shift2d]

void uv_pgrad_traditional(planetspec  *planet,
                          int          Kstart,
                          int          Kend,
                          EPIC_FLOAT **Buff2D)
{
  register int    
    K,J,I,
    kk,k_isen;
  register EPIC_FLOAT
    m_2jp1,n_2j;
  EPIC_FLOAT
   *fpara2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_pgrad_traditional";

  switch(grid.coord_type) {
    case COORD_ISOBARIC:
      k_isen = KLO-1;
    break;
    case COORD_ISENTROPIC:
      k_isen = grid.nk;
    break;
    case COORD_HYBRID:
      /*
       * For a change in algorithm at the seam, set k_isen = grid.k_sigma-1.
       *
       * To use the same algorithm for all the layers, set k_isen = grid.nk;
       * This is what Konor and Arakawa (1997) do.
       */
      k_isen = grid.nk;
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  if (var.fpara.on) {
    fpara2 = Buff2D[0];
  }

  /*
   * Check validity of layer range.
   */
  if (Kstart > Kend) {
    K      = Kstart;
    Kstart = Kend;
    Kend   = K;
  }
  if (Kstart < 1 || Kend > grid.nk) {
    sprintf(Message,"K=%d not implemented",K);
    epic_error(dbmsname,Message);
  }

  for (K = Kstart; K <= k_isen; K++) {
    kk = 2*K;

    if (var.fpara.on) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FPARA2(J,I) = .5*(FPARA(K,J,I)+FPARA(K-1,J,I));
        }
      }
      bc_lateral(fpara2,TWODIM);
    }

    for (J = JLO; J <= JHI; J++) {
      m_2jp1  = grid.m[kk][2*J+1];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= m_2jp1*(MONT2(K,J,I)-MONT2(K,J,I-1)-
                                              .5*(EXNER2(K,J,I)+EXNER2(K,J,I-1))*(THETA2(K,J,I)-THETA2(K,J,I-1)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB2(K,J,I)+FGIBB2(K,J,I-1))
                                           *(FPARA2(J,I)-FPARA2(J,I-1))*m_2jp1;
        }
      }
    }

    for (J = JFIRST; J <= JHI; J++) {
      n_2j = grid.n[kk][2*J];
      for (I = ILO; I <= IHI; I++) {

        DVDT(grid.it_uv_tend,K,J,I) -= n_2j*(MONT2(K,J,I)-MONT2(K,J-1,I)-
                                            .5*(EXNER2(K,J,I)+EXNER2(K,J-1,I))*(THETA2(K,J,I)-THETA2(K,J-1,I)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DVDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB2(K,J,I)+FGIBB2(K,J-1,I))
                                           *(FPARA2(J,I)-FPARA2(J-1,I))*n_2j;
        }
      }
    }
  }

  for (K = k_isen+1; K <= Kend; K++) {
    kk = 2*K;

    /*
     * NOTE: We were using eqns (2.22) and (3.22) from Arakawa and Konor (1996, MWR), but had problems with
     *       the derivative with respect to sigma across the seam in the hybrid case.
     *
     *       In practice, (1/rho) style p.g.f. terms seem to be more prone to numerical instability than
     *       those written in terms of the Exner function. 
     */

    if (var.fpara.on) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FPARA2(J,I) = .5*(FPARA(K,J,I)+FPARA(K-1,J,I));
        }
      }
      bc_lateral(fpara2,TWODIM);
    }

    for (J = JLO; J <= JHI; J++) {
      m_2jp1  = grid.m[kk][2*J+1];
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= m_2jp1*(     PHI2(  K,J,I)-PHI2(  K,J,I-1)
                                               +.5*(THETA2(K,J,I)+THETA2(K,J,I-1))
                                                  *(EXNER2(K,J,I)-EXNER2(K,J,I-1)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DUDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB2(K,J,I)+FGIBB2(K,J,I-1))
                                           *(FPARA2(J,I)-FPARA2(J,I-1))*m_2jp1;
        }
      }
    }

    for (J = JFIRST; J <= JHI; J++) {
      n_2j = grid.n[kk][2*J];
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= n_2j*(     PHI2(  K,J,I)-PHI2(  K,J-1,I)
                                             +.5*(THETA2(K,J,I)+THETA2(K,J-1,I))
                                                *(EXNER2(K,J,I)-EXNER2(K,J-1,I)));
      }
      if (var.fpara.on) {
        for (I = ILO; I <= IHI; I++) {
          DVDT(grid.it_uv_tend,K,J,I) -= .5*(FGIBB2(K,J,I)+FGIBB2(K,J-1,I))
                                           *(FPARA2(J,I)-FPARA2(J-1,I))*n_2j;
        }
      }
    }
  }

  return;
}

#undef fpara2

/*======================= end of uv_pgrad_traditional() ============================*/

/*======================= uv_pgrad_green_gauss() ===================================*/

/*
 * Calculate the horizontal pressure-gradient force (PGF) components and add to DUDT and DVDT,
 * for layers Kstart to Kend.
 *
 * Use a Green-Gauss, finite-volume approach with a 3D control volume (CV) around each
 * U point and each V point.  The 8 corners of each CV can have different PHI values,
 * and the surface forces and mass are computed accordingly.  This method expands the
 * 2D Green-Gauss approach of Lin (1997, QJRMS 123, 1749-1762) to 3D. See
 *   Bradley ME, Dowling TE, 2012, Using 3D finite volume for the pressure gradient force
 *     in atmospheric models, Q.J.R. Meteorol. Soc. 138, 2126-2135, doi:10.1002/qj.1929
 * The Mathematica notebooks used to derive the formulas are included in
 * epic/tools/mathematica.
 *
 * NOTE: If the map factors vary with altitude, then there will be a slight difference
 *       in the calculated values for a horizontal face that is shared between two
 *       stacked finite volumes, depending on whether it is treated as the top of the
 *       lower volume or the bottom of the upper one; this is a consequence of the
 *       fact that this algorithm has been derived using the plane-parallel assumption.
 */

#define  RHO_HAT_EAST(k,j,i)  rho_hat_east[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define RHO_HAT_NORTH(k,j,i) rho_hat_north[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void uv_pgrad_green_gauss(planetspec *planet,
                          int         Kstart,
                          int         Kend)
{ 
  register int    
    K,J,I,
    kk,jj;
  EPIC_FLOAT
    pbsw,pbse,pbnw,pbne,        /* Given pressures at the 8 corners */
    ptsw,ptse,ptnw,ptne,
    pbw,pbe,ptw,pte,            /* Calculated pressures, center of edges */
    pbs,pbn,pts,ptn,
    pi_bw,pi_be,pi_tw,pi_te,    /* Change in pressure from center to corner along edge */
    pi_bs,pi_bn,pi_ts,pi_tn,
    phi_bsw,phi_bse,phi_bnw,
    phi_bne,phi_bw,phi_be,      /* Geopotentials at the 8 corners */
    phi_tsw,phi_tse,phi_tnw,
    phi_tne,phi_tw,phi_te,      /* Calculated geop. at center of edges */
    phi_bs,phi_bn,phi_ts,phi_tn,
    dbw2,dbe2,dtw2,dte2,        /* Square of "deltas" for the geopotentials:  change from center to corner. (Lon)*/
    dtn,dbn,dts,dbs,dtw,dte,    /* Deltas for geopotentials:  Latitudinal case isn't squared */
    dbw,dbe,
    rlt,rlnn,rlns,              /* Map factors: lat = const., lon = linear from south to north */
    rbar,delta_r,               /* Average and difference between south & north map factors */
    lne,lnw,lts,ltn,            /* longitude and latitude limits on the control volume */
    D_u,N_u,                    /* Denominator and numerator terms, from zonal component, as in Bradley & Dowling (2011) */
    pgflon,pgflat;              /* Pressure gradient forces, longitudinal and latitudinal */
  EPIC_FLOAT
    S_v,N_v,E_v,W_v,
    T_v,B_v,M_v,                /*  South, north, east, west, top, bottom and mass terms in Meridional PGF */
    rho_hat_top,
    rho_hat_bot;
  static int
    j_periodic  = FALSE,
    initialized = FALSE;
  static EPIC_FLOAT
   *rho_hat_east,
   *rho_hat_north;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_pgrad_green_gauss";

  if (!initialized) {
    /* Allocate memory */
    rho_hat_east  = fvector(0,Nelem3d-1,dbmsname);
    rho_hat_north = fvector(0,Nelem3d-1,dbmsname);

    if (strcmp(grid.geometry,"f-plane") == 0 &&
        strcmp(grid.f_plane_map,"cartesian") == 0) {
      j_periodic = TRUE;
    }

    initialized = TRUE;
  }

  /*
   * Check validity of layer range.
   */
  if (Kstart > Kend) {
    K      = Kstart;
    Kstart = Kend;
    Kend   = K;
  }
  if (Kstart < 1 || Kend > grid.nk) {
    sprintf(Message,"K=%d not implemented",K);
    epic_error(dbmsname,Message);
  }

  /*------------------------------------*
   * Latitudinal (meridional) component *
   *------------------------------------*/

  /*
   * Compute RHO_HAT_NORTH(K,J,I) and RHO_HAT_EAST(K,J,I) on v-grid.
   */
  for (K = IMAX(KLO,Kstart-1); K <= IMIN(KHI,Kend+1); K++) {
    if (IS_SPOLE) {
      /*
       * NOTE: For J = 0 we need RHO_HAT_NORTH(K,J,I), but not RHO_HAT_EAST(K,J,I).
       */
      J = 0;

      /* 
       * Prime I loop with east values.
       */
      I = ILO-1;

      ptne = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I+1));
      pbne = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I+1));

      phi_tne = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I+1));
      phi_bne = .5*(PHI3(K,  J,I)+PHI3(K,  J,I+1));

      for (I = ILO; I <= IHI; I++) {
        /* 
         * Reuse shared-face values.
         */
        ptnw    = ptne;
        pbnw    = pbne;

        phi_tnw = phi_tne;
        phi_bnw = phi_bne;

        ptne    = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I+1));
        pbne    = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I+1));

        phi_tne = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I+1));
        phi_bne = .5*(PHI3(K,  J,I)+PHI3(K,  J,I+1));

        ptn    = .5*(ptne+ptnw);
        pbn    = .5*(pbne+pbnw);

        phi_tn = .5*(phi_tne+phi_tnw);
        phi_bn = .5*(phi_bne+phi_bnw);

        RHO_HAT_NORTH(K,J,I) = (pbn-ptn)/(phi_tn-phi_bn);
      }
    }

    for (J = JFIRST; J <= JHI; J++) {
      /* 
       * Prime I loop with east values.
       */
      I = ILO-1;

      ptne = .5*(EXNER3(K-1,J,  I)+EXNER3(K-1,J,  I+1));
      ptse = .5*(EXNER3(K-1,J-1,I)+EXNER3(K-1,J-1,I+1));  
      pbne = .5*(EXNER3(K,  J,  I)+EXNER3(K,  J,  I+1));
      pbse = .5*(EXNER3(K,  J-1,I)+EXNER3(K,  J-1,I+1));  

      phi_tne = .5*(PHI3(K-1,J,  I)+PHI3(K-1,J,  I+1));
      phi_tse = .5*(PHI3(K-1,J-1,I)+PHI3(K-1,J-1,I+1));  
      phi_bne = .5*(PHI3(K,  J,  I)+PHI3(K,  J,  I+1));
      phi_bse = .5*(PHI3(K,  J-1,I)+PHI3(K,  J-1,I+1));  

      for (I = ILO; I <= IHI; I++) {
        /* 
         * Reuse shared-face values.
         */
        ptnw    = ptne;
        pbnw    = pbne;
        phi_tnw = phi_tne;
        phi_bnw = phi_bne;

        ptne = .5*(EXNER3(K-1,J,  I)+EXNER3(K-1,J,  I+1));
        ptse = .5*(EXNER3(K-1,J-1,I)+EXNER3(K-1,J-1,I+1));
        pbne = .5*(EXNER3(K,  J,  I)+EXNER3(K,  J,  I+1));
        pbse = .5*(EXNER3(K,  J-1,I)+EXNER3(K,  J-1,I+1));

        phi_tne = .5*(PHI3(K-1,J,  I)+PHI3(K-1,J,  I+1));
        phi_tse = .5*(PHI3(K-1,J-1,I)+PHI3(K-1,J-1,I+1));
        phi_bne = .5*(PHI3(K,  J,  I)+PHI3(K,  J,  I+1));
        phi_bse = .5*(PHI3(K,  J-1,I)+PHI3(K,  J-1,I+1));

        ptn  = .5*(ptne+ptnw);
        pte  = .5*(ptne+ptse);
        pbn  = .5*(pbne+pbnw);
        pbe  = .5*(pbne+pbse);

        phi_tn  = .5*(phi_tne+phi_tnw);
        phi_te  = .5*(phi_tne+phi_tse);
        phi_bn  = .5*(phi_bne+phi_bnw);
        phi_be  = .5*(phi_bne+phi_bse);

        RHO_HAT_NORTH(K,J,I) = (pbn-ptn)/(phi_tn-phi_bn);
        RHO_HAT_EAST( K,J,I) = (pbe-pte)/(phi_te-phi_be);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(rho_hat_north+(K-Kshift)*Nelem2d,TWODIM);
    bc_lateral(rho_hat_east+ (K-Kshift)*Nelem2d,TWODIM);
  }

  /* Top and bottom values */
  if (Kstart == KLO) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        RHO_HAT_NORTH(KLO-1,J,I) = RHO_HAT_NORTH(KLO,J,I);
        RHO_HAT_EAST( KLO-1,J,I) = RHO_HAT_EAST( KLO,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  if (Kend == KHI) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        RHO_HAT_NORTH(KHI+1,J,I) = RHO_HAT_NORTH(KHI,J,I);
        RHO_HAT_EAST( KHI+1,J,I) = RHO_HAT_EAST( KHI,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  for (K = Kstart; K <= Kend; K++) {
    kk = 2*K;
    for (I = ILO; I <= IHI; I++) {
      /* 
       * Prime J loop with north values.
       */
      J  = JFIRST-1;
      jj = 2*J;

      /*
       * Scale rln with planet->re to reduce difference-of-large-numbers issues.
       */
      rlnn = grid.rln[kk][jj+1]/planet->re;
      ltn  = grid.lat[jj+1]*DEG;

      /*
       * Calculate 4 north corners for p and phi.
       */
      ptnw  = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I-1));
      ptne  = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I+1));
      pbnw  = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I-1));
      pbne  = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I+1));

      phi_tnw  = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I-1));
      phi_tne  = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I+1));
      phi_bnw  = .5*(PHI3(K,  J,I)+PHI3(K,  J,I-1));
      phi_bne  = .5*(PHI3(K,  J,I)+PHI3(K,  J,I+1));

      /*  Compute pressures at the center of the edges of the CV */
      ptn = .5*(ptne+ptnw);
      pbn = .5*(pbne+pbnw);

      /* NOTE: "delta pressures" are not needed for meridional component. */

      /*  Compute geopotentials at the center of the edges of the CV */
      phi_tn = .5*(phi_tne+phi_tnw);
      phi_bn = .5*(phi_bne+phi_bnw);

      /*  Compute "delta geopotentials" for the edges of the CV */
      dtn = .5*(phi_tne-phi_tnw);
      dbn = .5*(phi_bne-phi_bnw);
      N_v = rlnn*(3.*(phi_tn-phi_bn)*(ptn+pbn)-(dtn*dtn-dbn*dbn)*(RHO_HAT_NORTH(K,J,I)));

      for (J = JFIRST; J <= JHI; J++) {
        jj    = 2*J;
        rlt   = grid.rlt[kk][jj];

        /* 
         * Reuse shared-face values.
         */
        rlns    = rlnn;
        lts     = ltn;
        ptsw    = ptnw;
        ptse    = ptne;
        pbsw    = pbnw;
        pbse    = pbne;
        phi_tsw = phi_tnw;
        phi_tse = phi_tne;
        phi_bsw = phi_bnw;
        phi_bse = phi_bne;
        pts     = ptn;
        pbs     = pbn;
        phi_ts  = phi_tn;
        phi_bs  = phi_bn;
        dts     = dtn;
        dbs     = dbn;
        S_v     = N_v;

        /* Calculate new North values */

        /*
         * Scale rln with planet->re to reduce difference-of-large-number issues.
         */
        rlnn  = grid.rln[kk][jj+1]/planet->re;
        ltn   = grid.lat[jj+1]*DEG;

        /*
         * Calculate 4 north corners for p and phi.
         */
        ptnw  = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I-1));
        ptne  = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J,I+1));
        pbnw  = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I-1));
        pbne  = .5*(EXNER3(K,  J,I)+EXNER3(K,  J,I+1));

        phi_tnw  = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I-1));
        phi_tne  = .5*(PHI3(K-1,J,I)+PHI3(K-1,J,I+1));
        phi_bnw  = .5*(PHI3(K,  J,I)+PHI3(K,  J,I-1));
        phi_bne  = .5*(PHI3(K,  J,I)+PHI3(K,  J,I+1));

        /*  Compute pressures at the center of the edges of the CV */
        ptn = .5*(ptne+ptnw);
        ptw = .5*(ptsw+ptnw);
        pte = .5*(ptse+ptne);
        pbn = .5*(pbne+pbnw);
        pbw = .5*(pbsw+pbnw);
        pbe = .5*(pbse+pbne);

        /* NOTE: "delta pressures" are not needed for meridional component. */

        /*  Compute geopotentials at the center of the edges of the CV */
        phi_tn = .5*(phi_tne+phi_tnw);
        phi_tw = .5*(phi_tsw+phi_tnw);
        phi_te = .5*(phi_tse+phi_tne);
        phi_bn = .5*(phi_bne+phi_bnw);
        phi_bw = .5*(phi_bsw+phi_bnw);
        phi_be = .5*(phi_bse+phi_bne);

        /*  Compute "delta geopotentials" for the edges of the CV */
        dtn = .5*(phi_tne-phi_tnw);
        dtw = .5*(phi_tnw-phi_tsw);
        dte = .5*(phi_tne-phi_tse);
        dbn = .5*(phi_bne-phi_bnw);
        dbw = .5*(phi_bnw-phi_bsw);
        dbe = .5*(phi_bne-phi_bse);

        /*  Compute rho_hat_top and rho_hat_bot; see Bradley and Dowling (2011), eq.(12). */
        rho_hat_top = .25*(RHO_HAT_NORTH(K,J,I)+RHO_HAT_NORTH(K-1,J,I)+RHO_HAT_NORTH(K,J-1,I)+RHO_HAT_NORTH(K-1,J-1,I));
        rho_hat_bot = .25*(RHO_HAT_NORTH(K,J,I)+RHO_HAT_NORTH(K+1,J,I)+RHO_HAT_NORTH(K,J-1,I)+RHO_HAT_NORTH(K+1,J-1,I));

        delta_r = (rlns-rlnn);  

        /* Bradley & Dowling (2011), eqns.(17-22) */
        N_v = rlnn*(3.*(phi_tn-phi_bn)*(ptn+pbn)-(dtn*dtn-dbn*dbn)*(RHO_HAT_NORTH(K,J,I)));
        E_v = 0.5*delta_r*(3.*(phi_te-phi_be)*(pte+pbe)-(dte*dte-dbe*dbe)*(RHO_HAT_EAST(K,J,I  )));
        W_v = 0.5*delta_r*(3.*(phi_tw-phi_bw)*(ptw+pbw)-(dtw*dtw-dbw*dbw)*(RHO_HAT_EAST(K,J,I-1)));

        T_v = (phi_tn-phi_ts)*((2.*rlnn+rlns)*ptn+(rlnn+2.*rlns)*pts)-(dtn+dts)*(rlns*dtn-rlnn*dts)*(rho_hat_top);
        B_v = (phi_bn-phi_bs)*((2.*rlnn+rlns)*pbn+(rlnn+2.*rlns)*pbs)-(dbn+dbs)*(rlns*dbn-rlnn*dbs)*(rho_hat_bot);
        M_v = (2.*rlnn+rlns)*(pbn-ptn)+(rlnn+2.*rlns)*(pbs-pts);

        pgflat = (S_v-N_v-E_v-W_v+T_v-B_v)/(M_v*rlt*(ltn-lts));

        DVDT(grid.it_uv_tend,K,J,I) += pgflat;
      }
    }   
  }

  /*--------------------------------*
   * Longitudinal (zonal) component *
   *--------------------------------*/

  /*
   * Compute RHO_HAT_EAST(K,J,I) on u-grid.
   */
  for (K = IMAX(KLO,Kstart-1); K <= IMIN(KHI,Kend+1); K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /* 
         * Put appropriate monotonic function of p here; denoted pi(p)
         * in Bradley & Dowling (2012)
         */
        ptse = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J-1,I));
        ptne = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J+1,I));
        pbse = .5*(EXNER3(K,  J,I)+EXNER3(K,  J-1,I));
        pbne = .5*(EXNER3(K,  J,I)+EXNER3(K,  J+1,I));

        pte  = .5*(ptse+ptne);
        pbe  = .5*(pbse+pbne);

        phi_tse = .5*(PHI3(K-1,J,I)+PHI3(K-1,J-1,I));
        phi_tne = .5*(PHI3(K-1,J,I)+PHI3(K-1,J+1,I));
        phi_bse = .5*(PHI3(K,  J,I)+PHI3(K,  J-1,I));
        phi_bne = .5*(PHI3(K,  J,I)+PHI3(K,  J+1,I));

        phi_te  = .5*(phi_tse+phi_tne);
        phi_be  = .5*(phi_bse+phi_bne);

        RHO_HAT_EAST(K,J,I) = (pbe-pte)/(phi_te-phi_be);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(rho_hat_east+(K-Kshift)*Nelem2d,TWODIM);
  }

  /* Top and bottom values */
  if (Kstart == KLO) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        RHO_HAT_EAST(KLO-1,J,I) = RHO_HAT_EAST(KLO,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }
  if (Kend == KHI) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        RHO_HAT_EAST(KHI+1,J,I) = RHO_HAT_EAST(KHI,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  for (K = Kstart; K <= Kend; K++) {
    kk = 2*K;

    for (J = JLO; J <= JHI; J ++) {
      jj = 2*J+1;

      rlns    = grid.rln[kk][jj-1];
      rlnn    = grid.rln[kk][jj+1];
      /*  Compute rbar (average) and delta_r  */
      rbar    = .5*(rlns+rlnn);
      delta_r =    (rlns-rlnn);

      if (rbar == 0.) {
        /*
         * This can arise if grid.nj = 0, and can be skipped.
         */
        continue;
      }

      if (J == grid.jlo || J == grid.nj) {
        /*
         * Southern or northern edge of model; revert to Lin (1997) 2D algorithm,
         * to avoid 3D volumes that hang over the edge of model.
         */

        /* 
         * Prime I loop with east values.
         */
        I = ILO-1;

        lne = grid.lon[2*I+1]*DEG;

        /* 
         * Compute pressures at the center of the edges of the CV.
         */
        pte = EXNER3(K-1,J,I);
        pbe = EXNER3(K,  J,I);

        /* 
         * Compute geopotentials at the center of the edges of the CV.
         */
        phi_te = PHI3(K-1,J,I);
        phi_be = PHI3(K,  J,I);

        for (I = ILO; I <= IHI; I++) {
          /* 
           * Reuse shared-face values.
           */
          lnw    = lne;
          pbw    = pbe;
          ptw    = pte;
          phi_bw = phi_be;
          phi_tw = phi_te;

          /* 
           * Calculate east values.
           */
          lne = grid.lon[2*I+1]*DEG;

          /*
           * Compute pressures at the center of the edges of the CV.
           */
          pte = EXNER3(K-1,J,I);
          pbe = EXNER3(K,  J,I);

          /* 
           * Compute geopotentials at the center of the edges of the CV.
           */
          phi_te = PHI3(K-1,J,I);
          phi_be = PHI3(K,  J,I);

          pgflon = ((pbe-ptw)*(phi_bw-phi_te)+(pbw-pte)*(phi_tw-phi_be))
                   /(rbar*(lne-lnw)*((pbe-ptw)+(pbw-pte)));

          DUDT(grid.it_uv_tend,K,J,I) += pgflon;
        }
      }
      else {
        /*
         * Interior point; use Bradley & Dowling (2012) 3D algorithm.
         */

        /* 
         * Prime I loop with east values.
         */
        I = ILO-1;

        /*
         * Calculate 4 east corners for p and phi.
         */
        ptse = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J-1,I));
        ptne = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J+1,I));
        pbse = .5*(EXNER3(K,  J,I)+EXNER3(K,  J-1,I));
        pbne = .5*(EXNER3(K,  J,I)+EXNER3(K,  J+1,I));

        phi_tse = .5*(PHI3(K-1,J,I)+PHI3(K-1,J-1,I));
        phi_tne = .5*(PHI3(K-1,J,I)+PHI3(K-1,J+1,I));
        phi_bse = .5*(PHI3(K,  J,I)+PHI3(K,  J-1,I));
        phi_bne = .5*(PHI3(K,  J,I)+PHI3(K,  J+1,I));

        lne = grid.lon[2*I+1]*DEG;

        /* 
         * Compute pressures at the center of the edges of the CV.
         */
        pte = .5*(ptse+ptne);
        pbe = .5*(pbse+pbne);

        /* 
         * Compute "delta pressures" for the edges of the CV, after scaling.
         */
        pi_te = .5*(ptne-ptse);
        pi_be = .5*(pbne-pbse);

        /*  Compute geopotentials at the center of the edges of the CV */
        phi_te = .5*(phi_tse+phi_tne);
        phi_be = .5*(phi_bse+phi_bne);

        /*  Compute square "delta geopotentials" for the edges of the CV */
        dte2 = SQR(.5*(phi_tne-phi_tse));
        dbe2 = SQR(.5*(phi_bne-phi_bse));

        for (I = ILO; I <= IHI; I++) {
          /* 
           * Reuse shared-face values.
           */
          lnw    = lne;
          pbw    = pbe;
          ptw    = pte;
          pi_bw  = pi_be;
          pi_tw  = pi_te;
          phi_bw = phi_be;
          phi_tw = phi_te;
          dbw2   = dbe2;
          dtw2   = dte2;

          /*
           * Calculate 4 east corners for p and phi.
           */
          ptse = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J-1,I));
          ptne = .5*(EXNER3(K-1,J,I)+EXNER3(K-1,J+1,I));
          pbse = .5*(EXNER3(K,  J,I)+EXNER3(K,  J-1,I));
          pbne = .5*(EXNER3(K,  J,I)+EXNER3(K,  J+1,I));

          phi_tse = .5*(PHI3(K-1,J,I)+PHI3(K-1,J-1,I));
          phi_tne = .5*(PHI3(K-1,J,I)+PHI3(K-1,J+1,I));
          phi_bse = .5*(PHI3(K,  J,I)+PHI3(K,  J-1,I));
          phi_bne = .5*(PHI3(K,  J,I)+PHI3(K,  J+1,I));

          /* 
           * Calculate east values.
           */
          lne = grid.lon[2*I+1]*DEG;

          /*
           * Compute pressures at the center of the edges of the CV.
           */
          pte = .5*(ptse+ptne);
          pbe = .5*(pbse+pbne);

          /*
           * Compute "delta pressures" for the edges of the CV, after scaling.
           */
          pi_te = .5*(ptne-ptse);
          pi_be = .5*(pbne-pbse);

          /*  Compute geopotentials at the center of the edges of the CV */
          phi_te = .5*(phi_tse+phi_tne);
          phi_be = .5*(phi_bse+phi_bne);

          /*  Compute square "delta geopotentials" for the edges of the CV */
          dte2 = SQR(.5*(phi_tne-phi_tse));
          dbe2 = SQR(.5*(phi_bne-phi_bse));

          /*  Compute rho_hat_top and rho_hat_bot; see Bradley and Dowling (2011), eq.(6). */
          rho_hat_top = .25*(RHO_HAT_EAST(K,J,I)+RHO_HAT_EAST(K-1,J,I)+RHO_HAT_EAST(K,J,I-1)+RHO_HAT_EAST(K-1,J,I-1));
          rho_hat_bot = .25*(RHO_HAT_EAST(K,J,I)+RHO_HAT_EAST(K+1,J,I)+RHO_HAT_EAST(K,J,I-1)+RHO_HAT_EAST(K+1,J,I-1));

          /* 3D terms not present in the Lin (1997) 2D formulation. */
          N_u = -(dtw2-dbw2)*(RHO_HAT_EAST(K,J,I-1))+(dte2-dbe2)*(RHO_HAT_EAST(K,J,I))
                -(dte2-dtw2)*(rho_hat_top)          +(dbe2-dbw2)*(rho_hat_bot);

          D_u = -(delta_r/rbar)*(pi_be+pi_bw-pi_te-pi_tw);

          pgflon = ((pbe-ptw)*(phi_bw-phi_te)+(pbw-pte)*(phi_tw-phi_be)+N_u/3.)
                   /(rbar*(lne-lnw)*((pbe-ptw)+(pbw-pte)+D_u/6.));

          DUDT(grid.it_uv_tend,K,J,I) += pgflon;
        }
      }
    }
  }

  return;
}

#undef RHO_HAT_EAST
#undef RHO_HAT_NORTH

/*======================= end of uv_pgrad_green_gauss() =========================*/

/*======================= uv_sponge() ===========================================*/

/*
 * Add sponge layers using Rayleigh drag.
 *
 * For the Adams-Bashforth timestep, use the current state.
 * For the leapfrog timestep, use the previous state (IT_MINUS1), which is
 * called lagging (using the current state is numerically unstable).
 */

/*
 * Pick values for sponge stiffness, from 0. to 1. 
 */
#define TOP_SPONGE (0.1)
#define LATERAL_SPONGE (0.1)

void uv_sponge(planetspec *planet)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    tmp,
    t_sponge_inv;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
     max_nu_horizontal[2+1];
  wind_variable
    *u,*v;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_sponge";

  if (!initialized) {
    set_max_nu(max_nu_horizontal);

    initialized = TRUE;
  }

  /*
   * Sponge at the top of the model.
   *
   * Add Rayleigh friction to U and V, relative to no motion,
   * to dampen gravity-wave reflections at the model's top.
   */
  for (K = KLO; K <= grid.k_sponge; K++) {
    tmp          = (EPIC_FLOAT)(grid.k_sponge+1-K)/(grid.k_sponge);
    t_sponge_inv = TOP_SPONGE*max_nu_horizontal[0]*.5*(1.-cos(M_PI*tmp));
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DUDT(grid.it_uv_tend,K,J,I) -= t_sponge_inv*(U(grid.it_uv_dis,K,J,I)-0.);
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        DVDT(grid.it_uv_tend,K,J,I) -= t_sponge_inv*(V(grid.it_uv_dis,K,J,I)-0.);
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Sponges at the northern and southern (lateral) edges of the model.
   *
   * Add Rayleigh friction to V, relative to no motion,
   * to control numerical instabilities at channel walls.
   * Do not add Rayleigh friction to U.
   */
  if (grid.j_sponge > 0) {
    for (J = JFIRST; J <= JHI; J++) {
      /*
       * Southern sponge
       */
      tmp = (EPIC_FLOAT)(grid.j_sponge+1-J)/(grid.j_sponge);
      if (tmp > 0.) {
        t_sponge_inv = LATERAL_SPONGE*max_nu_horizontal[0]*.5*(1.-cos(M_PI*tmp));
        /* Avoid overlapping the top and lateral sponge */
        for (K = grid.k_sponge+1; K <= KHI; K++) {
          for (I = ILO; I <= IHI; I++) {
            DVDT(grid.it_uv_tend,K,J,I) -= t_sponge_inv*(V(grid.it_uv_dis,K,J,I)-0.);
          }
        }
      }
      /*
       * Northern sponge
       */
      tmp = (EPIC_FLOAT)(grid.j_sponge+1-(grid.nj-(J-grid.jfirst)))/(grid.j_sponge);
      if (tmp > 0.) {
        t_sponge_inv = LATERAL_SPONGE*max_nu_horizontal[0]*.5*(1.-cos(M_PI*tmp));
        /* Avoid overlapping the top and lateral sponge */
        for (K = grid.k_sponge+1; K <= KHI; K++) {
          for (I = ILO; I <= IHI; I++) {
            DVDT(grid.it_uv_tend,K,J,I) -= t_sponge_inv*(V(grid.it_uv_dis,K,J,I)-0.);
          }
        }
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  /*
   * For the Held-Suarez test case, add planetary boundary layer (PBL) Rayleigh drag.
   * Otherwise, the PBL is handled by epic_subgrid.c subroutines.
   */
  if (strcmp(planet->name,"Held_Suarez") == 0) {
    EPIC_FLOAT
      p2,pbot,amp,nu0;
    const EPIC_FLOAT
      held_suarez_pbl_nu0 = (1./(24.*60.*60.));

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pbot  = 0.5*(P3(grid.nk,J,I)+P3(grid.nk,J,I-1));
        for (K = KHI; K >= KLO; K--) {
          p2  = .5*(P2(K,J,I)+P2(K,J,I-1));
          amp = (p2/pbot-0.7)/(1.-0.7);
          if (amp < 0.) {
            break;
          }
          else {
            nu0                          = held_suarez_pbl_nu0*amp;
            DUDT(grid.it_uv_tend,K,J,I) -= nu0*WIND(&var.u,grid.it_uv_dis,K,J,I);
          }
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pbot  = 0.5*(P3(grid.nk,J,I)+P3(grid.nk,J-1,I));
        for (K = KHI; K >= KLO; K--) {
          p2  = .5*(P2(K,J,I)+P2(K,J-1,I));
          amp = (p2/pbot-0.7)/(1.-0.7);
          if (amp < 0.) {
            break;
          }
          else {
            nu0                          = held_suarez_pbl_nu0*amp;
            DVDT(grid.it_uv_tend,K,J,I) -= nu0*WIND(&var.v,grid.it_uv_dis,K,J,I);
          }
        }
      }
    }
  }

  return;
}

/*======================= end of uv_sponge() ===================================*/

/* * * * * * * * * * * * end of epic_timestep.c * * * * * * * * * * * * * * * * */
