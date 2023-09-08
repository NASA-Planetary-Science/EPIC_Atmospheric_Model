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

/* * * * * * * * * * epic_subgrid.c  * * * * * * * * * * * * * * * *
 *                                                                 *
 * TE Dowling, RP LeBeau, VK Parimi                                *
 *                                                                 *
 * Subgrid-scale closure subroutines for the EPIC model.           *
 *                                                                 *
 * This file contains the following subroutines:                   *
 *                                                                 *
 *  ____Called from outside___    ____Called from inside_____      *
 *                                                                 *
 *  make_arrays_subgrid()                                          *
 *                                                                 *
 *  free_arrays_subgrid()                                          *
 *                                                                 *
 *  init_subgrid()                                                 *
 *                                                                 *
 *  set_max_nu()                  max_nu_nondim()                  *
 *                                                                 *
 *  set_hyperviscosity()          set_max_nu()                     *
 *                                                                 *
 *  set_diffusion_coef()          dwall_SA()                       *
 *                                law_of_the_wall()                *
 *                                                                 *
 *  scalar_horizontal_subgrid()   scalar_horizontal_diffusion()    *
 *                                  laplacian_h()                  *
 *                                scalar_hyperviscosity()          *
 *                                  laplacian_h()                  *
 *                                                                 *
 *  scalar_vertical_subgrid()     scalar_vertical_diffusion()      *
 *                                  tau_surface()                  *
 *                                  stability_factor()             *
 *                                  crank_nicolson()               *
 *                                adiabatic_adjustment()           *
 *                                                                 *
 *  uv_horizontal_subgrid()       uv_horizontal_diffusion()        *
 *                                divergence_damping()             *
 *                                  divergence()                   *
 *                                                                 *
 *  uv_vertical_subgrid()         uv_vertical_diffusion()          *
 *                                  tau_surface()                  *
 *                                  crank_nicolson()               *
 *                                                                 *
 *  uv_hyperviscosity()           laplacian_uv()                   *
  *                                                                *
 *  source_sink_turb()            source_sink_SA()                 *
 *                                  dwall_SA()                     *
 *                                  delta_SA()                     *
 *                                  law_of_the_wall()              *
 *                                  invert_fv1()                   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Variables with scope of epic_subgrid.c.
 */
EPIC_FLOAT
  *d_wall;

/*=============== max_nu_nondim() =================================*/

EPIC_FLOAT max_nu_nondim(int order)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="max_nu_nondim";

  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {

    /*
     * Ray LeBeau's estimates of safe max nu's (Dowling et al, 1998) for
     * the 3rd-Order Adams-Bashforth timestep.
     */
    switch (order) {
      case 0: return 1./3.;
      case 2: return 1./30.;
      case 4: return 1./240.;
      case 6: return 1./800.;
      case 8: return 1./2400.;
      default:
        sprintf(Message,"hyperviscosity order %d not defined\n",order);
        epic_error(dbmsname,Message);
      break;
    }
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {

    /*
     * In the leapfrog-timestep case, we lag the viscosity terms with
     * a long forward step.
     *
     * order 0: 1/2 is mentioned in Durran (1991, eq(21)).
     * order 2: 1/8 is derived in Tannehill, Anderson, and Pletcher (1997, p. 137).
     *
     * We assume (1/8)^(order/2) in general.
     */
    switch (order) {
      case 0:  return .5;
      case 2:  return .125;
      default: return pow(.125,(double)(order/2));
    }
  }
  else {
    sprintf(Message,"unrecognized uv_timestep_scheme: %s \n",
                    grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }

  /* Should never get here.*/
  sprintf(Message,"should never get here");
  epic_error(dbmsname,Message);
  return 0.;
}

/*=============== end of max_nu_nondim() ==========================*/

/*=============== set_max_nu() ====================================*/

void set_max_nu(double *max_nu_horizontal)
{
  register int
    K,kk;
  register double
    dt,t2,
    dz,dz0,dz2,
    dx0,dx2;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_max_nu";

  dt = (double)(grid.dt);

  dz0 = FLOAT_MAX;
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    /*
     * Estimate layer temperature.
     */
    t2 = grid.t_ref[kk];
    /*
     * Estimate dz ~ H dlnp = (RT/g) dlnp and find minimum value.
     */
    dz  = log(grid.p_ref[kk+1]/grid.p_ref[kk-1])*planet->rgas*t2/grid.g[kk][2*grid.jtp+1];
    dz0 = MIN(dz0,dz);   
  }
  dz2 = dz0*dz0;

  if (grid.ni == 1 && grid.nj == 0) {
    /* Model is a 1D vertical column. */
    dx2 = dz2;
  }
  else {
    if (strcmp(grid.geometry,"globe") == 0) {
      dx0 = MIN(grid.dln*DEG*grid.re[grid.nk]/sqrt(1.+pow(grid.rp[grid.nk]/grid.re[grid.nk]*tan(LAT0*DEG),2.)),
                grid.dlt*DEG*SQR(grid.rp[grid.nk])/grid.re[grid.nk]);
    }
    else if (strcmp(grid.geometry,"f-plane")  == 0) { 
      if (strcmp(grid.f_plane_map,"polar") == 0) {
        dx0 = grid.f_plane_half_width/grid.nj;
      }
      else {
        dx0 = grid.f_plane_half_width/grid.ni;
      }
    }
    else {
      sprintf(Message,"unrecognized grid.geometry %s",grid.geometry);
      epic_error(dbmsname,Message);
    }
    dx2 = dx0*dx0;
  }

  max_nu_horizontal[0] = max_nu_nondim(0)/dt;
  max_nu_horizontal[2] = max_nu_nondim(2)*(dx2/dt);
  max_nu_horizontal[4] = max_nu_nondim(4)*(dx2/dt)*dx2;
  max_nu_horizontal[6] = max_nu_nondim(6)*(dx2/dt)*dx2*dx2;
  max_nu_horizontal[8] = max_nu_nondim(8)*(dx2/dt)*dx2*dx2*dx2;

  return;
}

/*============== end of set_max_nu() ==============================*/

/*============== set_hyperviscosity() =============================*/

/* 
 * Initialize hyperviscosity order and strength.
 *
 *  NOTE: Values of nu_nondim greater than (1./2.)^order are not useful, because the gridscale
 *        amplitude becomes negative. Values greater than 2.*(1./2)^order are numerically unstable.
 */

void set_hyperviscosity(void)
{
  char
    header[N_STR];
  register int
    ii;
  double  
    max_nu_horizontal[MAX_NU_ORDER+1];  /* NOTE: declared as double, not EPIC_FLOAT */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_hyperviscosity";

  set_max_nu(max_nu_horizontal);

  sprintf(Message,"Divergence damping coeff, fraction of max\n");
  grid.nudiv_nondim = input_float(Message,grid.nudiv_nondim);

  sprintf(header,"Hyperviscosity order [4");
  for (ii = 6; ii <= MAX_NU_ORDER; ii+=2) {
    sprintf(Message,",%d",ii);
    strcat(header,Message);
  }
  strcat(header,"; 0 => off]\n");
  grid.nu_order = input_int(header,grid.nu_order);

  if (grid.nu_order >= 4) {
    sprintf(Message,"nu%d, fraction of full strength\n",grid.nu_order);
    grid.nu_nondim = (grid.nu_nondim <= 0.) ?  0.5 : grid.nu_nondim;
    grid.nu_nondim = input_float(Message,grid.nu_nondim);
    grid.nu_hyper  = (double)grid.nu_nondim*max_nu_horizontal[grid.nu_order];
  }
  else {
    grid.nu_order  = 0;
    grid.nu_nondim = 0.;
    grid.nu_hyper  = 0.;
  }

  return;
}

/*============== end of set_hyperviscosity() ======================*/

/*============== scalar_hyperviscosity() ==========================*/

/*
 * Apply hyperviscosity to layer K of a scalar field on the h-grid.
 */

void scalar_hyperviscosity(int          nu_order,
                           double       nu_hyper,
                           EPIC_FLOAT **Buff2D,
                           int          kstart,
                           int          kend,
                           EPIC_FLOAT  *h)
{
  int
    K,J,I,kk,
    itmp,
    sign;
  register double
    tmp,rln,
    taper;
  EPIC_FLOAT
    *hh,
    *lphh,
    *diff_coef,
    *buff1,
    *buff2,
    *a,
    *ptmp;
  static double
    max_nu_horizontal[MAX_NU_ORDER+1];
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *m0;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_hyperviscosity";

  if (!initialized) {
    set_max_nu(max_nu_horizontal);

    /* Allocate memory */
    m0 = fvector(0,grid.nk,dbmsname);
    for (K = KLO; K <= KHI; K++) {
      rln   = grid.re[K]/sqrt(1.+SQR(grid.rp[K]/grid.re[K]*tan(LAT0*DEG)));
      m0[K] = 1./(rln*grid.dln*DEG);
    }

    initialized = TRUE;
  }

  hh        = Buff2D[0];
  lphh      = Buff2D[1];
  diff_coef = Buff2D[2];
  buff1     = Buff2D[3];
  buff2     = Buff2D[4];

  if (nu_order < 4) {
    return;
  }
  else if (nu_order > MAX_NU_ORDER) {
    sprintf(Message,"grid.nu_order=%d > MAX_NU_ORDER=%d\n",grid.nu_order,MAX_NU_ORDER);
    epic_error(dbmsname,Message);
  }
  if (nu_hyper < 0.) {
    sprintf(Message,"nu_hyper=%e < 0.",nu_hyper);
    epic_error(dbmsname,Message);
  }

  /*
   * Factor the hyperviscosity coefficient to
   * help prevent floating point overflow/underflow.
   */
  tmp = pow(nu_hyper,2./nu_order);

  for (K = kstart; K <= kend; K++) {
    kk = 2*K;

    for (itmp = 2; itmp <= grid.nu_order; itmp+=2) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          DIFF_COEF(J,I) = tmp;
        }
      }
      /* No need to call bc_lateral() here */
    }

    /* Point A(J,I) to input variable's K-layer */
    a = h+(K-Kshift)*Nelem2d;

    /* Copy A(J,I) into LPHH(J,I)  */
    memcpy(lphh,a,Nelem2d*sizeof(EPIC_FLOAT));

    sign = -1;
    for (itmp = 2; itmp <= nu_order; itmp+=2) {
      sign *= -1;
      ptmp  = hh;
      hh    = lphh;
      lphh  = ptmp;

      laplacian_h(kk,hh,diff_coef,lphh,buff1,buff2);
    }

    /* 
     * Apply hyperviscosity. Use a forward (Euler) step.
     */
    for (J = JLO; J <= JHI; J++) {
      /*
       * Taper viscosity coefficient to prevent numerical instability
       */
      taper = MIN(1.,(max_nu_horizontal[nu_order]/nu_hyper)*pow(m0[K]/grid.m[kk][2*J+1],nu_order));

      tmp   = DT*(EPIC_FLOAT)sign*taper;
      for (I = ILO; I <= IHI; I++) {
        A(J,I) += tmp*LPHH(J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(a,TWODIM);
  }

  return;
}

/*============== end of scalar_hyperviscosity() ===================*/

/*============== scalar_horizontal_subgrid() ======================*/

void scalar_horizontal_subgrid(planetspec  *planet,
                               EPIC_FLOAT **Buff2D)
{
  register int
    iq,
    kstart,kend;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_horizontal_subgrid";

  if (strcmp(grid.turbulence_scheme,"on") == 0) {
    scalar_horizontal_diffusion(planet,Buff2D);
  }
  else if (strcmp(grid.turbulence_scheme,"on_vertical_only") == 0 ||
           strcmp(grid.turbulence_scheme,"off")              == 0)  {
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

 /*
  * To avoid the sponge layers at the top of the model, set
  *   kstart = IMAX(KLO,grid.k_sponge+1);
  * otherwise, set
  *   kstart = KLO;
  */
  kstart = KLO;
  kend   = KHI;

  /*-----------------------*
   * Apply hyperviscosity. *
   *-----------------------*/

  /*
   * The iq loop is set up to reference only the species/phase fields that have been turned on.
   */
  if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
    for (iq = 0; iq < grid.nq; iq++) {
      scalar_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D,kstart,kend,var.species[grid.is[iq]].phase[grid.ip[iq]].q);
      restore_mass(planet,grid.is[iq],grid.ip[iq]);
    }
    sync_x_to_q(planet);
  }

  if (var.h.on) {
    switch(grid.coord_type) {
      case COORD_ISOBARIC:
        ;
      break;
      case COORD_ISENTROPIC:
      case COORD_HYBRID:
        /*
         * Apply hyperviscosity to H.
         */
        scalar_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D,kstart,kend,var.h.value);
        /* NOTE: restore_mass() is called in set_p2_etc() below */
      break;
      default:
        sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
        epic_error(dbmsname,Message);
      break;
    }
  }

  /*
   * Update P2, etc.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  if (var.theta.on) {
    switch(grid.coord_type) {
      case COORD_ISENTROPIC:
        ;
      break;
      case COORD_ISOBARIC:
      case COORD_HYBRID:
        scalar_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D,kstart,kend,var.theta.value);
      break;
      default:
        sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
        epic_error(dbmsname,Message);
      break;
    }
  }

  if (var.nu_turb.on) {
    scalar_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D,kstart,kend,var.nu_turb.value);
    restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
  }

  if (var.fpara.on) {
    scalar_hyperviscosity(grid.nu_order,grid.nu_hyper,Buff2D,kstart,kend,var.fpara.value);
  }

  return;
}

/*============== end of scalar_horizontal_subgrid() ===============*/

/*============== scalar_horizontal_diffusion ======================*/

void scalar_horizontal_diffusion(planetspec  *planet,
                                 EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,kk,
    iq,is,ip;
  register EPIC_FLOAT
    dt,
    rln,
    taper,
    tmp;
  const EPIC_FLOAT
    sigma_inv = 3./2.;
  EPIC_FLOAT
    *hh,
    *diff_coef,
    *laph,
    *a;
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *m0;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_horizontal_diffusion";

  if (!initialized) {
    /* Allocate memory */
    m0 = fvector(0,grid.nk,dbmsname);
    for (K = KLO; K <= KHI; K++) {
      rln   = grid.re[K]/sqrt(1.+SQR(grid.rp[K]/grid.re[K]*tan(LAT0*DEG)));
      m0[K] = 1./(rln*grid.dln*DEG);
    }

    initialized = TRUE;
  }

  hh        = Buff2D[0];
  laph      = Buff2D[1];
  diff_coef = Buff2D[2];

  dt = (EPIC_FLOAT)grid.dt;

  /*
   * Apply diffusion to THETA.
   *
   * NOTE: We have not included molecular diffusion for THETA.
   *       If it is added, may have to deal with
   *       density weighting and units of the transport coefficient.
   */
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      ;
    break;
    case COORD_ISOBARIC:
    case COORD_HYBRID:
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K+1;

        /* Copy DIFFUSION_COEF_THETA(K,J,I) into DIFF_COEF(J,I)  */
        a = var.diffusion_coef_theta.value+(K-Kshift)*Nelem2d;
        memcpy(diff_coef,a,Nelem2d*sizeof(EPIC_FLOAT));

        /* Copy THETA(K,J,I) into HH(J,I)  */
        a = var.theta.value+(K-Kshift)*Nelem2d;
        memcpy(hh,a,Nelem2d*sizeof(EPIC_FLOAT));

        laplacian_h(kk,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

        for (J = JLO; J <= JHI; J++) {
          /*
           * Old taper to help prevent numerical instability:
           *
           * taper = MIN(1.,pow(m0[K]/grid.m[kk][2*J+1],2.));
           */
          taper = 1.;

          tmp   = dt*taper;
          for (I = ILO; I <= IHI; I++) {
            THETA(K,J,I) += tmp*LAPH(J,I);
          }
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(var.theta.value,THREEDIM);
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * NOTE: Currently not applying diffusion to H.
   */

  /*
   * Apply diffusion to mixing ratios, Q, which are carried on the interfaces.
   */
  if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
    for (iq = 0; iq < grid.nq; iq++) {
      is = grid.is[iq];
      ip = grid.ip[iq];
      for (K = KLO; K <= KHI; K++) {
        kk = 2*K+1;

        /* Copy DIFFUSION_COEF_MASS(K,J,I) into DIFF_COEF(J,I)  */
        a = var.diffusion_coef_mass.value+(K-Kshift)*Nelem2d;
        memcpy(diff_coef,a,Nelem2d*sizeof(EPIC_FLOAT));

        if (grid.ip[iq] == VAPOR) {
          /*
           * Add molecular diffusion.
           */
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              DIFF_COEF(J,I) += mass_diffusivity(planet,is,T3(K,J,I),P3(K,J,I));
            }
          }
        }
        /* No need to apply bc_lateral() here. */

        /* Copy Q(is,ip,K,J,I) into HH(J,I)  */
        a = var.species[grid.is[iq]].phase[grid.ip[iq]].q+(K-Kshift)*Nelem2d;
        memcpy(hh,a,Nelem2d*sizeof(EPIC_FLOAT));
      
        laplacian_h(kk,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

        for (J = JLO; J <= JHI; J++) {
          /*
           * Old taper to help prevent numerical instability:
           *
           * taper = MIN(1.,pow(m0[K]/grid.m[kk][2*J+1],2.));
           */
          taper = 1.;

          tmp   = dt*taper;
          for (I = ILO; I <= IHI; I++) {
            Q(is,ip,K,J,I) += tmp*LAPH(J,I);
          }
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);

      /*
       * Clean up any negative mass introduced by diffusion truncation error.
       */
      restore_mass(planet,is,ip);
    }
    sync_x_to_q(planet);
  }

  /*
   * Update pressures, etc.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  /*
   * Apply diffusion to NU_TURB.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        DIFF_COEF(J,I) = sigma_inv*(planet->kinvisc+HH(J,I));
      }
    }
    /* No need to apply bc_lateral() here. */

    /* Copy NU_TURB(K,J,I) into HH(J,I)  */
    a = var.nu_turb.value+(K-Kshift)*Nelem2d;
    memcpy(hh,a,Nelem2d*sizeof(EPIC_FLOAT));

    laplacian_h(kk,hh,diff_coef,laph,Buff2D[3],Buff2D[4]);

    for (J = JLO; J <= JHI; J++) {
      /*
       * Old taper to help prevent numerical instability:
       *
       * taper = MIN(1.,pow(m0[K]/grid.m[kk][2*J+1],2.));
       */
      taper = 1.;

      tmp   = dt*taper;
      for (I = ILO; I <= IHI; I++) {
        NU_TURB(K,J,I) += tmp*LAPH(J,I);
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.nu_turb.value,THREEDIM);

  restore_mass(planet,NU_TURB_INDEX,NO_PHASE);

  return;
}

/*============== end of scalar_horizontal_diffusion ===============*/

/*============== laplacian_h() ====================================*/

/*
 *  Calculates the 2D Laplacian of an h-grid scalar field,
 *  multiplied by the input diffusion coefficient, DIFF_COEF(J,I).
 *
 *  The pointer hh should point to the appropriate JI plane.
 *
 *  Pointers to memory for two working JI-plane buffers
 *  are passed in as buff1 and buff2.
 */

void laplacian_h(int         kk,
                 EPIC_FLOAT *hh,
                 EPIC_FLOAT *diff_coef,
                 EPIC_FLOAT *lph,
                 EPIC_FLOAT *buff1,
                 EPIC_FLOAT *buff2)
{
  int 
    J,I;
  EPIC_FLOAT 
    h_edge,
    m_2j,n_2j,m_2jp1,n_2jp1,
    m_2j_inv,m_2jp2_inv;
  EPIC_FLOAT
    *gh1,*gh2;
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
    dbmsname[]="laplacian_h";

  /* 
   * Check that the pointer hh is not NULL.
   */
  if (!hh) {
    sprintf(Message,"input pointer hh=NULL");
    epic_error(dbmsname,Message);
  }

  /* Zero working buffers: */
  memset(buff1,0,Nelem2d*sizeof(EPIC_FLOAT));
  memset(buff2,0,Nelem2d*sizeof(EPIC_FLOAT));

  gh1 = buff1;
  gh2 = buff2;

  /* 
   * Compute grad h = (gh1,gh2).
   */
  for (J = JLO; J <= JHI; J++) {
    m_2jp1 = grid.m[kk][2*J+1];
    for (I = ILO; I <= IHI; I++) {
      GH1(J,I) = m_2jp1*(HH(J,I)-HH(J,I-1));    
    }
  }
  /* update gh1 edges below */
  for (J = JFIRST; J <= JHI; J++) {
    n_2j = grid.n[kk][2*J];
    for (I = ILO; I <= IHI; I++) {
      GH2(J,I) = n_2j*(HH(J,I)-HH(J-1,I));
    }
  }

  /* 
   * Fill in gh2 for top and bottom channel boundaries. 
   */
  if (JLO == grid.jlo && !IS_SPOLE) {
    /* southern edge */
    h_edge = 0.;
    for (I = ILO; I <= IHI; I++) {
      h_edge += HH(JLO,I);
    }

#if defined(EPIC_MPI)
    mpi_tmp = h_edge;
    MPI_Allreduce(&mpi_tmp,&h_edge,1,float_type,MPI_SUM,para.comm_JLO);
#endif

    h_edge /= grid.ni;

    n_2j = grid.n[kk][2*JLO];
    for (I = ILOPAD; I <= IHIPAD; I++) {
      GH2(JLO,I) = n_2j*(HH(JLO,I)-h_edge);
    }
  }
  if (JHI == grid.nj && !IS_NPOLE) {
    /* northern edge */
    h_edge = 0.;
    for (I = ILO; I <= IHI; I++) {
      h_edge += HH(JHI,I);
    }

#if defined(EPIC_MPI)
    mpi_tmp = h_edge;
    MPI_Allreduce(&mpi_tmp,&h_edge,1,float_type,MPI_SUM,para.comm_JLO);
#endif

    h_edge /= grid.ni;

    n_2j = grid.n[kk][2*(JHI+1)];
    for (I = ILOPAD; I <= IHIPAD; I++) {
      GH2(JHI+1,I) = n_2j*(h_edge-HH(JHI,I));
    }
  }
  /* update gh2 edges below */

  /* Update edges for gh1, gh2: */
  bc_lateral(gh1,TWODIM);
  bc_lateral(gh2,TWODIM);

  /*
   * Multiply by DIFF_COEF.
   */
  /* GH1 is on the U-grid. */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      GH1(J,I) *= .5*(DIFF_COEF(J,I)+DIFF_COEF(J,I-1));
    }
  }
  bc_lateral(gh1,TWODIM);

  /* GH2 is on the V-grid. */
  for (J = JFIRST; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      GH2(J,I) *= .5*(DIFF_COEF(J,I)+DIFF_COEF(J-1,I));
    }
  }
  bc_lateral(gh2,TWODIM);

  /* 
   * Compute laplacian.
   */
  for (J = JLO; J <= JHI; J++) {
    m_2jp1 = grid.m[kk][2*J+1];
    n_2jp1 = grid.n[kk][2*J+1];
    if (J == grid.jlo && IS_SPOLE) {
      m_2j_inv = 0.;
    }
    else {
      m_2j_inv = 1./grid.m[kk][2*J];
    }
    if (J == grid.nj && IS_NPOLE) {
      m_2jp2_inv = 0.;
    }
    else {
      m_2jp2_inv = 1./grid.m[kk][2*J+2];
    }
    for (I = ILO; I <= IHI; I++) {
      LPH(J,I) = m_2jp1*( (GH1(J,I+1)-GH1(J,I))
                  +n_2jp1*(GH2(J+1,I)*m_2jp2_inv-GH2(J,I)*m_2j_inv) );
    }
  }
  bc_lateral(lph,TWODIM);

  return;
}

/*============== end of laplacian_h() =============================*/

/*============== uv_hyperviscosity ================================*/

/*
 * Apply hyperviscosity to U,V to suppress computational modes.
 *
 * NOTE: The option to apply hyperviscosity to the velocity tendencies
 *       instead of their values was removed, since our experience has
 *       been that computational modes still grow when hyperviscosity
 *       is only applied to the tendencies.
 */

#undef  KLEN
#define KLEN (kend-kstart+1)

void uv_hyperviscosity(int          nu_order,
                       double       nu_hyper,
                       EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,kk,
    kstart,kend,
    itmp,
    sign;
  EPIC_FLOAT
    *uu,
    *vv,
    *lpuu,
    *lpvv,
    *buff1,
    *buff2,
    *ptmp;
  double
     visc_coef;
  register EPIC_FLOAT
     rln,taper,tmp;
  static double
    max_nu_horizontal[MAX_NU_ORDER+1];
  static int
    initialized = FALSE;
  static EPIC_FLOAT
    *m0;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_hyperviscosity";

  if (!initialized) {
    set_max_nu(max_nu_horizontal);

    /* Allocate memory */
    m0 = fvector(0,grid.nk,dbmsname);
    for (K = KLO; K <= KHI; K++) {
      rln   = grid.re[K]/sqrt(1.+SQR(grid.rp[K]/grid.re[K]*tan(LAT0*DEG)));
      m0[K] = 1./(rln*grid.dln*DEG);
    }

    initialized = TRUE;
  }

  if (nu_order < 4) {
    return;
  }
  else if (nu_order > MAX_NU_ORDER) {
    sprintf(Message,"grid.nu_order=%d > MAX_NU_ORDER=%d\n",grid.nu_order,MAX_NU_ORDER);
    epic_error(dbmsname,Message);
  }
  if (nu_hyper < 0.) {
    sprintf(Message,"nu_hyper=%e < 0.",nu_hyper);
    epic_error(dbmsname,Message);
  }

  uu    = Buff2D[0];
  vv    = Buff2D[1];
  lpuu  = Buff2D[2];
  lpvv  = Buff2D[3];
  buff1 = Buff2D[4];
  buff2 = Buff2D[5];

  /*
   * Factor hyperviscosity coefficient to avoid overflow/underflow.
   */
  visc_coef = pow(nu_hyper,2./nu_order);

  /*
   * To avoid the top sponge, set
   *   kstart = IMAX(KLO,grid.k_sponge+1);
   * otherwise, set
   *   kstart = KLO;
   */
  kstart = KLO;
  kend   = KHI;

  for (K = kstart; K <= kend; K++) {
    kk = 2*K;

    /*
     * Copy U and V into LPUU and LPVV.
     *
     * NOTE: Do not use grid.it_uv_dis here, which enables the lagged bookkeeping
     *       that yields numerical stability for the leapfrog scheme, since the hyperviscosity
     *       is applied directly to the variables, with a forward (Euler) step.
     */
    memcpy(lpuu,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    memcpy(lpvv,var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));

    sign = -1;
    for (itmp = 2; itmp <= nu_order; itmp+=2) {
      sign *= -1;
      ptmp  = uu;
      uu    = lpuu;
      lpuu  = ptmp;
      ptmp  = vv;
      vv    = lpvv;
      lpvv  = ptmp;

      laplacian_uv(K,uu,vv,visc_coef,lpuu,lpvv,buff1,buff2);
    }

    /*
     * Apply hyperviscosity. Use a forward (Euler) step.
     */
    for (J = JLO; J <= JHI; J++) {
      /*
       * Taper viscosity coefficient to prevent numerical instability
       */

      taper = MIN(1.,(max_nu_horizontal[nu_order]/nu_hyper)*pow(m0[K]/grid.m[kk][2*J+1],nu_order));

      tmp   = DT*(EPIC_FLOAT)sign*taper;
      for (I = ILO; I <= IHI; I++) {
        U(grid.it_uv,K,J,I) += tmp*LPUU(J,I);
      }
    }

    for (J = JFIRST; J <= JHI; J++) {
      /*
       * Taper viscosity coefficient to prevent numerical instability
       */
      taper = MIN(1.,(max_nu_horizontal[nu_order]/nu_hyper)*pow(m0[K]/grid.m[kk][2*J],nu_order));

      tmp   = DT*(EPIC_FLOAT)sign*taper;
      for (I = ILO; I <= IHI; I++) {
        V(grid.it_uv,K,J,I) += tmp*LPVV(J,I);
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

  return;
}

/*============== end of uv_hyperviscosity =========================*/

/*============== laplacian_uv() ===================================*/

/*
 *  Calculate the Laplacian of input 2D velocity vector (uu,vv)
 *  on the C-grid. See eq.(46) of Dowling et al. (1998).
 *
 *  The input viscosity coefficient is assumed to be a constant.
 *
 *  NOTE: If a spatially varying viscosity is needed, it is implemented 
 *        in uv_horizontal_diffusion() in terms of the stress tensor,
 *        following Tannehill, Anderson, and Pletcher (1997), p. 269,
 *        and that approach could be implemented.
 *
 *  The input pointers uu and vv should point to the appropriate JI plane.
 *
 *  Pointers to memory for two working JI-plane buffers
 *  are passed in as buff1 and buff2.
 */
void laplacian_uv(int         K,
                  EPIC_FLOAT *uu,
                  EPIC_FLOAT *vv,
                  EPIC_FLOAT  viscosity,
                  EPIC_FLOAT *lpuu,
                  EPIC_FLOAT *lpvv,
                  EPIC_FLOAT *buff1,
                  EPIC_FLOAT *buff2)
{
  int 
    J,I,
    kk = 2*K;
  EPIC_FLOAT 
    m_2j,n_2j,m_2jp1,n_2jp1;
  EPIC_FLOAT
    *ze,*di;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="laplacian_uv";

  /* 
   * Check that the input pointers are not NULL.
   */
  if (!uu) {
    sprintf(Message,"input pointer uu=NULL");
    epic_error(dbmsname,Message);
  }
  if (!vv) {
    sprintf(Message,"input pointer vv=NULL");
    epic_error(dbmsname,Message);
  }

  ze = buff1;
  di = buff2;

  /*
   * Calculate relative vorticity, ze, and horizontal divergence, di.
   */
  vorticity(ON_SIGMATHETA,RELATIVE,kk,uu,vv,NULL,ze);
  divergence(kk,uu,vv,di);

  /* Zero output arrays */
  memset(lpuu,0,Nelem2d*sizeof(EPIC_FLOAT));
  memset(lpvv,0,Nelem2d*sizeof(EPIC_FLOAT));

  /*
   * Compute zonal component of the Laplacian.
   */
  for (J = JLO; J <= JHI; J++) {
    m_2jp1 = viscosity*grid.m[kk][2*J+1];
    n_2jp1 = viscosity*grid.n[kk][2*J+1];
    for (I = ILO; I <= IHI; I++) {
      LPUU(J,I) = -n_2jp1*(ZE(J+1,I)-ZE(J,I))
                  +m_2jp1*(DI(J,I)-DI(J,I-1));
    }
  }
  bc_lateral(lpuu,TWODIM);

  /*
   * Compute meridional component of the Laplacian.
   */
  for (J = JFIRST; J <= JHI; J++) {
    m_2j = viscosity*grid.m[kk][2*J];
    n_2j = viscosity*grid.n[kk][2*J];
    for (I = ILO; I <= IHI; I++) {
      LPVV(J,I) = m_2j*(ZE(J,I+1)-ZE(J,I))
                 +n_2j*(DI(J,I)-DI(J-1,I));
    }
  }
  bc_lateral(lpvv,TWODIM);

  /*
   * NOTE: The domain's northern and southern boundaries of LPVV(J,I)
   *       are currently defaulted to zero.
   */

  return;
}

/*============== end of laplacian_uv() ============================*/


/*============== scalar_vertical_subgrid() ========================*/

#undef  KLEN
#define KLEN (kend-kstart+1)

void scalar_vertical_subgrid(planetspec  *planet,
                             EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    iq;
  int
    kstart,kend;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_vertical_subgrid";

  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    /*
     * NOTE: The vertical diffusion includes a provision for handling convectively unstable regions.
     */
    scalar_vertical_diffusion(planet,Buff2D);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    /*
     * Adjust convectively unstable regions to be neutrally stable.
     */
    adiabatic_adjustment(planet,Buff2D);
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============== end of scalar_vertical_subgrid() =================*/

/*============== scalar_vertical_diffusion() ======================*/

/*
 * In superadiabatic regions (for non-isentropic coordinates), 
 * convective adjustment is handled by increasing the
 * vertical turbulent diffusion for appropriate fields.
 */

#undef  STAB_MULT
#define STAB_MULT(k,j,i) stab_mult[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void scalar_vertical_diffusion(planetspec  *planet,
                               EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kay,kturb,delta_k_convect,
    iq;
  static int
    nnk,
    initialized = FALSE;
  EPIC_FLOAT
    diffusion_coeff,
    delta_z_convect,
    brunt2,
   *tau_wall;
  const EPIC_FLOAT
    sigma_inv = 3./2.;
  static EPIC_FLOAT
    *stab_factor,
    *nu_convect,
    *zee,
    *aaa,
    *dee,
    *ans,
    *stab_mult;
  unsigned long
    nbytes_2d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="scalar_vertical_diffusion";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  memset(Buff2D[0],0,nbytes_2d);
  tau_wall = Buff2D[0];

  if (!initialized) {
    nnk = KHI-KLO+1;

    /* Allocate memory: */
    stab_factor = fvector(0,2*KHI+1,  dbmsname);
    zee         = fvector(0,KHI+2,    dbmsname);
    aaa         = fvector(0,KHI+2,    dbmsname);
    dee         = fvector(0,KHI+2,    dbmsname);
    ans         = fvector(0,KHI+2,    dbmsname);
    nu_convect  = fvector(0,KHI+1,    dbmsname);
    stab_mult   = fvector(0,Nelem3d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Calculate STAB_MULT
   */
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      stability_factor(J,I,stab_factor);
      for (K = KLO; K <= KHI; K++) {
        STAB_MULT(K,J,I) = stab_factor[2*K];
      }
    }
  }
  /* No need to apply bc_lateral() here. */

  if (grid.coord_type == COORD_ISENTROPIC) {
    ;
  }
  else if (grid.coord_type == COORD_HYBRID) {
    sprintf(Message,"not yet implemented for grid.coord_type == COORD_HYBRID");
    epic_error(dbmsname,Message);
  }
  else if (grid.coord_type == COORD_ISOBARIC) {
    /*
     * Apply vertical diffusion to THETA, 
     * which is carried on the layer interfaces.
     *
     * NOTE: Have not yet included molecular diffusion for THETA.
     *
     * NOTE: This scheme is not yet working well for the hybrid-coordinate model.
     *       Crossing the seam between the hybrid and sigma regions, even just to load for the crank_nicolson() routine,
     *       can lead to a numerical instability at the seam.
     */

    /*
     * We handle convective adjustment via vertical turbulent diffusion of THETA.
     * The parameter kturb refers to K for the highest-altitude (lowest K) THETA that is changed
     * by vertical turbulent diffusion. Currently not implemented for the hybrid-coordinate case.
     *
     * To turn on:  set kturb < KHI
     * To turn off: set kturb = KHI
     */
    kturb = KLO;

    if (kturb < KHI) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          for (K = kturb; K <= KHI; K++) {
            kay      = KHI-K;
            dee[kay] = STAB_MULT(K,J,I)*.5*(DIFFUSION_COEF_THETA(K,J,I)+DIFFUSION_COEF_THETA(K-1,J,I));
          }

          /*
           * Use fixed boundary conditions at bottom and top.
           */
          for (K = kturb-1; K <= KHI; K++) {
            /*
             * Load in positive-z direction (which unfortunately
             * fights against the top-down K numbering of layers).
             */
            kay      = KHI-K;
            aaa[kay] = THETA(K,J,I);
            zee[kay] = Z3(K,J,I);
          }

          crank_nicolson(KHI-kturb,DT,zee,aaa,dee,NULL,ans);

          for (K = kturb; K < KHI; K++) {
            kay          = KHI-K;
            THETA(K,J,I) = ans[kay];
          }
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(var.theta.value,THREEDIM);
      /*
       * Clean up any negative potential temperature introduced by diffusion truncation errors.
       */
      restore_mass(planet,THETA_INDEX,NO_PHASE);
    }
  }
  else {
    sprintf(Message,"not implemented for grid.coord_type == %d",grid.coord_type);
    epic_error(dbmsname,Message);
  }

  /*
   * NOTE: Not applying diffusion to H.
   */

  /*
   * Apply vertical diffusion to mixing ratios, Q,
   * which are carried on the layer interfaces.
   *
   * Loop over all activated species/phase variables.
   */
  if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
    for (iq = 0; iq < grid.nq; iq++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          /*
           * The dee array is on the h-grid grid.
           */
          if (grid.ip[iq] == VAPOR) {
            /*
             * Include molecular mass diffusivity for vapor phase.
             */
            for (K = KLO; K <= KHI; K++) {
              kay      = KHI-K;
              dee[kay] = mass_diffusivity(planet,grid.is[iq],T2(K,J,I),P2(K,J,I))
                        +STAB_MULT(K,J,I)*DIFFUSION_COEF_MASS(K,J,I);
            }
          }
          else {
            for (K = KLO; K <= KHI; K++) {
              kay      = KHI-K;
              dee[kay] = STAB_MULT(K,J,I)*DIFFUSION_COEF_MASS(K,J,I);
            }
          }

          for (K = KLO-1; K <= KHI; K++) {
            kay      = KHI-K;
            aaa[kay] = Q(grid.is[iq],grid.ip[iq],K,J,I);
            zee[kay] = Z3(K,J,I);
          }

          crank_nicolson(KHI-KLO,DT,zee,aaa,dee,NULL,ans);

          for (K = KLO; K < KHI; K++) {
            kay                              = KHI-K;
            Q(grid.is[iq],grid.ip[iq],K,J,I) = ans[kay];
          }
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);
      /*
       * Clean up any negative mass introduced by diffusion truncation errors.
       */
      restore_mass(planet,grid.is[iq],grid.ip[iq]);
    }
    sync_x_to_q(planet);
  }

  /*
   * Apply vertical diffusion to FPARA,
   * which are carried on the layer interfaces.
   * Using DIFFUSION_COEF_MASS.
   */
  if (var.fpara.on) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /*
         * The dee array is on the h-grid grid.
         */
        for (K = KLO; K <= KHI; K++) {
          kay      = KHI-K;
          dee[kay] = STAB_MULT(K,J,I)*DIFFUSION_COEF_MASS(K,J,I);
        }

        for (K = KLO-1; K <= KHI; K++) {
          kay      = KHI-K;
          aaa[kay] = FPARA(K,J,I);
          zee[kay] = Z3(K,J,I);
        }

        crank_nicolson(KHI-KLO,DT,zee,aaa,dee,NULL,ans);

        for (K = KLO; K < KHI; K++) {
          kay          = KHI-K;
          FPARA(K,J,I) = ans[kay];
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(var.fpara.value,THREEDIM);
    /*
     * Clean up any negative fpara introduced by diffusion truncation errors.
     */
    restore_mass(planet,FPARA_INDEX,NO_PHASE);
  }

  /*
   * Update P2, etc.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);

  /*
   * Apply vertical diffusion to NU_TURB, which is carried in the layer.
   * Not modifying for superadiabatic regions.
   *
   * NOTE: The function tau_surface() currently sets TAU_WALL to zero for the gas-giant case.
   */
  if (var.nu_turb.on) {
    tau_surface(planet,NU_TURB_INDEX,tau_wall,Buff2D[1]);

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /*
         * The diffusion coefficient is on the p3-grid.
         */
        for (K = KLO; K < KHI; K++) {
          kay      = KHI-K;
          dee[kay] = sigma_inv*(planet->kinvisc+.5*(NU_TURB(K,J,I)+NU_TURB(K+1,J,I)));
        }
        dee[KHI  ] = dee[KHI-1];
        dee[KLO-1] = dee[KLO  ];

        /*
         * Load in positive-z direction (which fights
         * against the top-down K numbering of layers).
         */
        for (K = KLO; K <= KHI; K++) {
          kay      = KHI-K+1;
          aaa[kay] = NU_TURB(K,J,I);
          zee[kay] = Z2(K,J,I);
        }
        /*
         * Use no-flux boundary condition at top.
         * The bottom b.c. is set below, consistent with TAU_WALL.
         */
        aaa[KHI+1] = 1.e+20;
        zee[KHI+1] = Z3(KLO-1,J,I);

        /*
         * Set bottom boundary condition consistent with TAU_WALL.
         */
        zee[0] = Z3(KHI,J,I);
        aaa[0] = aaa[1]-TAU_WALL(J,I)*(zee[1]-zee[0])/dee[0];

        crank_nicolson(KHI-KLO+1,DT,zee,aaa,dee,NULL,ans);

        for (K = KLO; K <= KHI; K++) {
          kay            = KHI-K+1;
          NU_TURB(K,J,I) = ans[kay];
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(var.nu_turb.value,THREEDIM);
    /*
     * Clean up any negative nu_turb introduced by diffusion truncation errors.
     */
    restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
  }

  return;
}

/*============== end of scalar_vertical_diffusion() ===============*/

/*============== adiabatic_adjustment() ===========================*/

/*
 * In superadiabatic regions, N^2 < 0, relax the potential temperature towards N^2 = 0,
 * and relax any species towards the well-mixed state.
 *
 * The relaxtion rate is assumed to be proportional to 1/|N|, based on the local negative N^2,
 * but diluted as a subgrid process.
 *
 * NOTE: This algorithm is not applied to the isentropic coordinate or the isentropic part
 *       of the hybrid coordinate.
 */

void adiabatic_adjustment(planetspec  *planet,
                          EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    klo,iq,
    modification;
  EPIC_FLOAT
    avg,sum,
    fpara,mu,theta_ortho,theta_para,
    alpha,max_alpha,dilute;
  EPIC_FLOAT
    N2[2],d_enthalpy[3],
    orig_theta1,orig_theta3,
    avg_enthalpy,enthalpy1,enthalpy3,
    cp1,cp3,dp1,dp3,fpara1,fpara3,Temp1,Temp3;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="adiabatic_adjustment";
#ifdef EPIC_MPI
  int
    mpi_itmp;
#endif

  if (grid.coord_type == COORD_ISENTROPIC) {
    ;
  }
  else if (grid.coord_type == COORD_ISOBARIC || grid.coord_type == COORD_HYBRID) {
    if (grid.coord_type == COORD_ISOBARIC) {
      klo = KLO;
    }
    else {
      klo =  grid.k_sigma;
    }

    modification = FALSE;

    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /*
         * Move down the column starting at the first non-isentropic layer.
         */
        for (K = klo; K <= KHI; K++) {
          /*
           * get_brunt2() < 0 is a better indication of (hydro)static instability than dtheta/dz < 0. 
           * See the comments in get_brunt2() in epic/src/core/epic_funcs_diag.c.
           *
           * NOTE: No adjustment is made when THETA(K,J,I) <= THETA(K-1,J,I), since in practice this case
           *       has the potential for an unstable positive feedback.
           */
          N2[0] = get_brunt2(planet,2*K,J,I);
          if (N2[0] < 0 && THETA(K,J,I) > THETA(K-1,J,I)) {
            modification = TRUE;

            /*
             * Record original values of THETA to facilitate cases where no adjustment is made.
             */
            orig_theta1 = THETA(K-1,J,I);
            orig_theta3 = THETA(K,  J,I);

            /*
             * Adjust the size of alpha by modifying the parameters dilute and max_alpha.
             * Set the fraction of adjustment: alpha = dilute * |N| dt.
             */
            dilute    = .05;
            max_alpha = .002;

            alpha = MIN(dilute*sqrt(-N2[0])*(EPIC_FLOAT)grid.dt,max_alpha);

            if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
              /*
               * The mixing ratios, Q, are carried on the layer interfaces.
               * Relax towards the well-mixed state.
               *
               * NOTE: alpha*target+(1-alpha)*initial is practically the same as 
               *       (1-exp(-t/tau))*target+exp(-t/tau)*initial, since t/tau = dt/(1/|N|) = alpha << 1.
               *
               * Loop over all activated species/phase variables.
               */
              for (iq = 0; iq < grid.nq; iq++) {
                avg = .5*(Q(grid.is[iq],grid.ip[iq],K-1,J,I)+Q(grid.is[iq],grid.ip[iq],K,J,I));
                Q(grid.is[iq],grid.ip[iq],K-1,J,I) = alpha*avg+(1.-alpha)*Q(grid.is[iq],grid.ip[iq],K-1,J,I);
                Q(grid.is[iq],grid.ip[iq],K,  J,I) = alpha*avg+(1.-alpha)*Q(grid.is[iq],grid.ip[iq],K,  J,I);
              }
              /*
               * Sync X to Q, needed for avg_molar_mass();
               */
              sum = planet->rgas/R_GAS;
              for (iq = 0; iq < grid.nq; iq++) {
                sum += Q(grid.is[iq],grid.ip[iq],K,J,I)/var.species[grid.is[iq]].molar_mass;
              }
              for (iq = 0; iq < grid.nq; iq++) {
                X(grid.is[iq],grid.ip[iq],K,J,I)   = Q(grid.is[iq],grid.ip[iq],K,J,I)/(var.species[grid.is[iq]].molar_mass*sum);
                X(grid.is[iq],grid.ip[iq],K-1,J,I) = X(grid.is[iq],grid.ip[iq],K,J,I);
              }
            }

            if (var.fpara.on) {
              /*
               * Relax fpara, which is carried on the layer interfaces.
               */
              avg = .5*(FPARA(K-1,J,I)+FPARA(K,J,I));
              FPARA(K-1,J,I) = alpha*avg+(1.-alpha)*FPARA(K-1,J,I);
              FPARA(K,  J,I) = alpha*avg+(1.-alpha)*FPARA(K,  J,I);
            }

            /*
             * N^2 = 0 does not in general correspond to THETA = constant. The algorithm in get_brunt2() is quite
             * general; see the comments there.  The point is that molar mass gradients affect the neutral profile,
             * as does the ortho-para hydrogen mix for gas giants and ice giants.
             *
             * Number the top, middle and bottom of the layer by:   --- 1 ---
             *                                                          2 
             *                                                      --- 3 ---
             *
             * We want the adjustment to preserve the total enthalpy. Define
             *    avg_enthalpy = .5*(cp3*T3*dp3+cp1*T1*dp1)
             *    d_enthalpy   = .5*(cp3*T3*dp3-cp1*T1*dp1)
             *
             * To avoid weirdness on the bottom of the bottom layer, take the pressure differences, which
             * facilitate mass averaging, to be
             *   dp1 = P2(K)-P3(K-1)
             *   dp3 = P3(K)-P2(K  )  
             * We assume that pressure does not change.
             *
             * Model N2^2 = N2[] as a linear function of d_enthalpy[]. 
             * We start off with the d_enthalpy[0] that produces N2[0] < 0.
             * We guess a reasonable d_enthalpy[1] and find N2[1].  Then we assume a straight line to estimate 
             * the d_enthalpy[2] that produces N2[2] = 0, and update all the affected variables.
             */

            if (var.fpara.on) {
              fpara1 = FPARA(K-1,J,I);
              fpara3 = FPARA(K,  J,I);
            }
            else {
              fpara1 = return_fpe(T3(K-1,J,I));
              fpara3 = return_fpe(T3(K,  J,I));
            }

            cp1 = return_cp(planet,fpara1,P3(K-1,J,I),T3(K-1,J,I));
            cp3 = return_cp(planet,fpara3,P3(K,  J,I),T3(K,  J,I));
 
            dp1 = P2(K,J,I)-P3(K-1,J,I);
            dp3 = P3(K,J,I)-P2(K,  J,I);

            enthalpy1 = dp1*cp1*T3(K-1,J,I);
            enthalpy3 = dp3*cp3*T3(K,  J,I);

            /* avg_enthalpy is constant during the adiabatic adjustment */
            avg_enthalpy  = .5*(enthalpy3+enthalpy1);
            d_enthalpy[0] = .5*(enthalpy3-enthalpy1);

            /* Use theta = constant to get a good d_enthalpy[1]. */
            avg   = .5*(THETA(K-1,J,I)+THETA(K,J,I));
            Temp1 = return_temp(planet,fpara1,P3(K-1,J,I),avg);
            Temp3 = return_temp(planet,fpara3,P3(K,  J,I),avg);

            cp1 = return_cp(planet,fpara1,P3(K-1,J,I),Temp1);
            cp3 = return_cp(planet,fpara3,P3(K,  J,I),Temp3);

            enthalpy1 = dp1*cp1*Temp1;
            enthalpy3 = dp3*cp3*Temp3;

            d_enthalpy[1] = 0.5*(enthalpy3-enthalpy1);

            enthalpy1 = avg_enthalpy-d_enthalpy[1];
            enthalpy3 = avg_enthalpy+d_enthalpy[1];

            /*
             * Update T3.
             */
            T3(K-1,J,I) = enthalpy1/(dp1*cp1);
            T3(K,  J,I) = enthalpy3/(dp3*cp3);

            /*
             * Update the relevant variables for get_brunt2().
             */
            THETA(K-1,J,I) = return_theta(planet,fpara1,P3(K-1,J,I),T3(K-1,J,I),&theta_ortho,&theta_para);
            mu             = avg_molar_mass(planet,2*K-1,J,I);
            RHO3(K-1,J,I)  = return_density(planet,fpara1,P3(K-1,J,I),T3(K-1,J,I),mu,PASSING_T);

            THETA(K,J,I) = return_theta(planet,fpara3,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
            mu           = avg_molar_mass(planet,2*K+1,J,I);
            RHO3(K,J,I)  = return_density(planet,fpara3,P3(K,J,I),T3(K,J,I),mu,PASSING_T);

            if (var.fpara.on) {
              fpara = .5*(FPARA(K-1,J,I)+FPARA(K,J,I));
            }
            else {
              fpara = .5*(return_fpe(T3(K-1,J,I))+return_fpe(T3(K,J,I)));
            }
            THETA2(K,J,I) = .5*(THETA(K-1,J,I)+THETA(K,J,I));
            T2(    K,J,I) = return_temp(   planet,fpara,P2(K,J,I),THETA2(K,J,I));
            mu            = avg_molar_mass(planet,2*K,J,I);
            RHO2(K,J,I)   = return_density(planet,fpara,P2(K,J,I),T2(K,J,I),mu,PASSING_T);

            /*
             * Calculate N^2 for d_enthaphy[1].
             */
            N2[1] = get_brunt2(planet,2*K,J,I);

            /*
             * Assume the function N2(d_enthalpy) is linear. We have two data points:
             *   N2[0] at d_enthalpy[0]
             *   N2[1] at d_enthalpy[1]
             * plus N2[2] = 0. Use these to estimate d_enthalpy[2].
             */
            if (fcmp(N2[1],N2[0]) != 0) {
              d_enthalpy[2] = d_enthalpy[0]-N2[0]*(d_enthalpy[1]-d_enthalpy[0])/(N2[1]-N2[0]);
              /*
               * If the extrapolation step is too large, use d_enthalpy[0] or [1], whichever
               * yields |N2| closest to zero (to avoid positive feedback).
               */
              if (fabs(d_enthalpy[2]-d_enthalpy[1]) > fabs(d_enthalpy[1]-d_enthalpy[0])) {
                if (fabs(N2[1]) >= fabs(N2[0])) {
                  /* Make no adjustment to THETA */
                  THETA(K-1,J,I) = orig_theta1;
                  THETA(K,  J,I) = orig_theta3;
                  continue;
                }
                else {
                  d_enthalpy[2] = d_enthalpy[1];
                }
              }
            }
            else {
              /* Make no adjustment to THETA */
              THETA(K-1,J,I) = orig_theta1;
              THETA(K,  J,I) = orig_theta3;
              continue;
            }

            enthalpy1 = avg_enthalpy-d_enthalpy[2];
            enthalpy3 = avg_enthalpy+d_enthalpy[2];

            /*
             * Update T3.
             */
            T3(K-1,J,I) = enthalpy1/(dp1*cp1);
            T3(K,  J,I) = enthalpy3/(dp3*cp3);

            /*
             * Full THETA update.
             */
            THETA(K-1,J,I) = return_theta(planet,fpara1,P3(K-1,J,I),T3(K-1,J,I),&theta_ortho,&theta_para);
            THETA(K,  J,I) = return_theta(planet,fpara3,P3(K,  J,I),T3(K,  J,I),&theta_ortho,&theta_para);

            /*
             * Actual THETA update.
             * All the affected diagnostic variables are synchronized below.
             */
            if ((THETA(K-1,J,I)-orig_theta1) > 0. && (THETA(K,J,I)-orig_theta3) < 0.) {
              /*
               * The full update is a twist in the right direction, so relax towards it.
               */
              THETA(K-1,J,I) = alpha*THETA(K-1,J,I)+(1.-alpha)*orig_theta1;
              THETA(K,  J,I) = alpha*THETA(K,  J,I)+(1.-alpha)*orig_theta3;
            }
            else {
              /*
               * The full update is a drift without a twist, so discard and make no adjustment.
               */
              THETA(K-1,J,I) = orig_theta1;
              THETA(K,  J,I) = orig_theta3;
            }

          } /* If N^2 < 0 */
        } /* K loop */
      } /* I loop */
    } /* J loop */

#ifdef EPIC_MPI
  /*
   * Make status of modification global.
   */
  mpi_itmp = modification;
  MPI_Allreduce(&mpi_itmp,&modification,1,MPI_INT,MPI_SUM,para.comm);
#endif

    if (modification) {
      /*
       * Apply lateral boundary conditions.
       */
      if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
        for (iq = 0; iq < grid.nq; iq++) {
          bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);
        }
        sync_x_to_q(planet);
      }

      if (var.fpara.on) {
        bc_lateral(var.fpara.value,THREEDIM);
      }

      bc_lateral(var.theta.value, THREEDIM);

      /*
       * Update diagnoatic variables THETA2, HDRY2, HDRY3, etc.
       */
      set_p2_etc(planet,UPDATE_THETA,Buff2D);

      /*
       * Update diagnostic variables T2, T3, RHO2, RHO3, etc.
       */
      store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,CALC_PHI3NK);
    }
  }
  else {
    sprintf(Message,"not yet implemented for grid.coord_type == %d",grid.coord_type);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============== end of adiabatic_adjustment() ====================*/

/*============== uv_horizontal_subgrid() ==========================*/

void uv_horizontal_subgrid(planetspec  *planet,
                           EPIC_FLOAT **Buff2D)
{
  register int
    K;
  int
    kstart,kend;
  EPIC_FLOAT
     tmp,
    *pt_dudt,
    *pt_dvdt;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_horizontal_subgrid";

  if (strcmp(grid.turbulence_scheme,"on") == 0) {
    uv_horizontal_diffusion(planet,Buff2D);
  }
  else if (strcmp(grid.turbulence_scheme,"on_vertical_only") == 0 ||
           strcmp(grid.turbulence_scheme,"off")              == 0)  {
    ;
  }
  else {
    sprintf(Message,"Unrecognzied grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  if (grid.nudiv_nondim > 0.) {
    /*
     * Apply divergence damping.
     *
     * To avoid the top sponge, set
     *   kstart = IMAX(KLO,grid.k_sponge+1);
     * otherwise, set
     *   kstart = KLO;
     */
    kstart = KLO;
    kend   = KHI;

    for (K = kstart; K <= kend; K++) {
      divergence_damping(planet,K,grid.nudiv_nondim,Buff2D);
    }
  }

  /*
   * NOTE: Horizontal hyperviscosity may be applied directly to U,V elsewhere, but is not
   *       added into their tendencies here.
   */ 

  return;
}

/*============== end of uv_horizontal_subgrid() ===================*/

/*=================== uv_horizontal_diffusion() ===================*/

/*
 * We follow Tannehill, Anderson, and Pletcher (1997), p. 269
 * for the components of the viscosity in general coordinates.
 * Here,
 *   x1 = phi,         east longitude
 *   x2 = lambda,      planetographic latitude
 *   x3 = z,           geopotential height
 *
 *   h1 = r(lambda),
 *   h2 = R(lambda),
 *   h3 = 1,           plane-parallel approximation
 *
 *   u1 = u,
 *   u2 = v,
 *   u3 = 0,           shallow-atmosphere approximation
 *
 * We set all d/dz = 0 in this subroutine.  The vertical-gradient
 * terms are handled with an implicit timestep in uv_vertical_diffusion()
 * for numerical stability.
 *
 * NOTE: Currently not lagging the viscosity and density factors themselves in the 
 *       leapfrog case, which might possibly cause numerical stability problems, 
 *       but hasn't presented any obvious problems. We are lagging the (u,v) values 
 *       themselves via grid.it_uv_dis.
 */

#undef  COEFFD
#define COEFFD(j,i) coeffd[i+(j)*Iadim-Shift2d]

void uv_horizontal_diffusion(planetspec  *planet,
			     EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kk,jj;
  register EPIC_FLOAT
    rho,rho_inv,nu,
    rln,rln_inv,rlt_inv,
    e11,e12,e22,
    taper;
  EPIC_FLOAT
    *tau11,*tau12,*tau22,
    *uu,*vv,
    *coeffd,
     coeff;
  register EPIC_FLOAT
    m_2j,m_2jp1,
    n_2j,n_2jp1,n_2jp2;
  static int
   initialized = FALSE;
  static EPIC_FLOAT
   *m0;
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
    dbmsname[]="uv_horizontal_diffusion";

  if (!initialized) {
    /* Allocate memory */
    m0 = fvector(0,grid.nk,dbmsname);
    for (K = KLO; K <= KHI; K++) {
      rln   = grid.re[K]/sqrt(1.+SQR(grid.rp[K]/grid.re[K]*tan(LAT0*DEG)));
      m0[K] = 1./(rln*grid.dln*DEG);
    }

    initialized = TRUE;
  }

  tau11  = Buff2D[0];
  tau22  = Buff2D[1];
  tau12  = Buff2D[2];
  coeffd = Buff2D[3];
  uu     = Buff2D[4];
  vv     = Buff2D[5];

  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        COEFFD(J,I) = planet->dynvisc+RHO2(K,J,I)*DIFFUSION_COEF_UV(K,J,I);
      }
    }
    /* No need to apply bc_lateral() here. */

    /*
     * Copy U and V into UU and VV.
     */
    memcpy(uu,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));
    memcpy(vv,var.v.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,Nelem2d*sizeof(EPIC_FLOAT));

    /*
     * TAU11 and TAU22 adapt naturally to the h-grid.
     * Fortunately we do not have TAU33, because it is awkward in the Arakawa C-grid system.
     * We have e33 = 0., so do not include it.
     */
    for (J = JLO; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./grid.rln[kk][jj+1];
      m_2jp1  = grid.m[kk][jj+1];
      n_2jp1  = grid.n[kk][jj+1];

      for (I = ILO; I <= IHI; I++) {
        e11 = (UU(J,I+1)-UU(J,I))*m_2jp1
	      +.5*(VV(J,I)+VV(J+1,I))*rln_inv*(grid.rln[kk][jj+2]-grid.rln[kk][jj])*n_2jp1;
        e22 = (VV(J+1,I)-VV(J,I))*n_2jp1;

	/*
	 * TAU11 and TAU22 are calculated using e's.
	 */
        coeff      = (2./3.)*COEFFD(J,I);
	TAU11(J,I) = coeff*(2.*e11-e22);
	TAU22(J,I) = coeff*(2.*e22-e11);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(tau11,TWODIM);
    bc_lateral(tau22,TWODIM);

    /*
     * TAU12 falls naturally on the pv2-grid.
     * Call vorticity() to handle the northern and southern poles or edges.
     *
     * NOTE: the call to vorticity() must be by every node, even though we
     *       are only using the results at the model's northen and southern extremes.
     *       If the special edge cases were to be broken out into a separate subroutine,
     *       this would not be necessary. (This is not a pressing concern.)
     */
    vorticity(ON_SIGMATHETA,RELATIVE,kk,uu,vv,NULL,tau12);

    if (JHI == grid.nj) {
      /*
       * Northern pole or edge.
       */
      coeff = 0.;
      J     = grid.nj;
      for (I = ILO; I <= IHI; I++) {
        coeff += COEFFD(J,I);
      }

#if defined(EPIC_MPI)
      mpi_tmp = coeff;
      MPI_Allreduce(&mpi_tmp,&coeff,1,float_type,MPI_SUM,para.comm_JLO);
#endif

      coeff /= grid.ni;

      J = grid.nj+1;
      for (I = ILO; I <= IHI; I++) {
        /*
         * The minus sign comes from the fact that we need du/dy, not -du/dy.
         */
        TAU12(J,I) *= -coeff;
      }
    }
    if (JLO == grid.jlo) {
      /*
       * Southern pole or edge.
       */
      coeff = 0.;
      J     = 0;
      for (I = ILO; I <= IHI; I++) {
        coeff += COEFFD(J,I);
      }

#if defined(EPIC_MPI)
      mpi_tmp = coeff;
      MPI_Allreduce(&mpi_tmp,&coeff,1,float_type,MPI_SUM,para.comm_JLO);
#endif

      coeff /= grid.ni;

      for (I = ILO; I <= IHI; I++) {
        TAU12(J,I) *= -coeff;
      }
    }
    /*
     * Fill in interior of TAU12.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj   = 2*J;
      m_2j = grid.m[kk][jj];
      n_2j = grid.n[kk][jj];

      for (I = ILO; I <= IHI; I++) {
        e12 = (VV(J,I)-VV(J,I-1))*m_2j 
	     +grid.rln[kk][jj]*(UU(J,I)/grid.rln[kk][jj+1]-UU(J-1,I)/grid.rln[kk][jj-1])*n_2j;
	/*
	 * TAU12 is calculated using e12.
         * Average onto pv-grid.
	 */
        coeff      = .25*(COEFFD(J,I)+COEFFD(J,I-1)+COEFFD(J-1,I)+COEFFD(J-1,I-1));
	TAU12(J,I) = coeff*e12;
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(tau12,TWODIM);
    
    /*
     * Horizontal-gradient diffusion of U.
     */
    for (J = JLO; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./grid.rln[kk][jj+1];
      m_2jp1  = grid.m[kk][jj+1];
      n_2jp1  = grid.n[kk][jj+1];
      /*
       * Old taper to help prevent numerical instability:
       *
       * taper = MIN(1.,pow(m0[K]/m_2jp1,2.));
       */
      taper = 1.;

      for (I = ILO; I <= IHI; I++) {
        rho_inv                      = 2./(RHO2(K,J,I)+RHO2(K,J,I-1));
    	DUDT(grid.it_uv_tend,K,J,I) += taper*rho_inv*(m_2jp1*(TAU11(J,I)-TAU11(J,I-1))
                                                     +n_2jp1*rln_inv*(grid.rln[kk][jj+2]*TAU12(J+1,I)-grid.rln[kk][jj]*TAU12(J,I)
                                                           +.5*(TAU12(J+1,I)+TAU12(J,I))*(grid.rln[kk][jj+2]-grid.rln[kk][jj])));
      }
    }

    /*
     * Horizontal-gradient diffusion of V.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj      = 2*J;
      rln_inv = 1./grid.rln[kk][jj];
      m_2j    = grid.m[kk][jj];
      n_2j    = grid.n[kk][jj];
      /*
       * Old taper to help prevent numerical instability:
       *
       * taper = MIN(1.,pow(m0[K]/m_2j,2.));
       */
      taper = 1.;

      for (I = ILO; I <= IHI; I++) {
        rho_inv                      = 2./(RHO2(K,J,I)+RHO2(K,J-1,I));
    	DVDT(grid.it_uv_tend,K,J,I) += taper*rho_inv*(m_2j*(TAU12(J,I)-TAU12(J,I-1))
                                                     -n_2j*rln_inv*(grid.rln[kk][jj+1]*TAU22(J,I)-grid.rln[kk][jj-1]*TAU22(J-1,I)
                                                           -.5*(TAU11(J,I)+TAU11(J-1,I))*(grid.rln[kk][jj+1]-grid.rln[kk][jj-1])));
      }
    }
  }
 
  return;
}

/*================== end of uv_horizontal_diffusion() =============*/

/*================== divergence_damping() =========================*/

/*
 * Add artificial damping of horizontal divergence to control
 * numerical instabilities associated with gravity waves.
 * See Skamarock and Klemp (1992, Mon. Wea. Rev. 120, 2109-2127).
 */

void divergence_damping(planetspec  *planet,
                        int          K,
                        EPIC_FLOAT   nudiv_nondim,
                        EPIC_FLOAT **Buff2D)
{
  register int
    J,I,kay,
    kk = 2*K;
  register EPIC_FLOAT
    nudiv,
    coef,
    rln;
  static int
    initialized = FALSE;
  static double
    max_nu_horizontal[MAX_NU_ORDER+1];
  EPIC_FLOAT
   *div;
  static EPIC_FLOAT
   *m0;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="divergence_damping";

  if (!initialized) {
    set_max_nu(max_nu_horizontal);

    /* Allocate memory */
    m0 = fvector(0,grid.nk,dbmsname);
    for (kay = KLO; kay <= KHI; kay++) {
      rln   = grid.re[kay]/sqrt(1.+SQR(grid.rp[kay]/grid.re[kay]*tan(LAT0*DEG)));
      m0[kay] = 1./(rln*grid.dln*DEG);
    }

    initialized = TRUE;
  }

  if (nudiv_nondim < 0.) {
    sprintf(Message,"nudiv_nondim=%e < 0.",nudiv_nondim);
    epic_error(dbmsname,Message);
  }
  else if (nudiv_nondim == 0.) {
    return;
  }

  nudiv = nudiv_nondim*max_nu_horizontal[2];

  /*
   * Use grid.it_uv_dis for numerical stability (e.g. leapfrog timestep).
   */
  div = Buff2D[0];

  divergence(kk,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,
                var.v.value+(K-Kshift)*Nelem2d+grid.it_uv_dis*Nelem3d,div);

  /*
   * High-latitude, low-pass filter to prevent numerical instability.
   * This is done under the gradient operator to keep it irrotational.
   *
   * NOTE: We do not use a taper on the viscosity coefficient.  It does not
   *       work as well as zonal_filter() in practice.  Also, we cannot
   *       arrange to have the taper outside of the laplacian, since only
   *       the first derivative is actually applied.
   */
  zonal_filter(DIV_UV2_INDEX,div);

  for (J = JFIRST; J <= JHI; J++) {
    coef  = nudiv*grid.n[kk][2*J];
    for (I = ILO; I <= IHI; I++) {
      DVDT(grid.it_uv_tend,K,J,I) += coef*(DIV(J,I)-DIV(J-1,I));
    }
  }
  /* No need to call bc_lateral() here. */

  for (J = JLO; J <= JHI; J++) {
    coef  = nudiv*grid.m[kk][2*J+1];
    for (I = ILO; I <= IHI; I++) {
      DUDT(grid.it_uv_tend,K,J,I) += coef*(DIV(J,I)-DIV(J,I-1));
    }
  }
  /* No need to call bc_lateral() here. */

  return;
}

/*================== end of divergence_damping() ==================*/

/*============== uv_vertical_subgrid() ============================*/

void uv_vertical_subgrid(planetspec  *planet,
                         EPIC_FLOAT **Buff2D)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_vertical_subgrid";

  if (strcmp(grid.turbulence_scheme,"on")               == 0 ||
      strcmp(grid.turbulence_scheme,"on_vertical_only") == 0)  {
    uv_vertical_diffusion(planet,Buff2D);
  }
  else if (strcmp(grid.turbulence_scheme,"off") == 0) {
    /*
     * NOTE: convective adjustment is not currently applied to the horizonal momentum.
     */
    ;
  }
  else {
    sprintf(Message,"Unrecognized grid.turbulence_scheme=%s",grid.turbulence_scheme);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============== end of uv_vertical_subgrid() =====================*/

/*============== uv_vertical_diffusion() ==========================*/

/*
 * The vertical-gradient diffusion terms for (u,v) are handled here using 
 * an implicit timestep to avoid numerical instabilities arising from
 * relatively small dz values.  The horizontal-gradient terms are
 * handled in uv_horizontal_diffusion().
 *
 * The shallow-atmosphere approximation has been applied, which eliminates
 * terms involving the vertical velocity.
 *
 * NOTE: The timeplane for U and V is grid.it_uv rather than 
 *       grid.it_uv_dis, because the implicit timestep does not need
 *       to be lagged.
 *
 * NOTE: Killworth (1989, On the Parameterization of Deep Convection in Ocean Models)
 *       concludes that convective-adjustment mixing is not necessary for U, V.
 */

void uv_vertical_diffusion(planetspec  *planet,
                           EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kay;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    *tau_wall;
  static EPIC_FLOAT
    *zee,
    *aaa,
    *mu,
    *rho,
    *ans;
  unsigned long
    nbytes_2d;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_vertical_diffusion";

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);
  memset(Buff2D[0],0,nbytes_2d);
  tau_wall = Buff2D[0];

  if (!initialized) {
    /* Allocate memory: */
    zee          = fvector(0,KHI+1,  dbmsname);
    aaa          = fvector(0,KHI+1,  dbmsname);
    mu           = fvector(0,KHI+1,  dbmsname);
    rho          = fvector(0,KHI+1,  dbmsname);
    ans          = fvector(0,KHI+1,  dbmsname);

    initialized = TRUE;
  }

  /*
   * Apply vertical diffusion to U.
   */
  tau_surface(planet,U_INDEX,tau_wall,Buff2D[1]);
  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      /*
       * Use no-flux boundary condition at top.
       */
      K   = 0;
      kay = KHI-K+1;
      aaa[kay] = 1.e+20;
      zee[kay] = .5*(Z3(K,J,I)+Z3(K,J,I-1));
      for (K = KLO; K <= KHI; K++) {
        kay      = KHI-K+1;
        aaa[kay] = U(grid.it_uv,K,J,I);
        zee[kay] = .5*(  Z2(K,J,I)+  Z2(K,J,I-1));
        rho[kay] = .5*(RHO2(K,J,I)+RHO2(K,J,I-1));
      }

      /*
       * mu[kay] is staggered
       */
      for (K = KLO; K < KHI; K++) {
        kay     = KHI-K;
        mu[kay] = planet->dynvisc
                 +.5*(RHO3(K,J,I  )*.5*(DIFFUSION_COEF_UV(K,J,I  )+DIFFUSION_COEF_UV(K+1,J,I  ))
                     +RHO3(K,J,I-1)*.5*(DIFFUSION_COEF_UV(K,J,I-1)+DIFFUSION_COEF_UV(K+1,J,I-1)));
      }
      K       = KLO-1;
      kay     = KHI-K;
      mu[kay] = mu[kay-1];
      K       = KHI;
      kay     = KHI-K;
      mu[kay] = mu[kay+1];
      /*
       * Set bottom boundary condition consistent with TAU_WALL.
       * NOTE: TAU_WALL is zero for giant planets.
       */
      K   = KHI+1;
      kay = KHI-K+1;
      zee[kay] = .5*(Z3(K-1,J,I)+Z3(K-1,J,I-1));
      aaa[kay] = aaa[kay+1]-TAU_WALL(J,I)*(zee[kay+1]-zee[kay])/mu[kay];

      crank_nicolson(KHI,DT,zee,aaa,mu,rho,ans);

      for (K = KLO; K <= KHI; K++) {
        kay                 = KHI-K+1;
        U(grid.it_uv,K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  
  /*
   * Apply vertical diffusion to V.
   */
  tau_surface(planet,V_INDEX,tau_wall,Buff2D[1]);
  for (J = JFIRST; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      /*
       * Use no-flux boundary condition at top.
       */
      K   = 0;
      kay = KHI-K+1;
      aaa[kay] = 1.e+20;
      zee[kay] = .5*(Z3(K,J,I)+Z3(K,J-1,I));
      for (K = KLO; K <= KHI; K++) {
        kay      = KHI-K+1;
        aaa[kay] = V(grid.it_uv,K,J,I);
        zee[kay] = .5*(  Z2(K,J,I)+  Z2(K,J-1,I));
        rho[kay] = .5*(RHO2(K,J,I)+RHO2(K,J-1,I));
      }

      /*
       * mu[kay] is staggered
       */
      for (K = KLO; K < KHI; K++) {
        kay     = KHI-K;
        mu[kay] = planet->dynvisc
                 +.5*(RHO3(K,J,  I)*.5*(DIFFUSION_COEF_UV(K,J,  I)+DIFFUSION_COEF_UV(K+1,J,  I))
                     +RHO3(K,J-1,I)*.5*(DIFFUSION_COEF_UV(K,J-1,I)+DIFFUSION_COEF_UV(K+1,J-1,I)));
      }
      K       = KLO-1;
      kay     = KHI-K;
      mu[kay] = mu[kay-1];
      K       = KHI;
      kay     = KHI-K;
      mu[kay] = mu[kay+1];
      /*
       * Set bottom boundary condition consistent with TAU_WALL.
       * NOTE: TAU_WALL is zero for giant planets.
       */
      K   = KHI+1;
      kay = KHI-K+1;
      zee[kay] = .5*(Z3(K-1,J,I)+Z3(K-1,J-1,I));
      aaa[kay] = aaa[kay+1]-TAU_WALL(J,I)*(zee[kay+1]-zee[kay])/mu[kay];

      crank_nicolson(KHI,DT,zee,aaa,mu,rho,ans);

      for (K = KLO; K <= KHI; K++) {
        kay                 = KHI-K+1;
        V(grid.it_uv,K,J,I) = ans[kay];
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

  return;
}

/*============== end of uv_vertical_diffusion() ====================*/

/*======================= make_arrays_subgrid() ====================*/

void make_arrays_subgrid(void)
{
  register int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="make_arrays_subgrid";

  /*
   * Allocate memory.
   */
  d_wall = fvector(0,Nelem3d-1,dbmsname);

  return;
}

/*======================= end of make_arrays_subgrid() =============*/

/*======================= free_arrays_subgrid() ====================*/

void free_arrays_subgrid(void)
{
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_arrays_subgrid";

  /*
   * Free allocated memory.
   */
  free_fvector(d_wall,0,Nelem3d-1,dbmsname);

  return;
}

/*======================= end of free_arrays_subgrid() ============*/

/*======================= init_subgrid() ==========================*/

void init_subgrid(planetspec *planet)
{
  register int
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_subgrid";

  if (!var.nu_turb.on) {
    sprintf(Message,"var.nu_turb.on is off");
    epic_error(dbmsname,Message);
  }

  /*
   * Initialize arrays.
   */  
  for (K = KLOPAD; K <= KHIPAD; K++) {
    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        NU_TURB(K,J,I) = planet->kinvisc*20.;
      }
    }
    /* No need to apply bc_lateral() here. */
  }

  return;
}

/*======================= end of init_subgrid() ====================*/

/*======================= set_diffusion_coef() =====================*/

/*
  * Calculate turbulent diffusion coefficients.
  * h-grid:  DIFFUSION_COEF_MASS, DIFFUSION_COEF_UV
  * p3-grid: DIFFUSION_COEF_THETA
  *
  * Molecular diffusion should be accounted for elsewhere.
  */

void set_diffusion_coef(planetspec *planet)
{
  register int
    K,J,I,
    kk;
  register EPIC_FLOAT
    chi3,fv1,nu_turb,turb,
    u_tan,kin,
    tmp;
  EPIC_FLOAT
   *u2d,
   *v2d;
  const EPIC_FLOAT
    cv1   = 7.1,
    cv1_3 = cv1*cv1*cv1;
  static double
    max_nu_horizontal[MAX_NU_ORDER+1];
  static int
    initialized=FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="set_diffusion_coef";

  if (!initialized) {
    set_max_nu(max_nu_horizontal);

    initialized = TRUE;
  }

  if (strcmp(planet->type,"terrestrial") == 0) {
    /*
     * DIFFUSION_COEF_UV and DIFFUSION_COEF_MASS, both on the h-grid.
     * For terrestrial planets, the bottom layer, K = KHI,
     * is treated separately below.
     */
    for (K = KLO; K < KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          nu_turb  = NU_TURB(K,J,I);
          chi3     = nu_turb/planet->kinvisc;
          chi3    *= chi3*chi3;
          fv1      = chi3/(chi3+cv1_3);
          turb     = fv1*nu_turb;
          tmp      = turb;
          tmp      = LIMIT_RANGE(0.,tmp,max_nu_horizontal[2]);
          DIFFUSION_COEF_UV(K,J,I) = tmp;

          /*
           * Mass diffusivity.
           * Currently using the same value for turbulent mass diffusivity
           * as for temperature.
           */
          DIFFUSION_COEF_MASS(K,J,I) = tmp;
        }
      }
    }

    /*
     * DIFFUSION_COEF_UV and DIFFUSION_COEF_MASS, K = KHI
     */
    dwall_SA(planet,d_wall);

    K   = KHI;
    kk  = 2*K;
    u2d = var.u.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    v2d = var.v.value+(K-Kshift)*Nelem2d+(grid.it_uv_dis)*Nelem3d;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        nu_turb = NU_TURB(K,J,I);
        kin     = get_kin(planet,u2d,v2d,kk,J,I);
        u_tan   = sqrt(2.*kin);
        turb    = law_of_the_wall(planet,K,J,I,NU_TURB_INDEX,nu_turb,u_tan);
        tmp     = turb+planet->kinvisc;
        tmp     = LIMIT_RANGE(0.,tmp,max_nu_horizontal[2]);
        DIFFUSION_COEF_UV(KHI,J,I)   = tmp;
        DIFFUSION_COEF_MASS(KHI,J,I) = tmp;
      }
    }

    /*
     * DIFFUSION_COEF_THETA is on the p3-grid.
     */
    for (K = KLO; K < KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          /*
           * Take the turbulent Prandtl number to be unity,
           * following Collins et al. (2004, NCAR/TN-464+STR).
           * THETA is on the p3-grid.
           */
          nu_turb  = .5*(NU_TURB(K,J,I)+NU_TURB(K+1,J,I));
          chi3     = nu_turb/planet->kinvisc;
          chi3    *= chi3*chi3;
          fv1      = chi3/(chi3+cv1_3);
          turb     = fv1*nu_turb;
          tmp      = turb;
          tmp      = LIMIT_RANGE(0.,tmp,max_nu_horizontal[2]);
          DIFFUSION_COEF_THETA(K,J,I) = tmp;
        }
      }
    }
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /* THETA is on the p3-grid. */
        DIFFUSION_COEF_THETA(0,  J,I) = 0.;
        DIFFUSION_COEF_THETA(KHI,J,I) = 0.;
      }
    }
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * DIFFUSION_COEF_UV and DIFFUSION_COEF_MASS, both on the h-grid.
     * For gas-giant planets, the bottom layer is treated the same as
     * the other layers.
     */
    for (K = KLO; K <= KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          nu_turb  = NU_TURB(K,J,I);
          chi3     = nu_turb/planet->kinvisc;
          chi3    *= chi3*chi3;
          fv1      = chi3/(chi3+cv1_3);
          turb     = fv1*nu_turb;
          tmp      = turb;
          tmp      = LIMIT_RANGE(0.,tmp,max_nu_horizontal[2]);
          DIFFUSION_COEF_UV(K,J,I) = tmp;

          /*
           * Mass diffusivity.
           * Currently using the same value for turbulent mass diffusivity
           * as for temperature.
           */
          DIFFUSION_COEF_MASS(K,J,I) = tmp;
        }
      }
    }

    /*
     * DIFFUSION_COEF_THETA is on the p3-grid.
     */
    for (K = KLO; K < KHI; K++) {
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          /*
           * Take the turbulent Prandtl number to be unity,
           * following Collins et al. (2004, NCAR/TN-464+STR).
           */
          nu_turb  = .5*(NU_TURB(K,J,I)+NU_TURB(K+1,J,I));
          chi3     = nu_turb/planet->kinvisc;
          chi3    *= chi3*chi3;
          fv1      = chi3/(chi3+cv1_3);
          turb     = fv1*nu_turb;
          tmp      = turb;
          tmp      = LIMIT_RANGE(0.,tmp,max_nu_horizontal[2]);
          DIFFUSION_COEF_THETA(K,J,I) = tmp;
        }
      }
    }
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        /* THETA is on the p3-grid. */
        DIFFUSION_COEF_THETA(0,  J,I) = 0.;
        DIFFUSION_COEF_THETA(KHI,J,I) = 0.;
      }
    }
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }
 
  /* Need to apply bc_lateral() here, because get_kin() above doesn't work on pads. */
  bc_lateral(var.diffusion_coef_uv.value,THREEDIM);
  bc_lateral(var.diffusion_coef_theta.value,THREEDIM);
  bc_lateral(var.diffusion_coef_mass.value,THREEDIM);

  return;
}

/*======================= end of set_diffusion_coef() ==============*/

/*======================= stability_factor() =======================*/
/*
 * Stability factor for vertical turbulent diffusion.
 * This function fills in the entire kk column of stab_factor[],
 * which should include the range [1,2*nk+1] as in the declaration 
 * of Ri[] below.
 */
#define MAX_RI 100.

/*
 * Need to smooth Ri[kk].
 * Use a Savitzky-Golay filter of order m, width np.
 */
#undef  SAVITZKY_GOLAY_M
#define SAVITZKY_GOLAY_M   2
#undef  SAVITZKY_GOLAY_NP
#define SAVITZKY_GOLAY_NP (2*(SAVITZKY_GOLAY_M+2)+1)

void stability_factor(int         J,
                      int         I,
                      EPIC_FLOAT *stab_factor)
{
  static int
    nr,nl,
    nc = SAVITZKY_GOLAY_NP/2,
    initialized = FALSE;
  register int
    K,kk,kay,
    j,i;
  EPIC_FLOAT
    smooth_ri[2*grid.nk+2],
    extended_ri[2*grid.nk+2+2*nc],
    xleft[ SAVITZKY_GOLAY_M+1],
    xright[SAVITZKY_GOLAY_M+1],
   *padded_ri;
  EPIC_FLOAT
    dy;
  static EPIC_FLOAT
    c_sav[SAVITZKY_GOLAY_NP];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="stability_factor";

  if (!initialized) {
    /* Centered filter */
    nl = nr = nc;
    savitzky_golay(c_sav,SAVITZKY_GOLAY_NP,nl,nr,0,SAVITZKY_GOLAY_M);

    initialized = TRUE;
  }

  /*
   * Shift padded_ri.
   */
  padded_ri = extended_ri+nc;

  stab_factor[0] = 1.;

  /*
   * Smooth Ri to remove computational oscillations
   * in low N^2, low (du/dz)^2 regions.
   */
  for (kk = 1; kk <= 2*KHI+1; kk++) {
    padded_ri[kk] = get_richardson(planet,kk,J,I);
  }

  /*
   * Extend ends of padded_ri[kk] with a polynomial fit of order SAVITZKY_GOLAY_M,
   * fitted to the respective end-point data.
   */
  for (kay = 0; kay <= SAVITZKY_GOLAY_M; kay++) {
    xleft[ kay] = (EPIC_FLOAT)(kay+1);
    xright[kay] = (EPIC_FLOAT)(kay+2*grid.nk+1-SAVITZKY_GOLAY_M);
  }
  for (kay = -nc; kay < 0; kay++) {
    padded_ri[kay+1]           = poly_interp(SAVITZKY_GOLAY_M+1,
                                             xleft, padded_ri+1,                           (EPIC_FLOAT)(kay+1),          &dy);
    padded_ri[2*grid.nk+1-kay] = poly_interp(SAVITZKY_GOLAY_M+1,
                                             xright,padded_ri+2*grid.nk+1-SAVITZKY_GOLAY_M,(EPIC_FLOAT)(2*grid.nk+1-kay),&dy);
  }

  /* 
   * Limit range of |Ri| to MAX_RI.
   */
  for (kk = 1-nc; kk <= 2*KHI+1+nc; kk++) {
    padded_ri[kk] = LIMIT_RANGE(-MAX_RI,padded_ri[kk],MAX_RI);
  }

  /* Zero smooth array. */
  memset(smooth_ri,0,(2*grid.nk+2)*sizeof(EPIC_FLOAT));

  for (kk = 1; kk <= 2*KHI+1; kk++) {
    for (kay = kk-nl,i = 0; kay <= kk+nr; kay++,i++) {
      /* c_sav[] uses a wrap-around index */
      j             = (SAVITZKY_GOLAY_NP+nl-i)%SAVITZKY_GOLAY_NP;
      smooth_ri[kk] += padded_ri[kay]*c_sav[j];
    }
  }

  for (kk = 1; kk <= 2*KHI+1; kk++) {
    /*
     * From the NCAR Community Atmosphere Model (CAM 3.0),
     * as described by Collins et al. (2004, NCAR/TN-464+STR).
     */ 
    if (smooth_ri[kk] > 0.) {
      stab_factor[kk] = 1./(1.+10.*smooth_ri[kk]*(1.+8.*smooth_ri[kk]));
    }
    else {
      stab_factor[kk] = sqrt(1.-18.*smooth_ri[kk]);
    }
  }

  return;
}

#undef  MAX_RI

/*======================= end of stability_factor() ================*/

/*======================= source_sink_turb() =======================*/

/*
 * Wrapper for turbulence source-sink subroutine.
 * This allows flexibility in what turbulence model is used without
 * having to change the hook to the rest of the EPIC model.
 *
 * NOTE: There is a name conflict with "source_sink_subgrid" and 
 *       LAM MPI, in the file lam_config_file.h.
 */
void source_sink_turb(planetspec  *planet,
                      EPIC_FLOAT **Buff2D)
{

  source_sink_SA(planet,Buff2D);

  return;
}

/*======================= end of source_sink_turb() ================*/

/*======================= source_sink_SA() =========================*/

/*
 * See Dowling et al. (2006, Icarus 182, 259-273), Section 5.5 for details on the 
 * Spalart-Allmaras one-equation turbulence model implemented here.
 */

#define UUU(k,j,i) uuu[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define VVV(k,j,i) vvv[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define WWW(k,j,i) www[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void source_sink_SA(planetspec  *planet,
                    EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kk,jj;
  int
    itmp;
  const EPIC_FLOAT    
    cw2             = 0.3,
    cw3_6           = pow(2.,6.),
    cb1             = 0.1355,
    cb2             = 0.622,
    cv1             = 7.1,
    cv1_3           = cv1*cv1*cv1,
    sigma           = 2./3.,
    kappa           = 0.41,
    kappa_2         = kappa*kappa,
    cw1             = (cb1/kappa_2)+((1.+cb2)/sigma),
    C_DES           = 0.65;
  EPIC_FLOAT
    ptop,pbot,
    var1,var3,S,
    u_var1,u_var3,
    v_var1,v_var3,
    w_var1,w_var3,
    s11,s22,s33,s12,s13,s23,
    r_sa,g_sa,chi,fv1,fv2,fw,d_tilda,S_tilda,
    sig_z_conv,delta,
    tmp1,tmp2,dz_inv,
    nu_turb,u_tan,kin,chi3,
    dnudt;
  EPIC_FLOAT
    *u2d,*v2d;
  static EPIC_FLOAT
    *uuu,*vvv,*www;
  register EPIC_FLOAT
    m_2jp1,n_2j,n_2jp1,n_2jp2,
    m_2j_inv,n_2jp1_inv,
    mn_2jm1_inv,mn_2j_inv,
    mn_2jp1_inv,mn_2jp2_inv,
    mn_u,mn_v;
  static unsigned long
    nbytes_2d;
  static int
    solid_surface,
    initialized = FALSE;
  static double
    max_nu_horizontal[MAX_NU_ORDER+1];
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
    dbmsname[]="source_sink_SA";

  if (!initialized) {
    nbytes_2d = (unsigned long)(Nelem2d*sizeof(EPIC_FLOAT));

    /* Allocate memory */
    uuu = fvector(0,Nelem3d-1,dbmsname);
    vvv = fvector(0,Nelem3d-1,dbmsname);
    www = fvector(0,Nelem3d-1,dbmsname);

    if (strcmp(planet->type,"terrestrial") == 0) {
      solid_surface = TRUE;
    }
    else if (strcmp(planet->type,"gas-giant") == 0) {
      solid_surface = FALSE;
    }
    else {
      sprintf(Message,"unrecognized planet->type=%s",planet->type);
      epic_error(dbmsname,Message);
    }

    set_max_nu(max_nu_horizontal);

    initialized = TRUE;
  }

  /*
   * D_WALL(K,J,I) is the distance to the wall.
   * For gas giants (no wall), dwall_SA sets this to FLOAT_MAX.
   */
  dwall_SA(planet,d_wall);

  memcpy(uuu,var.u.value+grid.it_uv_dis*Nelem3d,Nelem3d*sizeof(EPIC_FLOAT));
  memcpy(vvv,var.v.value+grid.it_uv_dis*Nelem3d,Nelem3d*sizeof(EPIC_FLOAT));
  memcpy(www,var.dzdt2.value,                   Nelem3d*sizeof(EPIC_FLOAT));

  /*
   * The K = KHI case is done below for terrestrial planets.
   * For gas-giant planets, it is skipped to avoid issues with having the vertical differencing
   * stencil include the static abyssal wind.
   */

  for (K = KLO+1; K < KHI; K++) {
    kk = 2*K;

    /*
     * Assign zero'd memory to u2d and v2d:
     */
    memset(Buff2D[0],0,nbytes_2d);
    memset(Buff2D[1],0,nbytes_2d);
    u2d = Buff2D[0];
    v2d = Buff2D[1];
    /*
     * Averaging U values onto V grid and storing in U2D array.
     */
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;

      mn_2jp1_inv = 1./grid.mn[kk][jj+1];
      mn_2jm1_inv = 1./grid.mn[kk][jj-1];
      mn_v        = .5/(mn_2jm1_inv+mn_2jp1_inv);
      for (I = ILO; I <= IHI; I++) {
        U2D(J,I) = ((UUU(K,J,  I)+UUU(K,  J,I+1))*mn_2jp1_inv
                   +(UUU(K,J-1,I)+UUU(K,J-1,I+1))*mn_2jm1_inv)*mn_v;
      }
    }
    if (JLO == grid.jlo) {
      J      = grid.jlo;
      u_var1 = 0.;
      if (!IS_SPOLE) {
        /* 
         * Need U on southern-edge of channel.
         */
        for (I = ILO; I <= IHI; I++) {
          u_var1 += UUU(K,J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = u_var1;
        MPI_Allreduce(&mpi_tmp,&u_var1,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        u_var1 /= grid.ni;
      }
      for (I = ILO; I <= IHI; I++) {
        U2D(J,I) = u_var1;
      }
    }
    if (JHI == grid.nj) {
      J      = grid.nj;
      u_var3 = 0.;
      if (!IS_NPOLE) {
        /* 
         * Need U on northern-edge of channel.
         */
        for (I = ILO; I <= IHI; I++) {
          u_var3 += UUU(K,J,I);
        }

#if defined(EPIC_MPI)
        mpi_tmp = u_var3;
        MPI_Allreduce(&mpi_tmp,&u_var3,1,float_type,MPI_SUM,para.comm_JLO);
#endif

        u_var3 /= grid.ni;
      }
      for (I = ILO; I <= IHI; I++) {
        U2D(J+1,I) = u_var3;
      }
    }
    bc_lateral(u2d,TWODIM);

    /*
     * Averaging V values onto U grid and 
     * storing in V2D array.
     */
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J;

      mn_2j_inv   = 1./grid.mn[kk][jj  ];
      mn_2jp2_inv = 1./grid.mn[kk][jj+2];
      mn_u        = .5/(mn_2j_inv+mn_2jp2_inv);
      for (I = ILO; I <= IHI; I++) {
	    V2D(J,I) = ((VVV(K,J,  I)+VVV(K,J,  I-1))*mn_2j_inv
		       +(VVV(K,J+1,I)+VVV(K,J+1,I-1))*mn_2jp2_inv)*mn_u;
      }
    }
    bc_lateral(v2d,TWODIM);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J;

      m_2jp1 = grid.m[kk][jj+1];
      n_2jp1 = grid.n[kk][jj+1];
      n_2jp2 = grid.n[kk][jj+2];
      n_2j   = grid.n[kk][jj  ];
      for (I = ILO; I <= IHI; I++) {
        delta   = delta_SA(planet,K,J,I);
        d_tilda = MIN(D_WALL(K,J,I),C_DES*delta);
        chi     = NU_TURB(K,J,I)/planet->kinvisc;
        fv1     = pow(chi,3.)/(pow(chi,3.)+cv1_3);
        fv2     = 1.-(chi/(1.+chi*fv1));

        /*
         *  Turbulence eddy viscosity production
         */
        dz_inv = 1./(Z2(K-1,J,I)-Z2(K+1,J,I));

        s11    = (UUU(K,J,I+1)-UUU(K,J,I))*m_2jp1;
        s22    = (VVV(K,J+1,I)-VVV(K,J,I))*n_2jp1;
        s33    = (WWW(K-1,J,I)-WWW(K+1,J,I))*dz_inv;
        s12    = .5*((U2D(J+1,I  )-U2D(J,I))*n_2jp1
                    +(V2D(J,  I+1)-V2D(J,I))*m_2jp1);

        u_var1 = 0.5*(UUU(K-1,J,I)+UUU(K-1,J,I+1));
        u_var3 = 0.5*(UUU(K+1,J,I)+UUU(K+1,J,I+1));

        w_var1 = 0.5*(WWW(K,J,I)+WWW(K,J,I-1));
        w_var3 = 0.5*(WWW(K,J,I)+WWW(K,J,I+1));

	s13    = .5*((u_var1-u_var3)*dz_inv+(w_var3-w_var1)*m_2jp1);

        v_var1 = 0.5*(VVV(K-1,J,I)+VVV(K-1,J+1,I));
        v_var3 = 0.5*(VVV(K+1,J,I)+VVV(K+1,J+1,I));

        w_var1 = 0.5*(WWW(K,J,I)+WWW(K,J-1,I));
        w_var3 = 0.5*(WWW(K,J,I)+WWW(K,J+1,I));

	s23    = .5*((v_var1-v_var3)*dz_inv+(w_var3-w_var1)*n_2jp1);

        /*
         * 2.*s12 = s12+s21, etc. 
         */
	S = sqrt(2.*(s11*s11+s22*s22+s33*s33+2.*(s12*s12+s13*s13+s23*s23)));

        tmp1    = kappa_2*d_tilda*d_tilda;

        S_tilda = MAX(1.e-20,S+NU_TURB(K,J,I)*fv2/tmp1);

        tmp2    = S_tilda*tmp1;
        dnudt   = cb1*S_tilda*NU_TURB(K,J,I);

        dnudt += (cb2/sigma)*(SQR(    (NU_TURB(K,J,I+1)-NU_TURB(K,J,I-1))*m_2jp1*.5)
                             +SQR(.5*((NU_TURB(K,J+1,I)-NU_TURB(K,  J,I))*n_2jp2
                                     +(NU_TURB(K,J,  I)-NU_TURB(K,J-1,I))*n_2j  )  ));

        dnudt += (cb2/sigma)*SQR((NU_TURB(K-1,J,I)-NU_TURB(K+1,J,I))*dz_inv);

        /*
         *  SA-model turbulence destruction term.
         */
        if (tmp2 != 0.) {
          r_sa = NU_TURB(K,J,I)/tmp2;
          if (r_sa < 0.) {
            sprintf(Message,"r_sa=%g < 0., NU_TURB(%2d,%2d,%2d)=%g",r_sa,K,J,I,NU_TURB(K,J,I));
            epic_error(dbmsname,Message);
          }
          else if (r_sa < 10.) {
            g_sa = r_sa+cw2*(pow(r_sa,6.)-r_sa);
            fw   = g_sa*pow(((1.+cw3_6)/(pow(g_sa,6.)+cw3_6)),(1./6.));
          }
          else {
            fw = 1.;
          }
        }
        else {
          fw = 1.;
        }
        dnudt -= cw1*fw*SQR(NU_TURB(K,J,I)/d_tilda);

	NU_TURB(K,J,I) += DT*dnudt;

        NU_TURB(K,J,I) = LIMIT_RANGE(1.e-4*planet->kinvisc,NU_TURB(K,J,I),max_nu_horizontal[2]);
      }
    }
  }
  
  if (solid_surface) {
    K   = KHI;
    kk  = 2*K;
    u2d = uuu+(K-Kshift)*Nelem2d;
    v2d = vvv+(K-Kshift)*Nelem2d;
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        nu_turb          = NU_TURB(K,J,I);
        kin              = get_kin(planet,u2d,v2d,kk,J,I);
        u_tan            = sqrt(2.*kin);
        nu_turb          = law_of_the_wall(planet,K,J,I,NU_TURB_INDEX,nu_turb,u_tan);
        nu_turb          = invert_fv1(planet,nu_turb);
        NU_TURB(KHI,J,I) = nu_turb;
        NU_TURB(KHI,J,I) = LIMIT_RANGE(1.e-4*planet->kinvisc,NU_TURB(K,J,I),max_nu_horizontal[2]);
      }
    }
  }
  
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.nu_turb.value,THREEDIM);

  return;
}

#undef UUU
#undef VVV
#undef WWW

/*======================= end of source_sink_SA() ===========================*/

/*======================= dwall_SA() ========================================*/

void dwall_SA(planetspec *planet,
              EPIC_FLOAT *d_wall)
{
  register int
    K,J,I;
  register EPIC_FLOAT
    g0_inv;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="dwall_SA";

  if (strcmp(planet->type,"gas-giant") == 0) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD;I++) {
          /*
           * Set distance to "wall" to FLOAT_MAX.
           */
          D_WALL(K,J,I) = FLOAT_MAX;
        }
      }
    }
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    g0_inv = 1./grid.g[2*KHI+1][2*(grid.nj/2)+1];
    for (K = KLO-1; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD;I++) {
          /*
           * Calculating distance from the wall.
           */
          D_WALL(K,J,I) = (PHI2(K,J,I)-PHI_SURFACE(J,I))*g0_inv;
          /*
           * Screen for non-positive value.
           */
          if (D_WALL(K,J,I) <= 0.) {
            sprintf(Message,"KJI=%2d %2d %2d; PHI2=%g, PHI_SURFACE=%g, g0_inv=%g\n",
                             K,J,I,PHI2(K,J,I),PHI_SURFACE(J,I),g0_inv);
            epic_error(dbmsname,Message);
          }
        }
      }
    }
  }
  else {
    sprintf(Message,"unrecognized planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

  return;
}
      

/*======================= end of dwall_SA() ====================================*/

/*======================= delta_SA() ===========================================*/

/*
 * Compute delta in (5.9) of Dowling et al (2006) for the Spalart-Allmaras
 * turbulence model.
 *
 * NOTE: We have changed the definition from the Dowling et al (2006) paper.
 */

EPIC_FLOAT delta_SA(planetspec *planet,
                    int         K,
                    int         J,
                    int         I)
{
  register int
    kk = 2*K,
    jj = 2*J+1;
  EPIC_FLOAT
    dx,dy,dz,
    delta;
  
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="delta_SA";


  dx = 1./grid.m[kk][jj];
  dy = 1./grid.n[kk][jj];
  dz = (Z3(K-1,J,I)-Z3(K,J,I));

  delta = MIN(dx,dy);
  delta = MIN(delta,dz);

  return delta;
}

/*======================= end of delta_SA() =================================*/

/*======================= tau_surface() =====================================*/

/*
 * Calculate the shear stress at the surface.
 *
 * NOTE: Because of the call to get_kin, do not call this function
 *       from a pad position (like J=JLOPAD or I=IHIPAD).
 */

void tau_surface(planetspec  *planet,
                 int          index,
                 EPIC_FLOAT  *tau_wall,
                 EPIC_FLOAT  *buffji)
{
  register int
    K,J,I;
  EPIC_FLOAT
    chi3,fv1,
    rho,diffusion_coef,nu_turb,dz,
    kin,u_tan,theta0,
    *u2d,
    *v2d;
  EPIC_FLOAT
    *kie;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1;
  static int
    nbytes_2d,
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="tau_surface";

  if (!initialized) {
    nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

    initialized = TRUE;
  }

  K   = grid.nk;

  u2d = var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d;
  v2d = var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d;

  memset(buffji,0,nbytes_2d);
  kie = buffji;

  for (J = JLO; J <= JHI; J++) {
    for (I = ILO; I <= IHI; I++) {
      KIE(J,I) = get_kin(planet,u2d,v2d,2*K,J,I);
    }
  }
  bc_lateral(kie,TWODIM);

  if (strcmp(planet->type,"gas-giant") == 0) {
    switch(index) {
      case U_INDEX:
      case V_INDEX:
      case THETA_INDEX:
      case NU_TURB_INDEX:
        /*
         * NOTE: Placeholder, need better tau values for bottom of gas-giant atmosphere.
         */
        memset(tau_wall,0,sizeof(EPIC_FLOAT)*Nelem2d);
      break;
      default:
        sprintf(Message,"unrecognized index=%d",index);
        epic_error(dbmsname,Message);
      break;
    }
  }
  else if (strcmp(planet->type,"terrestrial") == 0) {
    dwall_SA(planet,d_wall);

    switch(index) {
      case U_INDEX:
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            rho            = .5*(RHO2(K,J,I)+RHO2(K,J,I-1));
            nu_turb        = .5*(NU_TURB(K,J,I)+NU_TURB(K,J,I-1));
            dz             = .5*(Z2(K,J,I)+Z2(K,J,I-1)-Z3(K,J,I)-Z3(K,J,I-1));
            kin            = .5*(KIE(J,I)+KIE(J,I-1));
            u_tan          = sqrt(2.*kin);
            diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
            TAU_WALL(J,I)  = (planet->dynvisc+rho*diffusion_coef)*(U(grid.it_uv,K,J,I)-0.)/dz;	
          }
        }
      break;
      case V_INDEX:
        for (J = JFIRST; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            rho            = .5*(RHO2(K,J,I)+RHO2(K,J-1,I));
            nu_turb        = .5*(NU_TURB(K,J,I)+NU_TURB(K,J-1,I));
            dz             = .5*(Z2(K,J,I)+Z2(K,J-1,I)-Z3(K,J,I)-Z3(K,J-1,I));
            kin            = .5*(KIE(J,I)+KIE(J-1,I));
            u_tan          = sqrt(2.*kin);
            diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
            TAU_WALL(J,I)  = (planet->dynvisc+rho*diffusion_coef)*(V(grid.it_uv,K,J,I)-0.)/dz;
          }
        }
      break;
      case NU_TURB_INDEX:
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            rho            = RHO2(K,J,I);
            nu_turb        = NU_TURB(K,J,I);
            dz             = (Z2(K,J,I)-Z3(K,J,I));
            kin            = get_kin(planet,u2d,v2d,2*K,J,I);
            u_tan          = sqrt(2.*kin);
            diffusion_coef = law_of_the_wall(planet,K,J,I,index,nu_turb,u_tan);
            nu_turb        = invert_fv1(planet,diffusion_coef);
            TAU_WALL(J,I)  = (planet->dynvisc+rho*nu_turb)*(u_tan-0.)/dz;
          }
        }
      break;
      default:
        sprintf(Message,"unrecognized index=%d",index);
        epic_error(dbmsname,Message);
      break;
    }
  }
  else {
    sprintf(Message,"unrecognzied planet->type=%s",planet->type);
    epic_error(dbmsname,Message);
  }

#if EPIC_CHECK == 1
  /*
   * Screen for nan.
   */
  for (J = JLO; J <= JHI; J++) {
     for (I = ILO; I <= IHI; I++) {
      if (!isfinite(TAU_WALL(J,I))) {
        sprintf(Message,"index=%d, TAU_WALL(%d,%d)=%g",index,J,I,TAU_WALL(J,I));
        epic_error(dbmsname,Message);
      }
    }
  }
#endif

  /* Need to apply bc_lateral() here. */
  bc_lateral(tau_wall,TWODIM);

  return;
}

/*======================= end of tau_surface() ====================*/

/*======================= law_of_the_wall() =======================*/

/*
 * Raymond P. LeBeau
 *
 * For the DES-type turbulence model there needs to be a boundary condition
 * placed on the value of the turbulent viscosity close to the wall.
 * Ideally, one gets close enough to the wall with the closest grid point
 * that the flow is effectively in the viscous sublayer, where the turbulent
 * viscosity is set at 0; however this requires a very small spacing.
 * Further out, one can use the law-of-the-wall correlation, which follows: 
 *   u+ = (1/.41)*log(y+)+5.0, u+ = u/u_tau, y+ = y*u_tau/nu,
 * u being the tagential velocity over the surface.
 * The routine returns the value for the turbulent visocity in 
 * t_vis_new on the bottom boundary.

 * Modified by Aaron Herrnstein on 03-13-06.  The extra argument "index"
 * is used to specify the h-grid, u-grid, or v-grid.
 *
 * Modified by Tim Dowling on 07-14-08. Added K argument.
 */

EPIC_FLOAT law_of_the_wall(planetspec *planet,
                           int         K,
                           int         J,
                           int         I,
                           int         index,
			   EPIC_FLOAT  t_vis,
			   EPIC_FLOAT  u_tan)
{
  EPIC_FLOAT
    t_vis_new,dwall,dwall_inv,
    u_tau,
    x1,x2,xl,dx,
    fl,f,
    rts,swap,
    chi3,fv1;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1,
    tol        = 1.e-6;
  int
    kk = 2*K,
    iter,
    itmax = 30;
  static int
    warned = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="law_of_the_wall";

  if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * For gas giants, set the turbulent viscosity on the boundary to be
     * negligible, for lack of better information.
     */
    return 1.0e-4*planet->kinvisc;
  }
  else {
    /*
     * Compute distance of grid point to surface.
     */
    if (index == U_INDEX) {
      dwall = 0.5*(D_WALL(K,J,I)+D_WALL(K,J,I-1));
    }
    else if (index == V_INDEX  ||  index == PV2_INDEX) {
      dwall = 0.5*grid.m[kk][2*J]*(D_WALL(K,J  ,I)/grid.m[kk][2*J+1] 
                                  +D_WALL(K,J-1,I)/grid.m[kk][2*J-1]);
    }
    else {
      /*
       * Default is the h-grid.
       */
      dwall = D_WALL(K,J,I);
    }

    /*
     * If the combination of distance and horizontal velocity is very small,
     * assume viscous sublayer such that turb visc << laminar visc.
     */
    if (dwall*u_tan < 0.00023) {
      return 0.01*planet->kinvisc;
    }

    /*
     * Check that the wall distance is positive.
     */
    if (dwall > 0.) {
      dwall_inv = 1./dwall;
    }
    else {
      sprintf(Message,"JI=%d %d, dwall = %g",J,I,dwall);
      epic_error(dbmsname,Message);
    }

    /*
     * Start secant method to solve for boundary condition.
     */
    x1 = sqrt(planet->kinvisc*(1.+ 0.01)*u_tan*dwall_inv);
    x2 = sqrt(planet->kinvisc*(1.+20.  )*u_tan*dwall_inv);
    fl = func_utau(x1,u_tan,dwall);
    f  = func_utau(x2,u_tan,dwall);

    /*
     * Essentially, we have two different definitions of u_tau:
     *   1) the correlation in func_utau()
     *   2) sqrt(t_vis*du/dy)
     * When they are equal, we have the necessary value of t_vis.
     */
    if(fabs(fl) < fabs(f)) {
      rts  = x1;
      xl   = x2;
      swap = fl;
      fl   = f;
      f    = swap;
    }
    else {
      xl  = x1;
      rts = x2;
    }

    /*
     * Iterative search to match different estimates of u_tau.
     */
    for (iter = 0; iter < itmax; iter++) {
      dx   = (xl-rts)*f/(f-fl);
      xl   = rts;
      fl   = f;
      rts += dx;
      f    = func_utau(rts,u_tan,dwall);
      /*
       * func_utau() returns the difference between the two u_tau calculations,
       * so when it is small, we have a solution.  Note that a relatively small
       * f will result in a small dx; see eqn. above.
       */
      if (f > .99e+20) {
        /*
         * A just-in-case boundary condition; hopefully it should not typically come up.
         */
        t_vis_new = 0.01*planet->kinvisc;

        return t_vis_new;
      }
      if (fabs(dx) < tol || f == 0.0) {
        u_tau     = rts;
        t_vis_new = MAX((u_tau*u_tau*dwall/u_tan)-planet->kinvisc,0.01*planet->kinvisc);

        return t_vis_new;
      }
    }
    /*
     * An emergency answer if the iteration above does not converge.
     */
    t_vis_new = 20.*planet->kinvisc;

    return t_vis_new;
  }

  /* Should never get here. */
  return 0.;
} 

/*=================== end of law_of_the_wall() ====================*/

/*=================== func_utau() =================================*/

EPIC_FLOAT func_utau(EPIC_FLOAT u_tau, 
                     EPIC_FLOAT u_tan, 
                     EPIC_FLOAT dwall)
{
  register EPIC_FLOAT
    res,u_plus,y_plus,
    yy;
  const EPIC_FLOAT
    C              = 5.,
    const1         = 0.127,
    const2         = 1./.41,
    const3         = 0.41*1.43e-3,
    one_third      = 1./3.,
    sqrt_one_third = sqrt(1./3.),
    pi_six         = M_PI/6.;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="func_utau";
  
  if (u_tau != 0.) {
    u_plus = u_tan/u_tau;
  }
  else {
    sprintf(Message,"u_tau=%g",u_tau);
    epic_error(dbmsname,Message);
  }
  y_plus = fabs(u_tau*dwall/planet->kinvisc);

  if (u_plus < 0.1) {
    /*
     * Very close to wall.
     */
    return 1.e+20;
  }
  
  if (y_plus < 5.) {
    /*
     * Close to wall, viscous sublayer.
     */
    res = y_plus-u_plus;
  }
  else if (y_plus < 70.) {
    /*
     * Between sublayer and log layer, smooth correlation.
     */
    yy  = const1*y_plus;
    res = (      one_third*log((1.+yy)/(sqrt(yy*(yy-1.)+1.)))
           +sqrt_one_third*(atan(sqrt_one_third*(2.*yy-1.))+pi_six) )/const1
         +.25*const2*log(1.+const3*pow(y_plus,4.))-u_plus;
  }
  else {
    /*
     * log layer
     */
    res = const2*log(y_plus)+C-u_plus;
  }

  /*
   * RP LeBeau:
   * Could add a correlation for the outer layer (y > 100s-1000s), but trickier
   * since dependent on pressure gradient and other things.
   */
  
  return res;
}

/*=================== end of func_utau() ==========================*/

/*======================= invert_fv1() ============================*/

/* 
 * Calculate nu_turb (nu_tilde) from DIFF_COEF (nu_t).
 */
EPIC_FLOAT invert_fv1(planetspec *planet,
                      EPIC_FLOAT  t_vis)

{
  EPIC_FLOAT
    nu_turb,
    chi3,fv1,chi;
  const EPIC_FLOAT
    cv1        = 7.1,
    cv1_3      = cv1*cv1*cv1,
    pi_factor  = M_PI*0.5/30.0;
  int
    iter,
    it_max=100;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="invert_fv1";

  /*
   * RP LeBeau:
   * The following fits are good to about +/-5%, will try to improve later.
   * Fits used to try to speed up the computation.
   * Have to go backwards from the nu_t we use to compute diffusion terms to the 
   * nu_tilde used in the computation of SA-DES model at times (BC and the like).
   * Note chi here is t_vis/planet->kinvisc, not the high Re version.
   */

  if (t_vis < 0.) {
    sprintf(Message,"t_vis=%g<0.",t_vis);
    epic_error(dbmsname,Message);
  }

  chi = t_vis/planet->kinvisc;
  if (chi > 30) {
    /* fvl ~ 1.0 */
    nu_turb = t_vis;
  }
  else if (chi < 0.4) {
    /* chi3 in denominator is small vs. 7.1^3  */
    nu_turb = planet->kinvisc*pow(chi*cv1_3,.25);
  }
  else {
    /* mid-range, fit */
    fv1 = cos(pi_factor*pow(chi,0.35));
    fv1 = 1.0-fv1*fv1;
    nu_turb = planet->kinvisc*chi/fv1;
  }

  return nu_turb;
} 

/*=================== end of invert_fv1() =========================*/

/* * * * * * * * * * end of epic_subgrid.c * * * * * * * * * * * * */




