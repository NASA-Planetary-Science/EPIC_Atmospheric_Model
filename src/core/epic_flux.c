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

/*
 * Functions contained in this file:
 *    flux_top()
 *    flux_bot()
 *    calc_w()
 *    advection()
 *      hsu_advection()
 *      upwind_3rd_order()
 *    uv_vertical_advection()
 *      uv_vert_upwind_3rd_order()
 */

#include <epic.h>

/*======================= flux_top() =======================================*/

/*
 * Hybrid mass flux for variable at the top of the model.
 *
 * For index == H_INDEX, this is WH.
 * For index == H_2O_INDEX, this is WH*(mixing ratio of H_2O), etc.
 *
 * Returns as an argument the hybrid vertical velocity at the
 * top of the model, wtop.
 */

EPIC_FLOAT flux_top(planetspec *planet,
                    int         index,
                    int         J,
                    int         I,
                    EPIC_FLOAT *wtop)
{
  EPIC_FLOAT
    flux;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="flux_top";

  flux  = 0.;
  if (wtop) {
   *wtop = 0.;
  }

  return flux;
}

/*======================= end of flux_top() ================================*/

/*======================= flux_bot() =======================================*/

/*
 * Hybrid mass flux for variable at the bottom of the model.
 *
 * For index == H_INDEX, this is WH.
 * For index == H_2O_INDEX, this is WH*(mixing ratio of phase), etc.
 *
 * Returns as an argument the hybrid vertical velocity at the
 * bottom of the model, wbot.
 */

EPIC_FLOAT flux_bot(planetspec *planet,
                    int         index,
                    int         J,
                    int         I,
                    EPIC_FLOAT *wbot)
{
  EPIC_FLOAT
    flux;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="flux_bot";

  flux  = 0.;
  if (wbot) {
    *wbot = 0.;
  }

  return flux;
}

/*======================= end of flux_bot() ================================*/

/*======================= calc_w() =========================================*/

/*
 * Calculate the vertical velocity, W3 = (D/Dt)sigmatheta,
 * at all the layer interfaces.
 *
 * NOTE: For action == STEP_PROGS_AND_SYNC_DIAGS, DZDT2 is started in
 *       store_pgrad_vars(), where the partial derivative with respect to time
 *       is calculated.
 *
 * Assumes the hybrid vertical coordinate, zeta (aka sigmatheta or sgth),
 * is given by zeta = F(p,pbot,theta) = f(sigma) + g(sigma)*theta.
 */

#define SUM_DIV3(k,j,i)       sum_div3[i+(j)*Iadim+(k)*Nelem2d-Shift3d]
#define HORIZONTAL_DIV(k,j,i) horizontal_div[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

void calc_w(planetspec *planet,
            int         action)
{
  register int
    K,J,I,
    kk,jj;
  EPIC_FLOAT
    pbot,sigma,sgdot,dFdsg,
    lnptoppbot_inv,dgdsg,
    dsumdiv,
    x,gsg,tmp,
    dx_inv3,dy_inv2,dy_inv4;
  const double
    coeff    = grid.hybrid_alpha/(tanh( grid.hybrid_alpha*(1.-grid.sigma_sigma))
                                 -tanh(-grid.hybrid_alpha*(   grid.sigma_sigma))),
    zeta_1_0 = grid.zeta1-grid.zeta0;
  static int
    initialized = FALSE;         
  static EPIC_FLOAT
    *sum_div3;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="calc_w";

  if (!initialized) {
    /* Allocate memory */
    sum_div3 = fvector(0,Nelem3d-1,dbmsname);

    initialized = TRUE;
  }

  /* 
   * Assign the coordinate vertical velocity, W3, for the
   * top and bottom of the model. The functions flux_top() and flux_bot() return
   * the coordinate vertical velocity as an argument.
   */
  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      flux_top(planet,H_INDEX,J,I,&W3(  0,J,I));
      flux_bot(planet,H_INDEX,J,I,&W3(KHI,J,I));
    }
  }
  /* No need to apply bc_lateral() here. */

  /*
   * Heating component of W3.
   */
  switch(grid.coord_type) {
    case COORD_ISOBARIC:
      ;
    break;
    case COORD_ISENTROPIC:
      /*
       * For isentropic vertical coordinates, the vertical velocity is just the heating, divided
       * by the Exner function.
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            W3(K,J,I) = HEAT3(K,J,I)/EXNER3(K,J,I);
          }
        }
      }
      /* No need to apply bc_lateral() here. */
    break;
    case COORD_HYBRID:
      /*
       * Calculate the heating contribution to W3 for the
       * model's interior interfaces (the dF/dtheta term).
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot = P3(KHI,J,I);
          for (K = KLO; K < KHI; K++) {
            sigma     = get_sigma(pbot,P3(K,J,I));
            W3(K,J,I) = g_sigma(sigma)*HEAT3(K,J,I)/EXNER3(K,J,I);
          }
        }
      }
      /* No need to apply bc_lateral() here. */
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * Calculate and store the summed horizontal-divergence terms.
   */

  /* Zero-out SUM_DIV array. */
  memset(sum_div3,0,Nelem3d*sizeof(EPIC_FLOAT));

  if (grid.coord_type != COORD_ISENTROPIC) {
    for (K = KLO; K <= KHI; K++) {
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          dsumdiv         = DIV_UV2(K,J,I)*(P3(K,J,I)-P3(K-1,J,I));
          SUM_DIV3(K,J,I) = SUM_DIV3(K-1,J,I)+dsumdiv;
        }
      }
    }
  }

  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
      ;
    break;
    case COORD_ISOBARIC:
      /*
       * In the isobaric case, internally W3 = D(zeta)/Dt [K/s], with zeta = f(sigma) [K],
       * sigma = log(p/pbot)/log(ptop/pbot), and constant pbot, ptop [Pa].
       */
      tmp = zeta_1_0/log(grid.ptop/grid.pbot);

      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            W3(K,J,I) = -SUM_DIV3(K,J,I)*tmp/P3(K,J,I);
          }
        }
      }
    break;
    case COORD_HYBRID:
      /*
       * Add dF/dsigma terms to W3.
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot           = P3(KHI,J,I);
          lnptoppbot_inv = 1./log(grid.ptop/pbot);

          /*
           * NOTE: Assumes g_sigma() is in tanh() form.
           */
          for (K = KLO; K < KHI; K++) {
            sigma      = get_sigma(pbot,P3(K,J,I));
            sgdot      = -lnptoppbot_inv*(SUM_DIV3(K,J,I)/P3(K,J,I)-SUM_DIV3(KHI,J,I)*(1.-sigma)/pbot);
            dgdsg      = coeff*sech2(grid.hybrid_alpha*(sigma-grid.sigma_sigma));
            dFdsg      = zeta_1_0*(1.-g_sigma(sigma))+dgdsg*(THETA(K,J,I)-(grid.zeta0+sigma*zeta_1_0));
            W3(K,J,I) += dFdsg*sgdot;
          }
        }
      }
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  if (action == STEP_PROGS_AND_SYNC_DIAGS) {
    /*
     * Finish DZDT2.
     *
     * NOTE: The calculation of DZDT2 starts in store_pgrad_vars(), where the
     *       partial derivative of Z2 with respect to time is calculated and
     *       stored in DZDT2.
     */
    for (K = KLO; K <= KHI; K++) {
      kk = 2*K;
      for (J = JLO; J <= JHI; J++) {
        jj = 2*J;
        dy_inv2 = grid.n[kk][jj  ];
        dy_inv4 = grid.n[kk][jj+2];
        dx_inv3 = grid.m[kk][jj+1];
        for (I = ILO; I <= IHI; I++) {
          DZDT2(K,J,I) += .5*(W3(          K,  J,  I  )*(Z2(K,  J,  I  )-Z2(K+1,J,  I  ))*grid.dsgth_inv[kk+1]
                             +W3(          K-1,J,  I  )*(Z2(K-1,J,  I  )-Z2(K,  J,  I  ))*grid.dsgth_inv[kk-1]
                             +V(grid.it_uv,K,  J,  I  )*(Z2(K,  J,  I  )-Z2(K,  J-1,I  ))*dy_inv2
                             +V(grid.it_uv,K,  J+1,I  )*(Z2(K,  J+1,I  )-Z2(K,  J,  I  ))*dy_inv4
                             +U(grid.it_uv,K,  J,  I  )*(Z2(K,  J,  I  )-Z2(K,  J,  I-1))*dx_inv3
                             +U(grid.it_uv,K,  J,  I+1)*(Z2(K,  J,  I+1)-Z2(K,  J,  I  ))*dx_inv3             );
        }
      }
    }
    /* Need to call bc_lateral() here. */
    bc_lateral(var.dzdt2.value,THREEDIM);
  }

  return;
}

#undef SUM_DIV3
#undef HORIZONTAL_DIV

/*======================= end of calc_w() ==================================*/

/*======================= advection() ======================================*/

/*
 * Advect scalar prognostic variables forward one timestep using 
 * the specified scheme(s). Ensure positive-definite status as necessary.
 *
 * NOTE: Non-flux form schemes should lead with the 13 characters 
 *       "Non-flux form" in their advection_scheme string.
 */

void advection(planetspec  *planet,
               EPIC_FLOAT **Buff2D)

{
  register int
    is,ip,iq,
    K,J,I;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="advection";

  /*
   * THETA
   */
  if (var.theta.on) {
    if (strncmp(var.theta.advection_scheme,"Non-flux form",13) == 0) {
      switch(grid.coord_type) {
        case COORD_ISENTROPIC:
          ;
        break;
        case COORD_ISOBARIC:
        case COORD_HYBRID:
          if (strcmp(var.theta.advection_scheme,"Non-flux form, 3rd-order upwind") == 0) {
            upwind_3rd_order(planet,THETA_INDEX,NO_PHASE,var.theta.value,HORIZONTAL_AND_VERTICAL,Buff2D);
          }
          else {
            sprintf(Message,"unrecognized var.theta.advection_scheme=%s",var.theta.advection_scheme);
            epic_error(dbmsname,Message);
          }
          restore_mass(planet,THETA_INDEX,NO_PHASE);
        break;
        default:
          sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
          epic_error(dbmsname,Message);
        break;
      }
    }
    else {
      switch(grid.coord_type) {
        case COORD_ISENTROPIC:
          ;
        break;
        case COORD_ISOBARIC:
        case COORD_HYBRID:
          /* Weight with H3. */
          for (K = KLOPAD; K <= KHIPAD; K++) {
            for (J = JLOPAD; J <= JHIPAD; J++) {
              for (I = ILOPAD; I <= IHIPAD; I++) {
                THETA(K,J,I) *= H3(K,J,I);
              }
            }
          }
          /*
           * Advect H3*THETA.
           */ 
          if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.theta.advection_scheme) == 0) {
            hsu_advection(planet,THETA_INDEX,NO_PHASE,var.theta.value,HORIZONTAL_AND_VERTICAL,Buff2D);
          }
          else {
            sprintf(Message,"unrecognized var.theta.advection_scheme=%s",var.theta.advection_scheme);
            epic_error(dbmsname,Message);
          }
        break;
        default:
          sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
          epic_error(dbmsname,Message);
        break;
      }
    }
  }

  /*
   * FPARA
   */
  if (var.fpara.on) {
    if (strncmp(var.fpara.advection_scheme,"Non-flux form",13) == 0) {
      if (strcmp(var.fpara.advection_scheme,"Non-flux form, 3rd-order upwind") == 0) {
        upwind_3rd_order(planet,FPARA_INDEX,NO_PHASE,var.fpara.value,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else {
        sprintf(Message,"unrecognized var.fpara.advection_scheme=%s",var.fpara.advection_scheme);
        epic_error(dbmsname,Message);
      }
      restore_mass(planet,FPARA_INDEX,NO_PHASE);
    }
    else {
      /* Weight with H3. */
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            FPARA(K,J,I) *= H3(K,J,I);
          }
        }
      }
      if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.fpara.advection_scheme) == 0) {
        hsu_advection(planet,FPARA_INDEX,NO_PHASE,var.fpara.value,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else {
        sprintf(Message,"unrecognized var.fpara.advection_scheme=%s",var.fpara.advection_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  /*
   * Optional Species, Qs.
   */
  if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
    for (iq = 0; iq < grid.nq; iq++) {
      if (strncmp(var.species[grid.is[iq]].advection_scheme,"Non-flux form",13) == 0) {
        if (strcmp(var.species[grid.is[iq]].advection_scheme,"Non-flux form, 3rd-order upwind") == 0) {
          upwind_3rd_order(planet,grid.is[iq],grid.ip[iq],
                           var.species[grid.is[iq]].phase[grid.ip[iq]].q,HORIZONTAL_AND_VERTICAL,Buff2D);
        }
        else {
          sprintf(Message,"unrecognized var.species[%d].advection_scheme=%s",grid.is[iq],var.species[grid.is[iq]].advection_scheme);
          epic_error(dbmsname,Message);
        }
        restore_mass(planet,grid.is[iq],grid.ip[iq]);
      }
      else {
        /* Weight with H3. */
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              /* Ensure Q >= Q_MIN */
              Q(grid.is[iq],grid.ip[iq],K,J,I) = MAX(Q(grid.is[iq],grid.ip[iq],K,J,I),Q_MIN)*H3(K,J,I);
            }
          }
        }
        if (strcmp(var.species[grid.is[iq]].advection_scheme,"Predictor-corrector (Hsu, Konor, and Arakawa)") == 0) {
          hsu_advection(planet,grid.is[iq],grid.ip[iq],
                        var.species[grid.is[iq]].phase[grid.ip[iq]].q,HORIZONTAL_AND_VERTICAL,Buff2D);
        }
        else {
          sprintf(Message,"unrecognized var.species[%d].advection_scheme=%s",
                          grid.is[iq],var.species[grid.is[iq]].advection_scheme);
          epic_error(dbmsname,Message);
        }
      }
    }

    /*
     * NOTE: call sync_x_to_q() below, not here, because of possible h-weighting.
     */
  }

  /*
   * NU_TURB
   */
  if (var.nu_turb.on) {
    if (strncmp(var.nu_turb.advection_scheme,"Non-flux form",13) == 0) {
      if (strcmp(var.nu_turb.advection_scheme,"Non-flux form, 3rd-order upwind") == 0) {
        upwind_3rd_order(planet,NU_TURB_INDEX,NO_PHASE,var.nu_turb.value,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else {
        sprintf(Message,"unrecognized var.nu_turb.advection_scheme=%s",var.nu_turb.advection_scheme);
        epic_error(dbmsname,Message);
      }
      restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
    }
    else {
      /* Weight with H. */
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            NU_TURB(K,J,I) *= H(K,J,I);
          }
        }
      }
      if (strcmp("Predictor-corrector (Hsu, Konor, and Arakawa)",var.nu_turb.advection_scheme) == 0) {
        hsu_advection(planet,NU_TURB_INDEX,NO_PHASE,var.nu_turb.value,HORIZONTAL_AND_VERTICAL,Buff2D);
      }
      else {
        sprintf(Message,"unrecognized var.nu_turb.advection_scheme=%s",var.nu_turb.advection_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  /*
   * Advect H. 
   * Update H3.
   */
  if (var.h.on) {
    switch(grid.coord_type) {
      case COORD_ISOBARIC:
        ;
      break;
      case COORD_ISENTROPIC:
      case COORD_HYBRID:
        if (strcmp(var.h.advection_scheme,"Predictor-corrector (Hsu, Konor, and Arakawa)") == 0) {
          hsu_advection(planet,H_INDEX,NO_PHASE,var.h.value,HORIZONTAL_AND_VERTICAL,Buff2D);
        }
        else {
          sprintf(Message,"unrecognized var.h.advection_scheme=%s",var.h.advection_scheme);
          epic_error(dbmsname,Message);
        }
        restore_mass(planet,H_INDEX,NO_PHASE);

        /*
         * Update H3, which depends on H.
         */
        for (J = JLOPAD; J <=JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            /*
             * NOTE: If this is changed, it should also be changed in set_p2_etc().
             */
            for (K = KLO; K < KHI; K++) {
              H3(K,J,I) = sqrt(H(K,J,I)*H(K+1,J,I));
            }
            K = KLO-1;
            H3(K,J,I) = SQR(H(K+1,J,I))/H3(K+1,J,I);
            K = KHI;
            H3(K,J,I) = SQR(H(K,J,I))/H3(K-1,J,I);
            K = KHI+1;
            H3(K,J,I) = SQR(H3(K-1,J,I))/H3(K-2,J,I);
          }
        }
      break;
      default:
        sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
        epic_error(dbmsname,Message);
      break;
    }
  }

  /*
   * Restore variables to non-H-weighted form, as necessary.
   */
  if (var.theta.on                       && 
      grid.coord_type != COORD_ISENTROPIC  ) {
    if (strncmp(var.theta.advection_scheme,"Non-flux form",13) != 0) {
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            THETA(K,J,I) /= H3(K,J,I);
          }
        }
      }
    }
    restore_mass(planet,THETA_INDEX,NO_PHASE);
  }

  if (var.fpara.on) {
    if (strncmp(var.fpara.advection_scheme,"Non-flux form",13) != 0) {
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            FPARA(K,J,I) /= H3(K,J,I);
          }
        }
      }
    }
    restore_mass(planet,FPARA_INDEX,NO_PHASE);
  }

  if (grid.cloud_microphysics != OFF && grid.cloud_microphysics != STEADY) {
    for (iq = 0; iq < grid.nq; iq++) {
      if (strncmp(var.species[grid.is[iq]].advection_scheme,"Non-flux form",13) != 0) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              Q(grid.is[iq],grid.ip[iq],K,J,I) /= H3(K,J,I);
            }
          }
        }
      }
      restore_mass(planet,grid.is[iq],grid.ip[iq]);
    }
    sync_x_to_q(planet);
  }

  if (var.nu_turb.on) {
    if (strncmp(var.nu_turb.advection_scheme,"Non-flux form",13) != 0) {
      for (K = KLOPAD; K <= KHIPAD; K++) {
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            NU_TURB(K,J,I) /= H(K,J,I);
          }
        }
      }
    }
    restore_mass(planet,NU_TURB_INDEX,NO_PHASE);
  }

  set_p2_etc(planet,DONT_UPDATE_THETA,Buff2D);

  /*
   * NOTE: For the hybrid-coordinate case, a call to set_p2_etc() with UPDATE_THETA should be made to update
   *       theta = theta_diag, which is a function of pressure, 
   *       as soon as possible after returning from this subroutine.   
   */

  return;
}

/*======================= end of advection() ===============================*/

/*======================= hsu_advection() ==================================*/

/*
 *  Hsu and Arakawa's predictor-corrector advection, as modified 
 *  by eqs. (B.5) and (B.6) of Konor and Arakawa (1997).
 *  We apply it to each direction sequentially.
 *
 *  This is a positive-definite, flux-form advection scheme
 *  that handles steep gradients.
 */

void hsu_advection(planetspec  *planet,
                   int          is,
                   int          ip,
                   EPIC_FLOAT  *buff3d,
                   int          direction,
                   EPIC_FLOAT **Buff2D)
{
  int   
    K,J,I,
    kk,k_first;
  static int
    initialized = FALSE;
  unsigned long
    nbytes_2d;
  EPIC_FLOAT
    uhatp,uhatm,
    vhatp,vhatm,
    whatp,whatm,
    mu,vertvel,
    al,gap,gam,bep,bem,behatp,behatm,g_hsu,
    mn_2jp1,mn_2j,m_2j_inv,n_2jp1_inv,tmp,
    ep,a_min,
    precip_density,
    rho_h,
    da;
  EPIC_FLOAT
    *a,*a_diff1,*a_diff2,
    *um,*up,*u2d,*vm,*vp,*v2d,
    *ff,*a_pred,*aa1,*aaa;
  static EPIC_FLOAT
    *wm,*wp,
    *w_a,*w_a_pred,
    *w_a_min,*w_a_diff1,*w_a_diff2,
    *w_ff,*w_aaa;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="hsu_advection";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    wm        = fvector(0,KHI,  dbmsname);
    wp        = fvector(0,KHI,  dbmsname);
    w_a       = fvector(0,KHI,  dbmsname);
    w_a_pred  = fvector(0,KHI,  dbmsname);
    w_a_min   = fvector(0,KHI,  dbmsname);
    w_ff      = fvector(0,KHI+1,dbmsname);
    w_a_diff1 = fvector(0,KHI,  dbmsname);
    w_a_diff2 = fvector(0,KHI,  dbmsname);
    w_aaa     = fvector(0,KHI,  dbmsname);

    initialized = TRUE;
  }

  /*
   * Small constant ep (epsilon) introduced in Konor and Arakawa's 
   * eqs. (B.5,B.6) as a modification to Hsu and Arakawa's eqs. (6.15) and (6.16):
   */
  ep = 1.e-10;

  nbytes_2d = Nelem2d*sizeof(EPIC_FLOAT);

  switch(is) {
    case H_INDEX:
      if (!var.h.on) return;
      k_first = KLO;
    break;
    case THETA_INDEX:
      if (!var.theta.on) return;
      k_first = KLO;
    break;
    case NU_TURB_INDEX:
      if (!var.nu_turb.on) return;
      k_first = KLO;
    break;
    case FPARA_INDEX:
      if (!var.fpara.on) return;
      k_first = KLO;
    break;
    default:
      if (is < FIRST_SPECIES || is > LAST_SPECIES) {
        sprintf(Message,"case is=%d not recognized",is);
        epic_error(dbmsname,Message);
      }
      if (!var.species[is].phase[ip].on) return;
      k_first = KLO;
    break;
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_VERTICAL) {

    /***********************
     * Vertical advection. *
     ***********************/

    switch(is) {
      case H_INDEX:
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = grid.h_min[K];
        }
      break;
      case THETA_INDEX:
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = 0.;
        }
      break;
      case NU_TURB_INDEX:
        for (K = k_first; K <= KHI; K++) {
          w_a_min[K] = NU_TURB_MIN;
        }
      break;
      case FPARA_INDEX:
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = Q_MIN;
        }
      break;
      default:
        for (K = k_first-1; K <= KHI; K++) {
          w_a_min[K] = Q_MIN;
        }
      break;
    }

    for (J = JLOPAD; J <= JHIPAD; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        for (K = k_first-1; K <= KHI; K++) {
          w_a[K] = BUFF3D(K,J,I);
        }
        if (is == H_INDEX || is == NU_TURB_INDEX) {
          /*
           * Variable is carried in layer.
           */
          for (K = k_first-1; K <= KHI; K++) {
            vertvel = W3(K,J,I);
            wm[K]   = .5*(vertvel-fabs(vertvel));
            wp[K]   = .5*(vertvel+fabs(vertvel));
          }
        }
        else {
          /*
           * Variable is carried on layer interfaces.
           */
          if (is >= FIRST_SPECIES && is <= LAST_SPECIES &&
              ip >= FIRST_PRECIP  && ip <= LAST_PRECIP    ) {
            for (K = k_first; K <= KHI; K++) {
              vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));
              /*
               * Include terminal velocity for precipitation.
               */
              precip_density  = sqrt(Q(is,ip,K,  J,I)*PDRY3(K,  J,I)/T3(K,  J,I)
                                    *Q(is,ip,K-1,J,I)*PDRY3(K-1,J,I)/T3(K-1,J,I))/planet->rgas;
              rho_h           = RHO2(K,J,I)/H(K,J,I);
              vertvel        -= rho_h*fabs(terminal_velocity(is,ip,P2(K,J,I),T2(K,J,I),precip_density));

              wm[K] = .5*(vertvel-fabs(vertvel));
              wp[K] = .5*(vertvel+fabs(vertvel));
            }
          }
          else {
            for (K = k_first; K <= KHI; K++) {
              vertvel = .5*(W3(K,J,I)+W3(K-1,J,I));

              wm[K] = .5*(vertvel-fabs(vertvel));
              wp[K] = .5*(vertvel+fabs(vertvel));
            }
          }
        }

        /* 
         * Predictor: 
         */
        if (is == H_INDEX || is == NU_TURB_INDEX) {
          /*
           *  Here, w_a[K] is located above wp[K],wm[K].
           */
          if (k_first == 1) {
            w_ff[k_first-1] = flux_top(planet,is,J,I,NULL);
          }
          else {
            w_ff[k_first-1] = 0.;
          }
          w_ff[KHI] = flux_bot(planet,is,J,I,NULL);
          for (K = k_first; K < KHI; K++) {
            w_ff[K] = wp[K]*w_a[K+1]+wm[K]*w_a[K];
          }
          for (K = k_first; K <= KHI; K++) {
            da          = (w_ff[K]-w_ff[K-1])*grid.dsgth_inv[2*K]*DT;
            w_a_pred[K] = w_a[K]+da;
          }
        }
        else {
          /*
           * Here, w_a[K] is located below wp[K],wm[K].
           */
          for (K = k_first; K <= KHI; K++) {
            w_ff[K] = wp[K]*w_a[K]+wm[K]*w_a[K-1];
          }
          w_ff[KHI+1] = 0.;
          w_a_pred[k_first-1] = w_a[k_first-1];
          for (K = k_first; K < KHI; K++) {
            da          = (w_ff[K+1]-w_ff[K])*grid.dsgth_inv[2*K+1]*DT;
            w_a_pred[K] = w_a[K]+da;
          }
        }

        /* 
         * Corrector: 
         */
        if (is == H_INDEX || is == NU_TURB_INDEX) {
          for (K = k_first; K < KHI; K++) {
            w_ff[K] = .5*(wp[K]*(w_a_pred[K  ]+w_a[K+1])
                         +wm[K]*(w_a_pred[K+1]+w_a[K  ]));
          }
          for (K = k_first; K <= KHI; K++) {
            if (K == KLO) {
              w_a_diff1[K] = w_a_pred[K]-w_a[K+1];
              w_a_diff2[K] = 0.;
              tmp          = fabs(w_a[K+1]-w_a[K])+ep;
              w_aaa[K]     = tmp*tmp;
            }
            else if (K == KHI) {
              w_a_diff1[K] = 0.;
              w_a_diff2[K] = w_a_pred[K]-w_a[K-1];
              tmp          = fabs(-w_a[K]+w_a[K-1])+ep;
              w_aaa[K]     = tmp*tmp;
            }
            else {
              w_a_diff1[K] = w_a_pred[K]-w_a[K+1];
              w_a_diff2[K] = w_a_pred[K]-w_a[K-1];
              tmp          = fabs(w_a[K+1]-2.*w_a[K]+w_a[K-1])+ep;
              w_aaa[K]     = tmp*tmp;
            }
          }

          for (K = k_first; K < KHI; K++) {
            whatp =  sqrt(wp[K]*wp[K+1]);
            whatm = -sqrt(wm[K]*wm[K-1]);
            /* Courant number, mu */
            mu     = (wp[K]-wm[K])*grid.dsgth_inv[2*K+1]*DT;
            al     = (1.+mu)/6.;
            gap    = w_aaa[K+1]/(w_aaa[K+1]+(w_a[K+1]-w_a_min[K+1])
                                           *(w_a[K  ]-w_a_min[K  ]));
            gap    = gap*gap;
            gam    = w_aaa[K  ]/(w_aaa[K  ]+(w_a[K+1]-w_a_min[K+1])
                                           *(w_a[K  ]-w_a_min[K  ]));
            gam    = gam*gam;
            bep    = 1.+(1./(2.*al)-1.)*gap;
            bem    = 1.+(1./(2.*al)-1.)*gam;
            behatp = 1.-gap;
            behatm = 1.-gam;
            g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                         -whatp*behatp*w_a_diff1[K+1]
                            +wm[K]*bem*w_a_diff2[K+1]
                         -whatm*behatm*w_a_diff2[K  ]);
            w_ff[K] += g_hsu;
          }
        }
        else {
          for (K = k_first; K <= KHI; K++) {
            w_ff[K] = .5*(wp[K]*(w_a_pred[K-1]+w_a[K  ])
                         +wm[K]*(w_a_pred[K  ]+w_a[K-1]));
          }

          for (K = k_first-1; K <= KHI; K++) {
            if (K == 0) {
              tmp      = fabs(w_a[K+1]-w_a[K])+ep;
              w_aaa[K] = tmp*tmp;
            }
            else if (K == KHI) {
              tmp      = fabs(-w_a[K]+w_a[K-1])+ep;
              w_aaa[K] = tmp*tmp;
            }
            else {
              tmp      = fabs(w_a[K+1]-2.*w_a[K]+w_a[K-1])+ep;
              w_aaa[K] = tmp*tmp;
            }
          }

          w_a_diff2[k_first-1] = 0.;
          for (K = k_first; K <= KHI; K++) {
            w_a_diff1[K] = w_a_pred[K-1]-w_a[K  ];
            w_a_diff2[K] = w_a_pred[K  ]-w_a[K-1];
          }

          for (K = k_first; K <= KHI; K++) {
            if (K == KLO) {
              whatm = 0.;
            }
            else {
              whatm = -sqrt(wm[K]*wm[K-1]);
            }
            if (K == KHI) {
              whatp = 0.;
            }
            else {
              whatp = sqrt(wp[K]*wp[K+1]);
            }
            /* Courant number, mu */
            mu     = (wp[K]-wm[K])*grid.dsgth_inv[2*K]*DT;
            al     = (1.+mu)/6.;
            gap    = w_aaa[K  ]/(w_aaa[K  ]+(w_a[K  ]-w_a_min[K  ])
                                           *(w_a[K-1]-w_a_min[K-1]));
            gap    = gap*gap;
            gam    = w_aaa[K-1]/(w_aaa[K-1]+(w_a[K  ]-w_a_min[K  ])
                                           *(w_a[K-1]-w_a_min[K-1]));
            gam    = gam*gam;
            bep    = 1.+(1./(2.*al)-1.)*gap;
            bem    = 1.+(1./(2.*al)-1.)*gam;
            behatp = 1.-gap;
            behatm = 1.-gam;
            if (K < KHI) {
              g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                           -whatp*behatp*w_a_diff1[K+1]
                              +wm[K]*bem*w_a_diff2[K  ]
                           -whatm*behatm*w_a_diff2[K-1]);
            }
            else {
              g_hsu  = -al*(   wp[K]*bep*w_a_diff1[K  ]
                              +wm[K]*bem*w_a_diff2[K  ]
                           -whatm*behatm*w_a_diff2[K-1]);
            }
            w_ff[K] += g_hsu;
          }
        }

        switch(is) {
          case H_INDEX:
          case NU_TURB_INDEX:
            for (K = k_first; K <= KHI; K++) {
              da             = (w_ff[K]-w_ff[K-1])*grid.dsgth_inv[2*K]*DT;
              BUFF3D(K,J,I) += da;
            }
          break;
          default:
            for (K = k_first; K < KHI; K++) {
              da             = (w_ff[K+1]-w_ff[K])*grid.dsgth_inv[2*K+1]*DT;
              BUFF3D(K,J,I) += da;
            }
          break;
        }
      } /* I loop */
    } /* J loop */
    /* No need to apply bc_lateral() here. */
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_HORIZONTAL) {

    /*************************
     * Meridional advection. *
     *************************/

    for (K = k_first; K <= KHI; K++) {
      switch(is) {
        case H_INDEX:
          if (!var.h.on) return;
          a_min = grid.h_min[K];
        break;
        case THETA_INDEX:
          if (!var.theta.on) return;
          a_min = 0.;
        break;
        case NU_TURB_INDEX:
          if (!var.nu_turb.on) return;
          a_min = NU_TURB_MIN;
        break;
        case FPARA_INDEX:
          if (!var.fpara.on) return;
          a_min = Q_MIN;
        break;
        default:
          if (!var.species[is].phase[ip].on) return;
          a_min = Q_MIN;
        break;
      }

      a = buff3d+(K-Kshift)*Nelem2d;

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);

      vm  = Buff2D[0];
      vp  = Buff2D[1];
      v2d = Buff2D[2];

      if (is == H_INDEX || is == NU_TURB_INDEX) {
        kk = 2*K;
        for (J = JLOPAD; J <= JHIPADPV; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            V2D(J,I) = V(grid.it_uv,K,J,I);
          }
        }
        /* No need to apply bc_lateral() here. */
      }
      else {
        kk = 2*K+1;
        for (J = JFIRST; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            V2D(J,I) = get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(v2d,TWODIM);
      }

      /*
       * Take range of J for VM and VP to be JLOPAD, JHIPADPV, thereby
       * eliminating the need to apply bc_lateral(). It is not necessary
       * to treat the I index in this manner for VM and VP.
       */
      for (J = JLOPAD; J <= JHIPADPV; J++) {
        if (fcmp(grid.m[kk][2*J],0.) == 0) {
          /* 
           * The corresponding V should be zero in this case, 
           * such that this value does not matter.
           */
          m_2j_inv = FLOAT_MAX;
        }
        else{
          m_2j_inv = 1./grid.m[kk][2*J];
        }
        for (I = ILO; I <= IHI; I++) {
          VM(J,I) = MIN(0.,V2D(J,I)*m_2j_inv);
          VP(J,I) = MAX(0.,V2D(J,I)*m_2j_inv);
        }
      }
      /* 
       * Do not need to apply bc_lateral() to VM and VP. 
       * Done with v2d memory.
       */

      /* 
       * Predictor: 
       */
      memset(Buff2D[4],0,nbytes_2d);
      ff = Buff2D[4];
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = VP(J,I)*A(J-1,I)+VM(J,I)*A(J,I);
        }
      }
      bc_lateral(ff,TWODIM);

      memset(Buff2D[5],0,nbytes_2d);
      a_pred = Buff2D[5];
      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = grid.mn[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da          = (FF(J,I)-FF(J+1,I))*mn_2jp1*DT;
          A_PRED(J,I) = A(J,I)+da;
        }
      }
      bc_lateral(a_pred,TWODIM);

      /*
       * Corrector:
       */
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = .5*(VP(J,I)*(A_PRED(J,I  )+A(J-1,I))
                       +VM(J,I)*(A_PRED(J-1,I)+A(J,I  )));
        }
      }
      /* No need to apply bc_lateral() here. */

      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);
      a_diff1 = Buff2D[2];
      a_diff2 = Buff2D[3];

      /*
       * NOTE: A_DIFF1, A_DIFF2 are on the v-grid.
       */
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          A_DIFF1(J,I) =  A_PRED(J,I)-A(J-1,I);
          A_DIFF2(J,I) = -A(J,I)+A_PRED(J-1,I);
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(a_diff1,TWODIM);
      bc_lateral(a_diff2,TWODIM);

      /* Finished with a_pred memory. */
      memset(Buff2D[5],0,nbytes_2d);
      aaa = Buff2D[5];
      memset(Buff2D[6],0,nbytes_2d);
      aa1 = Buff2D[6];

      /*
       * AAA is a curvature, so calculate in two sequential passes, including bc_lateral() calls.
       * This is how we properly take a second derivative with single-thickness pads.
       *
       * AA1 is on the v-grid, whereas AAA is on the h-grid.
       */
      for (J = JFIRST; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          AA1(J,I) = A(J,I)-A(J-1,I);
        }
      }
      bc_lateral(aa1,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          tmp      = fabs(AA1(J+1,I)-AA1(J,I))+ep;
          AAA(J,I) = tmp*tmp;
        }
      }
      bc_lateral(aaa,TWODIM);

      for (J = JFIRST; J <= JHI; J++) {
        mn_2j = grid.mn[kk][2*J];
        for (I = ILO; I <= IHI; I++) {
          vhatp  =  sqrt(VP(J,I)*VP(J-1,I));
          vhatm  = -sqrt(VM(J,I)*VM(J+1,I));
          /* Courant number, mu */
          mu     = (VP(J,I)-VM(J,I))*mn_2j*DT;
          al     = (1.+mu)/6.;
          gap    = AAA(J-1,I)/(AAA(J-1,I)+(A(J-1,I)-a_min)*(A(J,I)-a_min));
          gap    = gap*gap;
          gam    = AAA(J,I  )/(AAA(J,I  )+(A(J-1,I)-a_min)*(A(J,I)-a_min));
          gam    = gam*gam;

          bep    = 1.+(1./(2.*al)-1.)*gap;
          bem    = 1.+(1./(2.*al)-1.)*gam;
          behatp = 1.-gap;
          behatm = 1.-gam;

          g_hsu  = -al*( VP(J,I)*bep*A_DIFF1(J,  I)
                       -vhatp*behatp*A_DIFF1(J-1,I)
                        +VM(J,I)*bem*A_DIFF2(J,  I)
                       -vhatm*behatm*A_DIFF2(J+1,I));
          FF(J,I) += g_hsu;
        }  
      }
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = grid.mn[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da      = (FF(J,I)-FF(J+1,I))*mn_2jp1*DT;
          A(J,I) += da;
        }
      }
      bc_lateral(a,TWODIM);
    }

    /********************
     * Zonal advection. *
     ********************/

    for (K = k_first; K <= KHI; K++) {
      switch(is) {
        case H_INDEX:
          if (!var.h.on) return;
          a_min = grid.h_min[K];
        break;
        case THETA_INDEX:
          if (!var.theta.on) return;
          a_min = 0.;
        break;
        case NU_TURB_INDEX:
          if (!var.nu_turb.on) return;
          a_min = NU_TURB_MIN;
        break;
        case FPARA_INDEX:
          if (!var.fpara.on) return;
          a_min = Q_MIN;
        break;
        default:
          if (!var.species[is].phase[ip].on) return;
          a_min = Q_MIN;
        break;
      }

      /* Zero buffer memory. */
      memset(Buff2D[0],0,nbytes_2d);
      memset(Buff2D[1],0,nbytes_2d);
      memset(Buff2D[2],0,nbytes_2d);
      um  = Buff2D[0];
      up  = Buff2D[1];
      u2d = Buff2D[2];

      a = Buff2D[6];
      memcpy(a,buff3d+(K-Kshift)*Nelem2d,Nelem2d*sizeof(EPIC_FLOAT));

      if (is == H_INDEX || is == NU_TURB_INDEX) {
        kk = 2*K;
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            U2D(J,I) = U(grid.it_uv,K,J,I);
          }
        }
        /* Do not need to call bc_lateral() here. */
      }
      else {
        kk = 2*K+1;
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            U2D(J,I) = get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(u2d,TWODIM);
      }

      /*
       * Take range of I for UM and UP to be ILOPAD, IHIPAD, thereby
       * eliminating the need to apply bc_lateral(). It is not necessary
       * to treat the J index in this manner for UM and UP.
       */
      for (J = JLO; J <= JHI; J++) {
        n_2jp1_inv = 1./grid.n[kk][2*J+1];
        for (I = ILOPAD; I <= IHIPAD; I++) {
          UM(J,I) = MIN(0.,U2D(J,I)*n_2jp1_inv);
          UP(J,I) = MAX(0.,U2D(J,I)*n_2jp1_inv);
        }
      }
      /* No need to apply bc_lateral(). */

      /*
       * Done with u2d memory. 
       */

      /* 
       * Predictor: 
       */
      memset(Buff2D[4],0,nbytes_2d);
      ff = Buff2D[4];
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = UP(J,I)*A(J,I-1)+UM(J,I)*A(J,I);
        }
      }
      bc_lateral(ff,TWODIM);

      memset(Buff2D[5],0,nbytes_2d);
      a_pred = Buff2D[5];
      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = grid.mn[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da          = (FF(J,I)-FF(J,I+1))*mn_2jp1*DT;
          A_PRED(J,I) = A(J,I)+da;
        }
      }
      bc_lateral(a_pred,TWODIM);

      /* 
       * Corrector: 
       *
       * NOTE: A_DIFF1, A_DIFF2 are on the u-grid.
       */
      memset(Buff2D[2],0,nbytes_2d);
      memset(Buff2D[3],0,nbytes_2d);
      a_diff1 = Buff2D[2];
      a_diff2 = Buff2D[3];
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          FF(J,I) = .5*(UP(J,I)*(A_PRED(J,I  )+A(J,I-1))
                       +UM(J,I)*(A_PRED(J,I-1)+A(J,I  )));
          A_DIFF1(J,I) =  A_PRED(J,I)-A(J,I-1);
          A_DIFF2(J,I) = -A(J,I)+A_PRED(J,I-1);
        }
      }
      bc_lateral(a_diff1,TWODIM);
      bc_lateral(a_diff2,TWODIM);

      /* Finished with a_pred memory. */

      memset(Buff2D[5],0,nbytes_2d);
      aaa = Buff2D[5];
      memset(Buff2D[6],0,nbytes_2d);
      aa1 = Buff2D[7];

      /*
       * AAA is a curvature, so calculate in two sequential passes, including bc_lateral() calls.
       * This is how we properly take a second derivative with single-thickness pads.
       *
       * AA1 is on the u-grid, whereas AAA is on the h-grid.
       */
      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          AA1(J,I) = A(J,I)-A(J,I-1);
        }
      }
      bc_lateral(aa1,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          tmp      = fabs(AA1(J,I+1)-AA1(J,I))+ep;
          AAA(J,I) = tmp*tmp;
        }
      }
      bc_lateral(aaa,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = grid.mn[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          uhatp  =  sqrt(UP(J,I)*UP(J,I-1));
          uhatm  = -sqrt(UM(J,I)*UM(J,I+1));
          /* Courant number, mu */
          mu     = (UP(J,I)-UM(J,I))*mn_2jp1*DT;
          al     = (1.+mu)/6.;
          gap    = AAA(J,I-1)/(AAA(J,I-1)+(A(J,I-1)-a_min)*(A(J,I)-a_min));
          gap    = gap*gap;
          gam    = AAA(J,I  )/(AAA(J,I  )+(A(J,I-1)-a_min)*(A(J,I)-a_min));
          gam    = gam*gam;
          bep    = 1.+(1./(2.*al)-1.)*gap;
          bem    = 1.+(1./(2.*al)-1.)*gam;
          behatp = 1.-gap;
          behatm = 1.-gam;
          g_hsu  = -al*( UP(J,I)*bep*A_DIFF1(J,I  )
                       -uhatp*behatp*A_DIFF1(J,I-1)
                        +UM(J,I)*bem*A_DIFF2(J,I  )
                       -uhatm*behatm*A_DIFF2(J,I+1));
          FF(J,I) += g_hsu;
        }  
      }
      bc_lateral(ff,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        mn_2jp1 = grid.mn[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          da             = (FF(J,I)-FF(J,I+1))*mn_2jp1*DT;
          BUFF3D(K,J,I) += da;
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(buff3d,THREEDIM);
  }

  return;
}

/*======================= end of hsu_advection() ===========================*/

/*======================= upwind_3rd_order() ===============================*/

/*
 * See Tannehill, Anderson, Pletcher, 1997, Computational Fluid Mechanics
 *        and Heat Transfer, Taylor & Francis, p. 224.
 */

#undef  COEFP2
#define COEFP2(k,i) coefp2[i+(k-KLO)*3]
#undef  COEFM2
#define COEFM2(k,i) coefm2[i+(k-KLO)*3]

#undef  COEFP3
#define COEFP3(k,i) coefp3[i+(k-KLO)*3]
#undef  COEFM3
#define COEFM3(k,i) coefm3[i+(k-KLO)*3]

void upwind_3rd_order(planetspec  *planet,
                      int          is,
                      int          ip,
                      EPIC_FLOAT  *buff3d,
                      int          direction,
                      EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    k,klast,kk,kkshift;
  EPIC_FLOAT
    v,vp,vm,n,
    u,up,um,m,
    w,wp,wm,zp2,zp1,zm1,zm2,
    rho_h,
    precip_density;
  EPIC_FLOAT
   *aa1,
   *aaa,
   *buffji,
   *a,
   *u2d,
   *v2d;
  EPIC_FLOAT
    val0,tmp,
    buff[grid.nk+1];
  static EPIC_FLOAT
   *coefp3,*coefp2,
   *coefm3,*coefm2;
  static int
    j_periodic,
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="upwind_3rd_order";

  if (!initialized) {
    /* Allocate memory */
    coefp2 = fvector(0,3*(KHI-KLO+1)-1,dbmsname);
    coefm2 = fvector(0,3*(KHI-KLO+1)-1,dbmsname);

    coefp3 = fvector(0,3*(KHI-1-KLO+1)-1,dbmsname);
    coefm3 = fvector(0,3*(KHI-1-KLO+1)-1,dbmsname);

    if (strcmp(grid.geometry,"f-plane") == 0 &&
        strcmp(grid.f_plane_map,"cartesian") == 0) {
      j_periodic = TRUE;
    }
    else {
      j_periodic = FALSE;
    }

    /*
     * The vertical direction requires general coefficients, because the vertical grid
     * does not have uniform spacing. Used Mathematica to help derive the formulas
     * for the coefficients (Dowling April 2011).
     */

    /* 
     * Layer values
     */
    K           = KLO;
    kk          = 2*K;
    zp1         = (grid.sigmatheta[kk-1]-grid.sigmatheta[kk]); /* Fettered step */
    zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zm2         = (grid.sigmatheta[kk+4]-grid.sigmatheta[kk]);
    COEFP2(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
    COEFP2(K,1) = zp1*zm2/(zm1*(zm1-zm2)*(zm1-zp1));
    COEFP2(K,2) = zp1*zm1/(zm2*(zm2-zm1)*(zm2-zp1));
    for (K = KLO+1; K <= KHI-2; K++) {
      kk          = 2*K;
      zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      zm2         = (grid.sigmatheta[kk+4]-grid.sigmatheta[kk]);
      COEFP2(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
      COEFP2(K,1) = zp1*zm2/(zm1*(zm1-zm2)*(zm1-zp1));
      COEFP2(K,2) = zp1*zm1/(zm2*(zm2-zm1)*(zm2-zp1));
    }
    K           = KHI-1;
    kk          = 2*K;
    zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zm2         = 2.*zm1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFP2(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
    COEFP2(K,1) = zp1*(zm2*(zp1-zm2)+zm1*(zm1-zp1))/((zm1-zm2)*zm1*(zm1-zp1)*(zp1-zm2));

    K           = KLO+1;
    kk          = 2*K;
    zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zp2         = 2.*zp1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFM2(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
    COEFM2(K,1) = zm1*(zp2*(zm1-zp2)+zp1*(zp1-zm1))/((zp1-zp2)*zp1*(zp1-zm1)*(zm1-zp2));
    for (K = KLO+2; K <= KHI-1; K++) {
      kk          = 2*K;
      zp2         = (grid.sigmatheta[kk-4]-grid.sigmatheta[kk]);
      zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      COEFM2(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
      COEFM2(K,1) = zm1*zp2/(zp1*(zp1-zm1)*(zp1-zp2));
      COEFM2(K,2) = zm1*zp1/(zp2*(zp2-zm1)*(zp2-zp1));
    }
    K           = KHI;
    kk          = 2*K;
    zp2         = (grid.sigmatheta[kk-4]-grid.sigmatheta[kk]);
    zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1         = (grid.sigmatheta[kk+1]-grid.sigmatheta[kk]); /* Fettered step */
    COEFM2(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
    COEFM2(K,1) = zm1*zp2/(zp1*(zp1-zm1)*(zp1-zp2));
    COEFM2(K,2) = zm1*zp1/(zp2*(zp2-zm1)*(zp2-zp1));

    /* 
     * Interface values
     */
    for (K = KLO; K <= KHI-2; K++) {
      kk          = 2*K+1;
      zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      zm2         = (grid.sigmatheta[kk+4]-grid.sigmatheta[kk]);
      COEFP3(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
      COEFP3(K,1) = zp1*zm2/(zm1*(zm1-zm2)*(zm1-zp1));
      COEFP3(K,2) = zp1*zm1/(zm2*(zm2-zm1)*(zm2-zp1));
    }
    K           = KHI-1;
    kk          = 2*K+1;
    zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zm2         = 2.*zm1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFP3(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
    COEFP3(K,1) = zp1*(zm2*(zp1-zm2)+zm1*(zm1-zp1))/((zm1-zm2)*zm1*(zm1-zp1)*(zp1-zm2));

    for (K = KLO+1; K <= KHI-1; K++) {
      kk          = 2*K+1;
      zp2         = (grid.sigmatheta[kk-4]-grid.sigmatheta[kk]);
      zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      COEFM3(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
      COEFM3(K,1) = zm1*zp2/(zp1*(zp1-zm1)*(zp1-zp2));
      COEFM3(K,2) = zm1*zp1/(zp2*(zp2-zm1)*(zp2-zp1));
    }
    K           = KLO;
    kk          = 2*K+1;
    zp1         = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1         = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zp2         = 2.*zp1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFM3(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
    COEFM3(K,1) = zm1*(zp2*(zm1-zp2)+zp1*(zp1-zm1))/((zp1-zp2)*zp1*(zp1-zm1)*(zm1-zp2));

    initialized = TRUE;
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_VERTICAL) {

    /***********************
     * Vertical advection. *
     ***********************/
    switch(is) {
      case H_INDEX:
      case NU_TURB_INDEX:
        /*
         * Layer variable
         */
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            for (K = KLO+1; K <= KHI-1; K++) {
              w       = .5*(W3(K,J,I)+W3(K-1,J,I));
              wp      = .5*(w+fabs(w));
              wm      = .5*(w-fabs(w));
              val0    = BUFF3D(K,J,I);
              buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP2(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP2(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP2(K,2))
                            +wm*((BUFF3D(K+1,J,I)-val0)*COEFM2(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM2(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM2(K,2)));
            }
            /*
             * Assume variable is linear on ends.
             * Use a 1st-order difference on the fettered side.
             */
            K       = KLO;
            w       = .5*(W3(K,J,I)+W3(K-1,J,I));
            wp      = .5*(w+fabs(w));
            wm      = .5*(w-fabs(w));
            val0    = BUFF3D(K,J,I);
            tmp     = -.5*(BUFF3D(K+1,J,I)-val0);
            buff[K] = -DT*(wp*(tmp*COEFP2(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP2(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP2(K,2))
                          +wm*tmp/(grid.sigmatheta[2*K-1]-grid.sigmatheta[2*K]));

            K       = KHI;
            w       = .5*(W3(K,J,I)+W3(K-1,J,I));
            wp      = .5*(w+fabs(w));
            wm      = .5*(w-fabs(w));
            val0    = BUFF3D(K,J,I);
            tmp     = -.5*(BUFF3D(K-1,J,I)-val0);
            buff[K] = -DT*(wp*tmp/(grid.sigmatheta[2*K]-grid.sigmatheta[2*K+1])
                          +wm*(tmp*COEFM2(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM2(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM2(K,2)));

            for (K = KLO; K <= KHI; K++) {
              BUFF3D(K,J,I) += buff[K];
            }
          }
        }
        /* No need to apply bc_lateral() here. */
      break;
      default:
        /*
         * Interface variable
         */
        if (is >= FIRST_SPECIES && is <= LAST_SPECIES &&
            ip >= FIRST_PRECIP  && ip <= LAST_PRECIP    ) {
          /*
           * Include terminal velocity for precipitation.
           */
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              for (K = KLO+1; K <= KHI-2; K++) {
                /*
                 * Use rho_h to convert terminal_velocity() from [m/s] to W3-format [K/s].
                 */
                rho_h          = RHO3(K,J,I)/H3(K,J,I);
                precip_density = Q(is,ip,K,J,I)*PDRY3(K,J,I)/(planet->rgas*T3(K,J,I));;
                w              = W3(K,J,I)-rho_h*fabs(terminal_velocity(is,ip,P3(K,J,I),T3(K,J,I),precip_density));

                wp      = .5*(w+fabs(w));
                wm      = .5*(w-fabs(w));
                val0    =  BUFF3D(K,J,I);
                buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP3(K,2))
                              +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM3(K,2)));
              }
              K = KLO;
              rho_h          = RHO3(K,J,I)/H3(K,J,I);
              precip_density = Q(is,ip,K,J,I)*PDRY3(K,J,I)/(planet->rgas*T3(K,J,I));;
              w              = W3(K,J,I)-rho_h*fabs(terminal_velocity(is,ip,P3(K,J,I),T3(K,J,I),precip_density));

              wp      = .5*(w+fabs(w));
              wm      = .5*(w-fabs(w));
              val0    =  BUFF3D(K,J,I);
              buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP3(K,2))
                            +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)));
              K = KHI-1;
              rho_h          = RHO3(K,J,I)/H3(K,J,I);
              precip_density = Q(is,ip,K,J,I)*PDRY3(K,J,I)/(planet->rgas*T3(K,J,I));;
              w              = W3(K,J,I)-rho_h*fabs(terminal_velocity(is,ip,P3(K,J,I),T3(K,J,I),precip_density));

              wp      = .5*(w+fabs(w));
              wm      = .5*(w-fabs(w));
              val0    =  BUFF3D(K,J,I);
              buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1))
                            +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM3(K,2)));

              for (K = KLO; K < KHI; K++) {
                BUFF3D(K,J,I) += buff[K];
              }
            }
          }
          /* No need to apply bc_lateral() here. */
        }
        else {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              for (K = KLO+1; K <= KHI-2; K++) {
                w       =  W3(K,J,I);
                wp      = .5*(w+fabs(w));
                wm      = .5*(w-fabs(w));
                val0    =  BUFF3D(K,J,I);
                buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP3(K,2))
                              +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM3(K,2)));
              }
              K = KLO;
              w       =  W3(K,J,I);
              wp      = .5*(w+fabs(w));
              wm      = .5*(w-fabs(w));
              val0    =  BUFF3D(K,J,I);
              buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1)+(BUFF3D(K+2,J,I)-val0)*COEFP3(K,2))
                            +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)));
              K = KHI-1;
              w       =  W3(K,J,I);
              wp      = .5*(w+fabs(w));
              wm      = .5*(w-fabs(w));
              val0    =  BUFF3D(K,J,I);
              buff[K] = -DT*(wp*((BUFF3D(K-1,J,I)-val0)*COEFP3(K,0)+(BUFF3D(K+1,J,I)-val0)*COEFP3(K,1))
                            +wm*((BUFF3D(K+1,J,I)-val0)*COEFM3(K,0)+(BUFF3D(K-1,J,I)-val0)*COEFM3(K,1)+(BUFF3D(K-2,J,I)-val0)*COEFM3(K,2)));

              for (K = KLO; K < KHI; K++) {
                BUFF3D(K,J,I) += buff[K];
              }
            }
          }
          /* No need to apply bc_lateral() here. */
        }
      break;
    }
  }

  if (direction == HORIZONTAL_AND_VERTICAL || direction == JUST_HORIZONTAL) {

    aa1    = Buff2D[0];
    aaa    = Buff2D[1];
    buffji = Buff2D[2];
    a      = Buff2D[3];
    u2d    = Buff2D[4];
    v2d    = Buff2D[5];

    /* Zero V2D memory */
    memset(v2d,0,sizeof(EPIC_FLOAT)*Nelem2d);

    switch(is) {
      case H_INDEX:
      case NU_TURB_INDEX:
        klast   = KHI;
        kkshift = 0;
      break;
      default:
        klast   = KHI-1;
        kkshift = 1;
      break;
    }

    /*************************
     * Meridional advection. *
     *************************/

    if (j_periodic) {
      for (K = KLO; K <= klast; K++) {
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            AA1(J,I) = BUFF3D(K,J,I)-BUFF3D(K,J-1,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(aa1,TWODIM);

        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            AAA(J,I) = AA1(J+1,I)-AA1(J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(aaa,TWODIM);

        if (is == H_INDEX || is == NU_TURB_INDEX) {
          kk = 2*K;
          for (J = JLOPAD; J <= JHIPADPV; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              V2D(J,I) = V(grid.it_uv,K,J,I);
            }
          }
          /* No need to apply bc_lateral() here. */
        }
        else {
          kk = 2*K+1;
          for (J = JFIRST; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              V2D(J,I) = get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
            }
          }
          /* Need to apply bc_lateral() here. */
          bc_lateral(v2d,TWODIM);
        }

        for (J = JLO; J <= JHI; J++) {
          n = grid.n[kk][2*J+1];
          for (I = ILO; I <= IHI; I++) {
            v  = .5*(V2D(J,I)+V2D(J+1,I));
            vp = .5*(v+fabs(v));
            vm = .5*(v-fabs(v));

            BUFFJI(J,I) = -DT*n*( v*(BUFF3D(K,J+1,I)-BUFF3D(K,J-1,I))/2.
                                 -( vp*(AAA(J,  I)-AAA(J-1,I))
                                   +vm*(AAA(J+1,I)-AAA(J,  I)) )/6. );
          }
        }
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            BUFF3D(K,J,I) += BUFFJI(J,I);
          }
        }
      }
    }
    else {
      for (K = KLO; K <= klast; K++) {
        /* Clear memory to clear boundaries. */
        memset(aa1,0,Nelem2d*sizeof(EPIC_FLOAT));
        memset(aaa,0,Nelem2d*sizeof(EPIC_FLOAT));

        for (J = JFIRST; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            AA1(J,I) = BUFF3D(K,J,I)-BUFF3D(K,J-1,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(aa1,TWODIM);

        for (I = ILO; I <= IHI; I++) {
          for (J = JLO; J <= JHI; J++) {
            AAA(J,I) = AA1(J+1,I)-AA1(J,I);
          }
          /*
           * Assume zero curvature at northern and southern edges.
           */
          if (JLO == grid.jlo) {
            AAA(JLO-1,I) = AAA(JLO,I) = 0.;
          }
          if (JHI == grid.nj) {
            AAA(JHI+1,I) = AAA(JHI,I) = 0.;
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(aaa,TWODIM);

        if (is == H_INDEX || is == NU_TURB_INDEX) {
          kk = 2*K;
          for (J = JLOPAD; J <= JHIPADPV; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              V2D(J,I) = V(grid.it_uv,K,J,I);
            }
          }
          /* No need to apply bc_lateral() here. */
        }
        else {
          kk = 2*K+1;
          for (J = JFIRST; J <= JHI; J++) {
            for (I = ILO; I <= IHI; I++) {
              V2D(J,I) = get_var(planet,V_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
            }
          }
          /* Need to apply bc_lateral() here. */
          bc_lateral(v2d,TWODIM);
        }

        for (J = JLO; J <= JHI; J++) {
          if (J == grid.nj) {
            /* Northern edge */
            n = grid.n[kk][2*J];
            for (I = ILO; I <= IHI; I++) {
              v  = .5*(V2D(J,I)+V2D(J+1,I));
              vp = .5*(v+fabs(v));
              vm = .5*(v-fabs(v));

              BUFFJI(J,I) = -DT*n*( v*(BUFF3D(K,J,I)-BUFF3D(K,J-1,I))/1.
                                   -( vp*(AAA(J,  I)-AAA(J-1,I))
                                     +vm*(AAA(J+1,I)-AAA(J,  I)) )/6. );            }
          }
          else if (J == grid.jlo) {
            /* Southern edge */
            n = grid.n[kk][2*(J+1)];
            for (I = ILO; I <= IHI; I++) {
              v  = .5*(V2D(J,I)+V2D(J+1,I));
              vp = .5*(v+fabs(v));
              vm = .5*(v-fabs(v));

              BUFFJI(J,I) = -DT*n*( v*(BUFF3D(K,J+1,I)-BUFF3D(K,J,I))/1.
                                   -( vp*(AAA(J,  I)-AAA(J-1,I))
                                     +vm*(AAA(J+1,I)-AAA(J,  I)) )/6. );
            }
          }
          else {
            /* Interior point */
            n = grid.n[kk][2*J+1];
            for (I = ILO; I <= IHI; I++) {
              v  = .5*(V2D(J,I)+V2D(J+1,I));
              vp = .5*(v+fabs(v));
              vm = .5*(v-fabs(v));

              BUFFJI(J,I) = -DT*n*( v*(BUFF3D(K,J+1,I)-BUFF3D(K,J-1,I))/2.
                                   -( vp*(AAA(J,  I)-AAA(J-1,I))
                                     +vm*(AAA(J+1,I)-AAA(J,  I)) )/6. );
            }
          }
        }
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            BUFF3D(K,J,I) += BUFFJI(J,I);
          }
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(buff3d,THREEDIM);

    /********************
     * Zonal advection. *
     ********************/

    for (K = KLO; K <= klast; K++) {
      memcpy(a,buff3d+(K-Kshift)*Nelem2d,Nelem2d*sizeof(EPIC_FLOAT));

      if (is == H_INDEX || is == NU_TURB_INDEX) {
        kk = 2*K;
        for (J = JLOPAD; J <= JHIPAD; J++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            U2D(J,I) = U(grid.it_uv,K,J,I);
          }
        }
        /* Do not need to call bc_lateral() here. */
      }
      else {
        kk = 2*K+1;
        for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
            U2D(J,I) = get_var(planet,U_INDEX,NO_PHASE,grid.it_uv,kk,J,I);
          }
        }
        /* Need to apply bc_lateral() here. */
        bc_lateral(u2d,TWODIM);
      }

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          AA1(J,I) = A(J,I)-A(J,I-1);
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(aa1,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        for (I = ILO; I <= IHI; I++) {
          AAA(J,I) = AA1(J,I+1)-AA1(J,I);
        }
      }
      /* Need to apply bc_lateral() here. */
      bc_lateral(aaa,TWODIM);

      for (J = JLO; J <= JHI; J++) {
        m = grid.m[kk][2*J+1];
        for (I = ILO; I <= IHI; I++) {
          u  = .5*(U2D(J,I)+U2D(J,I+1));
          up = .5*(u+fabs(u));
          um = .5*(u-fabs(u));

          BUFF3D(K,J,I) += -DT*m*(  u*(A(  J,I+1)-A(  J,I-1))/2.
                                -( up*(AAA(J,I  )-AAA(J,I-1))
                                 +um*(AAA(J,I+1)-AAA(J,I)) )/6. );
        }
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(buff3d,THREEDIM);
  }

  return;
}

/*======================= end of upwind_3rd_order() ========================*/

/*======================= uv_vertical_advection() ==========================*/

/*
 * Wrapper function that allows us to change easily between algorithms.
 *
 * As background, note that horizontal advection of U and V is accomplished by
 * writing the momentum equations in the vector-invariant form with (zeta+f)u
 * and (zeta+f)v, and using the Arakawa C-grid.  Meanwhile, horizontal and vertical
 * advection of the scalar variables (mass variables like H, and also turbulence
 * variables like NU_TURB) are done with a positive-definite scheme. Because U and V
 * are signed quantities, they need a different approach.
 */
void uv_vertical_advection(planetspec *planet)
{
  uv_vert_upwind_3rd_order(planet);

  return;
}

/*======================= end of uv_vertical_advection() ===================*/

/*======================= uv_vert_upwind_3rd_order() =======================*/

#undef  COEFP
#define COEFP(k,i) coefp[i+(k-KLO)*3]
#undef  COEFM
#define COEFM(k,i) coefm[i+(k-KLO)*3]

void uv_vert_upwind_3rd_order(planetspec *planet)
{
  register int
    K,J,I,kk;
  EPIC_FLOAT
    w,wp,wm,
    zp2,zp1,zm1,zm2,
    val0;
  static EPIC_FLOAT
   *coefp,
   *coefm;
  static int
    initialized = FALSE;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_vert_upwind_3rd_order";

  if (!initialized) {
    /* Allocate memory */
    coefp  = fvector(0,3*(grid.nk)-1,dbmsname);
    coefm  = fvector(0,3*(grid.nk)-1,dbmsname);

    /*
     * The vertical direction requires general coefficients, because the vertical grid
     * does not have uniform spacing. Used Mathematica to help derive the formulas
     * for the coefficients (Dowling April 2011).
     */
    K          = KLO;
    kk         = 2*K;
    zp1        = (grid.sigmatheta[kk-1]-grid.sigmatheta[kk]); /* Fettered step */
    zm1        = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zm2        = (grid.sigmatheta[kk+4]-grid.sigmatheta[kk]);
    COEFP(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
    COEFP(K,1) = zp1*zm2/(zm1*(zm1-zm2)*(zm1-zp1));
    COEFP(K,2) = zp1*zm1/(zm2*(zm2-zm1)*(zm2-zp1));
    for (K = KLO+1; K <= KHI-2; K++) {
      kk         = 2*K;
      zp1        = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1        = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      zm2        = (grid.sigmatheta[kk+4]-grid.sigmatheta[kk]);
      COEFP(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
      COEFP(K,1) = zp1*zm2/(zm1*(zm1-zm2)*(zm1-zp1));
      COEFP(K,2) = zp1*zm1/(zm2*(zm2-zm1)*(zm2-zp1));
    }
    K          = KHI-1;
    kk         = 2*K;
    zp1        = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1        = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zm2        = 2.*zm1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFP(K,0) = zm1*zm2/(zp1*(zp1-zm2)*(zp1-zm1));
    COEFP(K,1) = zp1*(zm2*(zp1-zm2)+zm1*(zm1-zp1))/((zm1-zm2)*zm1*(zm1-zp1)*(zp1-zm2));

    K          = KLO+1;
    kk         = 2*K;
    zp1        = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1        = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
    zp2        = 2.*zp1;
    /* Assume 3-pt zero curvature for endpoint. */
    COEFM(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
    COEFM(K,1) = zm1*(zp2*(zm1-zp2)+zp1*(zp1-zm1))/((zp1-zp2)*zp1*(zp1-zm1)*(zm1-zp2));
    for (K = KLO+2; K <= KHI-1; K++) {
      kk         = 2*K;
      zp2        = (grid.sigmatheta[kk-4]-grid.sigmatheta[kk]);
      zp1        = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
      zm1        = (grid.sigmatheta[kk+2]-grid.sigmatheta[kk]);
      COEFM(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
      COEFM(K,1) = zm1*zp2/(zp1*(zp1-zm1)*(zp1-zp2));
      COEFM(K,2) = zm1*zp1/(zp2*(zp2-zm1)*(zp2-zp1));
    }
    K          = KHI;
    kk         = 2*K;
    zp2        = (grid.sigmatheta[kk-4]-grid.sigmatheta[kk]);
    zp1        = (grid.sigmatheta[kk-2]-grid.sigmatheta[kk]);
    zm1        = (grid.sigmatheta[kk+1]-grid.sigmatheta[kk]); /* Fettered step */
    COEFM(K,0) = zp1*zp2/(zm1*(zm1-zp1)*(zm1-zp2));
    COEFM(K,1) = zm1*zp2/(zp1*(zp1-zm1)*(zp1-zp2));
    COEFM(K,2) = zm1*zp1/(zp2*(zp2-zm1)*(zp2-zp1));

    initialized = TRUE;
  }

  for (I = ILO; I <= IHI; I++) {
    /*
     * U
     */
    for (J = JLO; J <= JHI; J++) {
      for (K = KLO+1; K <= KHI-1; K++) {
        w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J,I-1)+W3(K-1,J,I-1));
        wp   = .5*(w+fabs(w));
        wm   = .5*(w-fabs(w));
        val0 = U(grid.it_uv,K,J,I);

        DUDT(grid.it_uv_tend,K,J,I) -= 
               wp*((U(grid.it_uv,K-1,J,I)-val0)*COEFP(K,0)+(U(grid.it_uv,K+1,J,I)-val0)*COEFP(K,1)+(U(grid.it_uv,K+2,J,I)-val0)*COEFP(K,2))
              +wm*((U(grid.it_uv,K+1,J,I)-val0)*COEFM(K,0)+(U(grid.it_uv,K-1,J,I)-val0)*COEFM(K,1)+(U(grid.it_uv,K-2,J,I)-val0)*COEFM(K,2));
      }
      /*
       * Assume U = 0. at top of model. At the bottom, assume dU/dsgth = 0 in abyssal layer for gas giants, and U = 0 for terrestrial planets.
       * Use a 1st-order difference on the fettered side.
       */
      K    = KLO;
      w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J,I-1)+W3(K-1,J,I-1));
      wp   = .5*(w+fabs(w));
      wm   = .5*(w-fabs(w));
      val0 = U(grid.it_uv,K,J,I);
      DUDT(grid.it_uv_tend,K,J,I) -= 
             wp*((0.-val0)*COEFP(K,0)+(U(grid.it_uv,K+1,J,I)-val0)*COEFP(K,1)+(U(grid.it_uv,K+2,J,I)-val0)*COEFP(K,2))
            +wm*(0.-val0)/(grid.sigmatheta[2*K-1]-grid.sigmatheta[2*K]);

      K    = KHI;
      w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J,I-1)+W3(K-1,J,I-1));
      wp   = .5*(w+fabs(w));
      wm   = .5*(w-fabs(w));
      val0 = U(grid.it_uv,K,J,I);
      if (strcmp(planet->type,"gas-giant") == 0 || grid.coord_type == COORD_ISENTROPIC) {
        DUDT(grid.it_uv_tend,K,J,I) -=
               wp*0.
              +wm*((U(grid.it_uv,K+1,J,I)-val0)*COEFM(K,0)+(U(grid.it_uv,K-1,J,I)-val0)*COEFM(K,1)+(U(grid.it_uv,K-2,J,I)-val0)*COEFM(K,2));
      }
      else if (strcmp(planet->type,"terrestrial") == 0) {
        DUDT(grid.it_uv_tend,K,J,I) -= 
               wp*(val0-0.)/(grid.sigmatheta[2*K]-grid.sigmatheta[2*K+1])
              +wm*((0.-val0)*COEFM(K,0)+(U(grid.it_uv,K-1,J,I)-val0)*COEFM(K,1)+(U(grid.it_uv,K-2,J,I)-val0)*COEFM(K,2));
      }
      else {
        sprintf(Message,"unrecognized planet->type=%s",planet->type);
        epic_error(dbmsname,Message);
      }
    }

    /*
     * V
     */
    for (J = JFIRST; J <= JHI; J++) {
      for (K = KLO+1; K <= KHI-1; K++) {
        w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J-1,I)+W3(K-1,J-1,I));
        wp   = .5*(w+fabs(w));
        wm   = .5*(w-fabs(w));
        val0 = V(grid.it_uv,K,J,I);

        DVDT(grid.it_uv_tend,K,J,I) -= 
               wp*((V(grid.it_uv,K-1,J,I)-val0)*COEFP(K,0)+(V(grid.it_uv,K+1,J,I)-val0)*COEFP(K,1)+(V(grid.it_uv,K+2,J,I)-val0)*COEFP(K,2))
              +wm*((V(grid.it_uv,K+1,J,I)-val0)*COEFM(K,0)+(V(grid.it_uv,K-1,J,I)-val0)*COEFM(K,1)+(V(grid.it_uv,K-2,J,I)-val0)*COEFM(K,2));
      }
      /*
       * Assume V = 0. at top and bottom of model.
       * Use a 1st-order difference on the fettered side.
       */
      K    = KLO;
      w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J-1,I)+W3(K-1,J-1,I));
      wp   = .5*(w+fabs(w));
      wm   = .5*(w-fabs(w));
      val0 = V(grid.it_uv,K,J,I);
      DVDT(grid.it_uv_tend,K,J,I) -= 
             wp*((0.-val0)*COEFP(K,0)+(V(grid.it_uv,K+1,J,I)-val0)*COEFP(K,1)+(V(grid.it_uv,K+2,J,I)-val0)*COEFP(K,2))
            +wm*(0.-val0)/(grid.sigmatheta[2*K-1]-grid.sigmatheta[2*K]);

      K    = KHI;
      w    = .25*(W3(K,J,I)+W3(K-1,J,I)+W3(K,J-1,I)+W3(K-1,J-1,I));
      wp   = .5*(w+fabs(w));
      wm   = .5*(w-fabs(w));
      val0 = V(grid.it_uv,K,J,I);
      DVDT(grid.it_uv_tend,K,J,I) -= 
             wp*(val0-0.)/(grid.sigmatheta[2*K]-grid.sigmatheta[2*K+1])
            +wm*((0.-val0)*COEFM(K,0)+(V(grid.it_uv,K-1,J,I)-val0)*COEFM(K,1)+(V(grid.it_uv,K-2,J,I)-val0)*COEFM(K,2));
    }
  }
  /* No need to call bc_lateral() here */

  return;
}

/*======================= end of uv_vert_upwind_3rd_order() ================*/

/* * * * * * * * * * * * *  end of epic_flux.c  * * * * * * * * * * * * * * */
















