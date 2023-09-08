/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2002-2019 Csaba Palotai *A*                       *
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

/* * * * * * * * * *  epic_microphysics.c  * * * * * * * * * * * * *
 *                         v.2.1                                   *        
 *                                                                 *
 *       Functions governing the hydrological cycle.               *
 *       This file includes the following functions:               *
 *                                                                 *
 *           cloud_microphysics()                                  *
 *           instantaneous_processes()                             *
 *           finite_rate_processes()                               *
 *           terminal_velocity()                                   *
 *           restore_mass_min()                                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>


/*========================== cloud_microphysics() ============================*/

  /*
   * Add latent heating to HEAT and transfer mass between the phases of each species,
   * as appropriate.
   */

void cloud_microphysics(planetspec  *planet,
                        EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    is,ip;
  EPIC_FLOAT
    rh_ji,*rh;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0,
    non_precip = FALSE,    /* see line 100, 103 */
    no_heating = FALSE;
 static int
    warned_once = FALSE;
  static char
    dbmsname[]="cloud_microphysics";

#if EPIC_CHECK == 1 

    /*
     * Check that the phases are turned on for every species.
     */
    for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
      if (var.species[is].on) {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (!var.species[is].phase[ip].on) {
            sprintf(Message,"%s is not ON",var.species[is].phase[ip].info[0].name);
            epic_error(dbmsname,Message);
          }
        }
      }
    }

#endif 


  if (LAST_PHASE >= RAIN) {
    non_precip = FALSE;
  }

  rh = Buff2D[0];
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      for (K = KLO; K <= KHI; K++) {
        restore_mass_min(is,K);

	/*
         * Calculate relative humidity.
         */
        relative_humidity(planet,is,rh,K);

	for (J = JLO; J <= JHI; J++) {
          for (I = ILO; I <= IHI; I++) {
	    /* * * * * * * * * * * * * * * * * * * * * * * * * *
	     *   Instantaneous processes: melting, freezing    *
	     * * * * * * * * * * * * * * * * * * * * * * * * * */

            instantaneous_processes(is,K,J,I,non_precip,no_heating);

            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	     *   Finite-rate processes                                               *
	     *   Condensation/sublimation/evaporation, autoconversion, collection    *
	     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	    rh_ji = RH(J,I);
            finite_rate_processes(is,K,J,I,rh_ji,non_precip,no_heating);	    
          }
        }

        restore_mass_min(is,K);   

      }    /* end of K-loop */
      /*
       * Need to apply bc_lateral() here.
       */
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        bc_lateral(var.species[is].phase[ip].q,THREEDIM);
      }
    }
  }

  /*
   * Synchronize mole fractions, X, to mass mixing ratios, Q.
   */
  sync_x_to_q(planet);

  return;
}

/*======================== end of cloud_microphysics() ===========================*/

/*========================= instantaneous_processes() ============================*/

void instantaneous_processes(int    is,
                             int    K,
                             int    J,
                             int    I,
		             int    non_precip,
			     int    no_heating)
{
  double
    psmlti,psmlts;
  EPIC_FLOAT
    q_ice,q_liquid,q_snow,q_rain;  
  boolean
    warm;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="instantaneous_microphysics";

  psmlti   = 0.;
  psmlts   = 0.;
  q_ice    = Q(is,ICE,   K,J,I);
  q_liquid = Q(is,LIQUID,K,J,I);
  q_snow   = Q(is,SNOW,  K,J,I);
  q_rain   = Q(is,RAIN,  K,J,I);

  warm      = (fcmp(T3(K,J,I),T_triple_pt(is)) >= 0);
  if (warm && fcmp(q_ice,EPSILON) > 0) { 
    /*  Phase change: melting: ice cloud -> liquid cloud  */
    psmlti = q_ice/DT;
    if (fcmp(psmlti,EPSILON) != 1) {
      psmlti = 0.;
    }
  }
  else if (!warm && fcmp(q_liquid,EPSILON) > 0) {
    /* Phase change: freezing: liquid cloud -> ice cloud */
    psmlti = q_liquid/DT;
    if (fcmp(psmlti,EPSILON) != 1) {
      psmlti = 0.;
    }
    else {
      psmlti *= -1.;
    }  
  }

  if (!non_precip) {
    if (warm && fcmp(q_snow,EPSILON) > 0) {
      /*  Phase change: melting: snow -> rain  */
      psmlts = q_snow/DT;
      if (fcmp(psmlts,EPSILON) != 1) {
        psmlts = 0.;
      }
    }
    else if (!warm && fcmp(q_rain,EPSILON) > 0) {
      /* Phase change: freezing: rain -> snow */
      psmlts = q_rain/DT;
      if (fcmp(psmlts,EPSILON) != 1) {
        psmlts = 0.;
      }
      else {
        psmlts *= -1.;
      }	
    }
  }
  Q(is,LIQUID,K,J,I) += psmlti*DT;
  Q(is,   ICE,K,J,I) -= psmlti*DT;
  Q(is,  RAIN,K,J,I) += psmlts*DT;
  Q(is,  SNOW,K,J,I) -= psmlts*DT;

  if (!no_heating) {         
    HEAT3(K,J,I) -=(psmlti+psmlts)*Lf(is);
  }

  return;
}
/*==================== end of instantaneous_processes() =======================*/

/*========================= finite_rate_processes() ===========================*/
void finite_rate_processes(int        is,
                           int        K,
                           int        J,
                           int        I,
			   EPIC_FLOAT rh_ji,
			   int        non_precip,
			   int        no_heating)
{
  register EPIC_FLOAT
    tmp,rho,rho_inv,t3;
  static EPIC_FLOAT
    c_inv[LAST_SPECIES+1],d_inv[LAST_SPECIES+1],
    qi_crit[LAST_SPECIES+1];
  EPIC_FLOAT
    pcond,pint,pdepi,praut,psaut,prevap,psevap,
    psaci,psacw,pracw,psacr,
    required,d_mass,av_mass,
    subsat,supersat,av_ice,av_liquid,
    q_vapor,q_liquid,q_ice,q_rain,q_snow;           
  EPIC_FLOAT
    N_I0,Q_I0,Q_C0,E_c,Q_Icrit,E_SI,N_0S,N_I,dynvis,
    lambda_r,lambda_s,A_S,B_S,A_R,B_R,M_Imax;
  boolean
    saturated,warm;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  static int
    initialized = FALSE;
  int
temp,
    sp,
    idbms=0;
  static char
    dbmsname[]="finite_rate_microphysics";

  /*
   * Note: 3/25/2011 CJP
   * Variables subsat,supersat,av_ice,av_liquid are introduced so we wouldn't have to update the mixing ratios 
   * instantaneously after every finite-rate process, but we would still avoid using more than the available amount
   * of moisture by these processes.   
   *
   * These variables might get updated in the last process and never be used again during a function call. 
   * The reason for this is that should additional processes be added to the model, these moisture variables would give
   * the actual available amount and the user would not have to modify the old code.
   */
  
  if (!initialized) {
    for (sp = FIRST_SPECIES; sp <= LAST_SPECIES; sp++) {
      if (var.species[sp].on) {
        qi_crit[sp] = pow(COEFF_C(sp)*COEFF_M(sp)*pow(D_Icrit(is),COEFF_N(sp)),1./(1.-COEFF_D(sp))); 

	c_inv[sp]   = 1./COEFF_C(sp);
	d_inv[sp]   = 1./COEFF_D(sp);
      }
    }
    initialized = TRUE;
  }

  rho       = PDRY3(K,J,I)/(T3(K,J,I)*planet->rgas);
  rho_inv   = 1./rho;
  t3        = T3(K,J,I);
  q_vapor   = Q(is,VAPOR, K,J,I);
  q_liquid  = Q(is,LIQUID,K,J,I);
  q_ice     = Q(is,ICE,   K,J,I);
  q_rain    = Q(is,RAIN,  K,J,I);
  q_snow    = Q(is,SNOW,  K,J,I);
  pcond     = 0.;
  pint      = 0.;
  pdepi     = 0.;
  praut     = 0.;
  psaut     = 0.;
  pracw     = 0.;
  psaci     = 0.;
  prevap    = 0.;
  psevap    = 0.;
  av_ice    = q_ice; 
  av_liquid = q_liquid;

  saturated = (fcmp(rh_ji,1.) > 0);
  warm      = (fcmp(t3,T_triple_pt(is)) >= 0);
  
  supersat  = MAX(q_vapor*(1.-1./rh_ji),0.);
  subsat    = MAX(q_vapor*(1./rh_ji-1.),0.);

  /* PCOND liquid condensation/evaporation */
  if (warm) {    
    if (saturated) {             /* Phase change: condensation */
      pcond = supersat/(1.+ SQR(Lc(is,VAPOR,LIQUID,t3))*q_vapor/(rh_ji*planet->cp*GAS_R(is)*SQR(t3)))/DT;
      if (fcmp(pcond,EPSILON) <= 0) {
        pcond = 0.;
      }
      else { 
        supersat -= pcond*DT;  
      }
    }
    else if (!saturated && fcmp(av_liquid,EPSILON) > 0) {  
      /* Phase change: evaporation */
      pcond   = subsat/(1.+ SQR(Lc(is,VAPOR,LIQUID,t3))*q_vapor/(rh_ji*planet->cp*GAS_R(is)*SQR(t3)))/DT;
      if (fcmp(pcond,EPSILON) <= 0) {
        pcond = 0.;
      }
      else {
        pcond      = -MIN(pcond,av_liquid/DT);
        subsat    += pcond*DT;     /* pcond is now negative thus the + sign */
        av_liquid += pcond*DT; 
      }
    }
  }  /*====== End of PCOND ======*/


  /* PINIT + PDEPI: ice condensation/sublimation */
  A_S = Ls(is)*(Ls(is)/(GAS_R(is)*t3)-1.)/(conductivity(planet->name,t3)*t3);
  B_S = rh_ji/(mass_diffusivity(planet,is,T3(K,J,I),P3(K,J,I))*q_vapor*rho); 
  N_I = COEFF_C(is)*pow(rho*q_ice,COEFF_D(is)); 
  
  if (!warm) {     
    if (saturated) {   /* Phase change: vapor condensation to ice cloud */
      if (fcmp(q_ice,EPSILON) <= 0) {
        N_I0 = MIN(1.e+8,10000.*exp(0.1*(T_triple_pt(is)-t3)));
	Q_I0 = MIN(pow(c_inv[is]*N_I0,d_inv[is]),0.0015)*rho_inv;
	pint = MIN(Q_I0/DT,supersat/DT);

        if (fcmp(pint,EPSILON) <= 0) {
          pint = 0.;
        }
        else {
	  supersat -= pint*DT;
	}  
      }
      else if (fcmp(supersat,EPSILON) > 0) {    /* PDEPI: vapor deposition of a small ice crystal */
        tmp   = 1./COEFF_N(is);  
        pdepi = 4.0*pow(COEFF_M(is),-tmp)*(rh_ji-1.)*pow(rho*q_ice,tmp)*pow(N_I,1.-tmp)/(A_S+B_S);
        pdepi     = MIN(pdepi,supersat/DT);
 	if (fcmp(pdepi,EPSILON) <= 0) {
          pdepi = 0.;
        }
        else {
          supersat -= pdepi*DT;
        }
      }
    }
    else if (!saturated && fcmp(av_ice,EPSILON) > 0) {   
      /* Phase change: ice cloud evaporation */
      tmp   = 1./COEFF_N(is);  
      pdepi = 4.0*pow(COEFF_M(is),-tmp)*(1.-rh_ji)*pow(rho*q_ice,tmp)*pow(N_I,1.-tmp)/(A_S+B_S);
      pdepi = MIN(pdepi,av_ice/DT);
      if (fcmp(pdepi,EPSILON) <= 0) {
        pdepi = 0.;
      }
      subsat -= pdepi*DT;
      av_ice -= pdepi*DT; 
      pdepi   *= -1.;
    }
  }
  /*====== End of PINT + PDEPI ======*/
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * *
   *          Precipitation related processes          *
   * * * * * * * * * * * * * * * * * * * * * * * * * * */
  if (!non_precip) {         
    
    dynvis   = dynvisc(planet->name,t3);
    N_0S     = MIN(2.e+8,2.0e+6* exp(0.12*(T_triple_pt(is)-t3)));
    lambda_s = pow(M_PI*RHO_SNOW(is)*N_0S/(rho*q_snow),0.25);  
       
    /* PRAUT: Phase change: autoconversion: liquid cloud -> rain */
    if (fcmp(av_liquid,Q_LIQ_0(is)) > 0) {
      praut = MAX(ALPHA_RAUT(is)*(av_liquid-Q_LIQ_0(is)),0.); 
      if (fcmp(praut,EPSILON) <= 0) {
        praut = 0.;
      }
      else {
        av_liquid -= praut*DT;
      }
    }   /*============= End of PRAUT ============*/

    /* PSAUT: Phase change: autoconversion: ice crystals -> snow */
    Q_Icrit = qi_crit[is]*rho_inv; 
    if (fcmp(av_ice,Q_Icrit) > 0) {
      psaut   = MAX((av_ice-Q_Icrit)/DT,0.);   
      if (fcmp(psaut,EPSILON) <= 0) {
        psaut = 0.;
      }
      else {
        av_ice -= psaut*DT;
      }
    }   /*============= End of PSAUT ============*/

    /* PRACW: Accretion of cloud liquid by rain  */
    lambda_r = pow(M_PI*RHO_RAIN(is)*N_0R(is)/(rho*q_rain),0.25);  

    if (fcmp(q_rain,EPSILON) > 0 && fcmp(av_liquid,EPSILON) > 0) {
      pracw = M_PI*COEFF_XR(is)*q_liquid*E_R(is)*N_0R(is)*pow(P_REF/P3(K,J,I),P_EXP_LIQ(is))*gamma_nr(COEFF_YR(is)+3.)
             /(4.0*pow(lambda_r,COEFF_YR(is)+3.0));     
      pracw = MIN(av_liquid/DT,pracw);
       if (fcmp(pracw,EPSILON) <= 0) {
          pracw = 0.;
      }   /*============= End of PRACW ============*/
      else {
        av_liquid -= pracw*DT;
      }	    
    }
    /* PSACI: Accretion of cloud ice by snow */ 
    if (fcmp(q_snow,EPSILON) > 0 && fcmp(av_ice,EPSILON) > 0) {
      E_SI  = exp(0.05*(t3-T_triple_pt(is)));         
      psaci = M_PI*COEFF_XS(is)*av_ice*E_SI*N_0S*pow(P_REF/P3(K,J,I),P_EXP_ICE(is))*gamma_nr(COEFF_YS(is)+3.)
              /(4.0*pow(lambda_s,COEFF_YS(is)+3.0));
      psaci = MIN(av_ice/DT,psaci);
      if (fcmp(psaci,EPSILON) <= 0) {
        psaci = 0.;
      }
      else {
        av_ice -= psaci*DT;
      }	    
    }   /*============== End of PSACI ============*/

    /* Evaporation of precipitation takes place if gridbox is subsaturated after cloud evaporation */
    if (fcmp(q_rain,EPSILON) > 0 && warm) {
      A_R    = Lc(is,VAPOR,LIQUID,t3)*(Lc(is,VAPOR,LIQUID,t3)/(GAS_R(is)*t3)-1.)
             /(conductivity(planet->name,t3)*t3);
      B_R    = rh_ji/(mass_diffusivity(planet,is,T3(K,J,I),P3(K,J,I))*q_vapor*rho);
      tmp    = (f1r(is)/SQR(lambda_r))+f2r(is)*pow(SC(K,J,I),ONE_3)*sqrt(COEFF_XR(is)*rho/dynvis)*pow(P_REF/P3(K,J,I),P_EXP_LIQ(is)/2.)
               *gamma_nr((COEFF_YR(is)+5.)/2.)/pow(lambda_r,(COEFF_YR(is)+5.)/2.);
      if (fcmp(subsat,EPSILON) > 0) {         /* Evaporation of rain */ 
        prevap = 2.*M_PI*N_0R(is)*(1.-rh_ji)*tmp/(A_R+B_R);
	prevap = MIN(prevap,subsat/DT);
	prevap = MIN(prevap,q_rain/DT);
        if (fcmp(prevap, EPSILON) <= 0) {
          prevap = 0.;
        }
	else {
	  subsat -= prevap*DT;  /* reduce subsaturation by prevap */
	  prevap  = -prevap;
        }
      }
      else if (fcmp(supersat,EPSILON) > 0) {  /* Depositional growth of rain */ 
        prevap = 2*M_PI*N_0R(is)*(rh_ji-1.)*tmp/(A_R+B_R);
	prevap = MIN(prevap,supersat/DT);
        if (fcmp(prevap, EPSILON) <= 0) {
          prevap = 0.;
        }
        else {
          supersat -= prevap*DT;
	}
      }
    }   /*====== End of PREVAP ======*/
    else if (fcmp(q_snow,EPSILON) > 0 && !warm) { 
      A_S    = Ls(is)*rho*(Ls(is)-GAS_R(is)*t3)/(conductivity(planet->name,t3)*GAS_R(is)*SQR(t3));
      B_S    = rh_ji/(mass_diffusivity(planet,is,T3(K,J,I),P3(K,J,I))*q_vapor);
      tmp    = (f1s(is)/SQR(lambda_s))+f2s(is)*pow(SC(K,J,I),ONE_3)*sqrt(COEFF_XS(is)*rho/dynvis)*pow(P_REF/P3(K,J,I),P_EXP_ICE(is)/2.)
               *gamma_nr((COEFF_YS(is)+5.)/2.)/pow(lambda_s,(COEFF_YS(is)+5.)/2.);

      if (fcmp(subsat,EPSILON) > 0) {         /* Evaporation of snow */ 
        psevap  = 4.0*N_0S*(1.-rh_ji)*tmp/(A_S+B_S);
        if (fcmp(psevap, EPSILON) <= 0) {
          psevap = 0.;
        }
        else {
	  psevap  = MIN(psevap,subsat/DT);
	  psevap  = MIN(psevap,q_snow/DT);  
	  subsat -= psevap*DT;
	  psevap  = -psevap;
	} 
      }
      else if (fcmp(supersat,EPSILON) > 0) {  /* Depositional growth of snow */ 
        psevap = 4.0*N_0S*(rh_ji-1.)*tmp/(A_S+B_S);
        psevap = MIN(psevap,supersat/DT);
        if (fcmp(psevap, EPSILON) <= 0) {
          psevap = 0.;
        }
        else {
	  supersat -= psevap*DT;
        }
      }
    }  /*====== End of PSEVAP ======*/
  }    /* * * * * * * * * * * * * * * * * * * * * * * * * * *
        *      End of precipitation related processes       *
        * * * * * * * * * * * * * * * * * * * * * * * * * * */
	   
  Q(is, VAPOR,K,J,I) -= pcond*DT;
  Q(is,LIQUID,K,J,I) += pcond*DT;

  if (!no_heating) {         
    if (fcmp(fabs(pcond),EPSILON) > 0) {
      HEAT3(K,J,I) -= pcond*Lc(is,VAPOR,LIQUID,t3); 
    }
  }
  
  Q(is,VAPOR,K,J,I)  -= (pint+pdepi)*DT;
  Q(is,SOLID,K,J,I)  += (pint+pdepi)*DT; 

  Q(is,LIQUID,K,J,I) -= (praut+pracw)*DT;   /* autoconversion of rain */
  Q(is,  RAIN,K,J,I) += (praut+pracw)*DT;

  Q(is, ICE,K,J,I)   -= (psaut+psaci)*DT;   /* autoconversion of snow */
  Q(is,SNOW,K,J,I)   += (psaut+psaci)*DT;

  Q(is, RAIN,K,J,I)  += prevap*DT;
  Q(is,VAPOR,K,J,I)  -= prevap*DT;

  if (!no_heating) {         
    if (fcmp(fabs(prevap),EPSILON) > 0) {
      HEAT3(K,J,I) -= prevap*Lc(is,VAPOR,LIQUID,t3); 
    }
  }

  Q(is, SNOW,K,J,I)  += psevap*DT;   /* psevap is negative for evaporation */
  Q(is,VAPOR,K,J,I)  -= psevap*DT;

  if (!no_heating) {         
    HEAT3(K,J,I) += (pint+pdepi+psevap)*Ls(is);
  }
  return;
}
/*======================= end of finite_rate_processes() =========================*/

/*========================== terminal_velocity() =================================*/
/*
 * Cs. Palotai
 * Returns terminal velocity [m/s] as a positive number.
 *
 * Inputs:             is: species index
 *                     ip: phase index
 *               pressure: total pressure [Pa]
 *            temperature: physical temperature [K]
 *         precip_density: model precipitation density [kg/m^3]
 */
double terminal_velocity(int    is,
                         int    ip,
                         double pressure,
                         double temperature,
                         double precip_density)
{
  EPIC_FLOAT
    w,lambda_r,lambda_s,N_0S;
  int
    idbms=0;
  static char
    dbmsname[]="terminal_velocity";

  w = 0.;
  switch(ip) {
    case SNOW:
      N_0S     = MIN(2.e+8,2.e+6*exp(0.12*(T_triple_pt(is)-temperature)));
      lambda_s = pow(M_PI*RHO_SNOW(is)*N_0S/precip_density,0.25);
      w        = COEFF_XS(is)*gamma_nr(4.+COEFF_YS(is))*pow(lambda_s,-COEFF_YS(is))
                  *pow((P_REF/pressure),P_EXP_ICE(is))/6.;
    break;
    case RAIN:
      lambda_r = pow(M_PI*RHO_RAIN(is)*N_0R(is)/precip_density,0.25);
      w        = COEFF_XR(is)*gamma_nr(4.+COEFF_YR(is))*pow(lambda_r,-COEFF_YR(is))
                  *pow((P_REF/pressure),P_EXP_LIQ(is))/6.;
    break;
    default:
      sprintf(Message,"phase %s not yet implemented",var.species[is].phase[ip].info[0].name);
      epic_error(dbmsname,Message);
    break;
  }

  return w;
}

/*======================== end of terminal_velocity() ============================*/

/*=========================== restore_mass_min() =================================*/

void restore_mass_min(int  is,
                      int  K  )

{
  int
    ip,J,I;

  for (J = JLOPAD; J <= JHIPAD; J++) {
    for (I = ILOPAD; I <= IHIPAD; I++) {
      for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
        if (fcmp(Q(is,ip,K,J,I),EPSILON) <= 0) {
	  Q(is,ip,K,J,I) = EPSILON;
        }
      }
    }
  }

  return;
}

/*========================== end of restore_mass_min() ===========================*/

/* * * * * * * * * * * *  end of epic_microphysics.c  * * * * * * * * * * * * * * */
