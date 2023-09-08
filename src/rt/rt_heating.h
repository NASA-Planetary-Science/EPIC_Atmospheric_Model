/* * * * * * * rt_heating.h  * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2013  Thomas Greathouse, Timothy Dowling          *
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
 * Software Foundation, Inc., 59 Temple Place - Suite 330,         *
 * Boston, MA  02111-1307, USA.                                    *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef RT_HEATING_H
#define RT_HEATING_H
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Header file for rt_heating.c                                    *
 *                                                                 *
 * We try to follow the radiative-transfer notation used in        *
 *   A. Sanchez-Lavega 2011, An Introduction to Planetary          *
 *       Atmospheres, Taylor & Francis.                            *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Physical
 */
/* 
 * Loschmidt's number [molecules/m^3/amagat]
 *    Mohr PJ, Taylor BN, Newell DB, 2008, CODATA Recommended Values
 *    of the Fundamental Physical Constants: 2006", Rev. Mod. Phys.
 *    80 (2): 633--730, doi:10.1103/RevModPHYS.80.633
 */
#define N_STP 2.6867774e+25   

/*
 * Structures
 */

/* collision-induced absorption coeffs [1/(km amagat^2)] */
typedef struct {
  int
    on,
    nt;
  double
   *t,
   *h2_h2,
   *h2_he,
   *h2_ch4;
} rt_cia;

/* absorption coefficients [1/(km amagat)] */
typedef struct {
  int
    on,
    np,        /* number of pressures                                   */
    nt,        /* number of temperatures                                */
    ng;        /* number of g-ordinates                                 */
  double
   *t,         /* temperatures associated with data for an rt_band      */
   *p,         /* pressures                                             */
   *delg,      /* correlated-k statistical weights                      */
   *ch4,
   *c2h2,
   *c2h4,
   *c2h6,
   *nh3,
   *ph3;
} rt_k_nu;

/* scattering cross sections [cm^2] */
typedef struct {
  int
    on;
  double
   *ch4,
   *h2,
   *he;
} rt_sigma;

/* absorption cross sections [cm^2] */
typedef struct {
  int
    on;
  double
   *ch4,
   *c2h2,
   *c2h4,
   *c2h6;
} rt_k;

typedef struct {
  int
    on,
    nwn,       /* number of wavenumber bins spanning this band                          */
    planck,    /* true if thermal emission is turned on, false if not                   */
    scat_yes;  /* false implies single-scattering albedo, SSALB, is zero for all layers */
  double
    wn_start,  /* start of band's wavenumber interval [cm-1]                            */
    wn_end;    /* end of band's wavenumber interval [cm-1]                              */
  double
   *wn,        /* central wavenumber of bin [cm-1]                                      */
   *wnlo,      /* low wavenumber of bin [cm-1]                                          */
   *wnhi,      /* high wavenumber of bin [cm-1]                                         */
   *fluxtot,   /* total insolation in wavenumber bin                                    */
   *dtauc,     /* optical depth of computational layer                                  */
   *ssalb;     /* single scattering albedo of computational layer                       */
  rt_cia
   *cia;       /* collision-induced absorption coeffs [1/(km amagat^2)]                 */
  rt_k_nu
   *k_nu;      /* absorption coefficients [1/(km amagat)]                               */
  rt_k
   *k;         /* absorption cross sections [cm^2]                                      */
  rt_sigma
   *sigma;     /* scattering cross sections [cm^2]                                      */
} rt_band;

/*
 * Shift macros
 */
#define HEATING(k,j,i)                     heating[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/*
 * Shortwave bands
 */

/* UV */
#define UV_WN(iwn)                         uv->wn[iwn-1]
#define UV_WNLO(iwn)                       uv->wnlo[iwn-1]
#define UV_WNHI(iwn)                       uv->wnhi[iwn-1]
#define UV_K_CH4(iwn)                      uv->k->ch4[iwn-1]
#define UV_K_C2H2(iwn)                     uv->k->c2h2[iwn-1]
#define UV_K_C2H4(iwn)                     uv->k->c2h4[iwn-1]
#define UV_K_C2H6(iwn)                     uv->k->c2h6[iwn-1]
#define UV_SIGMA_H2(iwn)                   uv->sigma->h2[iwn-1]
#define UV_SIGMA_HE(iwn)                   uv->sigma->he[iwn-1]
#define UV_SIGMA_CH4(iwn)                  uv->sigma->ch4[iwn-1]
#define UV_FLUXTOT(iwn)                    uv->fluxtot[iwn-1]
#define UV_DTAUC(lc,iwn)                   uv->dtauc[lc-1+(iwn-1)*ds->nlyr]
#define UV_SSALB(lc,iwn)                   uv->ssalb[lc-1+(iwn-1)*ds->nlyr]

/* VIS */
#define VIS_WN(iwn)                        vis->wn[iwn-1]
#define VIS_WNLO(iwn)                      vis->wnlo[iwn-1]
#define VIS_WNHI(iwn)                      vis->wnhi[iwn-1]
#define VIS_K_NU_CH4(iwn)                  vis->k_nu->ch4[iwn-1]
#define VIS_CIA_T(itemp)                   vis->cia->t[itemp-1]
#define VIS_CIA_H2_H2(itemp,iwn)           vis->cia->h2_h2[itemp-1+(iwn-1)*vis->cia->nt]
#define VIS_SIGMA_H2(iwn)                  vis->sigma->h2[iwn-1]
#define VIS_SIGMA_HE(iwn)                  vis->sigma->he[iwn-1]
#define VIS_SIGMA_CH4(iwn)                 vis->sigma->ch4[iwn-1]
#define VIS_FLUXTOT(iwn)                   vis->fluxtot[iwn-1]
#define VIS_DTAUC(lc,iwn)                  vis->dtauc[lc-1+(iwn-1)*ds->nlyr]
#define VIS_SSALB(lc,iwn)                  vis->ssalb[lc-1+(iwn-1)*ds->nlyr]

/* NIR */
#define NIR_WN(iwn)                        nir->wn[iwn-1]
#define NIR_WNLO(iwn)                      nir->wnlo[iwn-1]
#define NIR_WNHI(iwn)                      nir->wnhi[iwn-1]
#define NIR_K_NU_T(itemp)                  nir->k_nu->t[itemp-1]
#define NIR_K_NU_P(ipress)                 nir->k_nu->p[ipress-1]
#define NIR_K_NU_DELG(ig)                  nir->k_nu->delg[ig-1]
#define NIR_K_NU_CH4(ig,ipress,itemp,iwn)  nir->k_nu->ch4[ig-1+(ipress-1+(itemp-1+(iwn-1)*nir->k_nu->nt)*nir->k_nu->np)*nir->k_nu->ng]
#define NIR_CIA_T(itemp)                   nir->cia->t[itemp-1]
#define NIR_CIA_H2_H2(itemp,iwn)           nir->cia->h2_h2[itemp-1+(iwn-1)*nir->cia->nt]
#define NIR_FLUXTOT(iwn)                   nir->fluxtot[iwn-1]
#define NIR_DTAUC(lc,ig,iwn)               nir->dtauc[lc-1+(ig-1+(iwn-1)*nir->k_nu->ng)*ds->nlyr]      

/*
 * Longwave bands
 */

/* MIDIR */
#define MIDIR_WN(iwn)                      midir->wn[iwn-1]
#define MIDIR_WNLO(iwn)                    midir->wnlo[iwn-1]
#define MIDIR_WNHI(iwn)                    midir->wnhi[iwn-1]
#define MIDIR_K_NU_T(itemp)                midir->k_nu->t[itemp-1]
#define MIDIR_K_NU_P(ipress)               midir->k_nu->p[ipress-1]
#define MIDIR_K_NU_CH4(ipress,itemp,iwn)   midir->k_nu->ch4[ ipress-1+(itemp-1+(iwn-1)*midir->k_nu->nt)*midir->k_nu->np]
#define MIDIR_K_NU_C2H2(ipress,itemp,iwn)  midir->k_nu->c2h2[ipress-1+(itemp-1+(iwn-1)*midir->k_nu->nt)*midir->k_nu->np]
#define MIDIR_K_NU_C2H6(ipress,itemp,iwn)  midir->k_nu->c2h6[ipress-1+(itemp-1+(iwn-1)*midir->k_nu->nt)*midir->k_nu->np]
#define MIDIR_K_NU_NH3(ipress,itemp,iwn)   midir->k_nu->nh3[ ipress-1+(itemp-1+(iwn-1)*midir->k_nu->nt)*midir->k_nu->np]
#define MIDIR_K_NU_PH3(ipress,itemp,iwn)   midir->k_nu->ph3[ ipress-1+(itemp-1+(iwn-1)*midir->k_nu->nt)*midir->k_nu->np]
#define MIDIR_CIA_T(itemp)                 midir->cia->t[itemp-1]
#define MIDIR_CIA_H2_H2(itemp,iwn)         midir->cia->h2_h2[ itemp-1+(iwn-1)*midir->cia->nt]
#define MIDIR_CIA_H2_HE(itemp,iwn)         midir->cia->h2_he[ itemp-1+(iwn-1)*midir->cia->nt]
#define MIDIR_CIA_H2_CH4(itemp,iwn)        midir->cia->h2_ch4[itemp-1+(iwn-1)*midir->cia->nt]
#define MIDIR_DTAUC(lc,iwn)                midir->dtauc[lc-1+(iwn-1)*ds->nlyr]

/* FIR */
#define FIR_WN(iwn)                        fir->wn[iwn-1]
#define FIR_WNLO(iwn)                      fir->wnlo[iwn-1]
#define FIR_WNHI(iwn)                      fir->wnhi[iwn-1]
#define FIR_CIA_T(itemp)                   fir->cia->t[itemp-1]
#define FIR_CIA_H2_H2(itemp,iwn)           fir->cia->h2_h2[ itemp-1+(iwn-1)*fir->cia->nt]
#define FIR_CIA_H2_HE(itemp,iwn)           fir->cia->h2_he[ itemp-1+(iwn-1)*fir->cia->nt]
#define FIR_CIA_H2_CH4(itemp,iwn)          fir->cia->h2_ch4[itemp-1+(iwn-1)*fir->cia->nt]
#define FIR_DTAUC(lc,iwn)                  fir->dtauc[lc-1+(iwn-1)*ds->nlyr]

/*
 * Function prototypes
 */

/*--------------------* 
 * Longwave functions *
 *--------------------*/

void rt_longwave(planetspec *planet,
                 EPIC_FLOAT *heating,
                 int         i_stride,
                 int         action);

void init_rt_longwave(rt_band      *fir,
                      rt_band      *midir,
                      disort_state *ds);

void free_rt_longwave(rt_band      *fir,
                      rt_band      *midir,
                      disort_state *ds);

void fir_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *fir,
                 disort_state *ds);

void midir_opacity(planetspec   *planet,
                   int           J,
                   int           I,
                   rt_band      *midir,
                   disort_state *ds);

/*---------------------*
 * Shortwave functions *
 *---------------------*/

void rt_shortwave(planetspec *planet,
                  EPIC_FLOAT *heating,
                  int         i_stride,
                  int         action);

void init_rt_shortwave(rt_band      *nir,
                       rt_band      *vis,
                       rt_band      *uv,
                       disort_state *ds);

void free_rt_shortwave(rt_band      *nir,
                       rt_band      *vis,
                       rt_band      *uv,
                       disort_state *ds);

void nir_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *nir,
                 disort_state *ds);

void vis_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *vis,
                 disort_state *ds);

void uv_opacity(planetspec   *planet,
                int           J,
                int           I,
                rt_band      *uv,
                disort_state *ds);

void cross_section_rayleigh(rt_band *band);

EPIC_FLOAT solar_irradiance(EPIC_FLOAT lambda);

void beer_law_only(disort_state  *ds,
                   disort_output *out,
                   double         radius);

/*---------------------* 
 * Planetary functions *
 *---------------------*/

double ring_shadow(planetspec *planet,
                   double      lat,
                   double      subsolar_lat);

/* * * * * * * * end of rt_heating.h * * * * * * * * * * * * * * * */
#endif
