/* * * * * * * rt_heating.c  * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 2013-2019 Thomas Greathouse, Timothy Dowling      *
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Radiative transfer code to drive seasonal forcing.              *
 *   Original Fortran by T. Greathouse, with inputs from           *
 *   J. Moses and others.                                          *
 *   Rewritten in C and adapted to EPIC by T. Dowling              *
 *                                                                 *
 * Unit-offset arrays are implemented using shift macros;          *
 * see the file rt_heating.h                                       *
 *                                                                 *
 *  Longwave functions:                                            *
 *    rt_longwave()                                                *
 *    init_rt_longwave()                                           *
 *    free_rt_longwave()                                           *
 *    fir_opacity()                                                *
 *    midir_opacity()                                              *
 *                                                                 *
 *  Shortwave functions:                                           *
 *    rt_shortwave()                                               *
 *    init_rt_shortwave()                                          *
 *    free_rt_shortwave()                                          *
 *    nir_opacity()                                                *
 *    vis_opacity()                                                *
 *    uv_opacity()                                                 *
 *    cross_section_rayleigh()                                     *
 *    solar_irradiance()                                           *
 *    beer_law_only()                                              *
 *                                                                 *
 *  Planetary functions:                                           *
 *    ring_shadow()                                                *
 *                                                                 *
 *  Vertical alignment of EPIC vs. DISORT [c_twostr()]:            *
 *                                                                 *
 *              EPIC                   DISORT                      *
 *                                                                 *
 *   ----kk = 1-------X,HEAT3----                                  *
 *        K = 1       T2         -----lev = 0-------               *
 *   -----------------X,HEAT3----     lc  = 1                      *
 *                    T2         -------------------               *
 *   -----------------X,HEAT3----     lc  = nlyr                   *
 *        K = nk      T2         -----lev = nlyr----               *
 *   ----kk = 2*nk+1--X,HEAT3----                                  *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cdisort.h>
#include <epic.h>
#include <rt_heating.h>

/*--------------------*
 * Longwave functions *
 *--------------------*/

/*============= rt_longwave() =====================================*/

/*
 * Partition of longwave wavelength regions:
 *
 *   FIR        0,    600 cm-1
 *   MIDIR    600,   1600 cm-1
 *
 * Input i_stride is used to subsample in the longitude direction,
 * to facilitate economical diurnal averaging.
 *
 * Fortran subroutine:
 *   heatingver12()
 */

void rt_longwave(planetspec *planet,
                 EPIC_FLOAT *heating,
                 int         i_stride,
                 int         action)
{
  int
    K,J,I,
    kk,jj,
    lc,lev,lu,ii,
    iwn;
  static int
    deltam,
    ierror[TWOSTR_NERR];
  double
    radius;
  EPIC_FLOAT
    x,dx;
  static double
   *gg,
   *tauc;
  static float_triplet
   *tau_table;
  static rt_band
    *fir,
    *midir;
  static disort_state
   *ds;
  static disort_output
   *out;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="rt_longwave";

  if (action == EPIC_APPLY) {
    /*
     * Zero longwave heating array.
     */
    memset(heating,0,Nelem3d*sizeof(EPIC_FLOAT));

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;

      /*
       * Start at grid.ilo so that i_stride is correct.
       */
      for (I = grid.ilo; I <= IHI; I += i_stride) {
        if (I < ILO) {
          continue;
        }

        /*
         * The Chapman function used to calculate spherical effects is
         * computationally expensive, so set ds->flag.spher to true only
         * where it is significant, as marked by cos(70.0*DEG) = 0.342.
         */
        ds->flag.spher = (ds->bc.umu0 < 0.342) ? TRUE : FALSE;

        if (ds->flag.spher == TRUE) {
          /*
           * The spherical correction needs the following radius information.
           */
          radius = grid.rlt[2*KHI+1][jj];
          for (K = KLO; K <= KHI; K++) {
            kk          = 2*K+1;
            lev         = K-1;
            ds->zd[lev] = grid.rlt[kk][jj]-radius;
          }
        }

        /*
         * Set independent variable for spline on tau data,
         * used to specify user-defined optical depth positions, UTAU.
         * Use log pressure.
         */
        for (lev = 0; lev <= ds->nlyr; lev++) {
          K                = lev+1;
          tau_table[lev].x = log(P2(K,J,I));
        }

        /*-----------------*
         * FIR: 0,600 cm-1 *
         *-----------------*/

        /*
         * SSALB refers to the single-scattering albedo.
         *
         * NOTE: We currently do not include single scattering for FIR.
         */
        fir->scat_yes = FALSE;
        memset(&SSALB(1),0,ds->nlyr*sizeof(double));

        fir_opacity(planet,J,I,fir,ds);

        /*
         * Loop over wavenumber bins.
         */
        for (iwn = 1; iwn <= fir->nwn; iwn++) {
          ds->wvnmhi = FIR_WNHI(iwn);
          ds->wvnmlo = FIR_WNLO(iwn);

          /*
           * ds->bc.fbeam is initialized to zero.
           * Single scattering albedo is initialized to zero.
           */

          /*
           * Set DISORT-layer optical depth, DTAUC.
           */
          memcpy(&DTAUC(1),&FIR_DTAUC(1,iwn),ds->nlyr*sizeof(double));

          /*
           * Construct TAUC(lev) in the same manner as is done
           * internally in DISORT, and then use it to specify user-defined
           * optical-depth positions, UTAU, for the output.
           *
           * NOTE: Inside DISORT, TAUC is vertically staggered with respect to DTAUC,
           *       hence it is aligned with EPIC layers.
           */
          TAUC(0) = 0.;
          for (lev = 1; lev <= ds->nlyr; lev++) {
            TAUC(lev) = TAUC(lev-1)+DTAUC(lev);
          }
 
          /*
           * EPIC needs heating on its interfaces, so spline TAUC
           * to locate UTAU at EPIC interfaces. Use log pressure as the independent variable.
           */
          for (lev = 0; lev <= ds->nlyr; lev++) {
            tau_table[lev].y = TAUC(lev);
          }
          spline_pchip(ds->nlyr+1,tau_table);

          ii = -2;
          for (lu = 1; lu <= ds->ntau; lu++) {
            K        = lu;
            x        = log(P3(K,J,I));
            ii       = hunt_place_in_table(ds->nlyr+1,tau_table,x,&dx,ii);
            UTAU(lu) = splint_pchip(x,tau_table+ii,dx);
          }

          /*
           * Calculate the radiative transfer.
           */
          if (!fir->planck && !fir->scat_yes) {
            beer_law_only(ds,out,radius);
          }
          else {
            c_twostr(ds,out,deltam,gg,ierror,radius);
          }

          /*
           * Sum flux divergence; convert to heating rate below.
           */
          for (K = KLO; K < KHI; K++) {
            lu              = K;
            HEATING(K,J,I) += DFDT(lu)*DTAUC(lu);
          }
        } /* iwn loop */

        /*----------------------*
         * MIDIR: 600,1600 cm-1 *
         *----------------------*/

        /*
         * SSALB refers to the single-scattering albedo.
         *
         * NOTE: We currently do not include single scattering for MIDIR.
         */
        midir->scat_yes = FALSE;
        memset(&SSALB(1),0,ds->nlyr*sizeof(double));

        midir_opacity(planet,J,I,midir,ds);

        /*
         * Loop over wavenumber bins.
         */
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          ds->wvnmhi = MIDIR_WNHI(iwn);
          ds->wvnmlo = MIDIR_WNLO(iwn);

          /*
           * ds->bc.fbeam is initialized to zero.
           * Single scattering albedo is initialized to zero.
           */

          /*
           * Set DISORT-layer optical depth, DTAUC.
           */
          memcpy(ds->dtauc,midir->dtauc+(iwn-1)*ds->nlyr,ds->nlyr*sizeof(double));

          /*
           * Construct TAUC(lev) in the same manner as is done
           * internally in DISORT, and then use it to specify user-defined
           * optical-depth positions, UTAU, for the output.
           *
           * NOTE: Inside DISORT, TAUC is vertically staggered with respect to DTAUC,
           *       hence it is aligned with EPIC layers.
           */
          TAUC(0) = 0.;
          for (lev = 1; lev <= ds->nlyr; lev++) {
            TAUC(lev) = TAUC(lev-1)+DTAUC(lev);
          }
 
          /*
           * EPIC needs heating on its interfaces, so spline TAUC
           * to locate UTAU at EPIC interfaces. Use log pressure as the independent variable.
           */
          for (lev = 0; lev <= ds->nlyr; lev++) {
            tau_table[lev].y = TAUC(lev);
          }
          spline_pchip(ds->nlyr+1,tau_table);

          ii = -2;
          for (lu = 1; lu <= ds->ntau; lu++) {
            K        = lu;
            x        = log(P3(K,J,I));
            ii       = hunt_place_in_table(ds->nlyr+1,tau_table,x,&dx,ii);
            UTAU(lu) = splint_pchip(x,tau_table+ii,dx);
          }

          /*
           * Calculate the radiative transfer.
           */
          if (!midir->planck && !midir->scat_yes) {
            beer_law_only(ds,out,radius);
          }
          else {
            c_twostr(ds,out,deltam,gg,ierror,radius);
          }

          /*
           * Sum flux divergence; convert to heating rate below.
           */
          for (K = KLO; K < KHI; K++) {
            lu              = K;
            HEATING(K,J,I) += DFDT(lu)*DTAUC(lu);
          }
        } /* iwn loop */

        /*------------------------------------------------------------------*
         * Convert flux divergence [W/m^2] to specific heating rate [W/kg]. *
         *------------------------------------------------------------------*/

        for (K = KLO; K < KHI; K++) {
          HEATING(K,J,I) /= H3(K,J,I)*grid.dsgth[2*K+1];
        }
      }  /* end I loop */
    }  /* end J loop */
  }
  else if (action == EPIC_ALLOC) {
    /*
     * Allocate longwave band memory and DISORT state and output memory.
     */
    fir = (rt_band *)calloc(1,sizeof(rt_band));
    if (!fir) {
      sprintf(Message,"error allocating fir");
      epic_error(dbmsname,Message);
    }

    midir = (rt_band *)calloc(1,sizeof(rt_band));
    if (!midir) {
      sprintf(Message,"error allocating midir");
      epic_error(dbmsname,Message);
    }

    ds = (disort_state *)calloc(1,sizeof(disort_state));
    if (!ds) {
      sprintf(Message,"error allocating ds");
      epic_error(dbmsname,Message);
    }

    out = (disort_output *)calloc(1,sizeof(disort_output));
    if (!out) {
      sprintf(Message,"error allocating out");
      epic_error(dbmsname,Message);
    }

    /*
     * DISORT layers correspond to EPIC interior interfaces
     */
    ds->nlyr = grid.nk-1;

    /*
     * Allocate opacity arrays and read in longwave data.
     */
    init_rt_longwave(fir,midir,ds);

    /*------------------*
     * c_twostr() setup *
     *------------------*/

    sprintf(ds->header,"EPIC Model (Version %4.2f): %s; %s",grid.epic_version,planet->name,dbmsname);

    /*
     * Limit to 2 streams for computational efficiency, and use c_twostr()
     * instead of c_disort (the former is optimized for two streams).
     *
     * NOTE: It is quick and easy to generalize this to multiple streams and c_disort().
     */
    ds->nstr = 2;

    /*
     * Set convergence criterion for azimuthal series, usually 0.
     */
    ds->accur = 0.;

    /*
     * Turn off diagnostic print statements and status messages.
     */
    ds->flag.prnt[0]=ds->flag.prnt[1]=ds->flag.prnt[2]=ds->flag.prnt[3]=ds->flag.prnt[4]=FALSE;
    ds->flag.quiet = TRUE;

    /*
     * Use original Nakajima-Tanaka intensity correction algorithm,
     *
     * NOTE: This needs to be reevaluated for the case of a strongly forward-scattering medium,
     *       see Buras, Dowling and Emde (2012).
     */
    ds->flag.old_intensity_correction = TRUE;

    /*
     * Indicate we do not want the special-case output of just albedo and
     * transmissivity of the entire medium.
     */
    ds->flag.ibcnd  = GENERAL_BC;

    /*
     * Indicate that we do not need intensity output, only flux-related output.
     */
    ds->flag.onlyfl = TRUE;

    /*
     * Indicate we may want to include thermal emission, so allocate relevant memory.
     */
    ds->flag.planck = TRUE;

    /*
     * Use internally calculated polar angles for returned radiant quantities,
     * not user-supplied ones.
     */
    ds->flag.usrang = FALSE;

    /*
     * Use an isotropically reflecting bottom boundary with zero albedo.
     */
    ds->flag.lamber = TRUE;
    ds->bc.albedo   = 0.;

    /*
     *  The primary inputs to DISORT are the differential optical depths at the DISORT computational
     *  layers, DTAUC(lc).  Since these depend on EPIC's mole fractions, X, which are carried
     *  on the EPIC layer interfaces, DISORT's computational layers are EPIC interfaces, running
     *  from lc = 1 at the kk=3 EPIC interface (bottom of top layer) to lc = nlyr = grid.nk-1
     *  at the kk=2*grid.nk-1 EPIC interface (top of the bottom layer).
     *
     *  EPIC requires heating on its layer interfaces, HEAT3.
     *  Hence, we indicate we want DISORT's dfdt = d(flux)/d(tau) output at those positions,
     *  by setting ds->flag.usrtau = true, ds->ntau = grid.nk-1, and assigning 
     *  UTAU(lu) to be EPIC's interior interface values.
     */
    ds->nlyr        = grid.nk-1;
    ds->flag.usrtau = TRUE;
    ds->ntau        = grid.nk-1;

    /*
     * Number of phase functions for scattering calculations.
     * Needs to be >= to ds->nstr; often set equal to ds->nstr.
     */
    ds->nmom = ds->nstr;

    /*
     * Set azimuth angle of incident beam to zero, since we are not sensitive to it.
     */
    ds->bc.phi0 = 0.;

    /*
     * Set intensity of top-boundary isotropic illumination to zero.
     */
    ds->bc.fisot = 0.; 

    /*
     * Set spher to TRUE if there is any chance it will be TRUE during the run,
     * so that ds->zd[lev] is allocated.
     */
    ds->flag.spher = TRUE;

    /*
     * Allocate DISORT arrays.
     */
    c_twostr_state_alloc(ds);
    c_twostr_out_alloc(ds,out);

    /*
     * Allocate related vectors not included in DISORT structures.
     */
    gg        = dvector( 0,ds->nlyr-1,dbmsname);
    tauc      = dvector( 0,ds->nlyr,  dbmsname);
    tau_table = ftriplet(0,ds->nlyr,  dbmsname);

    deltam = FALSE;
    for (lc = 1; lc <= ds->nlyr; lc++) {
      /*
       * GG is the asymmetry factor; 1: complete forward scattering
       *                             0: isotropic scattering, Rayleigh scattering
       *                            -1: complete backscattering
       */
      GG(lc) = 0.;

      if (GG(lc) >= 0.95) {
        deltam = TRUE;
      }
    }

    if (ds->nstr > 2) {
      for (lc = 1; lc <= ds->nlyr; lc++) {
        c_getmom(RAYLEIGH,GG(lc),ds->nmom,&PMOM(0,lc));
      }
    }

    /*
     * Neglect longwave component of solar insolation.
     */
    ds->bc.fbeam = 0.;

    /*----------------------*
     * End c_twostr() setup *
     *----------------------*/
  }
  else if (action == EPIC_FREE) {
    /*
     * Free memory.
     */
    if (gg)        free_dvector(gg,        0,ds->nlyr-1,dbmsname);
    if (tauc)      free_dvector(tauc,      0,ds->nlyr,  dbmsname);
    if (tau_table) free_ftriplet(tau_table,0,ds->nlyr,  dbmsname);

    if (fir || midir) free_rt_longwave(fir,midir,ds);

    if (fir)   free(fir);
    if (midir) free(midir);

    if (ds && out) c_twostr_out_free(ds,out);
    if (ds)        c_twostr_state_free(ds);

    if (ds)  free(ds);
    if (out) free(out);
  }
  else {
    sprintf(Message,"unrecoganized action=%d",action);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============= end of rt_longwave() ==============================*/

/*============ init_rt_longwave() =================================*/

/*
 * Allocate arrays and initialize absorption
 * as appropriate for far and mid infrared bands.
 *
 * Fortran subroutines:
 *   firver3()
 *   midpandtver3()
 *   inith2h2opacver2()
 *
 * NOTE: "Shortwave" and "longwave" are convenient partitions that
 *       are used in meteorological radiative transfer to be synonymous
 *       with "solar insolation" and "planetary emission", but this is a
 *       temperature-sensitive distinction.  Our current radiative-transfer
 *       partition, with shortwave = nir, vis, uv, and longwave = midir, fir,
 *       is appropriate for the outer solar system. However, hot planets like
 *       Venus and hot-Jupiters have significant emission at wavelengths
 *       shorter than midir, and for them we will need to revisit this
 *       partition.
 */

/*
 * If the band fraction for thermal emission is greater or equal to BAND_FRACTION_THRESHOLD
 * for any temperature in the reference profile, then planck is set to true, such
 * that thermal emission is included in the radiative transfer calculations.
 */
#define BAND_FRACTION_THRESHOLD 0.01

void init_rt_longwave(rt_band      *fir,
                      rt_band      *midir,
                      disort_state *ds)
{
  int
    iwn,itemp,ipress,
    K,kk,jj;
  int
    nc_err,nc_id,nc_dimid,nc_varid;
  size_t
    nc_len,
   *start,
   *count;
  double
    tmp,band_fraction,
    wavelength_start,wavelength_end;
  double
   *opac;
  char
    header[N_STR];
  FILE
   *infile;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_rt_longwave";

  /*--------------------------------------------------------------------*
   * FIR (far infrared, 0 to 600 cm-1)                                  *
   *--------------------------------------------------------------------*/

  /* Print setup progress to stdout */
  if (IAMNODE == NODE0) {
    fprintf(stdout," fir"); fflush(stdout);
  }

  fir->wn_start =   0.;
  fir->wn_end   = 600.;

  /*
   * Determine whether thermal radiation by the atmosphere itself 
   * should be included.
   */
  wavelength_start = 1./(100.*fir->wn_end);
  fir->planck      = FALSE;
  for (K = 0; K <= KHI; K++) {
    kk = 2*K+1;
    /*
     * Calculate the band-emission fraction for the reference profile
     * of temperature.
     *
     * NOTE: We have FIR defined to span all the way to infinite wavelength,
     *       which corresponds to a blackbody_fraction value of 1.
     */
    band_fraction = 1.-blackbody_fraction(wavelength_start,grid.t_ref[kk]);
    if (band_fraction >= BAND_FRACTION_THRESHOLD) {
      fir->planck = TRUE;
      break;
    }
  }

  /*
   * Input T. Greathouse's collision-induced absorption (cia) 
   * coefficients for H2-H2, H2-He, and H2-CH4.
   */

  /* Allocate memory for cia structure */
  fir->cia = (rt_cia *)calloc(1,sizeof(rt_cia));
  if (!fir->cia) {
    sprintf(Message,"error allocating fir->cia");
    epic_error(dbmsname,Message);
  }

  fir->cia->on = TRUE;

  /* 
   * Open FIR collision-induced absorption (cia) data file
   */
  sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/fir/fir_cia_greathouse.nc");
  nc_err = nc_open(header,NC_NOWRITE,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
    epic_error(dbmsname,Message);
  }

  /* Read number of temperatures */
  nc_err = nc_inq_dimid(nc_id,"temperature",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  fir->cia->nt = (int)nc_len;

  /* Allocate memory for temperatures */
  fir->cia->t = fvector(0,fir->cia->nt-1,dbmsname);

  /* Read in temperature array */
  nc_err = nc_inq_varid(nc_id,"temperature",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,fir->cia->t);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Read number of wavenumbers */
  nc_err = nc_inq_dimid(nc_id,"wavelength_center",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  fir->nwn = (int)nc_len;

  /* Allocate memory for wavenumbers */
  fir->wn   = fvector(0,fir->nwn-1,dbmsname);
  fir->wnlo = fvector(0,fir->nwn-1,dbmsname);
  fir->wnhi = fvector(0,fir->nwn-1,dbmsname);

  /* Allocate memory for cia data */
  fir->cia->h2_h2  = fvector(0,fir->nwn*fir->cia->nt-1, dbmsname);
  fir->cia->h2_he  = fvector(0,fir->nwn*fir->cia->nt-1, dbmsname);
  if (var.species[CH_4_INDEX].on == TRUE) {
    fir->cia->h2_ch4 = fvector(0,fir->nwn*fir->cia->nt-1,dbmsname);
  }

  /* 
   * Read in wavenumber arrays.
   */
  nc_err = nc_inq_varid(nc_id,"wavelength_center",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,fir->wn);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_inq_varid(nc_id,"wavelength_low",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,fir->wnlo);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_inq_varid(nc_id,"wavelength_high",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,fir->wnhi);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Read number of pairs of gases included for collision-induced absorption */
  nc_err = nc_inq_dimid(nc_id,"continuum_opacity_source",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Allocate memory for opac, start and count vectors */
  opac  = fvector(0,(int)nc_len-1,dbmsname);
  start = (size_t *)calloc(3,sizeof(size_t));
  if (!start) {
    sprintf(Message,"error allocating start vector");
    epic_error(dbmsname,Message);
  }
  count = (size_t *)calloc(3,sizeof(size_t));
  if (!count) {
    sprintf(Message,"error allocating count vector");
    epic_error(dbmsname,Message);
  }

  /* 
   * Read in cia coefficient data [1/(km amagat^2)]
   */
  nc_err = nc_inq_varid(nc_id,"continuum_opacity",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  start[2] = 0;
  count[0] = 1;
  count[1] = 1;
  count[2] = nc_len;
  for (iwn = 1; iwn <= fir->nwn; iwn++) {
    start[0] = iwn-1;
    for (itemp = 1; itemp <= fir->cia->nt; itemp++) {
      start[1] = itemp-1;
      nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,opac);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      FIR_CIA_H2_H2( itemp,iwn) = opac[0];
      FIR_CIA_H2_HE( itemp,iwn) = opac[1];
      if (var.species[CH_4_INDEX].on == TRUE) {
        FIR_CIA_H2_CH4(itemp,iwn) = opac[2];
      }
    }
  }
  /* Close data file */
  nc_err = nc_close(nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Allocate memory for dtauc */
  fir->dtauc = fvector(0,fir->nwn*ds->nlyr-1,dbmsname);

  /* Free allocated memory */
  free_fvector(opac,0,(int)nc_len-1,dbmsname);
  free(start);
  free(count);

  /*--------------------------------------------------------------------*
   * MIDIR (mid infrared, 600 to 1600 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Print setup progress to stdout */
  if (IAMNODE == NODE0) {
    fprintf(stdout,", midir"); fflush(stdout);
  }

  midir->wn_start =  600.;
  midir->wn_end   = 1600.;

  /*
   * Determine whether thermal radiation by the atmosphere itself 
   * should be included.
   */
  wavelength_start = 1./(100.*midir->wn_end  );
  wavelength_end   = 1./(100.*midir->wn_start);
  midir->planck    = FALSE;
  for (K = 0; K <= KHI; K++) {
    kk = 2*K+1;
    /*
     * Calculate the band-emission fraction for the reference profile
     * of temperature.
     */
    band_fraction = blackbody_fraction(wavelength_end,  grid.t_ref[kk])
                   -blackbody_fraction(wavelength_start,grid.t_ref[kk]);
    if (band_fraction >= BAND_FRACTION_THRESHOLD) {
      midir->planck = TRUE;
      break;
    }
  }

  /*
   * Input T. Greathouse's MIDIR correlated-k coefficients for single-species
   * absorption and pairwise collision-induced absorption (cia).
   *
   * NOTE: The correlated-k bookkeeping is slightly different here for the  
   *       MIDIR compared to the NIR below. For the NIR (Irwin et al data),
   *       a fourth dimension is included to handle the correlated-k approach,
   *       the other three dimensions being central wavenumber, temperature
   *       and pressure.  Instead, here a succession of data values are used with
   *       the same central wavenumber and varying (wnlo,wnhi) intervals that
   *       span each wavenumber neighborhood.      
   */

  /* Allocate memory for k_nu structure */
  midir->k_nu = (rt_k_nu *)calloc(1,sizeof(rt_k_nu));
  if (!midir->k_nu) {
    sprintf(Message,"error allocating midir->k_nu");
    epic_error(dbmsname,Message);
  }

  /* Open MIDIR CH4 absorption coefficient data file */
  sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_ch4_greathouse.nc");
  nc_err = nc_open(header,NC_NOWRITE,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
    epic_error(dbmsname,Message);
  }

  /* Read in number of temperatures */
  nc_err = nc_inq_dimid(nc_id,"temperature",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  midir->k_nu->nt = (int)nc_len;

  /* Allocate memory for temperatures */
  midir->k_nu->t = fvector(0,midir->k_nu->nt-1,dbmsname);

  /* Read in temperature array */
  nc_err = nc_inq_varid(nc_id,"temperature",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->k_nu->t);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Read in number of pressures */
  nc_err = nc_inq_dimid(nc_id,"pressure",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  midir->k_nu->np = (int)nc_len;

  /* Allocate memory for pressures */
  midir->k_nu->p = fvector(0,midir->k_nu->np-1,dbmsname);

  /* Read in pressure array */
  nc_err = nc_inq_varid(nc_id,"pressure",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
   epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->k_nu->p);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  /*
   * Convert input-file pressure units: [bar] -> [Pa]
   */
  for (ipress = 1; ipress <= midir->k_nu->np; ipress++) {
    MIDIR_K_NU_P(ipress) *= 1.e+5;
  }

  /* Read number of wavenumbers */
  nc_err = nc_inq_dimid(nc_id,"wavelength_center",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  midir->nwn = (int)nc_len;

  /* Allocate memory for wavenumbers */
  midir->wn   = fvector(0,midir->nwn-1,dbmsname);
  midir->wnlo = fvector(0,midir->nwn-1,dbmsname);
  midir->wnhi = fvector(0,midir->nwn-1,dbmsname);

  /* Read in wavenumber arrays */
  nc_err = nc_inq_varid(nc_id,"wavelength_center",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->wn);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_inq_varid(nc_id,"wavelength_low",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->wnlo);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_inq_varid(nc_id,"wavelength_high",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->wnhi);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* 
   * Start with k_nu off, turn on if any relevant species are invoked.
   */

  midir->k_nu->on = FALSE;

  /* Allocate memory for k_nu values */
  if (var.species[CH_4_INDEX].on == TRUE) {
    midir->k_nu->on   = TRUE;
    midir->k_nu->ch4 = fvector(0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
  }
  if (var.species[C_2H_2_INDEX].on == TRUE) {
    midir->k_nu->on    = TRUE;
    midir->k_nu->c2h2 = fvector(0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
  }
  if (var.species[C_2H_6_INDEX].on == TRUE) {
    midir->k_nu->on    = TRUE;
    midir->k_nu->c2h6 = fvector(0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
  }
  if (var.species[NH_3_INDEX].on == TRUE) {
    midir->k_nu->on   = TRUE;
    midir->k_nu->nh3 = fvector(0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
  }
  if (var.species[PH_3_INDEX].on == TRUE) {
    midir->k_nu->on   = TRUE;
    midir->k_nu->ph3 = fvector(0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
  }

  if (midir->k_nu->on) {
    /* Allocate bookkeeping memory */
    start = (size_t *)calloc(3,sizeof(size_t));
    if (!start) {
      sprintf(Message,"error allocating start vector");
      epic_error(dbmsname,Message);
    }
    count = (size_t *)calloc(3,sizeof(size_t));
    if (!count) {
      sprintf(Message,"error allocating count vector");
      epic_error(dbmsname,Message);
    }
 
    if (var.species[CH_4_INDEX].on == TRUE) {
      /* 
       * Read in CH4 absorption coefficient data
       */
      nc_err = nc_inq_varid(nc_id,"gaseous_opacity",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
      ipress   = 1;
      start[2] = ipress-1;
      count[0] = 1;
      count[1] = 1;
      count[2] = midir->k_nu->np;
      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        start[0] = iwn-1;
        for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
          start[1] = itemp-1;
          nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,&MIDIR_K_NU_CH4(ipress,itemp,iwn));
          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
      /* Close data file */
      nc_err = nc_close(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    if (var.species[C_2H_2_INDEX].on == TRUE) {
      /*
       * Read in C2H2 absorption coefficient data
       */
      sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_c2h2_greathouse.nc");
      nc_err = nc_open(header,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
        epic_error(dbmsname,Message);
      }

      nc_err = nc_inq_varid(nc_id,"gaseous_opacity",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
      ipress   = 1;
      start[2] = ipress-1;
      count[0] = 1;
      count[1] = 1;
      count[2] = midir->k_nu->np;
      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        start[0] = iwn-1;
        for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
          start[1] = itemp-1;
          nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,&MIDIR_K_NU_C2H2(ipress,itemp,iwn));
          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
      /* Close data file */
      nc_err = nc_close(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    if (var.species[C_2H_6_INDEX].on == TRUE) {
      /*
       * Read in C2H6 absorption coefficient data
       */
      sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_c2h6_greathouse.nc");
      nc_err = nc_open(header,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
        epic_error(dbmsname,Message);
      }

      nc_err = nc_inq_varid(nc_id,"gaseous_opacity",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
      ipress   = 1;
      start[2] = ipress-1;
      count[0] = 1;
      count[1] = 1;
      count[2] = midir->k_nu->np;
      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        start[0] = iwn-1;
        for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
          start[1] = itemp-1;
          nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,&MIDIR_K_NU_C2H6(ipress,itemp,iwn));
          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
      /* Close data file */
      nc_err = nc_close(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    if (var.species[NH_3_INDEX].on == TRUE) {
      /*
       * Read in NH3 absorption coefficient data
       */
      sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_nh3_greathouse.nc");
      nc_err = nc_open(header,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
        epic_error(dbmsname,Message);
      }

      nc_err = nc_inq_varid(nc_id,"gaseous_opacity",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
      ipress   = 1;
      start[2] = ipress-1;
      count[0] = 1;
      count[1] = 1;
      count[2] = midir->k_nu->np;
      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        start[0] = iwn-1;
        for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
          start[1] = itemp-1;
          nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,&MIDIR_K_NU_NH3(ipress,itemp,iwn));
          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
      /* Close data file */
      nc_err = nc_close(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    if (var.species[PH_3_INDEX].on == TRUE) {
      /*
       * Read in PH3 absorption coefficient data
       */
      sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_ph3_greathouse.nc");
      nc_err = nc_open(header,NC_NOWRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
        epic_error(dbmsname,Message);
      }

      nc_err = nc_inq_varid(nc_id,"gaseous_opacity",&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
      ipress   = 1;
      start[2] = ipress-1;
      count[0] = 1;
      count[1] = 1;
      count[2] = midir->k_nu->np;
      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        start[0] = iwn-1;
        for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
          start[1] = itemp-1;
          nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,&MIDIR_K_NU_PH3(ipress,itemp,iwn));
          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
      /* Close data file */
      nc_err = nc_close(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    /* Free allocated memory */
    free(start);
    free(count);
  }

  /*
   * Input MIDIR collision-induced absorption (cia) 
   * coefficients for H2-H2, H2-He, and H2-CH4.
   */

  /* Allocate memory for cia structure */
  midir->cia = (rt_cia *)calloc(1,sizeof(rt_cia));
  if (!midir->cia) {
    sprintf(Message,"error allocating midir->cia");
    epic_error(dbmsname,Message);
  }

  midir->cia->on = TRUE;

  /*
   * Open MIDIR collision-induced absorption (cia) data file.
   */
  sprintf(header,"%s",EPIC_PATH"/data/rt/absorption/midir/midir_cia_greathouse.nc");
  nc_err = nc_open(header,NC_NOWRITE,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),header);
    epic_error(dbmsname,Message);
  }

  /* Read number of temperatures */
  nc_err = nc_inq_dimid(nc_id,"temperature",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  midir->cia->nt = (int)nc_len;

  /* Allocate memory for temperatures */
  midir->cia->t = fvector(0,midir->cia->nt-1, dbmsname);

  /* Read in temperature array */
  nc_err = nc_inq_varid(nc_id,"temperature",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_get_var_double(nc_id,nc_varid,midir->cia->t);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* 
   * Read number of wavenumbers, to spot check that it is the same as above.
   */
  nc_err = nc_inq_dimid(nc_id,"wavelength_center",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /*
   * Spot check that nwn is the same as above.
   */
  if ((int)nc_len != midir->nwn) {
    sprintf(Message,"no. wavenumbers in %s is %d, does not match midir->nwn=%d",header,(int)nc_len,midir->nwn);
    epic_error(dbmsname,Message);
  }

  /* Allocate memory for cia data */
  midir->cia->h2_h2  = fvector(0,midir->nwn*midir->cia->nt-1,dbmsname);
  midir->cia->h2_he  = fvector(0,midir->nwn*midir->cia->nt-1,dbmsname);
  midir->cia->h2_ch4 = fvector(0,midir->nwn*midir->cia->nt-1,dbmsname);

  /* Read number of pairs of gases included for collision-induced absorption */
  nc_err = nc_inq_dimid(nc_id,"continuum_opacity_source",&nc_dimid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  nc_err = nc_inq_dimlen(nc_id,nc_dimid,&nc_len);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Allocate memory */
  opac  = fvector(0,(int)nc_len-1,dbmsname);
  start = (size_t *)calloc(3,sizeof(size_t));
  if (!start) {
    sprintf(Message,"error allocating start");
    epic_error(dbmsname,Message);
   }
  count = (size_t *)calloc(3,sizeof(size_t));
  if (!count) {
    sprintf(Message,"error allocating count");
    epic_error(dbmsname,Message);
  }

  /* Read in absorption coefficient data */
  nc_err = nc_inq_varid(nc_id,"continuum_opacity",&nc_varid);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }
  start[2] = 0;
  count[0] = 1;
  count[1] = 1;
  count[2] = nc_len;
  for (iwn = 1; iwn <= midir->nwn; iwn++) {
    start[0] = iwn-1;
    for (itemp = 1; itemp <= midir->cia->nt; itemp++) {
      start[1] = itemp-1;
      nc_err = nc_get_vara_double(nc_id,nc_varid,start,count,opac);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      MIDIR_CIA_H2_H2( itemp,iwn) = opac[0];
      MIDIR_CIA_H2_HE( itemp,iwn) = opac[1];
      if (var.species[CH_4_INDEX].on == TRUE) {
        MIDIR_CIA_H2_CH4(itemp,iwn) = opac[2];
      }
    }
  }
  /* Close data file */
  nc_err = nc_close(nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /* Allocate memory for dtauc */
  midir->dtauc = fvector(0,midir->nwn*ds->nlyr-1,dbmsname);

  /* Free allocated memory */
  free_fvector(opac,0,(int)nc_len-1,dbmsname);
  free(start);
  free(count);

  return;
}

#undef BAND_FRACTION_THRESHOLD

/*============ end of init_rt_longwave() ==========================*/

/*============ free_rt_longwave() =================================*/

/*
 * Free memory dynamically allocated in init_rt_longwave()
 */

void free_rt_longwave(rt_band      *fir,
                      rt_band      *midir,
                      disort_state *ds)
{
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_rt_longwave";

  /*--------------------------------------------------------------------*
   * FIR (far infrared, 0 to 600 cm-1)                                  *
   *--------------------------------------------------------------------*/

  /* Free wavenumbers */
  if (fir->wn)   free_fvector(fir->wn,  0,fir->nwn-1,dbmsname);
  if (fir->wnlo) free_fvector(fir->wnlo,0,fir->nwn-1,dbmsname);
  if (fir->wnhi) free_fvector(fir->wnhi,0,fir->nwn-1,dbmsname);

  /* Free dtauc */
  if (fir->dtauc) free_fvector(fir->dtauc,0,fir->nwn*ds->nlyr-1,dbmsname);

  /* Free cia */
  if (fir->cia) {
    if (fir->cia->t)      free_fvector(fir->cia->t,0,fir->cia->nt-1,dbmsname);
    if (fir->cia->h2_h2)  free_fvector(fir->cia->h2_h2, 0,fir->nwn*fir->cia->nt-1,dbmsname);
    if (fir->cia->h2_he)  free_fvector(fir->cia->h2_he, 0,fir->nwn*fir->cia->nt-1,dbmsname);
    if (fir->cia->h2_ch4) free_fvector(fir->cia->h2_ch4,0,fir->nwn*fir->cia->nt-1,dbmsname);
    free(fir->cia);
  }

  /*--------------------------------------------------------------------*
   * MIDIR (mid infrared, 600 to 1600 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Free wavenumbers */
  if (midir->wn)   free_fvector(midir->wn,  0,midir->nwn-1,dbmsname);
  if (midir->wnlo) free_fvector(midir->wnlo,0,midir->nwn-1,dbmsname);
  if (midir->wnhi) free_fvector(midir->wnhi,0,midir->nwn-1,dbmsname);

  /* Free dtauc */
  if (midir->dtauc) free_fvector(midir->dtauc,0,midir->nwn*ds->nlyr-1,dbmsname);

  /* Free k_nu */
  if (midir->k_nu) {
    if (midir->k_nu->t)    free_fvector(midir->k_nu->t,0,midir->k_nu->nt-1,dbmsname);
    if (midir->k_nu->p)    free_fvector(midir->k_nu->p,0,midir->k_nu->np-1,dbmsname);
    if (midir->k_nu->ch4)  free_fvector(midir->k_nu->ch4, 0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
    if (midir->k_nu->c2h2) free_fvector(midir->k_nu->c2h2,0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
    if (midir->k_nu->c2h6) free_fvector(midir->k_nu->c2h6,0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
    if (midir->k_nu->nh3)  free_fvector(midir->k_nu->nh3, 0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
    if (midir->k_nu->ph3)  free_fvector(midir->k_nu->ph3, 0,midir->nwn*midir->k_nu->nt*midir->k_nu->np-1,dbmsname);
    free(midir->k_nu);
  }

  /* Free cia */
  if (midir->cia) {
    if (midir->cia->t)      free_fvector(midir->cia->t,0,midir->cia->nt-1,dbmsname);
    if (midir->cia->h2_h2)  free_fvector(midir->cia->h2_h2, 0,midir->nwn*midir->cia->nt-1,dbmsname);
    if (midir->cia->h2_he)  free_fvector(midir->cia->h2_he, 0,midir->nwn*midir->cia->nt-1,dbmsname);
    if (midir->cia->h2_ch4) free_fvector(midir->cia->h2_ch4,0,midir->nwn*midir->cia->nt-1,dbmsname);
    free(midir->cia);
  }

  return;
}

/*============ end of free_rt_longwave() ==========================*/

/*============ fir_opacity() ======================================*/

#define TCIA(itemp) tcia_table[itemp-1].x

/*
 * Fortran subroutine:
 *   firver3()
 */
void fir_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *fir,
                 disort_state *ds)
{
  int
    itemp,iwn,lc,
    K,lev;
  double
    dtau,m,b,tau,
    dz_km,dt,temp,
    n_amagat;
  static int
    initialized = FALSE;
  static float_triplet
   *tcia_table;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="fir_opacity";

  if (!initialized) {
    /* Allocate Memory */
    tcia_table = ftriplet(0,fir->cia->nt,dbmsname);

    /* Assign values */
    for (itemp = 1; itemp <= fir->cia->nt; itemp++) {
      TCIA(itemp) = FIR_CIA_T(itemp);
    }

    initialized = TRUE;
  }

  ds->flag.planck = fir->planck;

  if (ds->flag.planck) {
    ds->bc.temis = 1.;   /* emissivity at the top boundary */
    ds->bc.ttemp = T3(0,J,I);
    for (K = KLO; K <= KHI; K++) {
      lev         = K-1;
      TEMPER(lev) = T2(K,J,I);
    }
    ds->bc.btemp = T3(grid.nk,J,I);
  }

  /* 
   * Zero fir->dtauc
   */
  memset(fir->dtauc,0,fir->nwn*ds->nlyr*sizeof(double));

  for (lc = 1; lc <= ds->nlyr; lc++) {
    K = lc;

    /*
     * Calculate number density [amagat]
     */
    n_amagat = (P3(K,J,I)*273.15)/(T3(K,J,I)*1.01325e+5);

    /*
     * DISORT layer thickness [km]
     */
    dz_km = (Z2(K,J,I)-Z2(K+1,J,I))/1.e+3;

    if (fir->cia->on) {
      temp  = LIMIT_RANGE(TCIA(1),T3(K,J,I),TCIA(fir->cia->nt));
      itemp = hunt_place_in_table(fir->cia->nt,tcia_table,temp,&dt,itemp);
      itemp++;  /* Shift to unit-based */

      for (iwn = 1; iwn <= fir->nwn; iwn++) {
        /* H_2-H_2 */
        m     = (FIR_CIA_H2_H2(itemp+1,iwn)-FIR_CIA_H2_H2(itemp,iwn))/dt;
        b     = FIR_CIA_H2_H2(itemp,iwn)-m*FIR_CIA_T(itemp);
        dtau  = m*temp+b;
        tau   = dtau*planet->x_h2;

        /* H_2-He */
        m     = (FIR_CIA_H2_HE(itemp+1,iwn)-FIR_CIA_H2_HE(itemp,iwn))/dt;
        b     = FIR_CIA_H2_HE(itemp,iwn)-m*FIR_CIA_T(itemp);
        dtau  = m*temp+b;
        tau  += dtau*planet->x_he;

        tau               *= planet->x_h2*n_amagat*n_amagat*dz_km;
        FIR_DTAUC(lc,iwn) += tau;
      }

      if (var.species[CH_4_INDEX].on == TRUE) {
        for (iwn = 1; iwn <= fir->nwn; iwn++) {
          /* H_2-CH_4 */
          m     = (FIR_CIA_H2_CH4(itemp+1,iwn)-FIR_CIA_H2_CH4(itemp,iwn))/dt;
          b     = FIR_CIA_H2_CH4(itemp,iwn)-m*FIR_CIA_T(itemp);
          dtau  = m*temp+b;
          tau   = dtau*X(CH_4_INDEX,VAPOR,K,J,I);

          tau               *= planet->x_h2*n_amagat*n_amagat*dz_km;
          FIR_DTAUC(lc,iwn) += tau;
        }
      }
    }
  }

  return;
}

#undef TCIA

/*============ end of fir_opacity() ===============================*/

/*============ midir_opacity() ====================================*/

/*
 * Calculate opacities in the mid-IR, between 600 and 1600 cm-1.
 *
 * Fortran subroutine:
 *   midpandtver3()
 */

#define LOGPTAU(ipress) logp_table[ipress-1].x
#define TTAU(itemp)     t_table[itemp-1].x
#define TCIA(itemp)     tcia_table[itemp-1].x

void midir_opacity(planetspec   *planet,
                   int           J,
                   int           I,
                   rt_band      *midir,
                   disort_state *ds)
{
  int
    lc,iwn,
    K,lev,
    itemp,ipress,
    test;
  double
    temp,dtau,m,b,
    dtau1,dtau2,tau,
    dz_km,dt,
    invt1,invt2,invt,
    ec,logp,dlogp,
    n_amagat;
  const double
    hc_k = 1.43879;  /* h*c/k_b [K/cm-1] */
  static int
    initialized = FALSE;
  static float_triplet
   *logp_table,
   *t_table,
   *tcia_table;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="midir_opacity";

  if (!initialized) {
    if (midir->k_nu->on) {
      /* Allocate memory */
      logp_table = ftriplet(0,midir->k_nu->np-1,dbmsname);
      t_table    = ftriplet(0,midir->k_nu->nt-1,dbmsname);
      tcia_table = ftriplet(0,midir->cia->nt,   dbmsname);

      /* Assign values */
      for (ipress = 1; ipress <= midir->k_nu->np; ipress++) {
        LOGPTAU(ipress) = log(MIDIR_K_NU_P(ipress));
      }
      for (itemp = 1; itemp <= midir->k_nu->nt; itemp++) {
        TTAU(itemp) = MIDIR_K_NU_T(itemp);
      }
      for (itemp = 1; itemp <= midir->cia->nt; itemp++) {
        TCIA(itemp) = MIDIR_CIA_T(itemp);
      }
    }

    initialized = TRUE;
  }

  ds->flag.planck = midir->planck;

  if (ds->flag.planck) {
    ds->bc.temis = 1.;   /* emissivity at the top boundary */
    ds->bc.ttemp = T3(0,J,I);
    for (K = KLO; K <= KHI; K++) {
      lev         = K-1;
      TEMPER(lev) = T2(K,J,I);
    }
    ds->bc.btemp = T3(grid.nk,J,I);
  }

  /* 
   * Zero midir->dtauc
   */
  memset(midir->dtauc,0,midir->nwn*ds->nlyr*sizeof(double));

  for (lc = 1; lc <= ds->nlyr; lc++) {
    K = lc;

    /*
     * Calculate number density [amagat]
     */
    n_amagat = (P3(K,J,I)*273.15)/(T3(K,J,I)*1.01325e+5);

    /*
     * DISORT layer thickness [km]
     */
    dz_km = (Z2(K,J,I)-Z2(K+1,J,I))/1.e+3;

    if (midir->cia->on){
      temp  = LIMIT_RANGE(TCIA(1),T3(K,J,I),TCIA(midir->cia->nt));
      itemp = hunt_place_in_table(midir->cia->nt,tcia_table,temp,&dt,itemp);
      itemp++;  /* Shift to unit-based */

      for (iwn = 1; iwn <= midir->nwn; iwn++) {
        /* H_2-H_2 */
        m     = (MIDIR_CIA_H2_H2(itemp+1,iwn)-MIDIR_CIA_H2_H2(itemp,iwn))/dt;
        b     = MIDIR_CIA_H2_H2(itemp,iwn)-m*MIDIR_CIA_T(itemp);
        dtau  = m*temp+b;
        tau   = dtau*planet->x_h2;

        /* H_2-He */
        m     = (MIDIR_CIA_H2_HE(itemp+1,iwn)-MIDIR_CIA_H2_HE(itemp,iwn))/dt;
        b     = MIDIR_CIA_H2_HE(itemp,iwn)-m*MIDIR_CIA_T(itemp);
        dtau  = m*temp+b;
        tau  += dtau*planet->x_he;

        tau                 *= planet->x_h2*n_amagat*n_amagat*dz_km;
        MIDIR_DTAUC(lc,iwn) += tau;
      }

      if (var.species[CH_4_INDEX].on == TRUE) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          /* H_2-CH_4 */
          m     = (MIDIR_CIA_H2_CH4(itemp+1,iwn)-MIDIR_CIA_H2_CH4(itemp,iwn))/dt;
          b     = MIDIR_CIA_H2_CH4(itemp,iwn)-m*MIDIR_CIA_T(itemp);
          dtau  = m*temp+b;
          tau   = dtau*X(CH_4_INDEX,VAPOR,K,J,I);

          tau                 *= planet->x_h2*n_amagat*n_amagat*dz_km;
          MIDIR_DTAUC(lc,iwn) += tau;
        }
      }
    }

    if (midir->k_nu->on) {
      logp   = LIMIT_RANGE(LOGPTAU(1),log(P3(K,J,I)),LOGPTAU(midir->k_nu->np));
      ipress = hunt_place_in_table(midir->k_nu->np,logp_table,logp,&dlogp,ipress);
      ipress++;   /* Shift to unit-based */

      temp  = LIMIT_RANGE(TTAU(1),T3(K,J,I),TTAU(midir->k_nu->nt));
      itemp = hunt_place_in_table(midir->k_nu->nt,t_table,temp,&dt,itemp);
      itemp++;  /* Shift to unit-based */

      invt1 = 1./MIDIR_K_NU_T(itemp  );
      invt2 = 1./MIDIR_K_NU_T(itemp+1);
      invt  = 1./temp;

      if (var.species[CH_4_INDEX].on) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          m     = (MIDIR_K_NU_CH4(ipress+1,itemp,iwn)-MIDIR_K_NU_CH4(ipress,itemp,iwn))/dlogp;
          b     = MIDIR_K_NU_CH4(ipress,itemp,iwn)-m*LOGPTAU(ipress);  
          dtau1 = m*logp+b;

          m     = (MIDIR_K_NU_CH4(ipress+1,itemp+1,iwn)-MIDIR_K_NU_CH4(ipress,itemp+1,iwn))/dlogp;
          b     = MIDIR_K_NU_CH4(ipress,itemp+1,iwn)-m*LOGPTAU(ipress);  
          dtau2 = m*logp+b;

          ec    = log(dtau1/dtau2)/(hc_k*(invt2-invt1));
          dtau  = dtau1*exp(hc_k*ec*(invt1-invt));

          MIDIR_DTAUC(lc,iwn) += dtau*X(CH_4_INDEX,VAPOR,K,J,I)*n_amagat*dz_km;
        }
      }

      if (var.species[C_2H_2_INDEX].on) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          m     = (MIDIR_K_NU_C2H2(ipress+1,itemp,iwn)-MIDIR_K_NU_C2H2(ipress,itemp,iwn))/dlogp;
          b     = MIDIR_K_NU_C2H2(ipress,itemp,iwn)-m*LOGPTAU(ipress);  
          dtau1 = m*logp+b;

          m     = (MIDIR_K_NU_C2H2(ipress+1,itemp+1,iwn)-MIDIR_K_NU_C2H2(ipress,itemp+1,iwn))/dlogp;
          b     = MIDIR_K_NU_C2H2(ipress,itemp+1,iwn)-m*LOGPTAU(ipress);  
          dtau2 = m*logp+b;

          ec    = log(dtau1/dtau2)/(hc_k*(invt2-invt1));
          dtau  = dtau1*exp(hc_k*ec*(invt1-invt));

          MIDIR_DTAUC(lc,iwn) += dtau*X(C_2H_2_INDEX,VAPOR,K,J,I)*n_amagat*dz_km;
        }
      }

      if (var.species[C_2H_6_INDEX].on) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          m     = (MIDIR_K_NU_C2H6(ipress+1,itemp,iwn)-MIDIR_K_NU_C2H6(ipress,itemp,iwn))/dlogp;
          b     = MIDIR_K_NU_C2H6(ipress,itemp,iwn)-m*LOGPTAU(ipress);  
          dtau1 = m*logp+b;

          m     = (MIDIR_K_NU_C2H6(ipress+1,itemp+1,iwn)-MIDIR_K_NU_C2H6(ipress,itemp+1,iwn))/dlogp;
          b     = MIDIR_K_NU_C2H6(ipress,itemp+1,iwn)-m*LOGPTAU(ipress);  
          dtau2 = m*logp+b;

          ec    = log(dtau1/dtau2)/(hc_k*(invt2-invt1));
          dtau  = dtau1*exp(hc_k*ec*(invt1-invt));

          MIDIR_DTAUC(lc,iwn) += dtau*X(C_2H_6_INDEX,VAPOR,K,J,I)*n_amagat*dz_km;
        }
      }

      if (var.species[NH_3_INDEX].on) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          m     = (MIDIR_K_NU_NH3(ipress+1,itemp,iwn)-MIDIR_K_NU_NH3(ipress,itemp,iwn))/dlogp;
          b     = MIDIR_K_NU_NH3(ipress,itemp,iwn)-m*LOGPTAU(ipress);  
          dtau1 = m*logp+b;

          m     = (MIDIR_K_NU_NH3(ipress+1,itemp+1,iwn)-MIDIR_K_NU_NH3(ipress,itemp+1,iwn))/dlogp;
          b     = MIDIR_K_NU_NH3(ipress,itemp+1,iwn)-m*LOGPTAU(ipress);  
          dtau2 = m*logp+b;

          ec    = log(dtau1/dtau2)/(hc_k*(invt2-invt1));
          dtau  = dtau1*exp(hc_k*ec*(invt1-invt));

          MIDIR_DTAUC(lc,iwn) += dtau*X(NH_3_INDEX,VAPOR,K,J,I)*n_amagat*dz_km;
        }
      }

      if (var.species[PH_3_INDEX].on) {
        for (iwn = 1; iwn <= midir->nwn; iwn++) {
          m     = (MIDIR_K_NU_PH3(ipress+1,itemp,iwn)-MIDIR_K_NU_PH3(ipress,itemp,iwn))/dlogp;
          b     = MIDIR_K_NU_PH3(ipress,itemp,iwn)-m*LOGPTAU(ipress);  
          dtau1 = m*logp+b;

          m     = (MIDIR_K_NU_PH3(ipress+1,itemp+1,iwn)-MIDIR_K_NU_PH3(ipress,itemp+1,iwn))/dlogp;
          b     = MIDIR_K_NU_PH3(ipress,itemp+1,iwn)-m*LOGPTAU(ipress);  
          dtau2 = m*logp+b;

          ec    = log(dtau1/dtau2)/(hc_k*(invt2-invt1));
          dtau  = dtau1*exp(hc_k*ec*(invt1-invt));

          MIDIR_DTAUC(lc,iwn) += dtau*X(PH_3_INDEX,VAPOR,K,J,I)*n_amagat*dz_km;
        }
      }
    }
  }

  return;
}

#undef LOGPTAU
#undef TTAU
#undef TCIA

/*============ end of midir_opacity() =============================*/

/*---------------------*
 * Shortwave functions *
 *---------------------*/

/*============= rt_shortwave() ====================================*/

/*
 * Partition of wavelength regions:
 *
 *   NIR     2000,   9500 cm-1
 *   VIS    10000,  40000 cm-1
 *   UV     40000, 100000 cm-1
 *
 * Input i_stride is used to subsample in the longitude direction,
 * to facilitate economical diurnal averaging.
 *
 * Fortran subroutine:
 *   heatingver12()
 */

void rt_shortwave(planetspec *planet,
                  EPIC_FLOAT *heating,
                  int         i_stride,
                  int         action)
{
  int
    K,J,I,
    kk,jj,
    lc,lev,lu,ii,
    iwn,ig;
  static int
    deltam,
    ierror[TWOSTR_NERR];
  double
    solar_flux_attenuation,fbeam0,
    subsolar_lat,alpha,
    glat,clat,
    radius;
  EPIC_FLOAT
    x,dx;
  static double
   *gg,
   *tauc;
  static float_triplet
   *tau_table;
  static rt_band
   *nir,
   *vis,
   *uv;
  static disort_state
   *ds;
  static disort_output
   *out;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="rt_shortwave";

  if (action == EPIC_APPLY) {
    /*
     * Zero heating array.
     */
    memset(heating,0,Nelem3d*sizeof(EPIC_FLOAT));

    /* 
     * Inverse-square attenuation of solar flux.
     * The function radius_vector() returns the planet's distance from the Sun [AU].
     */
    solar_flux_attenuation = 1./SQR(radius_vector(planet,var.model_time));

    /*
     * Subsolar latitude (planetocentric) [deg].
     */
    subsolar_lat = solar_declination(planet,L_s);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;

      /*
       * EPIC uses planetographic latitude, here denoted glat [DEG] and glatr [rad].
       * Planetocentric latitude is clat [deg].
       */
      glat = grid.lat[2*J+1];
      clat = lat_graphic_to_centric(glat,planet->re/planet->rp);

      if (strcmp(planet->name,"saturn") == 0) {
        /*
         * Account for the ring shadow on the atmosphere.
         *
         * NOTE: The function ring_shadow() returns the daytime-averaged ring shadow,
         *       not the local, instantaneous value.
         *
         * NOTE: The input latitudes to ring_shadow() are planetocentric.
         */
        solar_flux_attenuation *= ring_shadow(planet,clat,subsolar_lat); 
      }

      /*
       * Start at grid.ilo so that i_stride is correct.
       */
      for (I = grid.ilo; I <= IHI; I += i_stride) {
        if (I < ILO) {
          continue;
        }

        /*
         * Hour angle, alpha [deg]
         */
        alpha = (grid.lon[2*I+1]-east_longitude_of_solar_noon(planet,var.model_time,L_s));
        /*
         * Calculate cosine of solar zenith angle.
         * See Meeus (2005), eqn. (13.6)
         */
        ds->bc.umu0 = sin(glat*DEG)*sin(subsolar_lat*DEG)+cos(glat*DEG)*cos(subsolar_lat*DEG)*cos(alpha*DEG);

        /*
         * The Chapman function used to calculate spherical effects is
         * computationally expensive, so set ds->flag.spher to true only
         * where it is significant, as marked by cos(70.0*DEG) = 0.342.
         */
        ds->flag.spher = (ds->bc.umu0 < 0.342) ? TRUE : FALSE;

        if (ds->bc.umu0 <= 0. && ds->flag.spher == FALSE) {
         /*
          * Break when Sun is below horizon.
          */
          continue;
        }

        if (ds->flag.spher == TRUE) {
          /*
           * The spherical correction needs the following radius information.
           */
          radius = grid.rlt[2*KHI+1][jj];
          for (K = 1; K <= KHI; K++) {
            kk  = 2*K+1;
            lev = K-1;

            ds->zd[lev] = grid.rlt[kk][jj]-radius;
          }
        }

        /*
         * Set independent variable for spline on tau data,
         * used to specify user-defined optical depth positions, UTAU.
         * Use log pressure.
         */
        for (lev = 0; lev <= ds->nlyr; lev++) {
          K                = lev+1;
          tau_table[lev].x = log(P2(K,J,I));
        }

        /*----------------------*
         * NIR: 2000, 9500 cm-1 *
         *----------------------*/

        /*
         * SSALB refers to the single-scattering albedo.
         *
         * NOTE: We currently do not include single scattering for NIR.
         */
        nir->scat_yes = FALSE;
        memset(&SSALB(1),0,ds->nlyr*sizeof(double));

        nir_opacity(planet,J,I,nir,ds);

        /*
         * Loop over correlated-k bins.
         */
        for (ig = 1; ig <= nir->k_nu->ng; ig++) {
          /*
           * Loop over wavenumber bins.
           */
          for (iwn = 1; iwn <= nir->nwn; iwn++) {
            ds->wvnmhi = NIR_WNHI(iwn);
            ds->wvnmlo = NIR_WNLO(iwn);

            /*
             * Beam flux (assumed infinitely wide) at top of atmosphere.
             */
            fbeam0 = NIR_FLUXTOT(iwn)*solar_flux_attenuation;

            /*
             * Apply correlated-k weighting to incident beam flux.
             */
            ds->bc.fbeam = fbeam0*NIR_K_NU_DELG(ig);

            /*
             * Set DISORT-layer optical depth, DTAUC.
             */
            memcpy(&DTAUC(1),&NIR_DTAUC(1,ig,iwn),ds->nlyr*sizeof(double));

            /*
             * Construct TAUC(lev) in the same manner as is done
             * internally in DISORT, and then use it to specify user-defined
             * optical-depth positions, UTAU, for the output.
             *
             * NOTE: TAUC is vertically staggered with respect to DTAUC,
             *       hence it is aligned with EPIC layers.
             */
            TAUC(0) = 0.;
            for (lev = 1; lev <= ds->nlyr; lev++) {
              TAUC(lev) = TAUC(lev-1)+DTAUC(lev);
            }
 
            /*
             * EPIC needs heating on its interfaces, so spline TAUC
             * to locate UTAU at EPIC interfaces.
             */
            for (lev = 0; lev <= ds->nlyr; lev++) {
              tau_table[lev].y = TAUC(lev);
            }
            spline_pchip(ds->nlyr+1,tau_table);

            ii = -2;
            for (lu = 1; lu <= ds->ntau; lu++) {
              K        = lu;
              x        = log(P3(K,J,I));
              ii       = hunt_place_in_table(ds->nlyr+1,tau_table,x,&dx,ii);
              UTAU(lu) = splint_pchip(x,tau_table+ii,dx);
            }

            /*
             * Calculate the radiative transfer.
             */
            if (!nir->planck && !nir->scat_yes) {
              beer_law_only(ds,out,radius);
            }
            else {
              c_twostr(ds,out,deltam,gg,ierror,radius);
            }

            /*
             * Sum flux divergence; convert to heating rate below.
             */
            for (K = KLO; K < KHI; K++) {
              lu              = K;
              HEATING(K,J,I) += DFDT(lu)*DTAUC(lu);
            }
          }  /* iwn loop */
        } /* ig loop */

        /*------------------------*
         * VIS: 10000, 40000 cm-1 *
         *------------------------*/

        vis_opacity(planet,J,I,vis,ds);

        /*
         * Loop over wavenumber bins.
         */
        for (iwn = 1; iwn <= vis->nwn; iwn++) {
          ds->wvnmhi = VIS_WNHI(iwn);
          ds->wvnmlo = VIS_WNLO(iwn);

          /*
           * Beam flux (assumed infinitely wide) at top of atmosphere.
           */
          ds->bc.fbeam = VIS_FLUXTOT(iwn)*solar_flux_attenuation;

          /*
           * Set single-scattering albedo, SSALB, and DISORT-layer optical depth, DTAUC.
           */
          if (!vis->sigma->on) {
            vis->scat_yes = FALSE;
            memset(&SSALB(1),0,ds->nlyr*sizeof(double));
          }
          else {
            vis->scat_yes = TRUE;
            memcpy(&SSALB(1),&VIS_SSALB(1,iwn),ds->nlyr*sizeof(double));
          }
          memcpy(&DTAUC(1),&VIS_DTAUC(1,iwn),ds->nlyr*sizeof(double));

          /*
           * Construct TAUC(lev) in the same manner as is done
           * internally in DISORT, and then use it to specify user-defined
           * optical-depth positions, UTAU, for the output.
           *
           * NOTE: TAUC is vertically staggered with respect to DTAUC,
           *       hence it is aligned with EPIC layers.
           */
          TAUC(0) = 0.;
          for (lev = 1; lev <= ds->nlyr; lev++) {
            TAUC(lev) = TAUC(lev-1)+DTAUC(lev);
          }
 
          /*
           * EPIC needs heating on its interfaces, so spline TAUC
           * to locate UTAU at EPIC interfaces.
           */
          for (lev = 0; lev <= ds->nlyr; lev++) {
            tau_table[lev].y = TAUC(lev);
          }
          spline_pchip(ds->nlyr+1,tau_table);

          ii = -2;
          for (lu = 1; lu <= ds->ntau; lu++) {
            K        = lu;
            x        = log(P3(K,J,I));
            ii       = hunt_place_in_table(ds->nlyr+1,tau_table,x,&dx,ii);
            UTAU(lu) = splint_pchip(x,tau_table+ii,dx);
          }

          /*
           * Calculate the radiative transfer.
           */
          if (!vis->planck && !vis->scat_yes) {
            beer_law_only(ds,out,radius);
          }
          else {
            c_twostr(ds,out,deltam,gg,ierror,radius);
          }

          /*
           * Sum flux divergence; convert to heating rate below.
           */
          for (K = KLO; K < KHI; K++) {
            lu              = K;
            HEATING(K,J,I) += DFDT(lu)*DTAUC(lu);
          }
        } /* iwn loop */

        /*------------------------*
         * UV: 40000, 100000 cm-1 *
         *------------------------*/

        uv_opacity(planet,J,I,uv,ds);

        /*
         * Loop over wavenumber bins.
         */
        for (iwn = 1; iwn <= uv->nwn; iwn++) {
          ds->wvnmhi = UV_WNHI(iwn);
          ds->wvnmlo = UV_WNLO(iwn);

          /*
           * Beam flux (assumed infinitely wide) at top of atmosphere.
           */
          ds->bc.fbeam = UV_FLUXTOT(iwn)*solar_flux_attenuation;

          /*
           * Set single-scattering albedo, SSALB, and layer optical depth, DTAUC.
           */
          if (!uv->sigma->on) {
            uv->scat_yes = FALSE;
            memset(&SSALB(1),0,ds->nlyr*sizeof(double));
          }
          else {
            uv->scat_yes = TRUE;
            memcpy(&SSALB(1),&UV_SSALB(1,iwn),ds->nlyr*sizeof(double));
          }
          memcpy(&DTAUC(1),&UV_DTAUC(1,iwn),ds->nlyr*sizeof(double));

          /*
           * Construct TAUC(lev) in the same manner as is done
           * internally in DISORT, and then use it to specify user-defined
           * optical-depth positions, UTAU, for the output.
           *
           * NOTE: TAUC is vertically staggered with respect to DTAUC,
           *       hence it is aligned with EPIC layers.
           */
          TAUC(0) = 0.;
          for (lev = 1; lev <= ds->nlyr; lev++) {
            TAUC(lev) = TAUC(lev-1)+DTAUC(lev);
          }
 
          /*
           * EPIC needs heating on its interfaces, so spline TAUC
           * to locate UTAU at EPIC interfaces.
           */
          for (lev = 0; lev <= ds->nlyr; lev++) {
            tau_table[lev].y = TAUC(lev);
          }
          spline_pchip(ds->nlyr+1,tau_table);

          ii = -2;
          for (lu = 1; lu <= ds->ntau; lu++) {
            K        = lu;
            x        = log(P3(K,J,I));
            ii       = hunt_place_in_table(ds->nlyr+1,tau_table,x,&dx,ii);
            UTAU(lu) = splint_pchip(x,tau_table+ii,dx);
          }

          /*
           * Calculate the radiative transfer.
           */
          if (!uv->planck && !uv->scat_yes) {
            beer_law_only(ds,out,radius);
          }
          else {
            c_twostr(ds,out,deltam,gg,ierror,radius);
          }

          /*
           * Sum flux divergence; convert to heating rate below.
           */
          for (K = KLO; K < KHI; K++) {
            lu              = K;
            HEATING(K,J,I) += DFDT(lu)*DTAUC(lu);
          }
        } /* iwn loop */

        /*------------------------------------------------------------------*
         * Convert flux divergence [W/m^2] to specific heating rate [W/kg]. *
         *------------------------------------------------------------------*/

        for (K = KLO; K < KHI; K++) {
          HEATING(K,J,I) /= H3(K,J,I)*grid.dsgth[2*K+1];
        }

      }  /* I loop */
    }  /* J loop */
  }
  else if (action == EPIC_ALLOC) {
    /*
     * Allocate shortwave band memory and DISORT state and output memory.
     */
    nir = (rt_band *)calloc(1,sizeof(rt_band));
    if (!nir) {
      sprintf(Message,"error allocating nir");
      epic_error(dbmsname,Message);
    }

    vis = (rt_band *)calloc(1,sizeof(rt_band));
    if (!vis) {
      sprintf(Message,"error allocating vis");
      epic_error(dbmsname,Message);
    }

    uv  = (rt_band *)calloc(1,sizeof(rt_band));
    if (!uv) {
      sprintf(Message,"error allocating uv");
      epic_error(dbmsname,Message);
    }

    ds = (disort_state *)calloc(1,sizeof(disort_state));
    if (!ds) {
      sprintf(Message,"error allocating ds");
      epic_error(dbmsname,Message);
    }

    out = (disort_output *)calloc(1,sizeof(disort_output));
    if (!out) {
      sprintf(Message,"error allocating out");
      epic_error(dbmsname,Message);
    }

    /*
     * DISORT layers correspond to EPIC interior interfaces
     */
    ds->nlyr = grid.nk-1;

    /*
     * Allocate opacity arrays and read in shortwave data.
     */
    init_rt_shortwave(nir,vis,uv,ds);

    /*------------------*
     * c_twostr() setup *
     *------------------*/

    sprintf(ds->header,"EPIC Model (Version %4.2f): %s; %s",grid.epic_version,planet->name,dbmsname);

    /*
     * Limit to 2 streams for computational efficiency, and use c_twostr()
     * instead of c_disort (the former is optimized for two streams).
     *
     * NOTE: It is quick and easy to generalize this to multiple streams and c_disort().
     */
    ds->nstr = 2;

    /*
     * Set convergence criterion for azimuthal series, usually 0.
     */
    ds->accur = 0.;

    /*
     * Turn off diagnostic print statements and status messages.
     */
    ds->flag.prnt[0]=ds->flag.prnt[1]=ds->flag.prnt[2]=ds->flag.prnt[3]=ds->flag.prnt[4]=FALSE;
    ds->flag.quiet = TRUE;

    /*
     * Use original Nakajima-Tanaka intensity correction algorithm,
     *
     * NOTE: This needs to be reevaluated for the case of a strongly forward-scattering medium,
     *       see Buras, Dowling and Emde (2012).
     */
    ds->flag.old_intensity_correction = TRUE;

    /*
     * Indicate we do not want the special-case output of just albedo and
     * transmissivity of the entire medium.
     */
    ds->flag.ibcnd  = GENERAL_BC;

    /*
     * Indicate we may want to include thermal emission, so allocate relevant memory.
     */
    ds->flag.planck = TRUE;

    /*
     * Indicate that we do not need intensity output, only flux-related output.
     */
    ds->flag.onlyfl = TRUE;

    /*
     * Use internally calculated polar angles for returned radiant quantities,
     * not user-supplied ones.
     */
    ds->flag.usrang = FALSE;

    /*
     * Use an isotropically reflecting bottom boundary with zero albedo.
     */
    ds->flag.lamber = TRUE;
    ds->bc.albedo   = 0.;

    /*
     *  The primary inputs to DISORT are the differential optical depths at the DISORT computational
     *  layers, DTAUC(lc).  Since these depend on EPIC's mole fractions, X, which are carried
     *  on the EPIC layer interfaces, DISORT's computational layers are EPIC interfaces, running
     *  from lc = 1 at the kk=3 EPIC interface (bottom of top layer) to lc = nlyr = grid.nk-1
     *  at the kk=2*grid.nk-1 EPIC interface (top of the bottom layer).
     *
     *  EPIC requires heating on its layer interfaces, HEAT3.
     *  Hence, we indicate we want DISORT's dfdt = d(flux)/d(tau) output at those positions,
     *  by setting ds->flag.usrtau = true, ds->ntau = grid.nk-1, and assigning 
     *  UTAU(lu) to be EPIC's interior interface values.
     */
    ds->flag.usrtau = TRUE;
    ds->ntau        = grid.nk-1;

    /*
     * Number of phase functions for scattering calculations.
     * Needs to be >= to ds->nstr; often set equal to ds->nstr.
     */
    ds->nmom = ds->nstr;

    /*
     * Set azimuth angle of incident beam to zero, since we are not sensitive to it.
     */
    ds->bc.phi0 = 0.;

    /*
     * Set intensity of top-boundary isotropic illumination to zero.
     */
    ds->bc.fisot = 0.; 

    /*
     * Set spher to TRUE if there is any chance it will be TRUE during the run,
     * so that ds->zd[lev] is allocated.
     */
    ds->flag.spher = TRUE;

    /*
     * Allocate DISORT arrays.
     */
    c_twostr_state_alloc(ds);
    c_twostr_out_alloc(ds,out);

    /*
     * Allocate related vectors not included in DISORT structures.
     */
    gg        = dvector( 0,ds->nlyr-1,dbmsname);
    tauc      = dvector( 0,ds->nlyr,  dbmsname);
    tau_table = ftriplet(0,ds->nlyr,  dbmsname);

    /*
     * Set independent variable for spline on tau data,
     * used to specify user-defined optical depth positions, UTAU,
     * on EPIC interfaces.
     */
    for (lev = 0; lev <= ds->nlyr; lev++) {
      kk               = 2*(lev+1);
      tau_table[lev].x = (EPIC_FLOAT)kk;
    }

    deltam = FALSE;
    for (lc = 1; lc <= ds->nlyr; lc++) {
      /*
       * GG is the asymmetry factor; 1: complete forward scattering
       *                             0: isotropic scattering, Rayleigh scattering
       *                            -1: complete backscattering
       */
      GG(lc) = 0.;

      if (GG(lc) >= 0.95) {
        deltam = TRUE;
      }
    }

    if (ds->nstr > 2) {
      for (lc = 1; lc <= ds->nlyr; lc++) {
        c_getmom(RAYLEIGH,GG(lc),ds->nmom,&PMOM(0,lc));
      }
    }

    /*----------------------*
     * End c_twostr() setup *
     *----------------------*/
  }
  else if (action == EPIC_FREE) {
    /*
     * Free memory.
     */
    if (gg)        free_dvector(gg,        0,ds->nlyr-1,dbmsname);
    if (tauc)      free_dvector(tauc,      0,ds->nlyr,  dbmsname);
    if (tau_table) free_ftriplet(tau_table,0,ds->nlyr,  dbmsname);

    if (nir || vis || uv) free_rt_shortwave(nir,vis,uv,ds);

    if (nir) free(nir);
    if (vis) free(vis);
    if (uv)  free(uv);

    if (ds && out) c_twostr_out_free(ds,out);
    if (ds)        c_twostr_state_free(ds);

    if (ds)  free(ds);
    if (out) free(out);
  }
  else {
    sprintf(Message,"unrecognized action=%d",action);
    epic_error(dbmsname,Message);
  }

  return;
}

/*============ end of rt_shortwave() ==============================*/

/*============ init_rt_shortwave() ================================*/

/*
 * Allocate arrays and initialize absorption and scattering data, 
 * as appropriate for the near-infrared, visible and ultraviolet bands.
 *
 * Fortran subroutines:
 *   irwin_rewrite()
 *   initvisanduvver3()
 *   inith2h2opacver2()
 *   solarinitver3()
 *
 * NOTE: "Shortwave" and "longwave" are convenient partitions that
 *       are used in meteorological radiative transfer to be synonymous
 *       with "solar insolation" and "planetary emission", but this is a
 *       temperature-sensitive distinction.  Our current radiative-transfer
 *       partition, with shortwave = nir, vis, uv, and longwave = midir, fir,
 *       is appropriate for the outer solar system. However, hot planets like
 *       Venus and hot-Jupiters have significant emission at wavelengths
 *       shorter than midir, and for them we will need to revisit this
 *       partition.
 */

#define WAVE(icia)        wave[icia-1]
#define OPAC(itemp,icia)  opac[icia-1+(itemp-1)*ncia]
#define METHANE(iopac)    methane[iopac-1]
#define VAC(iopac)        vac[iopac-1]
/*
 * Specify NIR wn-bin merge index (1 => no effect)
 */
#define NIR_MERGE 8

#define BAND_FRACTION_THRESHOLD 0.01

void init_rt_shortwave(rt_band      *nir,
                       rt_band      *vis,
                       rt_band      *uv,
                       disort_state *ds)
{
  int
    iwn,itemp,ipress,ig,
    iopac,nopac,
    ibin,bin_size,
    icia,start_icia,ncia,
    iflux,nflux,
    imerge,itmp,icount,
    K,kk;
  int
    nc_err,nc_id,nc_dimid,nc_varid,
    full_nir_nwn;
  size_t
    nc_len,
   *start;
  double
    tmp,band_fraction,
    wavelength_start,wavelength_end;
  EPIC_FLOAT
    lambda_bot,lambda_top,tol;
  double
   *wave,
   *opac,
   *methane,
   *vac;
  char
    header[N_STR];
  FILE
   *infile;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="init_rt_shortwave";

  /*--------------------------------------------------------------------*
   * NIR (near infrared, 2000 to 9500 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Print setup progress to stdout */
  if (IAMNODE == NODE0) {
    fprintf(stdout,", nir"); fflush(stdout);
  }

  nir->wn_start = 2000.;
  nir->wn_end   = 9500.;

  /*
   * Determine whether thermal radiation by the atmosphere itself 
   * should be included.
   */
  wavelength_start = 1./(100.*nir->wn_end  );
  wavelength_end   = 1./(100.*nir->wn_start);
  nir->planck      = FALSE;
  for (K = 0; K <= KHI; K++) {
    kk = 2*K+1;
    /*
     * Calculate the band-emission fraction for the reference profile
     * of temperature.
     */
    band_fraction = blackbody_fraction(wavelength_end,  grid.t_ref[kk])
                   -blackbody_fraction(wavelength_start,grid.t_ref[kk]);
    if (band_fraction >= BAND_FRACTION_THRESHOLD) {
      nir->planck = TRUE;
      break;
    }
  }

  /*
   * Input Irwin et al file of near-infrared CH4 (methane) k-coefs.  The file consists
   * of nir->k_nu->nt sections, each assigned a unique temperature.
   * In each temperature section there are nir->k_nu->np pressures, and each pressure
   * is associated with nir->k_nu->ng k-coefficents.
   *
   * Irwin PGJ, Sromovsky LA, Strong EK, Sihra K, Teanby NA, Bowles N, Calcutt SB,
   *   Remedios JJ, 2006, Improved near-infrared methane band models and k-distribution
   *   parameters from 2000 to 9500 cm-1 and implications for interpretation of outer
   *   planet spectra.  Icarus 181, 309-319.
   */

  /* Allocate memory for k_nu */
  nir->k_nu = (rt_k_nu *)calloc(1,sizeof(rt_k_nu));
  if (!nir->k_nu) {
    sprintf(Message,"error allocating nir->k_nu");
    epic_error(dbmsname,Message);
  }

  nir->k_nu->on = FALSE;
  if (var.species[CH_4_INDEX].on == TRUE) {
    nir->k_nu->on = TRUE;
  }

  /* 
   * Open NIR CH4 absorption data file
   *
   * NOTE: We currently get NIR_WN from this file, so it needs to be
   *       read even if methane is not invoked by the user.
   */
  sprintf(header,EPIC_PATH"/data/rt/absorption/nir/irwinH_iii.par");

  infile = fopen(header,"r");
  if (!infile) {
    sprintf(Message,"unable to open %s\n",header);
    epic_error(dbmsname,Message);
  }

  fgets(header,N_STR,infile);

  fscanf(infile,"%*s %*s %d %*lf",&full_nir_nwn);
  /*
   * Downsize nir->nwn based on value of NIR_MERGE
   */
  nir->nwn = full_nir_nwn/NIR_MERGE;
  if (nir->nwn*NIR_MERGE < full_nir_nwn) {
    nir->nwn++;
  }

  fgets(header,N_STR,infile);

  /* Read number of g weights */
  fgets(header,N_STR,infile);
  fscanf(infile,"%*s %*s %*s %d",&nir->k_nu->ng);
  fgets(header,N_STR,infile);

  /* Allocate memory for g weights */
  nir->k_nu->delg = fvector(0,nir->k_nu->ng-1,dbmsname);

  /* Read in g weights */
  for (ig = 1; ig <= nir->k_nu->ng+1; ig++) {
    fgets(header,N_STR,infile);
  }
  for (ig = 1; ig <= nir->k_nu->ng; ig++) {
    /* g ordinates, the k-coeff statistical weights */
    fscanf(infile,"%lf",&NIR_K_NU_DELG(ig));
  }

  /* Read number of pressures */
  fscanf(infile,"%*s %*s %*s %*s %d",&nir->k_nu->np);
  fgets(header,N_STR,infile);

  /* Allocate memory for pressures */
  nir->k_nu->p = fvector(0,nir->k_nu->np-1,dbmsname);

  /* Read in pressures */
  for (ipress = 1; ipress <= nir->k_nu->np; ipress++) {
    fscanf(infile,"%lf",&tmp);
    /* Convert [atm] to [Pa] */
    NIR_K_NU_P(ipress) = tmp*1.01325e+5;
  }

  /* Read number of temperatures */
  fscanf(infile,"%*s %*s %*s %*s %d",&nir->k_nu->nt);
  fgets(header,N_STR,infile);
  fgets(header,N_STR,infile);

  /* Allocate memory for temperatures */
  nir->k_nu->t = fvector(0,nir->k_nu->nt-1,dbmsname);

  /* Read in temperatures */
  for (itemp = 1; itemp <= nir->k_nu->nt; itemp++) {
    fscanf(infile,"%lf",&NIR_K_NU_T(itemp));
  }
  for (itemp = 1; itemp <= 14; itemp++) {
    fgets(header,N_STR,infile);
  }

  /* Allocate memory for k_nu CH4 data */
  nir->k_nu->ch4 = fvector(0,nir->nwn
                            *nir->k_nu->ng
                            *nir->k_nu->np
                            *nir->k_nu->nt-1,dbmsname);

  /* Allocate memory for wavenumbers */
  nir->wn   = fvector(0,nir->nwn-1,dbmsname);
  nir->wnlo = fvector(0,nir->nwn-1,dbmsname);
  nir->wnhi = fvector(0,nir->nwn-1,dbmsname);

  icount = 0;
  for (iwn = 1; iwn <= nir->nwn; iwn++) {
    itmp = 0;
    for (imerge = 1; imerge <= NIR_MERGE; imerge++) {
      icount++;
      if (icount > full_nir_nwn) {
        /*
         * Do not go beyond the full input-file nwn data count.
         */
        break;
      }

      /*
       * The last merged wn bin can have fewer contributions than NIR_MERGE,
       * so we keep track via itmp.
       */
      itmp++;

      fscanf(infile,"%*s %*s %lf",&tmp);
      NIR_WN(iwn) += tmp;

      for (itemp = 1; itemp <= nir->k_nu->nt; itemp++) {
        for (ipress = 1; ipress <= 6; ipress++) {
          fgets(header,N_STR,infile);
        }
        for (ipress = 1; ipress <= nir->k_nu->np; ipress++) {
          fscanf(infile,"%*lf");
          for (ig = 1; ig <= nir->k_nu->ng/2; ig++) {
            fscanf(infile,"%lf",&tmp);
            NIR_K_NU_CH4(ig,ipress,itemp,iwn) += tmp;
          }
          fscanf(infile,"%*lf");
          for (ig = nir->k_nu->ng/2+1; ig <= nir->k_nu->ng; ig++) {
            fscanf(infile,"%lf",&tmp);
            NIR_K_NU_CH4(ig,ipress,itemp,iwn) += tmp;
          }
          fgets(header,N_STR,infile);
        }
        fgets(header,N_STR,infile);
      }
    }

    /* 
     * Assign averaged values in merged bin.
     */
    NIR_WN(iwn) /= itmp;
    for (itemp = 1; itemp <= nir->k_nu->nt; itemp++) {
      for (ipress = 1; ipress <= nir->k_nu->np; ipress++) {
        for (ig = 1; ig <= nir->k_nu->ng; ig++) {
          NIR_K_NU_CH4(ig,ipress,itemp,iwn) /= itmp;
        }
      }
    }
  }

  NIR_WNLO(1) = NIR_WN(1)-.5*(NIR_WN(2)-NIR_WN(1));
  NIR_WNHI(1) = NIR_WN(1)+.5*(NIR_WN(2)-NIR_WN(1));
  for (iwn = 2; iwn <= nir->nwn; iwn++) {
    NIR_WNLO(iwn) = NIR_WNHI(iwn-1);
    NIR_WNHI(iwn) = 2.*NIR_WN(iwn)-NIR_WNLO(iwn);
  }

  /* Close data file */
  fclose(infile);

  /*
   * NIR collision-induced absorption data (cia) for hydrogen (H_2-H_2).
   */

  /* Allocate memory for cia */
  nir->cia = (rt_cia *)calloc(1,sizeof(rt_cia));
  if (!nir->cia) {
    sprintf(Message,"error allocating nir->cia");
    epic_error(dbmsname,Message);
  }

  nir->cia->on = TRUE;

  /*
   * Open NIR cia data file.
   */
  sprintf(header,EPIC_PATH"/data/rt/absorption/nir/H2_H2_CIA_Borysow2002.dat");
  infile = fopen(header,"r");
  if (!infile) {
    sprintf(Message,"unable to open %s",header);
    epic_error(dbmsname,Message);
  }

  /* Read number of temperatures */
  for (itemp = 1; itemp <= 5; itemp++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s %*s %d",&nir->cia->nt);

  /* Read number of wavenumbers */
  fscanf(infile,"%*s %*s %d",&ncia);
  fgets(header,N_STR,infile);

  /* Allocate memory for temperatures */
  nir->cia->t = fvector(0,nir->cia->nt-1,dbmsname);

  /* Read in temperatures */
  for (itemp = 1; itemp <= 4; itemp++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s");
  for (itemp = 1; itemp <= nir->cia->nt; itemp++) {
    fscanf(infile,"%lfK",&NIR_CIA_T(itemp));
  }
  fgets(header,N_STR,infile);
  fgets(header,N_STR,infile);

  /* Allocate memory for input cia data */
  wave = fvector(0,ncia-1,dbmsname);
  opac = fvector(0,ncia*nir->cia->nt-1,dbmsname);

  /* Allocate memory for H2-H2 cia data */
  nir->cia->h2_h2 = fvector(0,nir->nwn*nir->cia->nt-1,dbmsname);

  for (icia = 1; icia <= ncia; icia++) {
    fscanf(infile,"%lf",&WAVE(icia));
    for (itemp = 1; itemp <= nir->cia->nt; itemp++) {
      fscanf(infile,"%lf",&OPAC(itemp,icia));
      /*
       * Convert [1/(cm amagat^2)] to [1/(km amagat^2)]
       */
      OPAC(itemp,icia) *= 1.e+5;
    }
  }

  start_icia = 1;
  for (iwn = 1; iwn <= nir->nwn; iwn++) {
    /*
     * Find place in table.
     * The two wavenumber columns are ordered in the same direction.
     */
    if (NIR_WN(iwn) < WAVE(1) || NIR_WN(iwn) > WAVE(ncia)) {
      /*
       * Assign zero opacity for wavenumbers that fall off data table.
       */
      for (itemp = 1; itemp <= nir->cia->nt; itemp++) {
        NIR_CIA_H2_H2(itemp,iwn) = 0.;
      }
    }
    else {
      for (icia = start_icia; icia <= ncia; icia++) {
        if (NIR_WN(iwn) <= WAVE(icia+1)) {
          /*
           * Use linear interpolation; all temperatures can be done
           * at the given wavelength. 
           */
          tmp = (NIR_WN(iwn)-WAVE(icia))/(WAVE(icia+1)-WAVE(icia));
          for (itemp = 1; itemp <= nir->cia->nt; itemp++) {
            NIR_CIA_H2_H2(itemp,iwn) = OPAC(itemp,icia)+tmp*(OPAC(itemp,icia+1)-OPAC(itemp,icia));
          }
          start_icia = icia;
          break;
        }
      }
    }
  }

  /* Free allocated memory */
  free_fvector(wave,0,ncia-1,            dbmsname);
  free_fvector(opac,0,ncia*nir->cia->nt-1,dbmsname);

  /* Close data file */
  fclose(infile);

  /*
   * Initialize NIR solar flux.
   * The same input file is used for all shortwave bands.
   */

  /* Allocate memory for NIR solar flux */
  nir->fluxtot = fvector(0,nir->nwn-1,dbmsname);

  /* Fractional tolerance for Romberg integration of solar flux. */
  tol = 1.e-8;
  for (iwn = 1; iwn <= nir->nwn; iwn++) {
    /*
     * The function solar_irradiance() inputs wavelength, lambda, in [micron].
     * Converting [cm-1] to [micron].
     */
    lambda_bot = 1.e+4/NIR_WNHI(iwn);
    lambda_top = 1.e+4/NIR_WNLO(iwn);
    NIR_FLUXTOT(iwn) = romberg_integral(solar_irradiance,lambda_bot,lambda_top,tol);
  }

  /*
   * Allocate memory for dtauc
   */
  nir->dtauc = fvector(0,nir->nwn*nir->k_nu->ng*ds->nlyr-1,dbmsname);

  /*--------------------------------------------------------------------*
   * VIS (visible, 10000 to 40000 cm-1)                                 *
   *--------------------------------------------------------------------*/

  /* Print setup progress to stdout */
  if (IAMNODE == NODE0) {
    fprintf(stdout,", vis"); fflush(stdout);
  }

  vis->wn_start = 10000.;
  vis->wn_end   = 40000.;

  /*
   * Determine whether thermal radiation by the atmosphere itself 
   * should be included.
   */
  wavelength_start = 1./(100.*vis->wn_end  );
  wavelength_end   = 1./(100.*vis->wn_start);
  vis->planck      = FALSE;
  for (K = 0; K <= KHI; K++) {
    kk = 2*K+1;
    /*
     * Calculate the band-emission fraction for the reference profile
     * of temperature.
     */
    band_fraction = blackbody_fraction(wavelength_end,  grid.t_ref[kk])
                   -blackbody_fraction(wavelength_start,grid.t_ref[kk]);
    if (band_fraction >= BAND_FRACTION_THRESHOLD) {
      vis->planck = TRUE;
      break;
    }
  }

  /*
   * Input VIS CH4 absorption coefficient data.
   */
  sprintf(header,EPIC_PATH"/data/rt/absorption/vis/1995lo_karko.tab");
  infile = fopen(header,"r");
  if (!infile) {
    sprintf(Message,"unable to open %s",header);
    epic_error(dbmsname,Message);
  }

  for (iopac = 1; iopac <= 13; iopac++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s %*s %d",&nopac);
  fgets(header,N_STR,infile);
  fgets(header,N_STR,infile);

  /* Allocate memory */
  vac     = fvector(0,nopac-1,dbmsname);
  methane = fvector(0,nopac-1,dbmsname);

  /*
   * Set bin size
   *
   * NOTE: If bin_size is set less than 4, then the 3 leading "CH4=.0000" entries
   *       in 1995lo_karko.tab should be deleted, to avoid any zero-opacity bins (if
   *       opacity is zero but Rayleigh scattering is on, c_twostr() fails due to
   *       an ill-conditioned matrix, stemming from SSALB = 1.0 throughout the column).
   */
  bin_size = 10;

  vis->nwn = nopac/bin_size;

  /* Allocate memory for wavenumbers */
  vis->wn   = fvector(0,vis->nwn-1,dbmsname);
  vis->wnlo = fvector(0,vis->nwn-1,dbmsname);
  vis->wnhi = fvector(0,vis->nwn-1,dbmsname);

  for (iopac = nopac; iopac >= 1; iopac--) {
    fscanf(infile,"%lf %*lf %lf %*lf %*lf %*lf %*lf %*lf",&VAC(iopac),&METHANE(iopac));
  }

  /* Allocate memory for methane opacity */
  vis->k_nu = (rt_k_nu *)calloc(1,sizeof(rt_k_nu));
  if (!vis->k_nu) {
    sprintf(Message,"error allocating vis->k_nu"); epic_error(dbmsname,Message);
  }

  vis->k_nu->on  = TRUE;

  vis->k_nu->ch4 = fvector(0,vis->nwn-1,dbmsname);

  /*
   * Read in methane opacity data
   */
  for (iwn = 1; iwn <= vis->nwn; iwn++) {
    VIS_K_NU_CH4(iwn) = 0.;
    for (ibin = 1; ibin <= bin_size; ibin++) {
      /* [1/(km amagat)] */
      VIS_K_NU_CH4(iwn) += METHANE((iwn-1)*bin_size+ibin);
    }
    VIS_K_NU_CH4(iwn) /= bin_size;

    /* 1.e-7 = [cm/nm] */
    VIS_WNLO(iwn) = 1./((VAC((iwn-1)*bin_size+1         )+0.2)*1.e-7);
    VIS_WNHI(iwn) = 1./((VAC((iwn-1)*bin_size+bin_size  )-0.2)*1.e-7);
    VIS_WN(iwn)   = 1./((VAC((iwn-1)*bin_size+bin_size/2)-0.2)*1.e-7);
  }

  /* Free allocated memory */
  free_fvector(vac,    0,nopac-1,dbmsname);
  free_fvector(methane,0,nopac-1,dbmsname);

  /* Close data file */
  fclose(infile);

  /*
   * Calculate VIS Rayleigh scattering cross sections.
   */

  /* Allocate memory for scattering cross sections */
  vis->sigma = (rt_sigma *)calloc(1,sizeof(rt_sigma));
  if (!vis->sigma) {
    sprintf(Message,"error allocating vis->sigma"); epic_error(dbmsname,Message);
  }

  /*
   * Turn Rayleigh scattering on or off (vis->sigma->on = TRUE or FALSE) for VIS band.
   *
   * NOTE: Currently turned off, because the case with zero opacity in a column plus scattering
   *       yields a column of SSALB = 1, which generates an ill-conditioned matrix
   *       in cdisort.c.
   */
  vis->sigma->on = FALSE;

  if (vis->sigma->on) {
    vis->sigma->h2  = fvector(0,vis->nwn-1,dbmsname);
    vis->sigma->he  = fvector(0,vis->nwn-1,dbmsname);
    vis->sigma->ch4 = fvector(0,vis->nwn-1,dbmsname);

    cross_section_rayleigh(vis);
  }

  /*
   * Collision-induced absorption (cia) for hydrogen (H2-H2).
   *
   * Data from
   *   Borysow A, 2002, Collision-induced absorption coefficients
   *     of H2 pairs at temperatures from 60 K to 1000 K,
   *     Astron. Astrophys. 390, 779-782, doi: 10.1051/0004-6361:20020555
   *
   * The data were obtained in August 2011 from the link
   *   http://www.astro.ku.dk/~aborysow/programs/final_CIA_LT.dat
   */

  /* Allocate memory for cia H2-H2 */
  vis->cia = (rt_cia *)calloc(1,sizeof(rt_cia));
  if (!vis->cia) {
    sprintf(Message,"error allocating vis->cia");
    epic_error(dbmsname,Message);
  }

  vis->cia->on = TRUE;

  sprintf(header,EPIC_PATH"/data/rt/absorption/vis/H2_H2_CIA_Borysow2002.dat");
  infile = fopen(header,"r");
  if (!infile) {
    sprintf(Message,"unable to open %s",header);
    epic_error(dbmsname,Message);
  }

  /* Read number of temperatures */
  for (itemp = 1; itemp <= 5; itemp++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s %*s %d",&vis->cia->nt);

  /* Read number of wavenumbers */
  fscanf(infile,"%*s %*s %d",&ncia);
  fgets(header,N_STR,infile);

  /* Allocate memory for temperatures */
  vis->cia->t = fvector(0,vis->cia->nt-1,dbmsname);

  /* Read in temperatures */
  for (itemp = 1; itemp <= 4; itemp++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s");
  for (itemp = 1; itemp <= vis->cia->nt; itemp++) {
    fscanf(infile,"%lfK",&VIS_CIA_T(itemp));
  }
  fgets(header,N_STR,infile);
  fgets(header,N_STR,infile);

  /* Allocate memory for input cia data */
  wave = fvector(0,ncia-1,dbmsname);
  opac = fvector(0,ncia*nir->cia->nt-1,dbmsname);

  /* Allocate memory for H2-H2 cia data */
  vis->cia->h2_h2 = fvector(0,vis->nwn*vis->cia->nt-1,dbmsname);

  for (icia = 1; icia <= ncia; icia++) {
    fscanf(infile,"%lf",&WAVE(icia));
    for (itemp = 1; itemp <= vis->cia->nt; itemp++) {
      fscanf(infile,"%lf",&OPAC(itemp,icia));
      /*
       * Convert [1/(cm amagat^2)] to [1/(km amagat^2)]
       */
      OPAC(itemp,icia) *= 1.e+5;
    }
  }

  start_icia = 1;
  for (iwn = 1; iwn <= vis->nwn; iwn++) {
    /*
     * Find place in table.
     * The two wavenumber columns are ordered in the same direction.
     */
    if (VIS_WN(iwn) < WAVE(1) || VIS_WN(iwn) > WAVE(ncia)) {
      /*
       * Assign zero opacity for wavenumbers that fall off data table.
       */
      for (itemp = 1; itemp <= vis->cia->nt; itemp++) {
        VIS_CIA_H2_H2(itemp,iwn) = 0.;
      }
    }
    else {
      for (icia = start_icia; icia <= ncia; icia++) {
        if (VIS_WN(iwn) <= WAVE(icia+1)) {
          /*
           * Use linear interpolation; all temperatures can be done
           * at the given wavelength. 
           */
          tmp = (VIS_WN(iwn)-WAVE(icia))/(WAVE(icia+1)-WAVE(icia));
          for (itemp = 1; itemp <= vis->cia->nt; itemp++) {
            VIS_CIA_H2_H2(itemp,iwn) = OPAC(itemp,icia)+tmp*(OPAC(itemp,icia+1)-OPAC(itemp,icia));
          }
          start_icia = icia;
          break;
        }
      }
    }
  }

  /* Free allocated memory */
  free_fvector(wave,0,ncia-1,            dbmsname);
  free_fvector(opac,0,ncia*vis->cia->nt-1,dbmsname);

  fclose(infile);

  /*
   * Initialize VIS solar flux.
   */

  /* Allocate memory for VIS solar flux */
  vis->fluxtot = fvector(0,vis->nwn-1,dbmsname);

  /* Fractional tolerance for Romberg integration of solar flux. */
  tol = 1.e-8;
  for (iwn = 1; iwn <= vis->nwn; iwn++) {
    /*
     * The function solar_irradiance() inputs wavelength, lambda, in [micron].
     * Converting [cm-1] to [micron].
     */
    lambda_bot = 1.e+4/VIS_WNHI(iwn);
    lambda_top = 1.e+4/VIS_WNLO(iwn);
    VIS_FLUXTOT(iwn) = romberg_integral(solar_irradiance,lambda_bot,lambda_top,tol);
  }

  /*
   * Allocate memory for dtauc, ssalb
   */
  vis->dtauc = fvector(0,vis->nwn*ds->nlyr-1,dbmsname);

  if (vis->sigma->on) {
    vis->ssalb = fvector(0,vis->nwn*ds->nlyr-1,dbmsname);
  }

  /*--------------------------------------------------------------------*
   * UV (ultraviolet, 40000 to 100000 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Print setup progress to stdout */
  if (IAMNODE == NODE0) {
    fprintf(stdout,", uv"); fflush(stdout);
  }

  uv->wn_start =  40000.;
  uv->wn_end   = 100000.;

  /*
   * Determine whether thermal radiation by the atmosphere itself 
   * should be included.
   */
  wavelength_start = 1./(100.*uv->wn_end  );
  wavelength_end   = 1./(100.*uv->wn_start);
  uv->planck       = FALSE;
  for (K = 0; K <= KHI; K++) {
    kk = 2*K+1;
    /*
     * Calculate the band-emission fraction for the reference profile
     * of temperature.
     */
    band_fraction = blackbody_fraction(wavelength_end,  grid.t_ref[kk])
                   -blackbody_fraction(wavelength_start,grid.t_ref[kk]);
    if (band_fraction >= BAND_FRACTION_THRESHOLD) {
      uv->planck = TRUE;
      break;
    }
  }

  /*
   * Ultraviolet photoabsorption cross sections, as used by
   *
   *   Moses JI, Bezard B, Lellouch E, Gladstone GR, Feuchtgruber H,
   *     Allen M, 2000, Photochemistry of Saturn's Atmosphere, Icarus 143, 244-298.
   *
   * See the header of the input file epic/data/rt/absorption/uv/uv_cross_sections.inp
   * for a discussion and complete references.
   */

  sprintf(header,EPIC_PATH"/data/rt/absorption/uv/uv_cross_sections.inp");
  infile = fopen(header,"r");
  if (!infile) {
    sprintf(Message,"unable to open %s",header);
    epic_error(dbmsname,Message);
  }

  /* Read number of wavenumbers */
  for (iwn = 1; iwn <= 58; iwn++) {
    fgets(header,N_STR,infile);
  }
  fscanf(infile,"%*s %*s %d",&uv->nwn);
  fgets(header,N_STR,infile);
  fgets(header,N_STR,infile);

  /* Allocate memory for wavenumbers */
  uv->wn   = fvector(0,uv->nwn-1,dbmsname);
  uv->wnlo = fvector(0,uv->nwn-1,dbmsname);
  uv->wnhi = fvector(0,uv->nwn-1,dbmsname);

  /* 
   * Allocate memory for opacity cross sections
   *
   * NOTE: Allocate these, even if EPIC user has not invoked them.
   */

  uv->k = (rt_k *)calloc(1,sizeof(rt_k));
  if (!uv->k) {
    sprintf(Message,"error allocating uv->k"); 
    epic_error(dbmsname,Message);
  }

  uv->k->on   = TRUE;

  uv->k->ch4  = fvector(0,uv->nwn-1,dbmsname);
  uv->k->c2h2 = fvector(0,uv->nwn-1,dbmsname);
  uv->k->c2h4 = fvector(0,uv->nwn-1,dbmsname);
  uv->k->c2h6 = fvector(0,uv->nwn-1,dbmsname);
  
  /* Read in opacity cross sections and wavenumbers */
  for (iwn = uv->nwn; iwn >= 1; iwn--) {
    fscanf(infile,"%lf %lf %lf %lf %lf",&tmp,&UV_K_CH4(iwn),&UV_K_C2H2(iwn),&UV_K_C2H4(iwn),&UV_K_C2H6(iwn));
    /* 1.e-8 = [cm/angstrom] */
    UV_WN(iwn) = 1./(tmp*1.e-8);
  }

  UV_WNLO(1) = UV_WN(1)-.5*(UV_WN(2)-UV_WN(1));
  UV_WNHI(1) = UV_WN(1)+.5*(UV_WN(2)-UV_WN(1));
  for (iwn = 2; iwn < uv->nwn; iwn++) {
    UV_WNLO(iwn) = .5*(UV_WN(iwn)+UV_WN(iwn-1));
    UV_WNHI(iwn) = .5*(UV_WN(iwn)+UV_WN(iwn+1));
  }
  UV_WNLO(uv->nwn) = UV_WN(uv->nwn)-.5*(UV_WN(uv->nwn)-UV_WN(uv->nwn-1));
  UV_WNHI(uv->nwn) = UV_WN(uv->nwn)+.5*(UV_WN(uv->nwn)-UV_WN(uv->nwn-1));

  fclose(infile);

  /*
   * Calculate UV Rayleigh scattering cross sections.
   */

  /* Allocate memory for scattering cross sections */
  uv->sigma = (rt_sigma *)calloc(1,sizeof(rt_sigma));
  if (!uv->sigma) {
    sprintf(Message,"error allocating uv->sigma");
    epic_error(dbmsname,Message);
  }

  /*
   * Turn Rayleigh scattering on or off (uv->sigma->on = TRUE or FALSE) for UV band.
   *
   * NOTE: Currently turned off, because the case with zero opacity in a column plus scattering
   *       yields a column of SSALB = 1, which generates an ill-conditioned matrix
   *       in cdisort.c.
   */

  uv->sigma->on =FALSE;

  if (uv->sigma->on) {
    uv->sigma->h2  = fvector(0,uv->nwn-1,dbmsname);
    uv->sigma->he  = fvector(0,uv->nwn-1,dbmsname);
    uv->sigma->ch4 = fvector(0,uv->nwn-1,dbmsname);

    cross_section_rayleigh(uv);
  }

  /*
   * Initialize UV solar flux
   */

  /* Allocate memory for UV solar flux */
  uv->fluxtot = fvector(0,uv->nwn-1, dbmsname);

  /* Fractional tolerance for Romberg integration of solar flux. */
  tol = 1.e-8;
  for (iwn = 1; iwn <= uv->nwn; iwn++) {
    /*
     * The function solar_irradiance() inputs wavelength, lambda, in [micron].
     * Converting [cm-1] to [micron].
     */
    lambda_bot = 1.e+4/UV_WNHI(iwn);
    lambda_top = 1.e+4/UV_WNLO(iwn);
    UV_FLUXTOT(iwn) = romberg_integral(solar_irradiance,lambda_bot,lambda_top,tol);
  }

  /*
   * Allocate memory for dtauc, ssalb
   */
  uv->dtauc = fvector(0,uv->nwn*ds->nlyr-1,dbmsname);

  if (uv->sigma->on) {
    uv->ssalb = fvector(0,uv->nwn*ds->nlyr-1,dbmsname);
  }

  return;
}

#undef WAVE
#undef OPAC
#undef METHANE
#undef VAC
#undef NIR_MERGE
#undef BAND_FRACTION_THRESHOLD

/*============ end of init_rt_shortwave() =========================*/

/*============ free_rt_shortwave() ================================*/

/*
 * Free memory dynamically allocated in init_rt_shortwave()
 */

void free_rt_shortwave(rt_band      *nir,
                       rt_band      *vis,
                       rt_band      *uv,
                       disort_state *ds)
{
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="free_rt_shortwave";

  /*--------------------------------------------------------------------*
   * NIR (near infrared, 2000 to 9500 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Free wavenumbers */
  if (nir->wn)   free_fvector(nir->wn,  0,nir->nwn-1,dbmsname);
  if (nir->wnlo) free_fvector(nir->wnlo,0,nir->nwn-1,dbmsname);
  if (nir->wnhi) free_fvector(nir->wnhi,0,nir->nwn-1,dbmsname);

  /* Free NIR solar flux */
  if (nir->fluxtot) free_fvector(nir->fluxtot,0,nir->nwn-1,dbmsname);

  /* Free dtauc */
  if (nir->dtauc) free_fvector(nir->dtauc,0,nir->nwn*nir->k_nu->ng*ds->nlyr-1,dbmsname);

  if (nir->k_nu) {
    /* Free k_nu */
    if (nir->k_nu->delg) free_fvector(nir->k_nu->delg,0,nir->k_nu->ng-1,dbmsname);
    if (nir->k_nu->p)    free_fvector(nir->k_nu->p,   0,nir->k_nu->np-1,dbmsname);
    if (nir->k_nu->t)    free_fvector(nir->k_nu->t,   0,nir->k_nu->nt-1,dbmsname);
    if (nir->k_nu->ch4)  free_fvector(nir->k_nu->ch4, 0,nir->nwn
                                     *nir->k_nu->ng
                                     *nir->k_nu->np
                                     *nir->k_nu->nt-1,dbmsname);
    free(nir->k_nu);
  }

  /* Free cia */
  if (nir->cia) {
    if (nir->cia->t)     free_fvector(nir->cia->t,    0,nir->cia->nt-1,dbmsname);
    if (nir->cia->h2_h2) free_fvector(nir->cia->h2_h2,0,nir->nwn*nir->cia->nt-1,dbmsname);
    free(nir->cia);
  }

  /*--------------------------------------------------------------------*
   * VIS (visible, 10000 to 40000 cm-1)                                 *
   *--------------------------------------------------------------------*/

  /* Free wavenumbers */
  if (vis->wn)   free_fvector(vis->wn,  0,vis->nwn-1,dbmsname);
  if (vis->wnlo) free_fvector(vis->wnlo,0,vis->nwn-1,dbmsname);
  if (vis->wnhi) free_fvector(vis->wnhi,0,vis->nwn-1,dbmsname);

  /* Free VIS solar flux */
  if (vis->fluxtot) free_fvector(vis->fluxtot,0,vis->nwn-1,dbmsname);

  /* Free dtauc, ssalb */
  if (vis->dtauc) free_fvector(vis->dtauc,0,vis->nwn*ds->nlyr-1,dbmsname);
  if (vis->ssalb) free_fvector(vis->ssalb,0,vis->nwn*ds->nlyr-1,dbmsname);

  /* Free k_nu */
  if (vis->k_nu) {
    if (vis->k_nu->ch4) free_fvector(vis->k_nu->ch4,0,vis->nwn-1,dbmsname);
    free(vis->k_nu);
  }

  /* Free sigma */
  if (vis->sigma) {
    if (vis->sigma->h2)  free_fvector(vis->sigma->h2, 0,vis->nwn-1,dbmsname);
    if (vis->sigma->he)  free_fvector(vis->sigma->he, 0,vis->nwn-1,dbmsname);
    if (vis->sigma->ch4) free_fvector(vis->sigma->ch4,0,vis->nwn-1,dbmsname);
    free(vis->sigma);
  }

  /* Free cia */
  if (vis->cia) {
    if (vis->cia->t)     free_fvector(vis->cia->t,    0,vis->cia->nt-1,dbmsname);
    if (vis->cia->h2_h2) free_fvector(vis->cia->h2_h2,0,vis->nwn*vis->cia->nt-1,dbmsname);
    free(vis->cia);
  }

  /*--------------------------------------------------------------------*
   * UV (ultraviolet, 40000 to 100000 cm-1)                             *
   *--------------------------------------------------------------------*/

  /* Free wavenumbers */
  if (uv->wn)   free_fvector(uv->wn,  0,uv->nwn-1,dbmsname);
  if (uv->wnlo) free_fvector(uv->wnlo,0,uv->nwn-1,dbmsname);
  if (uv->wnhi) free_fvector(uv->wnhi,0,uv->nwn-1,dbmsname);

  /* Free UV solar flux */
  if (uv->fluxtot) free_fvector(uv->fluxtot,0,uv->nwn-1,dbmsname);

  /* Free dtauc, ssalb */
  if (uv->dtauc) free_fvector(uv->dtauc,0,uv->nwn*ds->nlyr-1,dbmsname);
  if (uv->ssalb) free_fvector(uv->ssalb,0,uv->nwn*ds->nlyr-1,dbmsname);

  /* Free k */
  if (uv->k) {
    if (uv->k->ch4)  free_fvector(uv->k->ch4, 0,uv->nwn-1,dbmsname);
    if (uv->k->c2h2) free_fvector(uv->k->c2h2,0,uv->nwn-1,dbmsname);
    if (uv->k->c2h4) free_fvector(uv->k->c2h4,0,uv->nwn-1,dbmsname);
    if (uv->k->c2h6) free_fvector(uv->k->c2h6,0,uv->nwn-1,dbmsname);
    free(uv->k);
  }

  /* Free sigma */
  if (uv->sigma) {
    if (uv->sigma->h2)  free_fvector(uv->sigma->h2, 0,uv->nwn-1,dbmsname);
    if (uv->sigma->he)  free_fvector(uv->sigma->he, 0,uv->nwn-1,dbmsname);
    if (uv->sigma->ch4) free_fvector(uv->sigma->ch4,0,uv->nwn-1,dbmsname);
    free(uv->sigma);
  }

  return;
}

/*============ end of free_rt_shortwave() =========================*/

/*============= nir_opacity() =====================================*/

/*
 * Interpolate in pressure and temperature to obtain near-infrared
 * opacity due to methane, using NIR_K_NU_CH4.
 *
 * Fortran subroutines:
 *   nir_irver4()
 *   heatingver12()
 */

#define LOGPTAU(ipress) logp_table[ipress-1].x
#define TTAU(itemp)     t_table[itemp-1].x
#define TCIA(itemp)     tcia_table[itemp-1].x

void nir_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *nir,
                 disort_state *ds)
{
  int
    lc,
    ipress,itemp,ig,iwn,
    K,lev,
    test;
  double
    dt,dlogp,logp,temp,
    m,b,
    mu,dz_km,
    logdtau1,logdtau2,
    dtau,h2_h2_factor,
    n_amagat;
  static int
    initialized = FALSE;
  static float_triplet
   *t_table,
   *logp_table,
   *tcia_table;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="nir_opacity";

  if (!initialized) {
    if (nir->k_nu->on) {
      /* Allocate memory */
      t_table    = ftriplet(0,nir->k_nu->nt-1,dbmsname);
      logp_table = ftriplet(0,nir->k_nu->np-1,dbmsname);
      tcia_table = ftriplet(0,nir->cia->nt-1, dbmsname);

      /* Assign values */
      for (ipress = 1; ipress <= nir->k_nu->np; ipress++) {
        LOGPTAU(ipress) = log(NIR_K_NU_P(ipress));
      }
      for (itemp = 1; itemp <= nir->k_nu->nt; itemp++) {
        TTAU(itemp) = NIR_K_NU_T(itemp);
      }
      for (itemp = 1; itemp <= nir->cia->nt; itemp++) {
        TCIA(itemp) = NIR_CIA_T(itemp);
      }
    }

    initialized = TRUE;
  }

  ds->flag.planck = nir->planck;

  if (ds->flag.planck) {
    ds->bc.temis = 1.;   /* emissivity at the top boundary */
    ds->bc.ttemp = T3(0,J,I);
    for (K = KLO; K <= KHI; K++) {
      lev         = K-1;
      TEMPER(lev) = T2(K,J,I);
    }
    ds->bc.btemp = T3(grid.nk,J,I);
  }

  /* 
   * Zero nir->dtauc
   */
  memset(nir->dtauc,0,nir->nwn*nir->k_nu->ng*ds->nlyr*sizeof(double));

  for (lc = 1; lc <= ds->nlyr; lc++) {
    K = lc;

    /*
     * Calculate number density [amagat]
     */
    n_amagat = (P3(K,J,I)*273.15)/(T3(K,J,I)*1.01325e+5);

    /*
     * DISORT layer thickness [km]
     */
    dz_km = (Z2(K,J,I)-Z2(K+1,J,I))/1.e+3;

    if (nir->k_nu->on) {
      logp   = LIMIT_RANGE(LOGPTAU(1),log(P3(K,J,I)),LOGPTAU(nir->k_nu->np));
      ipress = hunt_place_in_table(nir->k_nu->np,logp_table,logp,&dlogp,ipress);
      ipress++;  /* Shift to unit-based */

      temp  = LIMIT_RANGE(TTAU(1),T3(K,J,I),TTAU(nir->k_nu->nt));
      itemp = hunt_place_in_table(nir->k_nu->nt,t_table,temp,&dt,itemp);
      itemp++;   /* Shift to unit-based */

      /*
       * Layout of T and P in NIR_K_NU_CH4
       *
       *   NIR_K_NU_CH4(ig,ipress,  itemp,  iwn)   lo temp, lo P
       *   NIR_K_NU_CH4(ig,ipress,  itemp+1,iwn)   hi temp, lo P
       *   NIR_K_NU_CH4(ig,ipress+1,itemp+1,iwn)   hi temp, hi P
       *   NIR_K_NU_CH4(ig,ipress+1,itemp,  iwn)   lo temp, hi P
       */

      for (ig = 1; ig <= nir->k_nu->ng; ig++) {
        for (iwn = 1; iwn <= nir->nwn; iwn++) {
          /* 
           * Low-pressure temperature interpolation.
           * Opacities vary approximately linearly with temperature.
           */
          m        = (NIR_K_NU_CH4(ig,ipress,itemp+1,iwn)-NIR_K_NU_CH4(ig,ipress,itemp,iwn))/dt;
          b        = NIR_K_NU_CH4(ig,ipress,itemp,iwn)-m*NIR_K_NU_T(itemp);
          logdtau1 = log(m*temp+b);

          /* 
           * High-pressure temperature interpolation.
           */
          m        = (NIR_K_NU_CH4(ig,ipress+1,itemp+1,iwn)-NIR_K_NU_CH4(ig,ipress+1,itemp,iwn))/dt;
          b        = NIR_K_NU_CH4(ig,ipress+1,itemp,iwn)-m*NIR_K_NU_T(itemp);
          logdtau2 = log(m*temp+b);

          /*
           * Pressure interpolation.
           * Log opacity varies approximately linearly with log pressure.
           */
          m    = (logdtau2-logdtau1)/dlogp;
          b    = logdtau1-m*LOGPTAU(ipress);
          dtau = exp(m*logp+b);

          /*
           * Ensure non-negative dtauc.
           */
          NIR_DTAUC(lc,ig,iwn) += MAX(dtau*X(CH_4_INDEX,VAPOR,K,J,I)*n_amagat*dz_km,0.);
        }
      }
    }

    if (nir->cia->on) {
      /*
       * H_2-H_2 collisionally induced absorption.
       *
       * NOTE: For Uranus and Neptune, may need to extend data table to colder temperatures.
       * NOTE: For deep (warm) layers, may need to extend data table to warmer temperatures.
       */
      temp  = LIMIT_RANGE(TCIA(1),T3(K,J,I),TCIA(nir->cia->nt));
      itemp = hunt_place_in_table(nir->cia->nt,tcia_table,temp,&dt,itemp);
      itemp++;  /* Shift to unit-based */

      /*
       * Density times column density [km almagat^2]
       */
      h2_h2_factor = SQR(planet->x_h2*n_amagat)*dz_km;

      for (iwn = 1; iwn <= nir->nwn; iwn++) {
        m    = (NIR_CIA_H2_H2(itemp+1,iwn)-NIR_CIA_H2_H2(itemp,iwn))/dt;
        b    = NIR_CIA_H2_H2(itemp,iwn)-m*NIR_CIA_T(itemp);
        dtau = (m*temp+b)*h2_h2_factor;
        for (ig = 1; ig <= nir->k_nu->ng; ig++) {
          NIR_DTAUC(lc,ig,iwn) += dtau;
        }
      }
    }
  }

  return;
}

#undef LOGPTAU
#undef TTAU

/*============= end of nir_opacity() ==============================*/

/*============ vis_opacity() ======================================*/

/*
 * NOTE: Greathouse now has temperature-dependent CH_4 opacities
 *       from Karkoschka (2010), which we could use to upgrade here.
 */

#define TCIA(itemp) tcia_table[itemp-1].x

void vis_opacity(planetspec   *planet,
                 int           J,
                 int           I,
                 rt_band      *vis,
                 disort_state *ds)
{
  int
    K,lev,
    lc,iwn,itemp;
  double
    tot_scat,dtau,
    m,b,dt,temp,
    cd_xh2,cd_xhe,cd_xch4,
    n_amagat,h2_h2_factor,
    dz;
  static int
    initialized = FALSE;
  static float_triplet
   *tcia_table;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="vis_opacity";

  if (!initialized) {
    /* Allocate memory */
    tcia_table = ftriplet(0,vis->cia->nt-1,dbmsname);

    /* Assign values */
    for (itemp = 1; itemp <= vis->cia->nt; itemp++) {
      TCIA(itemp) = VIS_CIA_T(itemp);
    }

    initialized = TRUE;
  }

  ds->flag.planck = vis->planck;

  if (ds->flag.planck) {
    ds->bc.temis = 1.;   /* emissivity at the top boundary */
    ds->bc.ttemp = T3(0,J,I);
    for (K = KLO; K <= KHI; K++) {
      lev         = K-1;
      TEMPER(lev) = T2(K,J,I);
    }
    ds->bc.btemp = T3(grid.nk,J,I);
  }

  /* 
   * Zero vis->dtauc
   */
  memset(vis->dtauc,0,vis->nwn*ds->nlyr*sizeof(double));

  for (lc = 1; lc <= ds->nlyr; lc++) {
    K = lc;

    /*
     * Calculate number density [amagat]
     */
    n_amagat = (P3(K,J,I)*273.15)/(T3(K,J,I)*1.01325e+5);

    /*
     * DISORT layer thickness [m]
     */
    dz = Z2(K,J,I)-Z2(K+1,J,I);

    if (vis->cia->on) {
      /*
       * H_2-H_2 collisionally induced absorption.
       *
       * NOTE: For Uranus and Neptune, may need to extend data table to colder temperatures.
       * NOTE: For deep (warm) layers, may need to extend data table to warmer temperatures.
       */

      /*
       * Density times column density [km almagat^2]
       */
      h2_h2_factor = SQR(planet->x_h2*n_amagat)*(dz/1000.);

      temp  = LIMIT_RANGE(TCIA(1),T3(K,J,I),TCIA(vis->cia->nt));
      itemp = hunt_place_in_table(vis->cia->nt,tcia_table,temp,&dt,itemp);
      itemp++;  /* Shift to unit-based */

      for (iwn = 1; iwn <= vis->nwn; iwn++) {
        m    = (VIS_CIA_H2_H2(itemp+1,iwn)-VIS_CIA_H2_H2(itemp,iwn))/dt;
        b    = VIS_CIA_H2_H2(itemp,iwn)-m*VIS_CIA_T(itemp);
        dtau = (m*temp+b)*h2_h2_factor;

        VIS_DTAUC(lc,iwn) += dtau;
      }
    }

    if (vis->k_nu->on) {
      /*
       * Methane absorption.
       */
      if (var.species[CH_4_INDEX].on == TRUE) {
        cd_xch4 = (dz/1000.)*n_amagat*X(CH_4_INDEX,VAPOR,K,J,I);

        for (iwn = 1; iwn <= vis->nwn; iwn++) {
          VIS_DTAUC(lc,iwn) += cd_xch4*VIS_K_NU_CH4(iwn);
        }
      }
    }

    if (vis->sigma->on) {
      /*
       * Rayleigh scattering.
       *
       * 1.e-4 m2 = 1 cm2
       */
      cd_xh2  = dz*n_amagat*(N_STP*1.e-4)*planet->x_h2;
      cd_xhe  = dz*n_amagat*(N_STP*1.e-4)*planet->x_he;
      if (var.species[CH_4_INDEX].on == TRUE) {
        cd_xch4 = dz*n_amagat*(N_STP*1.e-4)*X(CH_4_INDEX,VAPOR,K,J,I);
      }
      else {
        cd_xch4 = 0.;
      }

      for (iwn = 1; iwn <= vis->nwn; iwn++) {
        tot_scat = VIS_SIGMA_H2( iwn)*cd_xh2
                  +VIS_SIGMA_HE( iwn)*cd_xhe
                  +VIS_SIGMA_CH4(iwn)*cd_xch4;
        VIS_DTAUC(lc,iwn) += tot_scat;
        VIS_SSALB(lc,iwn)  = tot_scat/VIS_DTAUC(lc,iwn);
      }
    }
  }

  return;
}

#undef TCIA

/*============ end of vis_opacity() ===============================*/

/*============ uv_opacity() =======================================*/

void uv_opacity(planetspec   *planet,
                int           J,
                int           I,
                rt_band      *uv,
                disort_state *ds)
{
  int
    K,lev,
    lc,iwn;
  double
    dz,n_amagat,
    cd,cd_xh2,cd_xhe,
    cd_xch4,cd_xc2h2,cd_xc2h4,cd_xc2h6,
    tot_scat;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="uv_opacity";

  ds->flag.planck = uv->planck;

  if (ds->flag.planck) {
    ds->bc.temis = 1.;   /* emissivity at the top boundary */
    ds->bc.ttemp = T3(0,J,I);
    for (K = KLO; K <= KHI; K++) {
      lev         = K-1;
      TEMPER(lev) = T2(K,J,I);
    }
    ds->bc.btemp = T3(grid.nk,J,I);
  }

  /* 
   * Zero uv->dtauc
   */
  memset(uv->dtauc,0,uv->nwn*ds->nlyr*sizeof(double));

  for (lc = 1; lc <= ds->nlyr; lc++) {
    K = lc;

    /*
     * Calculate number density [amagat]
     */
    n_amagat = (P3(K,J,I)*273.15)/(T3(K,J,I)*1.01325e+5);

    /*
     * DISORT layer thickness [m]
     */
    dz = Z2(K,J,I)-Z2(K+1,J,I);

    /* 
     * cd is column density [molecules/cm^2]
     * (1.e-4) = [m/cm]^2
     */
    cd = (1.e-4)*dz*n_amagat*N_STP;

    /*
     * Calculate partial column densities for each species [molecules/cm^2]
     */
    cd_xh2  = cd*planet->x_h2;
    cd_xhe  = cd*planet->x_he;

    if (var.species[CH_4_INDEX].on == TRUE) {
      cd_xch4 = cd*X(CH_4_INDEX,VAPOR,K,J,I);
    }
    else {
      cd_xch4 = 0.;
    }

    if (var.species[C_2H_2_INDEX].on == TRUE) {
      cd_xc2h2 = cd*X(C_2H_2_INDEX,VAPOR,K,J,I);
    }
    else {
      cd_xc2h2 = 0.;
    }

    if (var.species[C_2H_4_INDEX].on == TRUE) {
      cd_xc2h4 = cd*X(C_2H_4_INDEX,VAPOR,K,J,I);
    }
    else {
      cd_xc2h4 = 0.;
    }

    if (var.species[C_2H_6_INDEX].on == TRUE) {
      cd_xc2h6 = cd*X(C_2H_6_INDEX,VAPOR,K,J,I);
    }
    else {
      cd_xc2h6 = 0.;
    }

    if (uv->k->on) {
      /* 
       * UV_K [cm^2/molecule], cd [molecules/cm^2]
       */
      for (iwn = 1; iwn <= uv->nwn; iwn++) {
        UV_DTAUC(lc,iwn) += UV_K_CH4( iwn)*cd_xch4
                           +UV_K_C2H2(iwn)*cd_xc2h2
                           +UV_K_C2H4(iwn)*cd_xc2h4
                           +UV_K_C2H6(iwn)*cd_xc2h6;
      }
    }

    if (uv->sigma->on) {
      for (iwn = 1; iwn <= uv->nwn; iwn++) {
        tot_scat = UV_SIGMA_H2( iwn)*cd_xh2
                  +UV_SIGMA_HE( iwn)*cd_xhe
                  +UV_SIGMA_CH4(iwn)*cd_xch4;

        UV_DTAUC(lc,iwn) += tot_scat;
        UV_SSALB(lc,iwn)  = tot_scat/UV_DTAUC(lc,iwn);
      }
    }
  }

  return;
}

/*============ end of uv_opacity() ================================*/

/*============ cross_section_rayleigh() ===========================*/

/*
 * cross section: cm^2
 *    wavenumber: cm-1
 *
 * Fortran subroutine: crosssectionrayleigh()
 */

void cross_section_rayleigh(rt_band *band)
{
  register int
    iwn;
  register double
    w2;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="cross_section_rayleigh";


  for (iwn = 0; iwn < band->nwn; iwn++) {
    /* 1.e-8 = [cm/angstrom] */
    w2 = SQR(band->wn[iwn]*1.e-8);

    /* J. Moses fit to Ford & Browne (1973) */
    band->sigma->h2[iwn] = ((((((1.38e+19)*w2
                                +1.20e+12)*w2
                                +1.00e+06)*w2
                                +1.53e+00)*w2
                                +1.25e-06)*w2
                                +8.70e-13)*w2*w2;

    /* Edgington et al. (1998) */
    band->sigma->he[iwn] = (((3.06e-03)*w2
                             +2.66e-08)*w2
                             +5.78e-14)*w2*w2;

    /* Rages <NEED REFERENCE> */
    band->sigma->ch4[iwn] = (1.230866e-17/2.687e+19)*SQR(((3.408e+02)*w2
                                                          +4.318e-04)*w2);
  }

  return;
}


/*============ end of cross_section_rayleigh() ====================*/

/*============ solar_irradiance() =================================*/

/*
 * Input: lambda, wavelength [micron]
 * Output: solar irradiance [W/m2/micron]
 *
 * Using the Gueymard (2004, Solar Energy 76, 423-453) solar spectra.
 */
EPIC_FLOAT solar_irradiance(EPIC_FLOAT lambda)
{
  char
    header[N_STR];
  static int
    i,
    ndat,
    initialized = FALSE;
  EPIC_FLOAT
    lambda_d;
  static float_triplet
    *table;
  FILE
    *infile;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="cross_section_rayleigh";

  if (!initialized) {

    sprintf(header,EPIC_PATH"/data/rt/solar/Gueymard2004.dat");
    infile = fopen(header,"r");
    if (!infile) {
      sprintf(Message,"unable to open %s",header);
      epic_error(dbmsname,Message);
    }

    /*
     * Read in size of data table.
     */
    for (i = 1; i <= 9; i++) {
      fgets(header,N_STR,infile);
    }
    fscanf(infile,"%*s %*s %*s %d",&ndat);
    fgets(header,N_STR,infile);
    fgets(header,N_STR,infile);

    /* Allocate memory */
    table = ftriplet(0,ndat-1,dbmsname);

    /* Read in data */
    for (i = 0; i < ndat; i++) {
      fscanf(infile,"%lf %lf",&table[i].x,&table[i].y);
      /* Convert [nm] to [micron] and [W/m2/nm] to [W/m2/micron] */
      table[i].x /= 1000.;
      table[i].y *= 1000.;
    }
    fclose(infile);

    /* Set up spline */
    spline_pchip(ndat,table);

    initialized = TRUE;
  }

  /*
   * Return zero if outside table range.
   */
  if (lambda < table[     0].x ||
      lambda > table[ndat-1].x   ) {
    return 0.;
  }

  /* Hunt place in table */
  i = hunt_place_in_table(ndat,table,lambda,&lambda_d,i);

  /* Piecewise monotonic cubic spline */
  return splint_pchip(lambda,table+i,lambda_d);
}

/*============ end of solar_irradiance() ==========================*/

/*============ beer_law_only() ====================================*/

/*
 * The radiative-transfer special case of an incoming direct beam, 
 * ds->bc.fbeam, that is attenuated by opacity via the Beer-Lambert-Bougher
 * law, but otherwise has no scattering and no Planck thermal emission.
 *
 * The option of the pseudo-spherical correction is included.  Uses
 * the same inputs as for c_twostr().
 *
 * This function yields identical results to c_twostr() for this special
 * case, but runs about twice as fast; the speed-up is the result of eliminating
 * unnecessary lines of code.
 *
 * NOTE: Use c_twostr() and not this function if there is either thermal
 *       emission by the atmosphere itself (ds->flag.planck == TRUE), or
 *       scattering (SSALB(lu) != 0.0 for any user-defined position, lu). 
 *
 * The code is an adaptation of the relevant lines from c_twostr().
 * The function c_twostr() was written by A. Kylling and adapted into C
 * by T. Dowling.  The complete cdisort source code is in
 * $EPIC_PATH/src/rt/cdisort.
 */

#define FDNTOT(lu) fdntot[lu-1]
#define FNET(lu)   fnet[lu-1]

/*
 * ABSCUT is the cumulative optical-depth cut off
 */
#define ABSCUT 10.

void beer_law_only(disort_state  *ds,
                   disort_output *out,
                   double         radius)
{
  int
    lc;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="beer_law_only";

  /*
   * Check validity of inputs, and exit with a message if there is an error.
   */
  if (ds->flag.planck) {
    sprintf(Message,"called with ds->flag.planck == TRUE; call c_twostr() instead");
    epic_error(dbmsname,Message);
  }
  else if (ds->bc.fisot != 0.) {
    sprintf(Message,"called with ds->bc.fisot=%g != 0; call c_twostr() instead",ds->bc.fisot);
    epic_error(dbmsname,Message);
  }
  else if (ds->bc.fbeam < 0.) {
    sprintf(Message,"called with ds->bc.fbeam=%g < 0",ds->bc.fbeam);
    epic_error(dbmsname,Message);
  }

  if (!ds->flag.prnt[0] && !ds->flag.prnt[1]) {
    /*
     * Not printing inputs or outputs.
     * Fewer arrays are necessary in this case.
     */
    if ((ds->bc.fbeam == 0.                  ) || 
        (ds->bc.umu0 <= 0. && !ds->flag.spher)   ) {
      /*
       * Cases with no incoming beam to process.
       * Set UTAU if necessary, and zero output arrays.
       */
      if (!ds->flag.usrtau) {
        /*
         * Set output levels at computational layer boundaries
         */
        register double
          the_tauc;

        ds->ntau = ds->nlyr+1;
        the_tauc = 0.;
        for (lc = 1; lc <= ds->nlyr; lc++) {
          UTAU(lc)  = the_tauc;
          the_tauc += DTAUC(lc);
        }
        UTAU(lc) = the_tauc;
      }

      /*
       * Zero the disort_radiant-structure output array.
       */
      memset(out->rad,0,ds->ntau*sizeof(disort_radiant));
    }
    else {
      /*
       * The Sun is in the sky; process ds->bc.fbeam 
       */
      int
        lu,lyrcut,ncut;
      int
       *layru;
      double
        dirint,
        zenang,taup,chtau_tmp;
      double
       *tauc,
       *ch;

      /*
       * Allocate zeroed memory.
       *
       * NOTE: Avoid allocating unnecessary arrays like FDNTOT and FNET.
       */
      tauc  = dvector(0,ds->nlyr,dbmsname);
      ch    = dvector(0,ds->nlyr-1,dbmsname);
      layru = ivector(0,ds->ntau-1,dbmsname);

      /*
       * Calculate cumulative optical depth
       */
      for (lc = 1; lc <= ds->nlyr; lc++) {
        TAUC(lc) = TAUC(lc-1)+DTAUC(lc);
      }

      if (!ds->flag.usrtau) {
        /*
         * Set output levels at computational layer boundaries
         */
        ds->ntau = ds->nlyr+1;
        for (lc = 0; lc <= ds->ntau-1; lc++) {
          UTAU(lc+1) = TAUC(lc);
        }
      }

      /*
       * Zero the disort_radiant-structure output array.
       */
      memset(out->rad,0,ds->ntau*sizeof(disort_radiant));

      /*
       * Set LAYRU(lu), the layer containing UTAU(lu).
       */
      for (lu = 1; lu <= ds->ntau; lu++) {
        for (lc = 1; lc <= ds->nlyr-1; lc++) {
          if (UTAU(lu) >= TAUC(lc-1) && UTAU(lu) <= TAUC(lc)) {
            break;
          }
        }
        LAYRU(lu) = lc;
      }

      ncut = 0;
      for (lc = 1; lc <= ds->nlyr; lc++) {
        if (TAUC(lc) < ABSCUT) {
          ncut = lc;
        }
        else {
          break;
        }
      }

      /*
       * Cut off medium below absorption optical depth ABSCUT. 
       * Not worth the trouble for one-layer problems, though.
       */
      lyrcut = FALSE;
      if (TAUC(ds->nlyr) >= ABSCUT && ds->nlyr > 1) {
        lyrcut = TRUE;
      }
      if (!lyrcut) {
        ncut = ds->nlyr;
      }

      zenang = acos(ds->bc.umu0)/DEG;
    
      if (ds->flag.spher == TRUE) {
        /*
         * Calculate Chapman function to handle spherical geometry.
         */
        for (lc = 1; lc <= ncut; lc++) {
          taup = TAUC(lc-1)+DTAUC(lc)/2.;
          if (taup == 0.) {
            CH(lc) = ds->bc.umu0;
          }
          else {
            chtau_tmp = c_chapman(lc,0.5,tauc,ds->nlyr,ds->zd,ds->dtauc,zenang,radius);
            CH(lc)    = taup/chtau_tmp;
          }
        }
      }
      else {
        /*
         * Plane-parallel case.
         */
        for (lc = 1; lc <= ncut; lc++) {
          CH(lc) = ds->bc.umu0;
        }
      }

      /*
       * Loop over user-output levels.
       */
      for (lu = 1; lu <= ds->ntau; lu++) {
        if (lyrcut && LAYRU(lu) > ncut) {
          /*
           * No radiation reaches this region, or it is being neglected.
           */
          break;
        }
        else {
          /*
           * Apply Beer-Lambert-Bougher Law.
           */
          dirint     = ds->bc.fbeam*exp(-UTAU(lu)/CH(LAYRU(lu)));
          RFLDIR(lu) = dirint*fabs(ds->bc.umu0);

          /*
           * RFLDN(lu), FLUP(lu) and UAVG(lu) are zeroed out as part of out->rad above.
           */

          /*
           * Direct-beam intensity is the only source of total mean intensity.
           */
          UAVG(lu) = dirint/(4.*M_PI);

          /*
           * Flux divergence, not calculated as a finite difference but analytically as in KST(12),
           * here with single-scattering albedo and planck_source set to zero.
           */
          DFDT(lu) = dirint;
        }
      }

      /*
       * Free allocated memory
       */
      free(tauc);
      free(ch);
      free(layru);
    }
  }
  else {
    /*
     * Printing of information is requested, which requires additional arrays.
     */
    int
      lu,lyrcut,ncut;
    int
     *layru;
    double
      dirint,
      zenang,taup,chtau_tmp;
    double
     *tauc,
     *ch,
     *fdntot,
     *fnet;

    /*
     * Allocate zeroed memory.
     */
    tauc   = dvector(0,ds->nlyr,dbmsname);
    ch     = dvector(0,ds->nlyr-1,dbmsname);
    fdntot = dvector(0,ds->ntau-1,dbmsname);
    fnet   = dvector(0,ds->ntau-1,dbmsname);
    layru  = ivector(0,ds->ntau-1,dbmsname);

    /*
     * Calculate cumulative optical depth
     */
    for (lc = 1; lc <= ds->nlyr; lc++) {
      TAUC(lc) = TAUC(lc-1)+DTAUC(lc);
    }

    if (!ds->flag.usrtau) {
      /*
       * Set output levels at computational layer boundaries
       */
      ds->ntau = ds->nlyr+1;
      for (lc = 0; lc <= ds->ntau-1; lc++) {
        UTAU(lc+1) = TAUC(lc);
      }
    }

    /*
     * Zero the disort_radiant-structure output array.
     */
    memset(out->rad,0,ds->ntau*sizeof(disort_radiant));

    /*
     * Set LAYRU(lu), the layer containing UTAU(lu).
     */
    for (lu = 1; lu <= ds->ntau; lu++) {
      for (lc = 1; lc <= ds->nlyr-1; lc++) {
        if (UTAU(lu) >= TAUC(lc-1) && UTAU(lu) <= TAUC(lc)) {
          break;
        }
      }
      LAYRU(lu) = lc;
    }

    ncut = 0;
    for (lc = 1; lc <= ds->nlyr; lc++) {
      if (TAUC(lc) < ABSCUT) {
        ncut = lc;
      }
      else {
        break;
      }
    }

    /*
     * Cut off medium below absorption optical depth ABSCUT. 
     * Not worth the trouble for one-layer problems, though.
     */
    lyrcut = FALSE;
    if (TAUC(ds->nlyr) >= ABSCUT && ds->nlyr > 1) {
      lyrcut = TRUE;
    }
    if (!lyrcut) {
      ncut = ds->nlyr;
    }

    if (ds->flag.prnt[0]) {
      /*
       * Print input information
       */
      fprintf(stdout,"\n\n"
                     " ****************************************************************************************************\n"
                     " %s\n"
                     " ****************************************************************************************************\n",
                     ds->header);
      fprintf(stdout,"\n No. streams = %4d     No. computational layers =%4d\n",ds->nstr,ds->nlyr);
      fprintf(stdout,"%4d User optical depths :",ds->ntau);
      for (lu = 1; lu <= ds->ntau; lu++) {
        fprintf(stdout,"%10.4f",UTAU(lu));
        if (lu%10 == 0) {
          fprintf(stdout,"\n                          ");
        }
      }
      fprintf(stdout,"\n");

      if (ds->flag.spher) {
        fprintf(stdout," Pseudo-spherical geometry invoked\n");
      }
      fprintf(stdout," No thermal emission\n");
      fprintf(stdout,"    Incident beam with intensity =%11.3e and polar angle cosine = %8.5f\n"
                     "    plus isotropic incident intensity =%11.3e\n",
                     ds->bc.fbeam,ds->bc.umu0,ds->bc.fisot);
      fprintf(stdout,"    Bottom albedo (lambertian) =%8.4f\n",ds->bc.albedo);
      fprintf(stdout," Non-scattering case\n");
      if(lyrcut) {
        fprintf(stdout," Sets radiation = 0 below absorption optical depth %g\n",ABSCUT);
      }
      fprintf(stdout,"\n"
                     "\n                   total    single"
                     "\n       optical   optical   scatter"
                     "\n         depth     depth    albedo\n");
      for (lc = 1; lc <= ds->nlyr; lc++) {
        fprintf(stdout,"%4d%10.4f%10.4f%10.5f\n",
                       lc,DTAUC(lc),TAUC(lc),SSALB(lc));
      }
    }

    if (ds->bc.fbeam > 0.) {
      zenang = acos(ds->bc.umu0)/DEG;
    
      if (ds->flag.spher == TRUE) {
        /*
         * Calculate Chapman function to handle spherical geometry.
         */
        for (lc = 1; lc <= ncut; lc++) {
          taup = TAUC(lc-1)+DTAUC(lc)/2.;
          if (taup == 0.) {
            CH(lc) = ds->bc.umu0;
          }
          else {
            chtau_tmp = c_chapman(lc,0.5,tauc,ds->nlyr,ds->zd,ds->dtauc,zenang,radius);
            CH(lc)    = taup/chtau_tmp;
          }
        }
      }
      else {
        /*
         * Plane-parallel case.
         */
        for (lc = 1; lc <= ncut; lc++) {
          CH(lc) = ds->bc.umu0;
        }
      }
    }

    /*
     * Loop over user-output levels
     */
    for (lu = 1; lu <= ds->ntau; lu++) {
      if (lyrcut && LAYRU(lu) > ncut) {
        /*
         * No radiation reaches this region, or it is being neglected.
         */
        break;
      }
      else {
        /*
         * Radiatively active region.
         */
        if (ds->bc.fbeam > 0.) {
          if (ds->bc.umu0 > 0. || ds->flag.spher) {
            /* 
             * Incoming, external-source (direct) beam is above the horizon (typically the
             * shortwave, solar insolation), or is in play via the pseudo-spherical correction.
             *
             * Apply Beer-Lambert-Bougher Law.
             */
            dirint     = ds->bc.fbeam*exp(-UTAU(lu)/CH(LAYRU(lu)));
            RFLDIR(lu) = dirint*fabs(ds->bc.umu0);
          }
          else {
            /*
             * Incoming beam is below the effective horizon.
             */
            dirint = 0.;
          }
        }
        else {
          /* 
           * No external-source beam.
           */
          dirint = 0.;
        }

        /*
         * RFLDN(lu), FLUP(lu) and UAVG(lu) are zeroed out as part of out->rad above.
         */
        FDNTOT(lu) = RFLDN(lu)+RFLDIR(lu);
        FNET(lu)   = FDNTOT(lu)-FLUP(lu);

        /*
         * Direct-beam intensity is the only source of total mean intensity.
         */
        UAVG(lu) = dirint/(4.*M_PI);

        /*
         * Flux divergence, not calculated as a finite difference but analytically as in KST(12),
         * here with single-scattering albedo and planck_source set to zero.
         */
        DFDT(lu) = dirint;
      }
    }

    if (ds->flag.prnt[1]) {
      fprintf(stdout,"\n\n                     <----------------------- Fluxes ----------------------->\n"
                     "   optical  compu    downward    downward    downward      upward                    mean      Planck   d(net flux)\n"
                     "     depth  layer      direct     diffuse       total     diffuse         net   intensity      source   / d(op dep)\n");
      for (lu = 1; lu <= ds->ntau; lu++) {
        fprintf(stdout,"%10.4f%7d%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%12.3e%14.3e\n",
                       UTAU(lu),LAYRU(lu),RFLDIR(lu),RFLDN(lu),FDNTOT(lu),FLUP(lu),FNET(lu),UAVG(lu),0.,DFDT(lu));
      }
    }

    /*
     * Free allocated memory.
     */
    free(tauc);
    free(ch);
    free(fdntot);
    free(fnet);
    free(layru);
  }

  return;
}

#undef FDNTOT
#undef FNET
#undef ABSCUT
 
/*============ end of beers_law_only() ============================*/

/*---------------------*
 * Planetary functions *
 *---------------------*/

/*============= ring_shadow() =====================================*/

/*
 * Returns attenuation factor of insolation due to planet's rings for the given latitude,
 * averaged over the day side for the given planet latitude.
 *
 * Inputs:
 *            lat: planetocentric [deg]
 *   subsolar_lat: planetocentric [deg]
 *
 * NOTE: Currently implemented for Saturn only.
 *
 * Bruno Bezard derived the diurnally averaged equations:
 *     Bezard B, 1986, Variations saisonniere de la structure thermique
 *       et composition chimique de Jupiter, Saturn, et Uranus, Ph.D. thesis, Univ. Paris 7, Paris, France
 * A scan of the relevant pages may be found in the file
 *     $EPIC_PATH/data/saturn/archive/bezard_thesis_shadow.pdf.
 *
 * J. Moses and T. Greathouse implemented Bezard's algorithm as the Fortran shadow(), first used in  
 *   Moses JI, Greathouse TK, 2005, Latitudinal and seasonal models of stratospheric photochemistry on
 *       Saturn: Comparison with infrared data from IRTF/TEXES, J. Geophys. Res., 110, E09007,
 *       doi:10.1029/2005JE002450
 *
 * T. Dowling converted the Fortran code to C, and removed redundancies in conditionals, June 2010
 *
 * Fortran subroutine:
 *   shadow()
 */

double ring_shadow(planetspec *planet,
                   double      lat,
                   double      subsolar_lat)
{
  const double         /* Radii [km] for Saturn's main rings and the Cassini Division */
    rcint  =  74655.,
    rcext  =  91975.,
    rb1int =   rcext,
    rb1ext =  98500.,
    rb2int =  rb1ext,
    rb2ext = 117510.,
    rcdint =  rb2ext,
    rcdext = 122340.,
    raint  =  rcdext,
    raext  = 136780.;
  const double        /* Optical depths for Saturn's main rings and the Cassini Division */
    tauc  = 0.099,
    taub1 = 0.920,
    taub2 = 2.070,
    taucd = 0.140,
    taua  = 0.510;
  const double
    dec    = subsolar_lat*DEG,
    cdec   = cos(dec),
    sdec   = sin(dec),
    re_rp2 = SQR(planet->re/planet->rp);
  double
    latr,clat,slat,c2lat,s2lat,
    tlat,tdec,t2dec,alpha,phie,
    ec,eb1,eb2,ecd,ea,
    ringfrac,
    f0,f1,f2,tmp,
    aintc,aextc,aintb1,aintb2,aextb1,aextb2,aintcd,aextcd,ainta,aexta,
    phiinta,phiintb1,phiintb2,phiintcd,phiintc,phiexta,phiextb1,phiextb2,phiextcd,phiextc;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int
    idbms=0;
  static char
    dbmsname[]="ring_shadow";

  /*
   * Check validity of inputs
   */
  if (strcmp(planet->name,"saturn") != 0) {
    fprintf(stderr,"ERROR: rt_.c, ring_shadow: planet->name=%s not yet implemented\n",planet->name);
    exit(1);
  }

  latr  = lat*DEG;
  clat  = cos(latr);
  c2lat = clat*clat;
  slat  = sin(latr);
  s2lat = slat*slat;

  if (latr >= M_PI_2) {
    tlat  =  1.e+19;
  }
  else if (latr <= -M_PI_2) {
    tlat  = -1.e+19;
  }
  else {
    tlat  = tan(latr);
  }

  if (dec >= M_PI_2) {
    tdec  =  1.e+19;
    t2dec =  1.e+38;
  }
  else if (dec <= -M_PI_2) {
    tdec  = -1.e+19;
    t2dec =  1.e+38;
  }
  else {
    tdec  = tan(dec);
    t2dec = tdec*tdec;
  }

  alpha = -(re_rp2*tdec)*tlat;
  phie  = acos(LIMIT_RANGE(-1.,alpha,1.));

  tmp = slat*clat;
  if (tmp != 0.) {
    f0 = .5*tdec/tmp;
  }
  else {
    f0 = (tdec >= 0.) ? 1.e+38 : -1.e+38;
  }

  if (t2dec != 0.) {
    f1 = f0*(c2lat+s2lat/t2dec);
  }
  else {
    f1 = (f0*s2lat >= 0.) ? 1.e+38 : -1.e+38;
  }

  f2 = f0*(c2lat+s2lat*re_rp2);

  aintc  = f1-f2*SQR(rcint /planet->re);
  aintb1 = f1-f2*SQR(rb1int/planet->re);
  aintb2 = f1-f2*SQR(rb2int/planet->re);
  aintcd = f1-f2*SQR(rcdint/planet->re);
  ainta  = f1-f2*SQR(raint /planet->re);

  aextc  = f1-f2*SQR(rcext /planet->re);
  aextb1 = f1-f2*SQR(rb1ext/planet->re);
  aextb2 = f1-f2*SQR(rb2ext/planet->re);
  aextcd = f1-f2*SQR(rcdext/planet->re);
  aexta  = f1-f2*SQR(raext /planet->re);

  phiintc  = acos(LIMIT_RANGE(-1.,aintc, 1.));
  phiintb1 = acos(LIMIT_RANGE(-1.,aintb1,1.));
  phiintb2 = acos(LIMIT_RANGE(-1.,aintb2,1.));
  phiintcd = acos(LIMIT_RANGE(-1.,aintcd,1.));
  phiinta  = acos(LIMIT_RANGE(-1.,ainta, 1.));

  phiextc  = acos(LIMIT_RANGE(-1.,aextc, 1.));
  phiextb1 = acos(LIMIT_RANGE(-1.,aextb1,1.));
  phiextb2 = acos(LIMIT_RANGE(-1.,aextb2,1.));
  phiextcd = acos(LIMIT_RANGE(-1.,aextcd,1.));
  phiexta  = acos(LIMIT_RANGE(-1.,aexta, 1.));

  if (phie >= phiexta && phie <= phiintc) {
    if (phie <= phiinta) {
      if (cdec != 0.) {
        ea = exp(-taua/cdec);
      }
      else {
        ea = 0.;
      }
      /*--A ring--*/
      ringfrac  = clat*cdec*(sin(phiexta)+ea*(sin(phie)-sin(phiexta))) 
          +re_rp2*slat*sdec*(    phiexta +ea*(    phie -    phiexta ));
    }
    else {
      if (phie <= phiintcd) {
        if (cdec != 0.) {
          ecd = exp(-taucd/cdec);
        }
        else {
          ecd = 0.;
        }
        if (ainta > 1.) {
         /*--Cassini Division--*/
          ringfrac = clat*cdec*(sin(phiextcd)+ecd*(sin(phie)- sin(phiextcd))) 
             +re_rp2*slat*sdec*(    phiextcd +ecd*(    phie -     phiextcd ));
        }
        else {
          if (cdec != 0.) {
            ea  = exp(-taua /cdec);
          }
          else {
            ea  = 0.;
          }
          /*--A ring and Cassini Division--*/
          ringfrac = clat*cdec*(sin(phiexta)+ ea*(sin(phiinta)-sin(phiexta ))
                                            +ecd*(sin(phie   )-sin(phiextcd))) 
             +re_rp2*slat*sdec*(    phiexta + ea*(    phiinta -    phiexta  )  
                                            +ecd*(    phie    -    phiextcd ));
        }
      }
      else {
        if (phie <= phiintb2) {
          if (cdec != 0.) {
            eb2 = exp(-taub2/cdec);
          }
          else {
            eb2 = 0.;
          }
          if (ainta > 1.) {
            if (aintcd > 1.) {
              /*--B2 ring--*/
              ringfrac = clat*cdec*(sin(phiextb2)+eb2*(sin(phie)-sin(phiextb2))) 
                 +re_rp2*slat*sdec*(    phiextb2 +eb2*(    phie -    phiextb2 ));
            }
            else {
              if (cdec != 0.) {
                ecd = exp(-taucd/cdec);
              }
              else {
                ecd = 0.;
              }
              /*--Cassini Division and B2 ring--*/
              ringfrac = clat*cdec*(sin(phiextcd)+ecd*(sin(phiintcd)-sin(phiextcd)) 
                                                 +eb2*(sin(phie    )-sin(phiextb2))) 
                 +re_rp2*slat*sdec*(    phiextcd +ecd*(    phiintcd -    phiextcd )
                                                 +eb2*(    phie     -    phiextb2 ));
            }
          }
          else {
            /*--A ring, Cassini Division, and B2 ring--*/
            if (cdec != 0.) {
              ea  = exp(-taua /cdec);
              ecd = exp(-taucd/cdec);
            }
            else {
              ea  = 0.;
              ecd = 0.;
            }
            ringfrac = clat*cdec*(sin(phiexta)+ ea*(sin(phiinta )-sin(phiexta ))
                                              +ecd*(sin(phiintcd)-sin(phiextcd)) 
                                              +eb2*(sin(phie    )-sin(phiextb2))) 
               +re_rp2*slat*sdec*(    phiexta + ea*(    phiinta  -    phiexta  )
                                              +ecd*(    phiintcd -    phiextcd )
                                              +eb2*(    phie     -    phiextb2 ));
          }
        }
        else {
          if (phie <= phiintb1) {
            if (cdec != 0.) {
              eb1 = exp(-taub1/cdec);
            }
            else {
              eb1 = 0.;
            }
            if (ainta > 1.) {
              if (aintcd > 1.) {
                if (aintb2 > 1.) {
                  /*--B1 ring--*/
                  ringfrac = clat*cdec*(sin(phiextb1)+eb1*(sin(phie)-sin(phiextb1))) 
                     +re_rp2*slat*sdec*(    phiextb1 +eb1*(    phie -    phiextb1 ));
                }
                else {
                  /*--B2 and B1 rings--*/
                  if (cdec != 0.) {
                    eb2 = exp(-taub2/cdec);
                  }
                  else {
                    eb2 = 0.;
                  }
                  ringfrac = clat*cdec*(sin(phiextb2)+eb2*(sin(phiintb2)-sin(phiextb2))  
                                                     +eb1*(sin(phie    )-sin(phiextb1))) 
                     +re_rp2*slat*sdec*(    phiextb2 +eb2*(    phiintb2 -    phiextb2 )
                                                     +eb1*(    phie     -    phiextb1 ));

                }
              }
              else {
                /*--Cassini Division, B2 and B1 rings--*/
                if (cdec != 0.) {
                  ecd = exp(-taucd/cdec);
                  eb2 = exp(-taub2/cdec);
                }
                else {
                  ecd = 0.;
                  eb2 = 0.;
                }
                ringfrac = clat*cdec*(sin(phiextcd)+ecd*(sin(phiintcd)-sin(phiextcd)) 
                                                   +eb2*(sin(phiintb2)-sin(phiextb2))  
                                                   +eb1*(sin(phie    )-sin(phiextb1))) 
                   +re_rp2*slat*sdec*(    phiextcd +ecd*(    phiintcd -    phiextcd )
                                                   +eb2*(    phiintb2 -    phiextb2 )
                                                   +eb1*(    phie     -    phiextb1 ));
              }
            }
            else {
              /*--A ring, Cassini Division, B2 and B1 rings--*/
              if (cdec != 0.) {
                ea  = exp(-taua /cdec);
                ecd = exp(-taucd/cdec);
                eb2 = exp(-taub2/cdec);
              }
              else {
                ea  = 0.;
                ecd = 0.;
                eb2 = 0.;
              }
              ringfrac = clat*cdec*(sin(phiexta)+ ea*(sin(phiinta )-sin(phiexta ))
                                                +ecd*(sin(phiintcd)-sin(phiextcd)) 
                                                +eb2*(sin(phiintb2)-sin(phiextb2))  
                                                +eb1*(sin(phie    )-sin(phiextb1))) 
                 +re_rp2*slat*sdec*(    phiexta + ea*(    phiinta  -    phiexta  )
                                                +ecd*(    phiintcd -    phiextcd )
                                                +eb2*(    phiintb2 -    phiextb2 )
                                                +eb1*(    phie     -    phiextb1 ));
            }
          }
          else {
            if (cdec != 0.) {
              ec = exp(-tauc /cdec);
            }
            else {
              ec = 0.;
            }
            if (ainta > 1.) {
              if (aintcd > 1.) {
                if (aintb2 > 1.) {
                  if (aintb1 > 1.) {
                    /*--C ring--*/
                    ringfrac = clat*cdec*(sin(phiextc)+ec*(sin(phie)-sin(phiextc))) 
                       +re_rp2*slat*sdec*(    phiextc +ec*(    phie -    phiextc ));
                  }
                  else {
                    /*--B1 and C rings--*/
                    if (cdec != 0.) {
                      eb1 = exp(-taub1/cdec);
                    }
                    else {
                      eb1 = 0.;
                    }
                    ringfrac = clat*cdec*(sin(phiextb1)+eb1*(sin(phiintb1)-sin(phiextb1)) 
                                                       + ec*(sin(phie    )-sin(phiextc ))) 
                       +re_rp2*slat*sdec*(    phiextb1 +eb1*(    phiintb1 -    phiextb1 )
                                                       + ec*(    phie     -    phiextc  ));
                  }
                }
                else {
                  /*--B2, B1, and C rings--*/
                  if (cdec != 0.) {
                    eb2 = exp(-taub2/cdec);
                    eb1 = exp(-taub1/cdec);
                  }
                  else {
                    eb2 = 0.;
                    eb1 = 0.;
                  }
                  ringfrac = clat*cdec*(sin(phiextb2)+eb2*(sin(phiintb2)-sin(phiextb2))
                                                     +eb1*(sin(phiintb1)-sin(phiextb1)) 
                                                     + ec*(sin(phie    )-sin(phiextc ))) 
                     +re_rp2*slat*sdec*(    phiextb2 +eb2*(    phiintb2 -    phiextb2 )
                                                     +eb1*(    phiintb1 -    phiextb1 )
                                                     + ec*(    phie     -    phiextc  ));
                }
              }
              else {
                /*--Cassini Division, B2, B1, and C rings--*/
                if (cdec != 0.) {
                  ecd = exp(-taucd/cdec);
                  eb2 = exp(-taub2/cdec);
                  eb1 = exp(-taub1/cdec);
                }
                else {
                  ecd = 0.;
                  eb2 = 0.;
                  eb1 = 0.;
                }
                ringfrac = clat*cdec*(sin(phiextcd)+ecd*(sin(phiintcd)-sin(phiextcd)) 
                                                   +eb2*(sin(phiintb2)-sin(phiextb2))  
                                                   +eb1*(sin(phiintb1)-sin(phiextb1)) 
                                                   + ec*(sin(phie    )-sin(phiextc ))) 
                   +re_rp2*slat*sdec*(    phiextcd +ecd*(    phiintcd -    phiextcd )
                                                   +eb2*(    phiintb2 -    phiextb2 )
                                                   +eb1*(    phiintb1 -    phiextb1 )
                                                   + ec*(    phie     -    phiextc  ));
              }
            }
            else {
              /*--A ring, Cassini Division, B2, B1, and C rings--*/
              if (cdec != 0.) {
                ea  = exp(-taua /cdec);
                ecd = exp(-taucd/cdec);
                eb2 = exp(-taub2/cdec);
                eb1 = exp(-taub1/cdec);
              }
              else {
                ea  = 0.;
                ecd = 0.;
                eb2 = 0.;
                eb1 = 0.;
              }
              ringfrac = clat*cdec*(sin(phiexta)+ ea*(sin(phiinta )-sin(phiexta )) 
                                                +ecd*(sin(phiintcd)-sin(phiextcd)) 
                                                +eb2*(sin(phiintb2)-sin(phiextb2))  
                                                +eb1*(sin(phiintb1)-sin(phiextb1)) 
                                                + ec*(sin(phie    )-sin(phiextc ))) 
                 +re_rp2*slat*sdec*(    phiexta + ea*(    phiinta  -    phiexta  )
                                                +ecd*(    phiintcd -    phiextcd )
                                                +eb2*(    phiintb2 -    phiextb2 )
                                                +eb1*(    phiintb1 -    phiextb1 )
                                                + ec*(    phie     -    phiextc  ));
            }
          }
        }
      }
    }
    ringfrac /= clat*cdec*sin(phie)+re_rp2*slat*sdec*phie;
  }
  else {
    ringfrac = 1.;
  }

  return ringfrac;
}

/*============= end of ring_shadow() ==============================*/

/* * * * * * * end of rt_heating.c * * * * * * * * * * * * * * * * */


