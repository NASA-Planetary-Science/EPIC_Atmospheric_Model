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

/* * * * * * * * * * epic_change.c * * * * * * * * * * * * * * * * * 
 *                                                                 *
 *  Makes changes to parameters in epic.nc                         *
 *                                                                 *
 *  Options for introducing perturbations:                         *
 *                                                                 *
 *    Use -spots spots.dat to add vortices                         *
 *        -waves waves.dat to add waves                            *
 *    These can both be done at the same time.                     *
 *                                                                 *
 *  Options for converting data to isentropic coordinates:         *
 *                                                                 *
 *    Use -openmars to process an OpenMARS file                    *
 *    Use -emars to process an EMARS file                          *
 *    Use -weizmann to process a gas-giant Weizmann Institute file *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*
 * Local function prototypes:
 */
void add_spots(planetspec  *planet,
               char        *spots_file,
               EPIC_FLOAT  *pert,
               EPIC_FLOAT **Buff2D);

void add_waves(planetspec  *planet,
               char        *waves_file,
               EPIC_FLOAT  *pert,
               EPIC_FLOAT **Buff2D);

void modify_progs_from_pert(planetspec *planet,
                            EPIC_FLOAT *pert);

int number_objects_in_file(char *objects_file);

void read_spots_file(char       *spots_file,
                     EPIC_FLOAT *ampspot,
                     EPIC_FLOAT *lonspot,
                     EPIC_FLOAT *latspot,
                     EPIC_FLOAT *pspot,
                     EPIC_FLOAT *aspot,
                     EPIC_FLOAT *bspot,
                     EPIC_FLOAT *cspot_up,
                     EPIC_FLOAT *cspot_down,
                     int         adjust_amplitude);

void read_waves_file(char       *waves_file,
                     EPIC_FLOAT *latwave,
                     EPIC_FLOAT *ampwave,
                     EPIC_FLOAT *wnwave,
                     EPIC_FLOAT *pwave,
                     EPIC_FLOAT *cwave,
                     EPIC_FLOAT *fwhmwave);

void read_defaults(change_defaultspec  *def);

void write_defaults(change_defaultspec *def);

/*
 * Conversion of non-EPIC input into isentropic-coordinate output.
 *
 *   Input options:
 *     OpenMars Mars data
 *     EMARS Mars data (emars_v1.0_back_*.nc format)
 *     Weizmann-format gas-giant data
 *
 *   The customized data-processing functions are in epic_funcs_init.c,
 *   the associated function prototypes and shift macros are in epic.h,
 *   and the associated data types are in epic_datatypes.h.
 */

#undef  PERT
#define PERT(k,j,i) pert[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

/*======================= main() =====================================*/

/*
 *  NOTE: structures planet, grid, and var are declared globally in epic.h.
 */

int main(int   argc,
         char *argv[])
{
  char   
    spots_file[FILE_STR], /*  added-spot locations and sizes         */
    waves_file[FILE_STR], /*  wave perturbation parameters           */
    sflag[80],            /*  string to hold command-line flags      */
    infile[ FILE_STR],
    outfile[FILE_STR],
    openmars_infile[FILE_STR],
    openmars_outfile_qb[FILE_STR],
    openmars_outfile_uvpt[FILE_STR],
    emars_infile[FILE_STR],
    emars_outfile_qb[FILE_STR],
    emars_outfile_uvpt[FILE_STR],
    weizmann_infile[FILE_STR],
    weizmann_outfile_qb[FILE_STR],
    weizmann_outfile_uvpt[FILE_STR],
    buffer[16];
  int    
    time_index,
    openmars_itime,
    emars_itime,
    weizmann_itime,
    K,J,I,
    kk,jj,
    is,itmp,
    count,index,ii;
  int
    spots      = FALSE,
    waves      = FALSE,
    stretch_ni = FALSE,
    openmars   = FALSE,
    emars      = FALSE,
    weizmann   = FALSE;
  openmars_gridspec
    *openmars_grid;
  emars_gridspec
    *emars_grid;
  weizmann_gridspec
    *weizmann_grid;
  EPIC_FLOAT  
    dx0,dt;
  EPIC_FLOAT
    *p,*h;
  static EPIC_FLOAT
    *Buff2D[NUM_WORKING_BUFFERS];
  change_defaultspec
    defaults;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="epic_change";

#if defined(EPIC_MPI)
  sprintf(Message,"not designed to be run on multiple processors");
  epic_error(dbmsname,Message);
#endif

  declare_copyright();

  /* 
   * Interpret command-line arguments: 
   */
  /* Start with defaults: */
  sprintf(spots_file,"none");
  sprintf(waves_file,"none");
  if (argc > 1) {
    /* Read flags: */
    for (count = 1; count < argc; count++) {
      sscanf(argv[count],"%s",sflag);
      if (strcmp(sflag,"-spots") == 0) {
        sscanf(argv[++count],"%s",spots_file);
        spots = TRUE;
      }
      else if (strcmp(sflag,"-waves") == 0) {
        sscanf(argv[++count],"%s",waves_file);
        waves = TRUE;
      }
      else if (strcmp(sflag,"-stretch_ni") == 0) {
        sscanf(argv[++count],"%d",&stretch_ni);
        /*
         * Verify that stretch_ni is a power of 2.
         */
        if (frexp((double)stretch_ni,&itmp) != 0.5) {
          sprintf(Message,"-stretch_ni %d is not a power of 2",stretch_ni);
          epic_error(dbmsname,Message);
        }
      }
      else if (strcmp(sflag,"-openmars") == 0) {
        /*
         * OpenMARS (formerly MACDA)
         * Input an OpenMARS .nc file and convert to isentropic-coordinate QB* and UVPT* files.
         */
        openmars_grid = (openmars_gridspec *)calloc(1,sizeof(openmars_gridspec));

        openmars = TRUE;
      }
      else if (strcmp(sflag,"-emars") == 0) {
        /*
         * EMARS (Ensemble Mars Atmosphere Reanalysis System)
         * Input an EMARS .nc file and convert to isentropic-coordinate QB* and UVPT* files.
         */
        emars_grid = (emars_gridspec *)calloc(1,sizeof(emars_gridspec));

        emars = TRUE;
      }
      else if (strcmp(sflag,"-weizmann") == 0) {
        /*
         * Weizmann Institute gas-giant data.
         * Input a Weizmann .nc file and convert to isentropic-coordinate QB* and UVPT* files.
         */
        weizmann_grid = (weizmann_gridspec *)calloc(1,sizeof(weizmann_gridspec));

        weizmann = TRUE;
      }
      else if (strcmp(sflag,"-help") == 0 ||
               strcmp(sflag,"-h")    == 0) {
        /* Print help, exit: */
        system("more "EPIC_PATH"/help/epic_change.help");
        exit(1);
      }
      else {
        sprintf(Message,"Unrecognized change command-line flag: %s \n",sflag);
        epic_error(dbmsname,Message);
      }
    }
  }

  if (openmars || emars) {
    planet = &mars;
  }
  else if (weizmann) {
    planet = &jupiter;
  }
  else {
    /* Allocate memory */
    if((planet=( planetspec *)malloc(sizeof(planetspec))) == 0) {
      sprintf(Message,"allocating space for planetspec \n");
      epic_error(dbmsname,Message);
    }
  }

#if defined(EPIC_MPI)
  MPI_Init(&argc,&argv);
  para.comm = MPI_COMM_WORLD;
  MPI_Comm_set_errhandler(para.comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(para.comm,&para.iamnode);
  MPI_Comm_size(para.comm,&para.nproc);
  para.ndim = NINT(log((EPIC_FLOAT)para.nproc)/log(2.));
  /*
   * epic_change.c is not set up to be run across multiple processors.
   */
  if (para.nproc > 1) {
    sprintf(Message,"Not set up to be run across multiple processors");
    epic_error(dbmsname,Message);
  }
#endif

  /*
   *  Read in default parameter settings:
   */
  read_defaults(&defaults);

  /* Time-plane index. */
  time_index = 0;

  /*
   * Determine model size and allocate memory for arrays.
   */
  if (openmars) {
    int
      nc_id,nc_err;

    input_string("Input OpenMARS data file [netCDF format]\n",defaults.openmars_infile,openmars_infile);
    openmars_itime = 0;
    openmars_var_read(planet,openmars_grid,openmars_infile,SIZE_DATA,openmars_itime);
    openmars_make_arrays(planet,openmars_grid);
    /*
     * NOTE: For OpenMARS, the portion POST_SIZE_DATA refers to the 1D constant
     *       arrays for the dimensions, but not the 2D and 3D variable fields.  The latter
     *       portion is referred to by VAR_DATA. The reason for the two-part size-data
     *       partition for OpenMARS is that grid.nj depends on the latitude spacing in OpenMARS,
     *       which comes from the POST_SIZE_DATA segment.
     */
    openmars_var_read(planet,openmars_grid,openmars_infile,POST_SIZE_DATA,openmars_itime);

     /*
      * Use an EPIC file, openmars_epic.nc, to establish an appropriate EPIC Mars environment.
      * Check whether openmars_epic.nc already exists, and if not, create it.
      */
    nc_err = nc_open("./openmars_epic.nc",NC_NOWRITE,&nc_id);
    if (nc_err) {
      openmars_epic_nc(planet);
    }
    else {
      /*
       * The file openmars_epic.nc already exists.
       */
      nc_close(nc_id);
    }

    sprintf(defaults.infile,"./openmars_epic.nc");
    sprintf(infile,"%s",defaults.infile);
  }
  else if (emars) {
    int
      nc_id,nc_err;

    input_string("Input EMARS data file [netCDF format]\n",defaults.emars_infile,emars_infile);
    emars_itime = 0;
    emars_var_read(planet,emars_grid,emars_infile,SIZE_DATA,emars_itime);
    emars_make_arrays(planet,emars_grid);
    /*
     * NOTE: For EMARS, the portion POST_SIZE_DATA refers to the 1D constant
     *       arrays for the dimensions, but not the 2D and 3D variable fields.  The latter
     *       portion is referred to by VAR_DATA. The reason for the two-part size-data
     *       partition for EMARS is that grid.nj depends on the latitude spacing in EMARS,
     *       which comes from the POST_SIZE_DATA segment.
     */
    emars_var_read(planet,emars_grid,emars_infile,POST_SIZE_DATA,emars_itime);

     /*
      * Use an EPIC file, emars_epic.nc, to establish an appropriate EPIC Mars environment.
      * Check whether emars_epic.nc already exists, and if not, create it.
      */
    nc_err = nc_open("./emars_epic.nc",NC_NOWRITE,&nc_id);
    if (nc_err) {
      emars_epic_nc(planet);
    }
    else {
      /*
       * The file emars_epic.nc already exists.
       */
      nc_close(nc_id);
    }

    sprintf(defaults.infile,"./emars_epic.nc");
    sprintf(infile,"%s",defaults.infile);
  }
  else if (weizmann) {
    int
      nc_id,nc_err;

    input_string("Input Weizmann (gas giant) data file [netCDF format]\n",defaults.weizmann_infile,weizmann_infile);
    weizmann_itime = 0;
    weizmann_var_read(planet,weizmann_grid,weizmann_infile,SIZE_DATA,weizmann_itime);
    weizmann_make_arrays(planet,weizmann_grid);
    /*
     * NOTE: For Weizmann data, the portion POST_SIZE_DATA refers to the 1D constant
     *       arrays for the dimensions, but not the 2D variable fields.  The latter
     *       portion is referred to by VAR_DATA. The reason for the two-part size-data
     *       partition for Weizmann is that grid.nj depends on the latitude spacing,
     *       which comes from the POST_SIZE_DATA segment.
     */
    weizmann_var_read(planet,weizmann_grid,weizmann_infile,POST_SIZE_DATA,weizmann_itime);

     /*
      * Use an EPIC file, weizmann_epic.nc, to establish an appropriate EPIC gas-giant environment.
      * Check whether weizmann_epic.nc already exists, and if not, create it.
      */
    nc_err = nc_open("./weizmann_epic.nc",NC_NOWRITE,&nc_id);
    if (nc_err) {
      weizmann_epic_nc(planet);
    }
    else {
      /*
       * The file weizmann_epic.nc already exists.
       */
      nc_close(nc_id);
    }

    sprintf(defaults.infile,"./weizmann_epic.nc");
    sprintf(infile,"%s",defaults.infile);
  }
  else {
    input_string("Input EPIC file [netCDF format]\n",defaults.infile,infile);
  }

  /* NOTE: time_index is not used for SIZE_DATA */
  var_read(planet,infile,SIZE_DATA,time_index);

  set_var_props(planet);
  make_arrays(planet);

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    Buff2D[I] = fvector(0,Nelem2d-1,dbmsname);
  }

  /*
   * Read in rest of input data.
   */
  var_read(planet,infile,POST_SIZE_DATA,time_index);

  /* timeplane_bookkeeping() must come after reading in variables. */
  timeplane_bookkeeping();

  /* 
   * Set lon, lat, etc. 
   */
  set_lonlat();
  set_fmn(planet);
  set_gravity(planet);
  set_dsgth();

  /*
   * Set up thermodynamics subroutines.
   *
   * NOTE: The return value cpr is a low-temperature reference value, and
   *       should not be used otherwise.  Use return_cp() for a given
   *       thermodynamical state.
   */
  thermo_setup(planet,&planet->cpr);

  /* 
   * Store diagnostic variables. 
   */
  if (!openmars && !emars && !weizmann) {
    fprintf(stdout,"\nCalculating and storing diagnostic variables...");
  }

  /*
   * Reconstitute the prognostic variable H from input P.
   */

  /* Allocate memory */
  p = fvector(0,2*grid.nk+1,dbmsname);
  h = fvector(0,2*grid.nk+1,dbmsname);

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

  if (openmars) {
    openmars_conversion(planet,openmars_grid,openmars_infile,openmars_outfile_qb,openmars_outfile_uvpt,Buff2D);

    goto cleanup;
  }
  else if (emars) {
    emars_conversion(planet,emars_grid,emars_infile,emars_outfile_qb,emars_outfile_uvpt,Buff2D);

    goto cleanup;
  }
  else if (weizmann) {
    weizmann_conversion(planet,weizmann_grid,weizmann_infile,weizmann_outfile_qb,weizmann_outfile_uvpt,Buff2D);

    goto cleanup;
  }
  else {
    /*
     * No long jump (goto).
     */
    ;
  }

  /*
   * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
   *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
   *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
   */
  set_p2_etc(planet,UPDATE_THETA,Buff2D);
  store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,CALC_PHI3NK);
  store_diag(planet);

  /*
   * Print out a listing of important model parameters.
   */
  print_model_description(planet);

  /* 
   *  Print out vertical information:
   */
  print_vertical_column(planet,JLO,ILO,"vertical.dat");

  /*
   *  Change parameters as instructed.
   */
  grid.dt = input_int("\nInput timestep\n", grid.dt);

  /*
   * Inquire about radiation scheme.
   */
  inquire_radiation_scheme(planet);

  if (var.fpara.on) {
    /*
     * Inquire about fpara_rate_scaling.
     */
    sprintf(Message,"Input ortho-para H2 conversion rate scaling [nominal is 1.0]:\n");
    var.fpara_rate_scaling = input_float(Message,var.fpara_rate_scaling);
  }

  /*
   * Inquire whether to change the status of cloud microphysics.
   */
  if (grid.cloud_microphysics != OFF) {
    sprintf(Message,"Cloud microphysics: %2d => active  (latent heat, phase changes, precipitation), or \n"
                    "                    %2d => passive (advection only)\n"
                    "                    %2d => steady  (maintains starting condition)\n",ACTIVE,PASSIVE,STEADY);
    grid.cloud_microphysics = input_int(Message,grid.cloud_microphysics);
  }

  /*
   * Set sponges:
   */
  grid.k_sponge = input_int("Input k_sponge (-1 = no effect):\n",grid.k_sponge);
  grid.j_sponge = input_int("Input j_sponge (-1 = no effect; 3 is typical):\n",grid.j_sponge);

  /*
   * Revisit hyperviscosity coefficients (since dt may have changed).
   */
  set_hyperviscosity();
  
  /*
   * Add perturbations if requested, in the form of spots (vortices) and/or waves.
   */
  if (spots || waves) {
    EPIC_FLOAT
      *pert;

    /* 
     * Allocate memory for perturbation streamfunction
     */
    pert = fvector(0,Nelem3d-1,dbmsname);

    if (spots) {
      add_spots(planet,spots_file,pert,Buff2D);
    }

    if (waves) {
      add_waves(planet,waves_file,pert,Buff2D);
    }

    modify_progs_from_pert(planet,pert);

    /*
     * Free allocated memory.
     */
    free_fvector(pert,0,Nelem3d-1,dbmsname);

    /*
     * Update most commonly used diagnostic variables, in case they are needed.
     *
     * NOTE: For planet->type "terrestrial" and grid.coord_type == COORD_ISENTROPIC,
     *       need to calculate PHI3(KHI,J,I) on the grid.thetabot isentropic surface
     *       and call store_pgrad_vars with PASSING_PHI3NK; this is not yet implemented.
     */
    set_p2_etc(planet,UPDATE_THETA,Buff2D);
    store_pgrad_vars(planet,Buff2D,SYNC_DIAGS_ONLY,CALC_PHI3NK);
    store_diag(planet);
  }

  /*
   * Prompt for which variables to write to extract.nc.
   */
  defaults.extract_str[0] = '\0';
  for (index = FIRST_INDEX; index <= LAST_INDEX; index++) {
    if (var.extract_on_list[index] == LISTED_AND_ON) {
      sprintf(buffer," %d",index);
      strcat(defaults.extract_str,buffer);
    }
  }
  prompt_extract_on(defaults.extract_str,&grid.extract_species_fraction_type);

  /*
   * Write epic.nc file.
   */
  input_string("Name of output file\n",defaults.outfile,outfile);
  var_write(planet,outfile,ALL_DATA,time_index,stretch_ni);

  /*---------------*
   * Cleanup block *
   *---------------*/
  cleanup:

  /* Write defaults file: */
  write_defaults(&defaults);

  /*
   * Free dynamically allocated memory.
   */
  free_fvector(p,0,2*grid.nk+1,dbmsname);
  free_fvector(h,0,2*grid.nk+1,dbmsname);

  free_arrays(planet);
  free_var_props(planet);

  if (openmars) {
    openmars_free_arrays(planet,openmars_grid);
    free(openmars_grid);
  }
  else if (emars) {
    emars_free_arrays(planet,emars_grid);
    free(emars_grid);
  }
  else if (weizmann) {
    weizmann_free_arrays(planet,weizmann_grid);
    free(weizmann_grid);
  }

  if (openmars || emars || weizmann) {
    /*
     * No need to free planet structure, since it
     * was not dynamically allocated.
     */
    ;
  }
  else {
    free(planet);
  }

  for (I = 0; I < NUM_WORKING_BUFFERS; I++) {
    free_fvector(Buff2D[I],0,Nelem2d-1,dbmsname);
  }

  return 0;
}

/*======================= end of main() =====================================*/

/*======================= add_spots() =======================================*/

/*
 * Add vortices via a perturbation streamfunction.
 *
 * Example spots.dat file:
 * ---------------------------------------------------------------------------------------
 *  Number of vortices: 2
 *  lon[deg] lat[deg] press[hPa]  a[deg] b[deg] c_up[scale_hts] c_down[scale_hts] amp[m/s]
 *   30.      -33.     680.       3.0    2.5        2.5               3.0          100.
 *   60.      -33.5    680.       3.0    2.5        2.5               3.0          100.
 * ---------------------------------------------------------------------------------------
 *
 * NOTE: Not MPI ready. 
 */

/*
 * Implemented styles for the vortex perturbation streamfunction:
 */
#define GAUSSIAN_ELLIPSOID              0
#define POLYNOMIAL_GAUSSIAN_ELLIPSOID   1
#define POLYNOMIAL_GAUSSIAN_ORDER       2.0

/*
 * Choose streamfunction style from those listed above:
 */
#define SPOT_PERT  GAUSSIAN_ELLIPSOID

void add_spots(planetspec  *planet,
               char        *spots_file,
               EPIC_FLOAT  *pert,
               EPIC_FLOAT **Buff2D)
{
  register int
    K,J,I,
    kk,jj,
    ispot;
  int
    nspots = 0;
  EPIC_FLOAT
    rr,xspot,yspot,zspot,
    lon_width,lon_half_width,
    pressure;
  EPIC_FLOAT
    *lonspot,*latspot,*pspot,
    *aspot,*bspot,*cspot_up,*cspot_down,*ampspot;
  char
    buffer[FILE_STR];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_spots";

  if (strcmp(spots_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }

  /* Read vortex description file: */
  lon_width      = grid.globe_lontop-grid.globe_lonbot;
  lon_half_width = .5*lon_width;

  nspots = number_objects_in_file(spots_file);

  if (nspots == 1) {
    fprintf(stdout,"Generating spot..."); fflush(stdout);
  }
  else {
    fprintf(stdout,"Generating spots..."); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  lonspot      = fvector(0,nspots-1,dbmsname);
  latspot      = fvector(0,nspots-1,dbmsname);
  pspot        = fvector(0,nspots-1,dbmsname);
  aspot        = fvector(0,nspots-1,dbmsname);
  bspot        = fvector(0,nspots-1,dbmsname);
  cspot_up     = fvector(0,nspots-1,dbmsname);
  cspot_down   = fvector(0,nspots-1,dbmsname);
  ampspot      = fvector(0,nspots-1,dbmsname);

  /* 
   * Read in vortex information: 
   */
  read_spots_file(spots_file,ampspot,lonspot,latspot,pspot,
                  aspot,bspot,cspot_up,cspot_down,ADJUST_AMPLITUDE);

  /* 
   * Calculate streamfunction for vortices, and store in PERT.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        pressure = P2(K,J,I);
        for (ispot = 0; ispot < nspots; ispot++) {
          /* Account for periodicity in x-direction: */
          xspot  = (grid.lon[2*I+1]-lonspot[ispot]); 
          if (xspot > lon_half_width) { 
            xspot -= lon_width; 
          } else if (xspot < -lon_half_width) { 
            xspot += lon_width; 
          } 
          xspot /= aspot[ispot]; 

          yspot  = (grid.lat[2*J+1]-latspot[ispot])/bspot[ispot];

	  rr     = xspot*xspot+yspot*yspot;
          if (pressure <= pspot[ispot]){
            zspot = -log(pressure/pspot[ispot])/cspot_up[ispot];
	  }
          else{
            zspot =  log(pressure/pspot[ispot])/cspot_down[ispot];
          }
          rr += zspot*zspot; 
      
#if (SPOT_PERT == GAUSSIAN_ELLIPSOID)
          PERT(K,J,I) += ampspot[ispot]*exp(-rr);
#elif (SPOT_PERT == POLYNOMIAL_GAUSSIAN_ELLIPSOID)
          rr = sqrt(rr);
          if ( rr <= pow(2.0,1./POLYNOMIAL_GAUSSIAN_ORDER) ) {
            rr = pow(rr,POLYNOMIAL_GAUSSIAN_ORDER);
            PERT(K,J,I) += ampspot[ispot]*(1.+(1.-rr)*exp(-rr+2.))/(1.+exp(2.));
          }
#endif
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(pert,THREEDIM);

  fprintf(stdout,"done.\n"); fflush(stdout);

  /* Free allocated memory: */
  free_fvector(ampspot,     0,nspots-1, dbmsname);
  free_fvector(cspot_down,  0,nspots-1, dbmsname);
  free_fvector(cspot_up,    0,nspots-1, dbmsname);
  free_fvector(bspot,       0,nspots-1, dbmsname);
  free_fvector(aspot,       0,nspots-1, dbmsname);
  free_fvector(pspot,       0,nspots-1, dbmsname);
  free_fvector(latspot,     0,nspots-1, dbmsname);
  free_fvector(lonspot,     0,nspots-1, dbmsname);

  return;
}

#undef GAUSSIAN_ELLIPSOID
#undef POLYNOMIAL_GAUSSIAN_ELLIPSOID
#undef POLYNOMIAL_GAUSSIAN_ORDER

/*======================= end of add_spots() ================================*/

/*======================= add_waves() =======================================*/

/*
 * Add sinusoidal perturbation based on parameters in waves_file.
 * The vertical size c is in scale heights [sc_ht].
 * The amplitude [m2s-2] refers to the perturbation streamfunction, i.e., 
 * PERT = MONT for theta coordinates or PERT = PHI for pressure coordinates.
 *
 * Example waves_file:
 * ------------------------------------------------------------------
 * Number of waves: 2
 * lat[deg]   amp[m2s-2]  wn[]  p[hPa]  c[sc_ht]   fwhm[deg]
 *   30.0      0.02       1.0   700.0      1.5        1.50
 *   40.0      0.10       5.0   700.0      1.5        3.00
 * ------------------------------------------------------------------
 *
 * NOTE: wn[] is domain wavenumber, not planetary wavenumber.
 *
 * Tim Dowling and Raul Morales-Juberias, June 2022
 *
 * NOTE: Not MPI ready. 
 */

void add_waves(planetspec  *planet,
               char        *waves_file,
               EPIC_FLOAT  *pert,
               EPIC_FLOAT **Buff2D)
{
  int
    K,J,I,
    kkbot,jj,
    iwave;
  int
    nwaves = 0;
  long
    seed = -1;
  double
    xwave,ywave,zwave,kwave,
    rlt,rln,rr,
    phase_offset,
    coef;
  EPIC_FLOAT
   *latwave,
   *ampwave,
   *wnwave,
   *pwave,
   *cwave,
   *fwhmwave;
  char
    buffer[FILE_STR];
  FILE
   *waves;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="add_waves";

  if (strcmp(waves_file,"none") == 0) {
    /* Return if there is nothing to do: */
    return;
  }
		
  /* 
   * Read number of waves:
   */
  nwaves = number_objects_in_file(waves_file);

  if (nwaves == 1) {
    fprintf(stdout,"Generating wave..."); fflush(stdout);
  }
  else {
    fprintf(stdout,"Generating waves..."); fflush(stdout);
  }

  /* 
   * Allocate memory:
   */
  latwave  = fvector(0,nwaves-1,dbmsname);
  ampwave  = fvector(0,nwaves-1,dbmsname);
  wnwave   = fvector(0,nwaves-1,dbmsname);
  pwave    = fvector(0,nwaves-1,dbmsname);
  cwave    = fvector(0,nwaves-1,dbmsname);
  fwhmwave = fvector(0,nwaves-1,dbmsname);

  /*
   * Read wave information:
   */
  read_waves_file(waves_file,latwave,ampwave,wnwave,pwave,cwave,fwhmwave);

  /* 
   * Calculate streamfunction for waves, and store in PERT.
   */

  if (strcmp(grid.geometry,"globe") != 0) {
    sprintf(Message,"grid.geometry = %s not yet implemented",grid.geometry);
    epic_error(dbmsname,Message);
  }

  /* Set kk to be kkbot = 2*grid.nk */
  kkbot = 2*grid.nk;

  /* Set meridional radius of curvature to be mid-channel value; note: nj = 2*nj/2 */
  rlt = grid.rlt[kkbot][grid.nj];

  for (iwave = 0; iwave < nwaves; iwave++) {
    phase_offset = 2.*M_PI*random_number(&seed);

    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;

      /* Normalized y coordinate, centered on wave */
      ywave  = rlt*(grid.lat[jj]-latwave[iwave])*DEG;
      ywave /= fwhmwave[iwave];

      /* Zonal radius of curvature; varies with latitude */
      rln = grid.rln[kkbot][jj];

      /* Wavenumber k = 2pi/wavelength [1/m] */
      kwave = 2.*M_PI/( rln*(grid.globe_lontop-grid.globe_lonbot)*DEG/wnwave[iwave] );

      for (I = ILO; I <= IHI; I++) {
        /* x coordinate [m], origin at west end of channel */
        xwave = rln*(grid.lon[2*I+1]-grid.globe_lonbot)*DEG;
        for (K = KLO; K < KHI; K++) {
          /* NOTE: bottom layer is not perturbed. */
          /* Normalized z coordinate, centered on wave */
          zwave        = -log(P2(K,J,I)/pwave[iwave])/cwave[iwave];
          rr           = ywave*ywave+zwave*zwave;
          PERT(K,J,I) += ampwave[iwave]*exp(-rr)*sin(kwave*xwave+phase_offset);
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(pert,THREEDIM);

  fprintf(stdout,"done.\n"); fflush(stdout);
	
  /* 
   * Free allocated memory:
   */
  free_fvector(latwave, 0,nwaves-1, dbmsname);
  free_fvector(ampwave, 0,nwaves-1, dbmsname);
  free_fvector(wnwave,  0,nwaves-1, dbmsname);
  free_fvector(pwave,   0,nwaves-1, dbmsname);
  free_fvector(cwave,   0,nwaves-1, dbmsname);
  free_fvector(fwhmwave,0,nwaves-1, dbmsname);
	
  return;
}

/*======================= end of add_waves() ================================*/

/*======================= modify_progs_from_pert() ==========================*/

/*
 * Use gradient-wind balance to modify the prognostic variables
 * given a streamfunction perturbation.
 */

#undef  UG
#define UG(k,j,i) ug[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  VG
#define VG(k,j,i) vg[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  ZETAG
#define ZETAG(j,i) zetag[i+(j)*Iadim-Shift2d]

#undef  KING
#define KING(j,i) king[i+(j)*Iadim-Shift2d]

#undef  DTEMP
#define DTEMP(k,j,i) dtemp[i+(j)*Iadim+(k)*Nelem2d-Shift3d]

#undef  LAT_MIN
#define LAT_MIN 5.0

void modify_progs_from_pert(planetspec *planet,
                            EPIC_FLOAT *pert)
{
  EPIC_FLOAT
   *ug,*vg,
   *zetag,*king,
   *dtemp,
   *p_hybrid,*theta_hybrid,*theta_sigma;
  EPIC_FLOAT
    fpara,rgas,sigma,
    pbot,gsg,xi2,xi4,
    theta_ortho,theta_para;
  int
    K,J,I,
    kk,jj;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="modify_progs_from_pert";

  /* Allocate memory */
  ug    = fvector(0,Nelem3d-1,dbmsname);
  vg    = fvector(0,Nelem3d-1,dbmsname);
  dtemp = fvector(0,Nelem3d-1,dbmsname);

  zetag = fvector(0,Nelem2d-1,dbmsname);
  king  = fvector(0,Nelem2d-1,dbmsname);

  p_hybrid     = fvector(0,KHI,dbmsname);
  theta_hybrid = fvector(0,KHI,dbmsname);
  theta_sigma  = fvector(0,KHI,dbmsname);

  /*
   * Estimate the geostrophic wind (UG,VG) of the perturbed system.
   */
  for (K = KLO; K <= KHI; K++) {
    for (J = JLOPAD; J <= JHIPADPV; J++) {
      for (I = ILOPAD; I <= IHIPAD; I++) {
        UG(K,J,I) = U(grid.it_uv,K,J,I);
        VG(K,J,I) = V(grid.it_uv,K,J,I);
      }
    }
  }
  for (J = JLO+1; J <= JHI-1; J++) {
    jj = 2*J+1;
    if (fabs(grid.lat[jj]) > LAT_MIN) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          UG(K,J,I) -= (PERT(K,J+1,I-1)+PERT(K,J+1,I)-PERT(K,J-1,I-1)-PERT(K,J-1,I))*.25*grid.n[2*K][jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(ug,THREEDIM);

  for (J = JFIRST; J <= JHI; J++) {
    jj = 2*J;
    if (fabs(grid.lat[jj]) > LAT_MIN) {
      for (I = ILO; I <= IHI; I++) {
        for (K = KLO; K <= KHI; K++) {
          VG(K,J,I) += (PERT(K,J,I+1)+PERT(K,J-1,I+1)-PERT(K,J,I-1)-PERT(K,J-1,I-1))*.25*grid.m[2*K][jj]/grid.f[jj];
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(vg,THREEDIM);

  /* 
   * Modify P, THETA, FPARA as appropriate.
   */

  switch(grid.coord_type) {
    case COORD_HYBRID:
      /*
       * For the hybrid coordinate, this is a challenging task, and the current approach is
       * approximate (see Dowling's notes, 12/2/10).
       *
       * Note that hydrostatic balance in isentropic and isobaric coordinates can be written:
       *
       *      ____isentropic____         ____isobaric____
       *
       *     dM/dlog(theta) = cp T       dPhi/dlog(p) = - R T
       *
       * These are quite similar in form.  Let <T> be the unperturbed temperature.  Then,
       *
       *     T = <T> + dPERT/dxi ,
       *
       * where PERT is the perturbation to either the Montgomery potential or the geopotential, as set above,
       * and xi is either cp log(theta) or -R log(p), respectively.  Our strategy is to define xi to be a blend,
       * using the transition function g(sigma):
       *
       *     xi = g cp log(theta) + (1.-g)(-R log(p)) ,
       *
       * where sigma, theta, and p are all the unperturbed values.
       */
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot = P3(KHI,J,I);
          for (K = KLO; K < KHI; K++) {
            /*
             * Form the independent variable, xi, by blending from the pressure region to the theta region.
             */
            rgas         = R_GAS/avg_molar_mass(planet,2*K+1,J,I);

            sigma        = get_sigma(pbot,P2(K,J,I));
            gsg          = (double)g_sigma(sigma);
            xi2          = gsg*planet->cp*log(THETA2(K,J,I))+(1.-gsg)*(-rgas*log(P2(K,J,I)));

            sigma        = get_sigma(pbot,P2(K+1,J,I));
            gsg          = (double)g_sigma(sigma);
            xi4          = gsg*planet->cp*log(THETA2(K+1,J,I))+(1.-gsg)*(-rgas*log(P2(K+1,J,I)));

            /* 
             * DTEMP(K,J,I) is on the p3-grid.
             */
            DTEMP(K,J,I) = (PERT(K,J,I)-PERT(K+1,J,I))/(xi2-xi4);
          }
        }
      }

      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          pbot = P3(KHI,J,I);
          for (K = KLO; K < grid.k_sigma; K++) {
            T3(K,J,I) += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(planet,T3(K,J,I),RHO3(K,J,I),avg_molar_mass(planet,2*K+1,J,I));

            /*
             * Use diagnostic value for THETA.
             */
            sigma        = get_sigma(pbot,P3(K,J,I));
            gsg          = (double)g_sigma(sigma);
            THETA(K,J,I) = (EPIC_FLOAT)(((double)grid.sigmatheta[2*K+1]-(double)f_sigma(sigma))/gsg);
          }
          for (K = grid.k_sigma; K < KHI; K++) {
            T3(K,J,I) += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(planet,T3(K,J,I),RHO3(K,J,I),avg_molar_mass(planet,2*K+1,J,I));

            /*
             * Use thermodynamic value for THETA.
             */
            THETA(K,J,I) = return_theta(planet,fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
          }
        }
      }
    break;
    case COORD_ISOBARIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            rgas          = R_GAS/avg_molar_mass(planet,2*K+1,J,I);

            DTEMP(K,J,I)  = (PERT(K,J,I)-PERT(K+1,J,I))/(rgas*log(P2(K+1,J,I)/P2(K,J,I)));
            T3(K,J,I)    += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = return_fpe(T3(K,J,I));
            }

            THETA(K,J,I)  = return_theta(planet,fpara,P3(K,J,I),T3(K,J,I),&theta_ortho,&theta_para);
          }
        }
      }
    break;
    case COORD_ISENTROPIC:
      for (J = JLOPAD; J <= JHIPAD; J++) {
        for (I = ILOPAD; I <= IHIPAD; I++) {
          for (K = KLO; K < KHI; K++) {
            DTEMP(K,J,I)  = (PERT(K,J,I)-PERT(K+1,J,I))/(planet->cp*log(THETA(K,J,I)/THETA(K+1,J,I)));
            T3(K,J,I)    += DTEMP(K,J,I);

            /*
             * Estimate para-hydrogen fraction with equilibrium value.
             */
            fpara = return_fpe(T3(K,J,I));
            if (var.fpara.on) {
              FPARA(K,J,I) = fpara;
            }

            /*
             * NOTE: Assuming the density doesn't change.
             */
            P3(K,J,I) = p_from_t_rho_mu(planet,T3(K,J,I),RHO3(K,J,I),avg_molar_mass(planet,2*K+1,J,I));
          }
        }
      }
    break;
    default:
      sprintf(Message,"grid.coord_type=%d not yet implemented",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  for (K = KLO; K < KHI; K++) {
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          U(grid.it_uv,K,J,I) = UG(K,J,I);
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          V(grid.it_uv,K,J,I) = VG(K,J,I);
        }
      }
    }
  }
  /* Need to apply bc_lateral() here. */
  bc_lateral(var.u.value+grid.it_uv*Nelem3d,THREEDIM);
  bc_lateral(var.v.value+grid.it_uv*Nelem3d,THREEDIM);

  /*
   * The gradient-balance correction is based on the paper:
   *     McIntyre and Roulstone, 2002, Large-Scale Atmosphere-Ocean Dynamics. II Geometric Methods
   *         and Models. Cambridge Univ. Press, Ch. 8.
   */
  for (K = KLO; K < KHI; K++) {
    kk = 2*K;
    /*
     * Calculate geostrophic relative vorticity and kinetic energy per mass.
     */
    vorticity(ON_SIGMATHETA,RELATIVE,2*K,ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,NULL,zetag);
    for (J = JLO; J <= JHI; J++) {
      for (I = ILO; I <= IHI; I++) {
        KING(J,I) = get_kin(planet,ug+(K-Kshift)*Nelem2d,vg+(K-Kshift)*Nelem2d,kk,J,I);
      }
    }
    /* Need to apply bc_lateral() here. */
    bc_lateral(king,TWODIM);
  
    for (J = JLO; J <= JHI; J++) {
      jj = 2*J+1;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          UG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J+1,I))*UG(K,J,I)
                        -.25*grid.n[kk][jj]*(KING(J+1,I)+KING(J+1,I-1)-KING(J-1,I)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
    for (J = JFIRST; J <= JHI; J++) {
      jj = 2*J;
      if (fabs(grid.lat[jj]) > LAT_MIN) {
        for (I = ILO; I <= IHI; I++) {
          VG(K,J,I) += (-.5*(ZETAG(J,I)+ZETAG(J,I+1))*VG(K,J,I)
                        +.25*grid.m[kk][jj]*(KING(J,I+1)+KING(J-1,I+1)-KING(J,I-1)-KING(J-1,I-1)))/grid.f[jj];
        }
      }
    }
  }

  /* Free allocated memory. */
  free_fvector(ug,   0,Nelem3d-1,dbmsname);
  free_fvector(vg,   0,Nelem3d-1,dbmsname);
  free_fvector(dtemp,0,Nelem3d-1,dbmsname);

  free_fvector(zetag,0,Nelem2d-1,dbmsname);
  free_fvector(king, 0,Nelem2d-1,dbmsname);

  free_fvector(p_hybrid,    0,KHI,dbmsname);
  free_fvector(theta_hybrid,0,KHI,dbmsname);
  free_fvector(theta_sigma, 0,KHI,dbmsname);

  return;
}

/*======================= end of modify_progs_from_pert() ===================*/

/*====================== number_objects_in_file() ===========================*/

/*
 * Retrieve the integer after the first colon, ':', which
 * is interpreted as the number of objects (spots, waves, etc).
 * in the given parameter file.
 */

int number_objects_in_file(char *objects_file)
{
  int 
    n_objects = 0;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *objects;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="number_objects_in_file";

  objects = fopen(objects_file,"r");
  if (!objects) {
    sprintf(Message,"unable to open %s \n",objects_file);
    epic_error(dbmsname,Message);
  }
  while (n_objects == 0) {
    fgets(buffer,FILE_STR,objects);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&n_objects);
    }
  }
  fclose(objects);

  return n_objects;
}

/*====================== end of number_objects_in_file() ====================*/

/*====================== read_spots_file() ==================================*/

void read_spots_file(char       *spots_file,
                     EPIC_FLOAT *ampspot,
                     EPIC_FLOAT *lonspot,
                     EPIC_FLOAT *latspot,
                     EPIC_FLOAT *pspot,
                     EPIC_FLOAT *aspot,
                     EPIC_FLOAT *bspot,
                     EPIC_FLOAT *cspot_up,
                     EPIC_FLOAT *cspot_down,
                     int         adjust_amplitude)
{
  register int 
    ispot;
  int
    nspots=0;
  EPIC_FLOAT
    fspot,
    factor;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *spots;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_spots_file";

  spots = fopen(spots_file,"r");
  while (nspots == 0) {
    fgets(buffer,FILE_STR,spots);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nspots);
    }
  }

  fgets(buffer,FILE_STR,spots);
  for (ispot = 0; ispot < nspots; ispot++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(spots,"%lf %lf %lf %lf %lf %lf %lf %lf",
           lonspot+ispot,latspot+ispot,pspot+ispot,
           aspot+ispot,bspot+ispot,cspot_up+ispot,cspot_down+ispot,ampspot+ispot);
#else
    fscanf(spots,"%f %f %f %f %f %f %f %f",
           lonspot+ispot,latspot+ispot,pspot+ispot,
           aspot+ispot,bspot+ispot,cspot_up+ispot,cspot_down+ispot,ampspot+ispot);
#endif
    /* Convert pspot from hPa to Pa: */
    pspot[ispot] *= 100.;


    /* 
     * Convert ampspot to amp for mont.  The parameter 'factor' 
     * should make the maximum spot velocity close to the input 
     * ampspot[m/s] for a gaussian spot perturbation streamfunction.
     */
    if (adjust_amplitude == ADJUST_AMPLITUDE) {
      fspot           = 2.*planet->omega_sidereal*sin(latspot[ispot]*DEG);
      factor          = 1.166;
      ampspot[ispot] *= factor*bspot[ispot]*DEG*planet->re*fabs(fspot);
    }

    fgets(buffer,FILE_STR,spots);
  }
  fclose(spots);

  return;

}

/*====================== end of read_spots_file() ===========================*/

/*====================== read_waves_file() ==================================*/

void read_waves_file(char       *waves_file,
                     EPIC_FLOAT *latwave,
                     EPIC_FLOAT *ampwave,
                     EPIC_FLOAT *wnwave,
                     EPIC_FLOAT *pwave,
                     EPIC_FLOAT *cwave,
                     EPIC_FLOAT *fwhmwave)
{
  register int 
    iwave;
  int
    kkbot,
    nwaves = 0;
  EPIC_FLOAT
    coef,rlt;
  char
    *char_pt,
    buffer[FILE_STR];
  FILE
    *waves;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_waves_file";

  /* 
   * Coefficient for using FWHM (full-width at half-maximum) to specify the width of a 1D gaussian.
   */
  coef = 2.*sqrt(log(2.));

  /* Set kk to be kkbot = 2*grid.nk */
  kkbot = 2*grid.nk;

  /* Set meridional radius of curvature to be mid-channel value; note: nj = 2*nj/2 */
  rlt = grid.rlt[kkbot][grid.nj];

  waves = fopen(waves_file,"r");
  while (nwaves == 0) {
    fgets(buffer,FILE_STR,waves);
    char_pt = strchr(buffer,':');
    if (char_pt) {
      sscanf(char_pt+1,"%d",&nwaves);
    }
  }

  fgets(buffer,FILE_STR,waves);
  for (iwave = 0; iwave < nwaves; iwave++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(waves,"%lf %lf %lf %lf %lf %lf",
                 latwave+iwave,ampwave+iwave,wnwave+iwave,pwave+iwave,cwave+iwave,fwhmwave+iwave);
#else
    fscanf(waves,"%f %f %f %f %f %f",
                 latwave+iwave,ampwave+iwave,wnwave+iwave,pwave+iwave,cwave+iwave,fwhmwave+iwave);
#endif
		
    /* Convert pwave from hPa (mbar) to Pa: */
    pwave[iwave] *= 100.;

    /* Convert FWHM for use in 1D gaussian. */
    fwhmwave[iwave] *= rlt*DEG/coef;

    fgets(buffer,FILE_STR,waves);
  }

  fclose(waves);

  return;
}

/*====================== end of read_waves_file() ===========================*/

/*======================= read_defaults() ===================================*/

void read_defaults(change_defaultspec *def) 
{
  int
    nc_id,
    nc_err,
    index;
  char
    min_element[4];
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_defaults";

  nc_err = lookup_netcdf("change_defaults.nc",
                         &nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err == NC_NOERR) {
    READI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
    READC(def->infile,def_infile,N_STR);
    READC(def->outfile,def_outfile,N_STR);
    READC(def->openmars_infile,def_openmars_infile,N_STR);
    READC(def->emars_infile,def_emars_infile,N_STR);
    READC(def->weizmann_infile,def_weizmann_infile,N_STR);
    READC(def->extract_str,def_extract_str,N_STR);
  }
  else {
    /*
     * If the file is not readable, use standard defaults.
     */
    def->uv_timestep_scheme = 0;

    strcpy(def->infile,"epic.nc");
    strcpy(def->outfile,"epic.nc");
    strcpy(def->openmars_infile,"");
    strcpy(def->emars_infile,"");
    strcpy(def->weizmann_infile,"");
    strcpy(def->extract_str,"");
  }

  return;
}

/*======================= end of read_defaults() ============================*/

/*======================= write_defaults() ==================================*/

void write_defaults(change_defaultspec *def)
{
  int
    nc_id,
    nc_err;
  nc_type
    the_nc_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="write_defaults";

  nc_err = nc_create("change_defaults.nc",NC_CLOBBER,&nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  WRITEI(&def->uv_timestep_scheme,def_uv_timestep_scheme,1);
  WRITEC(def->infile,def_infile,N_STR);
  WRITEC(def->outfile,def_outfile,N_STR);
  WRITEC(def->openmars_infile,def_openmars_infile,N_STR);
  WRITEC(def->emars_infile,def_emars_infile,N_STR);
  WRITEC(def->weizmann_infile,def_weizmann_infile,N_STR);
  WRITEC(def->extract_str,def_extract_str,N_STR);

  nc_close(nc_id);

  return;
}

/*======================= end of write_defaults() ===========================*/

/* * * * * * * * * * * * end of epic_change.c * * * * * * * * * * * * * * * * */

