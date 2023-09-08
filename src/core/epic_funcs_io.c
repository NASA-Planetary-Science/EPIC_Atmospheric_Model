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

/* * * * * * * * *  epic_funcs_io.c  * * * * * * * * * * * * * * * * 
 *                                                                 *
 *       Functions for EPIC model input and output.                *
 *       This file includes the following:                         *
 *                                                                 *
 *           setup_read_array()                                    *
 *           setup_write_array()                                   *
 *           read_array()                                          *
 *           write_array()                                         *
 *           var_read()                                            *
 *           var_write()                                           *
 *           lookup_netcdf()                                       *
 *           define_netcdf()                                       *
 *           get_jlohi()                                           *
 *           get_ilohi()                                           *
 *           prompt_extract_on()                                   *
 *           prompt_species_on()                                   *
 *           bcast_char(),bcast_int(),                             *
 *           bcast_float(),bcast_double()                          *
 *           read_spacing_file()                                   *
 *           get_sounding()                                        *
 *           read_t_vs_p()                                         *
 *           read_meridional_plane()                               *
 *           inquire_radiation_scheme()                            *
 *           input_float()                                         *
 *           input_int(),input_string()                            *
 *           print_model_description()                             *
 *           print_zonal_info()                                    *
 *           print_vertical_column()                               *
 *           node0_barrier_print()                                 *
 *           scdswap()                                             *
 *           epic_error()                                          *
 *           epic_warning()                                        *
 *           declare_copyright()                                   *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*======================= setup_read_array() ================================*/
/*
 * A function that uses read_array() must first call this setup function.
 *
 * Returns the number of processors running the model.
 */
int setup_read_array(void)
{
  int
    node,
   *jlo,
   *jhi,
   *ilo,
   *ihi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_read_array";

  /* 
   * Allocate memory.
   */
  jlo = ivector(0,grid.we_num_nodes-1,dbmsname);
  jhi = ivector(0,grid.we_num_nodes-1,dbmsname);
  ilo = ivector(0,grid.we_num_nodes-1,dbmsname);
  ihi = ivector(0,grid.we_num_nodes-1,dbmsname);

  /*
   * Exchange jlo, jhi and ilo, ihi information.
   */
  for (node = 0; node < grid.we_num_nodes; node++) {
    if (node == IAMNODE) {
      jlo[node] = JLO;
      jhi[node] = JHI;
      ilo[node] = ILO;
      ihi[node] = IHI;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(jlo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(jhi+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(ilo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(ihi+node,1,MPI_INT,node,para.comm);
#endif

  }

  /*
   * Setup get_jlohi() and get_ilohi().
   */
  get_jlohi(SETUP_GET_JLOHI,grid.we_num_nodes,jlo,jhi);
  get_ilohi(SETUP_GET_ILOHI,grid.we_num_nodes,ilo,ihi);

  /*
   * Free allocated memory.
   */
  free_ivector(jlo,0,grid.we_num_nodes-1,dbmsname);
  free_ivector(jhi,0,grid.we_num_nodes-1,dbmsname);
  free_ivector(ilo,0,grid.we_num_nodes-1,dbmsname);
  free_ivector(ihi,0,grid.we_num_nodes-1,dbmsname);

  return grid.we_num_nodes;
}

/*======================= end of setup_read_array() =========================*/

/*======================= setup_write_array() ===============================*/
/*
 * A function that calls write_array() must first call this setup function.
 *
 * Returns the number of processors running the model.
 */

int setup_write_array(void)
{
  int
    node;
  static int
    initialized = FALSE,
   *jlo,
   *jhi,
   *ilo,
   *ihi;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="setup_write_array";

  if (!initialized) {

    /* 
     * Allocate memory.
     */
    jlo = ivector(0,grid.we_num_nodes-1,dbmsname);
    jhi = ivector(0,grid.we_num_nodes-1,dbmsname);
    ilo = ivector(0,grid.we_num_nodes-1,dbmsname);
    ihi = ivector(0,grid.we_num_nodes-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Exchange jlo, jhi and ilo, ihi information.
   */
  for (node = 0; node < grid.we_num_nodes; node++) {
    if (node == IAMNODE) {
      jlo[node] = JLO;
      jhi[node] = JHI;
      ilo[node] = ILO;
      ihi[node] = IHI;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(jlo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(jhi+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(ilo+node,1,MPI_INT,node,para.comm);
    MPI_Bcast(ihi+node,1,MPI_INT,node,para.comm);
#endif

  }

  /*
   * Setup get_jlohi() and get_ilohi().
   */
  get_jlohi(SETUP_GET_JLOHI,grid.we_num_nodes,jlo,jhi);
  get_ilohi(SETUP_GET_ILOHI,grid.we_num_nodes,ilo,ihi);
 
  return grid.we_num_nodes;
}

/*======================= end of setup_write_array() ========================*/

/*======================= read_array() ======================================*/

/*
 *   Basic steps:
 *
 *   1) NODE0 reads subarray data from the file into its buff_subarray. 
 *
 *   2) NODE0 sends data from its buff_subarray to the target node's buff_subarray.
 *
 *   3) The target node transfers data from its buff_subarray to the appropriate place.
 *
 * A function that uses read_array() must first call setup_read_array().
 * Call from all nodes. 
 *
 * For array types that have striped data (array_type != EPIC_FLOAT_ARRAY), pass
 * the name of the first float array that will be read in, and the naming
 * convention will be reconstructed. 
 *
 * NOTE: Does not apply boundary conditions.
 */

void read_array(int   node,
                int   dim,
                int  *start,
                int  *end,
                char *name,
                int   index,
                void *array,
                int   array_type,
                int   nc_id)
{
  register int
    n,
    K,J,I,
    nlo,nhi,n_len,nkji_len,
    klo,khi,k_len,kji_len,
    jlo,jhi,j_len,ji_len,
    ilo,ihi,i_len,
    shift_buff,offset,
    i_bytes,
    al=FOURDIM-dim;
  int
    nc_varid,
    nc_err;
  char
    the_name[VAR_NM_SZ];
  size_t
    nc_start[FOURDIM], 
    nc_count[FOURDIM];
  EPIC_FLOAT
    *buff_subarray,
    *epic_float_array;

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
    dbmsname[]="read_array";

  if (IAMNODE != node  &&
      IAMNODE != NODE0) {
    /*
     * Return if not a participating node.
     */
    return;
  }

  /*
   * Cast pointers based on array_type.
   */
  nlo  = 0;
  if (array_type == EPIC_FLOAT_ARRAY) {
    epic_float_array = (EPIC_FLOAT *)array;
    nhi              = 1-1;
  }
  else {
    sprintf(Message,"unrecognized array_type=%d",array_type);
    epic_error(dbmsname,Message);
  }

  khi = klo = 0;
  jhi = jlo = 0;
  ihi = ilo = 0;

  if (dim >= TWODIM) {
    ilo = start[0];
    ihi = end[  0]; 

    jlo = start[1];
    jhi = end[  1];
   
    if (dim >= THREEDIM) {
      klo = start[2];
      khi = end[  2];
    }
  }
  else {
    sprintf(Message,"dim = %d not recognized",dim);
    epic_error(dbmsname,Message);
  }

  n_len      = nhi-nlo+1;
  k_len      = khi-klo+1;
  j_len      = jhi-jlo+1;
  i_len      = ihi-ilo+1;
  i_bytes    = i_len*sizeof(EPIC_FLOAT);
  ji_len     = j_len*i_len;
  kji_len    = k_len*ji_len;
  nkji_len   = n_len*kji_len;
  shift_buff = ilo+(jlo)*i_len+(klo)*ji_len+(nlo)*kji_len;

  /*
   * Allocate memory for buff_subarray:
   */
  buff_subarray = fvector(0,nkji_len-1,dbmsname);

 /*
  * 1) NODE0 reads subarray data and stores it in buff_subarray.
  */
  if (IAMNODE == NODE0) {
    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;  
    nc_count[NETCDF_I_INDEX] = i_len; 

    nc_start[NETCDF_T_INDEX] = start[3];
    nc_start[NETCDF_I_INDEX] = ilo-grid.ilo;

    for (n = nlo; n <= nhi; n++) {
      /*
       * Get the netCDF variable ID.
       * First, reconstruct its name.
       */
      if (array_type == EPIC_FLOAT_ARRAY) {
        strcpy(the_name,name);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
      nc_err = nc_inq_varid(nc_id,the_name,&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      for (K = klo; K <= khi; K++) {
        nc_start[NETCDF_K_INDEX] = K-klo;
        for (J = jlo; J <= jhi; J++) {
          nc_start[NETCDF_J_INDEX] = J-grid.jlo;
 
#if EPIC_PRECISION == DOUBLE_PRECISION 
          nc_err = nc_get_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,
                                      &BUFF_SUBARRAY(n,K,J,ilo));
#else
          nc_err = nc_get_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,
                                      &BUFF_SUBARRAY(n,K,J,ilo));
#endif

          if (nc_err != NC_NOERR) {
            sprintf(Message,"%s, K,J=%d,%d; %s",name,K,J,nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
        }
      }
    }
  }

 /*
  * 2) NODE0 sends data from its buff_subarray to the target node's buff_subarray.
  */

  if (node != NODE0) {
    /*
     * Send data from NODE0 to node.
     */

#if defined(EPIC_MPI)
    if (IAMNODE == NODE0) {
      MPI_Send(buff_subarray,nkji_len,float_type,node,index,para.comm);
    }
    else if (IAMNODE == node) {
      int
        count;
      MPI_Status
        status;

      MPI_Recv(buff_subarray,nkji_len,float_type,NODE0,index,para.comm,&status);

      /* 
       * Verify number of items received.
       */
      MPI_Get_count(&status,float_type,&count);
      if (count != nkji_len) {
        sprintf(Message,"count=%d != nkji_len=%d",count,nkji_len);
        epic_error(dbmsname,Message);
      }
    }
#endif

  }
  /*
   * Free allocated memory and return if not the target node.
   */
  if (IAMNODE != node) {
    free_fvector(buff_subarray,0,nkji_len-1,dbmsname);
    return;
  }

 /*
  * 3) The target node transfers data from its buff_subarray to the appropriate place.
  */
  for (K = klo; K <= khi; K++) {
    for (J = jlo; J <= jhi; J++) {
      if (array_type == EPIC_FLOAT_ARRAY) {
        /*
         * Use a fast string copy, since the data are contiguous.
         */
        offset = ilo+(J)*Iadim+(K)*Nelem2d-Shift3d;
        memcpy(epic_float_array+offset,&BUFF_SUBARRAY(0,K,J,ilo),i_bytes);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
    }
  }

  /*
   * Free allocated memory.
   */
  free_fvector(buff_subarray,0,nkji_len-1,dbmsname);

  return;
}

/*======================= end of read_array() ===============================*/

/*======================= write_array() =====================================*/

/*
 *   Basic steps:
 *
 *   1) Source node transfers data from the appropriate place to its buff_subarray.
 *
 *   2) Source node sends data from its buff_subarray to NODE0's buff_subarray.
 *
 *   3) NODE0 writes source node's data from its buff_subarray to the file.
 *
 * A function that uses write_array() must first call setup_write_array().
 * Call from all nodes.
 *
 * If stretch_ni > 0, stretch KJ-plane zonal average into a KJI cube.
 * If stretch_ni < 0, suppress the I (zonal) dimension (as in a KJ plane).
 *
 * For array types that have striped data (array_type != EPIC_FLOAT_ARRAY), pass
 * the name of the first float array that will be read in, and the naming
 * convention will be deduced from this.
 *
 */

void write_array(int    node,
                 int    dim,
                 int   *start,
                 int   *end,
                 int    stretch_ni,
                 char  *name,
                 int    index,
                 void  *array,
                 int    array_type,
                 int    nc_id)
{
  register int
    n,K,J,I,
    nlo,nhi,n_len,nkji_len,
    klo,khi,k_len,kji_len,
    jlo,jhi,j_len,ji_len,
    ilo,ihi,i_len,
    shift_buff,offset,
    i_bytes,
    al=FOURDIM-dim;
  int
    nc_varid,
    nc_err;
  char
    the_name[VAR_NM_SZ];
  size_t
    nc_start[FOURDIM], 
    nc_count[FOURDIM];
  EPIC_FLOAT
    *buff_subarray,
    *buffer,
    *epic_float_array,
     tmp;

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
    dbmsname[]="write_array";

  /*
   * Return if not a participating node.
   */
  if (IAMNODE != node && 
      IAMNODE != NODE0) {
    return;
  }

  /*
   * Cast pointers based on array_type.
   */
  nlo  = 0;
  if (array_type == EPIC_FLOAT_ARRAY) {
    epic_float_array = (EPIC_FLOAT *)array;
    nhi              = 1-1;
  }
  else {
    sprintf(Message,"unrecognized array_type=%d",array_type);
    epic_error(dbmsname,Message);
  }

  khi = klo = 0;
  jhi = jlo = 0;
  ihi = ilo = 0;

  if (dim >= TWODIM) {
    ilo = start[0];
    ihi = end[  0]; 
    jlo = start[1];
    jhi = end[  1];
    if (dim >= THREEDIM) {
      klo = start[2];
      khi = end[  2];
    }
  }
  else {
    sprintf(Message,"dim = %d not recognized",dim);
    epic_error(dbmsname,Message);
  }

  n_len      = nhi-nlo+1;
  k_len      = khi-klo+1;
  j_len      = jhi-jlo+1;
  if (stretch_ni < 0) {
    /* I (zonal) dimension suppressed */
    i_len = 1;
  }
  else {
    i_len = ihi-ilo+1;
  }
  i_bytes    = i_len*sizeof(EPIC_FLOAT);
  ji_len     = j_len*i_len;
  kji_len    = k_len*ji_len;
  nkji_len   = n_len*kji_len;
  shift_buff = ilo+(jlo)*i_len+(klo)*ji_len+(nlo)*kji_len;

  /*
   * Allocate memory.
   */
  buff_subarray = fvector(0,nkji_len-1,dbmsname);
  if (stretch_ni > 0) {
    buffer = fvector(0,stretch_ni-ilo,dbmsname);
  }

  /*
   * 1) node obtains subarray data and stores it in buff_subarray.
   */
  if (IAMNODE == node) {
    for (K = klo; K <= khi; K++) {
      for (J = jlo; J <= jhi; J++) {
        if (array_type == EPIC_FLOAT_ARRAY) {
          if (stretch_ni < 0) {
            /*
             * I (zonal) dimension is suppressed.
             */
            offset = J+(K)*Jadim-Shiftkj;
          }
          else {
            /*
             * Use a fast string copy, since the data are contiguous.
             */
            offset = ilo+(J)*Iadim+(K)*Nelem2d-Shift3d;
          }
          memcpy(&BUFF_SUBARRAY(0,K,J,ilo),epic_float_array+offset,i_bytes);
        }
        else {
          sprintf(Message,"unrecognized array_type=%d",array_type);
          epic_error(dbmsname,Message);
        }
      }
    }
  }

  /*
   * 2) node sends data from its buff_subarray to NODE0's buff_subarray.
   */

#if defined(EPIC_MPI)
  if (node != NODE0) {
    if (IAMNODE == node) {
      MPI_Send(buff_subarray,nkji_len,float_type,NODE0,index,para.comm);
    }
    else if (IAMNODE == NODE0) {
      int
        count;
      MPI_Status
        status;

      MPI_Recv(buff_subarray,nkji_len,float_type,node,index,para.comm,&status);

      /* 
       * Verify number of items received.
       */
      MPI_Get_count(&status,float_type,&count);
      if (count != nkji_len) {
        sprintf(Message,"count=%d != nkji_len=%d",count,nkji_len);
        epic_error(dbmsname,Message);
      }
    }
  }
#endif

  /*
   * 3) NODE0 writes data from its buff_subarray.
   */
  if (IAMNODE == NODE0) {
    nc_count[NETCDF_T_INDEX] = 1;
    nc_count[NETCDF_K_INDEX] = 1; 
    nc_count[NETCDF_J_INDEX] = 1;

    if (stretch_ni > 0) {
      nc_count[NETCDF_I_INDEX] = stretch_ni-grid.ilo+1;
    }
    else {
      nc_count[NETCDF_I_INDEX] = ihi-ilo+1;
    }

    nc_start[NETCDF_T_INDEX] = start[3];
    nc_start[NETCDF_I_INDEX] = ilo-grid.ilo;

    for (n = nlo; n <= nhi; n++) {
      /*
       * Get the netCDF variable ID.
       * First, reconstruct its name.
       */
      if (array_type == EPIC_FLOAT_ARRAY) {
        strcpy(the_name,name);
      }
      else {
        sprintf(Message,"unrecognized array_type=%d",array_type);
        epic_error(dbmsname,Message);
      }
      nc_err = nc_inq_varid(nc_id,the_name,&nc_varid);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s, %s",the_name,nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }

      for (K = klo; K <= khi; K++) {
        nc_start[NETCDF_K_INDEX] = K-klo;
        for (J = jlo; J <= jhi; J++) {
          nc_start[NETCDF_J_INDEX] = J-grid.jlo;

#if EPIC_PRECISION == DOUBLE_PRECISION
          if (stretch_ni > 0) {
            if (ilo == grid.ilo) {
              /*
               * Assume zonal symmetry and copy the grid.ilo value into an I buffer.
               */
              tmp = BUFF_SUBARRAY(n,K,J,ilo);
              for (I = ilo; I <= stretch_ni; I++) {
                buffer[I-ilo] = tmp;
              }
              nc_err = nc_put_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,buffer);
            }
            else{
              nc_err = NC_NOERR;
            }
          }
          else {
            nc_err = nc_put_vara_double(nc_id,nc_varid,nc_start+al,nc_count+al,&BUFF_SUBARRAY(n,K,J,ilo));
          }
          if (nc_err != NC_NOERR) {
            sprintf(Message,"nc_put_vara_double(),%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
#else
          if (stretch_ni > 0) {
            if (ilo == grid.ilo) {
              /*
               * Assume zonal symmetry and copy the grid.ilo value into an I buffer.
               */
              tmp = BUFF_SUBARRAY(n,K,J,ilo);
              for (I = ilo; I <= stretch_ni; I++) {
                buffer[I-ilo] = tmp;
              }
              nc_err = nc_put_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,buffer);
            }
            else {
              nc_err = NC_NOERR;
            }
          }
          else {
            nc_err = nc_put_vara_float(nc_id,nc_varid,nc_start+al,nc_count+al,&BUFF_SUBARRAY(n,K,J,ilo));
          }
          if (nc_err != NC_NOERR) {
            sprintf(Message,"nc_put_vara_float(),%s",nc_strerror(nc_err));
            epic_error(dbmsname,Message);
          }
#endif
        }
      }
    }
  }

  /*
   * Free allocated memory for buff_subarray:
   */
  free_fvector(buff_subarray,0,nkji_len-1,dbmsname);
  if (stretch_ni > 0) {
    free_fvector(buffer,0,stretch_ni-ilo,dbmsname);
  }

  return;
}

/*======================= end of write_array() ==============================*/

/*======================= var_read() ========================================*/

/*
 * Read in the variables and applies boundary conditions.
 *
 * NOTE: When adding new variables, remember to apply boundary conditions
 *       here.
 */

void var_read(planetspec   *planet,
              char         *infile,
              int           portion,
              unsigned int  time_index)
{
  int 
    K,J,I,
    it,itlo,ithi,
    nk,
    node,
    is,iq,
    i;
  int
    start[FOURDIM],
    end[FOURDIM],
    num_nodes;
  static char
    **gattname=NULL,
    **varname =NULL;
  static int
    ngatts    =0,
    num_progs =0;
  int
    nc_err,nc_id;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="var_read";

  /*
   * NOTE: Call lookup_netcdf() from all nodes so that
   *       its calls to MPI_Bcast() will work.
   */
  nc_err = lookup_netcdf(infile,&nc_id,&ngatts,&gattname,&num_progs,&varname);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),infile);
    epic_error(dbmsname,Message);
  }

  /*
   * Open NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    switch(portion) {
      case SIZE_DATA:
        fprintf(stdout,"Reading SIZE_DATA from %s \n",infile);
      break;
      case POST_SIZE_DATA:
        fprintf(stdout,"Reading POST_SIZE_DATA from %s \n",infile);
      break;
      case HEADER_DATA:
        fprintf(stdout,"Reading HEADER_DATA from %s \n",infile);
      break;
      case EXTRACT_HEADER_DATA:
        fprintf(stdout,"Reading EXTRACT_HEADER_DATA from %s \n",infile);
      break;
      case POST_HEADER_DATA:
        fprintf(stdout,"Reading POST_HEADER_DATA from %s \n",infile);
      break;
      case EXTRACT_DATA:
        fprintf(stdout,"Reading EXTRACT_DATA from %s \n",infile);
      break;
      case ALL_DATA:
        fprintf(stdout,"Reading ALL_DATA from %s \n",infile);
      break;
      default:
        sprintf(Message,"unrecognized portion = %d",portion);
        epic_error(dbmsname,Message);
      break;
    }
    fflush(stderr);
  }

  /*
   * Read in size of model and other data needed by make_arrays().
   */
  if (portion == SIZE_DATA           || 
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    READF(&grid.epic_version,grid_epic_version,1);

    READTIME(&var.start_time,var_start_time);

    READI(&planet->index,planet_index,1);
    READC(planet->name,planet_name,32);
    READC(planet->type,planet_type,16);
    READC(planet->orbital_epoch,planet_orbital_epoch,8);
    READF(&planet->re,planet_re,1);
    READF(&planet->rp,planet_rp,1);
    READF(&planet->obliquity,planet_obliquity,1);
    READF(&planet->omega_sidereal,planet_omega_sidereal,1);
    READF(&planet->omega_synodic,planet_omega_synodic,1);
    READF(&planet->cp,planet_cp,1);
    READF(&planet->rgas,planet_rgas,1);
    READF(&planet->p0,planet_p0,1);

    planet->cpr   = planet->cp/planet->rgas;
    planet->kappa = 1./planet->cpr;

    READF(&planet->GM,planet_GM,1);
    READF(&planet->J2,planet_J2,1);
    READF(&planet->x_he,planet_x_he,1);
    READF(&planet->x_h2,planet_x_h2,1);
    READF(&planet->x_3,planet_x_3,1);
    READF(&planet->a,planet_a,1);
    READF(&planet->e,planet_e,1);
    READF(&planet->i,planet_i,1);
    READF(&planet->lon_ascending_node,planet_lon_ascending_node,1);
    READF(&planet->lon_perihelion,planet_lon_perihelion,1);
    READF(&planet->mean_lon,planet_mean_lon,1);
    READF(&planet->orbit_period,planet_orbit_period,1);
    READF(&planet->vernal_equinox_anomaly,planet_vernal_equinox_anomaly,1);
    READF(&planet->kinvisc,planet_kinvisc,1);
    READF(&planet->dynvisc,planet_dynvisc,1);
    READF(&planet->k_a,planet_k_a,1);

    READI(&grid.nk,grid_nk,1);
    READI(&grid.nj,grid_nj,1);
    READI(&grid.ni,grid_ni,1);

    READI(&grid.jtp,grid_jtp,1);

    READC(grid.geometry,grid_geometry,GEOM_STR);
    READC(grid.vertical_coordinate,grid_vertical_coordinate,N_STR);
    READI(&grid.coord_type,grid_coord_type,1);
    READI(&grid.radiation_index,grid_radiation_index,1);

    READC(grid.uv_timestep_scheme,grid_uv_timestep_scheme,N_STR);
    READC(grid.radiation_scheme,grid_radiation_scheme,N_STR);
    READC(grid.turbulence_scheme,grid_turbulence_scheme,N_STR);

    READC(var.h.advection_scheme,var_h_advection_scheme,N_STR);
    READC(var.theta.advection_scheme,var_theta_advection_scheme,N_STR);
    READC(var.fpara.advection_scheme,var_fpara_advection_scheme,N_STR);
    READC(var.nu_turb.advection_scheme,var_nu_turb_advection_scheme,N_STR);
    READC(var.species[FIRST_SPECIES].advection_scheme,var_species_advection_scheme,N_STR);
    for (is = FIRST_SPECIES+1; is <= LAST_SPECIES; is++) {
      sprintf(var.species[is].advection_scheme,"%s",var.species[FIRST_SPECIES].advection_scheme);
    }

    READF(&grid.globe_lonbot,grid_globe_lonbot,1);
    READF(&grid.globe_lontop,grid_globe_lontop,1);
    READF(&grid.globe_latbot,grid_globe_latbot,1);
    READF(&grid.globe_lattop,grid_globe_lattop,1);
    READC(grid.f_plane_map,grid_f_plane_map,GEOM_STR);
    READF(&grid.f_plane_lat0,grid_f_plane_lat0,1);
    READF(&grid.f_plane_half_width,grid_f_plane_half_width,1);

    READI(grid.wrap,grid_wrap,TOPDIM);
    READI(grid.pad,grid_pad,TOPDIM);
    READI(&grid.jlo,grid_jlo,1);
    READI(&grid.jfirst,grid_jfirst,1);
    READI(&grid.jlast,grid_jlast,1);
    READI(&grid.ilo,grid_ilo,1);
    READI(&grid.k_sponge,grid_k_sponge,1);
    READI(&grid.j_sponge,grid_j_sponge,1);
    READI(&grid.extract_species_fraction_type,grid_extract_species_fraction_type,1);
    READI(&grid.k_sigma,grid_k_sigma,1);

    READF(&grid.du_vert,grid_du_vert,1);

    READI(&grid.cloud_microphysics,grid_cloud_microphysics,1);
    READI(&grid.include_nontrad_accel,grid_include_nontrad_accel,1);
    READI(&grid.zonal_average_rt,grid_zonal_average_rt,1);
    READI(var.on_list,var_on_list,LAST_INDEX-FIRST_INDEX+1);
    READI(&var.extract_on,var_extract_on,1);
    READI(var.extract_on_list,var_extract_on_list,LAST_INDEX-FIRST_INDEX+1);

    READI(&var.ntp,var_ntp,1);
  }

  if (portion == SIZE_DATA) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /********************************************** 
   *                                            *
   * Seam between SIZE_DATA and POST_SIZE_DATA. *
   *                                            * 
   **********************************************/

  if (portion == POST_SIZE_DATA      ||
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    EPIC_FLOAT
      *buffer;

    nk = grid.nk;

    READF(&grid.dlt,grid_dlt,1);
    READF(&grid.dln,grid_dln,1);
    READI(&grid.dt,grid_dt,1);

    READF(&var.fpara_rate_scaling,var_fpara_rate_scaling,1);

    if (var.ntp > 0) {
      /*
       * Read the reference temperature sounding profile.
       */
      READF(var.pdat,var_pdat,var.ntp);
      READF(var.tdat,var_tdat,var.ntp);
      READF(var.dtdat,var_dtdat,var.ntp);
      /*
       * Convert input pressure from hPa to Pa.
       */
      for (i = 0; i < var.ntp; i++) {
        var.pdat[i] *= 100.;
      }
    }

    READF(grid.re,grid_re,nk+1);
    READF(grid.rp,grid_rp,nk+1);

    READD(grid.sigmatheta,grid_sigmatheta,2*(nk+1)+2);

    READF(grid.p_ref,grid_p_ref,2*(nk+1)+2);
    READF(grid.t_ref,grid_t_ref,2*(nk+1)+2);
    READF(grid.rho_ref,grid_rho_ref,2*(nk+1)+2);
    READF(grid.theta_ref,grid_theta_ref,2*(nk+1)+2);
    READF(grid.h_min,grid_h_min,nk+1);

    READF(&grid.sgth_bot,grid_sgth_bot,1);
    READF(&grid.sgth_top,grid_sgth_top,1);
    /*
     * Double precision to improve calculation of diagnostic theta.
     */
    READD(&grid.zeta0,grid_zeta0,1);
    READD(&grid.zeta1,grid_zeta1,1);
    READD(&grid.hybrid_alpha,grid_hybrid_alpha,1);
    READD(&grid.sigma_sigma,grid_sigma_sigma,1);

    READI(&grid.newt_cool_adjust,grid_newt_cool_adjust,1);
    READC(grid.eos,grid_eos,8);

    READF(&grid.ptop,grid_ptop,1);
    READF(&grid.pbot,grid_pbot,1);
    READF(&grid.thetatop,grid_thetatop,1);
    READF(&grid.thetabot,grid_thetabot,1);
    READF(&grid.phi0,grid_phi0,1);

    READI(&grid.nu_order,grid_nu_order,1);
    READF(&grid.nudiv_nondim,grid_nudiv_nondim,1);
    READF(&grid.nu_nondim,grid_nu_nondim,1);
    READF(&grid.nu_hyper,grid_nu_hyper,1);
  }

  if (portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA   ) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /********************************************************************** 
   *                                                                    *
   * Seam between HEADER_DATA/EXTRACT_HEADER_DATA and POST_HEADER_DATA. *
   *                                                                    * 
   **********************************************************************/

  /*
   * Read non-array data that is time dependent.
   */
  READTIME(&var.model_time,var_model_time);

  READI(&grid.cfl_dt,grid_cfl_dt,1);

  READI(&grid.aux_a,grid_aux_a,1);
  READI(&grid.aux_b,grid_aux_b,1);
  READI(&grid.aux_c,grid_aux_c,1);
  READF(&grid.aux_fa,grid_aux_fa,1);
  READF(&grid.aux_fb,grid_aux_fb,1);
  READF(&grid.aux_fc,grid_aux_fc,1);

  if (IAMNODE == NODE0) {
    fprintf(stdout,"  0%%");
    fflush(stdout);
  }

  /*
   * Setup read_array().
   */
  num_nodes = setup_read_array();

  /*
   * Load start and end vectors.
   */

  /*
   * Read layers 0 to grid.nk+1 for vertical dimension (grid.nk+1 is used for abyssal-layer
   * values for gas giants).
   */
  start[2]  = 0;
  end[  2]  = grid.nk+1;
  start[3]  = time_index;
  end[  3]  = time_index;

  for (node = 0; node < num_nodes; node++) {
    if (IAMNODE == NODE0) {
      fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)node/num_nodes));
      fflush(stdout);
    }

    get_jlohi(node,num_nodes,start+1,end+1);
    get_ilohi(node,num_nodes,start+0,end+0);

    /*
     * Read parameter and key diagnostic arrays.
     */
    if (var.phi_surface.on) {
      /*
       * Read surface geopotential.
       */
      read_array(node,TWODIM,start,end,var.phi_surface.info[0].name,
                 var.phi_surface.info[0].index,var.phi_surface.value,EPIC_FLOAT_ARRAY,nc_id);
    }

    if (var.pbot.on) {
      /*
       * Read pressure bottom boundary condition.
       */
      read_array(node,TWODIM,start,end,var.pbot.info[0].name,
                 var.pbot.info[0].index,var.pbot.value,EPIC_FLOAT_ARRAY,nc_id);
    }

    if (var.dzdt2.on) {
      /*
       * Read standard vertical velocity [m/s], carried in the layer. 
       */
      read_array(node,FOURDIM,start,end,var.dzdt2.info[0].name,
                 var.dzdt2.info[0].index,var.dzdt2.value,EPIC_FLOAT_ARRAY,nc_id);
    }

    if (portion == EXTRACT_DATA) {
      sprintf(Message,"portion == EXTRACT_DATA not yet implemented");
      epic_error(dbmsname,Message);
    }

    /* 
     * Read variables, and their associated tendency data as appropriate.
     */
    if (var.u.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.u.info[0].name,
                     var.u.info[0].index,var.u.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
        if (portion != VAR_DATA) {
          read_array(node,THREEDIM,start,end,var.u.info_tend[0].name,
                     var.u.info_tend[0].index,var.u.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,THREEDIM,start,end,var.u.info_tend[1].name,
                     var.u.info_tend[1].index,var.u.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.u.info[0].name,
                     var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,FOURDIM,start,end,var.u.info[1].name,
                     var.u.info[1].index,var.u.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.v.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.v.info[0].name,
                     var.v.info[0].index,var.v.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
        if (portion != VAR_DATA) {
          read_array(node,THREEDIM,start,end,var.v.info_tend[0].name,
                     var.v.info_tend[0].index,var.v.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,THREEDIM,start,end,var.v.info_tend[1].name,
                     var.v.info_tend[1].index,var.v.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion != EXTRACT_DATA) {
          read_array(node,FOURDIM,start,end,var.v.info[0].name,
                     var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          read_array(node,FOURDIM,start,end,var.v.info[1].name,
                     var.v.info[1].index,var.v.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"unrecognized portion=%d",portion);
          epic_error(dbmsname,Message);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.p3.on) {
      if ((portion == EXTRACT_DATA && var.p3.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.p3.info[0].name,
                   var.p3.info[0].index,var.p3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.theta.on) {
      if ((portion == EXTRACT_DATA && var.theta.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.theta.info[0].name,
                   var.theta.info[0].index,var.theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.fpara.on) {
      if ((portion == EXTRACT_DATA && var.fpara.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.fpara.info[0].name,
                   var.fpara.info[0].index,var.fpara.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    for (iq = 0; iq < grid.nq; iq++) {
      if ((portion == EXTRACT_DATA && var.species[grid.is[iq]].phase[grid.ip[iq]].extract_on) ||
          (portion != EXTRACT_DATA                                        )   ) {
        read_array(node,FOURDIM,start,end,var.species[grid.is[iq]].phase[grid.ip[iq]].info[MASS].name,
                   var.species[grid.is[iq]].phase[grid.ip[iq]].info[MASS].index,
                   var.species[grid.is[iq]].phase[grid.ip[iq]].q,
                   EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.nu_turb.on) {
      if ((portion == EXTRACT_DATA && var.nu_turb.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        read_array(node,FOURDIM,start,end,var.nu_turb.info[0].name,
                   var.nu_turb.info[0].index,var.nu_turb.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }
  }

  /*
   * Close NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    nc_close(nc_id);
    fprintf(stdout,"\b\b\b\b%3d%%\n",100);
    fflush(stdout);
  }

  /*
   * Set solar longitude, L_s [deg], which is a function of time.
   */
  L_s = solar_longitude(planet,var.model_time);

  /*
   * Apply lateral boundary conditions.
   *
   * NOTE: bc_lateral() must be called from all nodes.
   */

  if (var.dzdt2.on) {
    bc_lateral(var.dzdt2.value,THREEDIM);
  }

  if (var.u.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      bc_lateral(var.u.value+0*Nelem3d,THREEDIM);
      if (portion != VAR_DATA) {
        bc_lateral(var.u.tendency+IT_MINUS2*Nelem3d,THREEDIM);
        bc_lateral(var.u.tendency+IT_MINUS1*Nelem3d,THREEDIM);
      }
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      bc_lateral(var.u.value+IT_MINUS1*Nelem3d,THREEDIM);
      bc_lateral(var.u.value+IT_ZERO*Nelem3d,  THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }
  }

  if (var.v.on) {
    if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
      bc_lateral(var.v.value+0*Nelem3d,THREEDIM);
      if (portion != VAR_DATA) {
        bc_lateral(var.v.tendency+IT_MINUS2*Nelem3d,THREEDIM);
        bc_lateral(var.v.tendency+IT_MINUS1*Nelem3d,THREEDIM);
      }
    }
    else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
      bc_lateral(var.v.value+IT_MINUS1*Nelem3d,THREEDIM);
      bc_lateral(var.v.value+IT_ZERO*Nelem3d,  THREEDIM);
    }
    else {
      sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
      epic_error(dbmsname,Message);
    }
  }

  if (var.p3.on) {
    bc_lateral(var.p3.value,THREEDIM);
  }

  if (var.theta.on) {
    bc_lateral(var.theta.value,THREEDIM);
  }

  if (var.fpara.on) {
    bc_lateral(var.fpara.value,THREEDIM);
  }

  for (iq = 0; iq < grid.nq; iq++) {
    bc_lateral(var.species[grid.is[iq]].phase[grid.ip[iq]].q,THREEDIM);
  }
  /*
   * Synchronize mole fractions, X, to mass mixing ratios, Q.
   */
  sync_x_to_q(planet);

  if (var.nu_turb.on) {
    bc_lateral(var.nu_turb.value,THREEDIM);
  }

  if (var.phi_surface.on) {
    bc_lateral(var.phi_surface.value,TWODIM);
  }

  if (var.pbot.on) {
    bc_lateral(var.pbot.value,TWODIM);
  }

  /*
   * Make sure v = 0 at poles and channel walls.
   */
  if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
    itlo = 0;
    ithi = 0;
  }
  else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
    itlo = 0;
    ithi = 1;
  }
  else {
    sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
    epic_error(dbmsname,Message);
  }
  if (!grid.wrap[1]) {
    J = JLO;
    if (J == grid.jlo) {
      for (it = itlo; it <= ithi; it++) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            if (V(it,K,J,I) != 0.) {
              V(it,K,J,I) = 0.;
            }
          }
        }
      }
    }
    J = JHI+1;
    if (J == grid.nj+1) {
      for (it = itlo; it <= ithi; it++) {
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (I = ILOPAD; I <= IHIPAD; I++) {
            if (V(it,K,J,I) != 0.) {
              V(it,K,J,I) = 0.;
            }
          }
        }
      }
    }
  }

  return;
}

/*======================= end of var_read() =================================*/

/*======================= var_write() =======================================*/

/*
 * Write the variables u, v, p, theta, etc.
 *
 * If ni == 1 and stretch_ni > 1, stretch the model from 2D to 3D on output,
 * with the new ni being stretch_ni.
 */

void var_write(planetspec   *planet,
               char         *outfile,
               int           portion,
               unsigned int  time_index,
               int           stretch_ni)
{
  int 
    K,J,I,kk,
    nk,
    node,
    num_nodes,
    iq,
    i;
  int
    start[FOURDIM],
    end[FOURDIM],
    nc_err,
    nc_id;
  static int
    initialized = FALSE;
  size_t
    t_index[1];
  EPIC_FLOAT
    the_time[1],
    the_L_s[1],
    avg;
  register EPIC_FLOAT
    tmp;
  static EPIC_FLOAT
    *buff3d;
  nc_type
    the_nc_type;     /* NOTE: Used in i/o macros. */
  char 
    history[N_STR];
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
    dbmsname[]="var_write";

  if (!initialized) {
    /* Allocate memory. */
    buff3d = fvector(0,Nelem3d-1,dbmsname);

    initialized = TRUE;
  }

  /*
   * Open NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    if (portion == SIZE_DATA           ||
        portion == HEADER_DATA         ||
        portion == EXTRACT_HEADER_DATA ||
        portion == ALL_DATA              ) {
      /*
       * Create netCDF file and define variables and attributes.
       */
      handle_file_compression(outfile);
      nc_err = nc_create(outfile,NC_CLOBBER,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),outfile);
        epic_error(dbmsname,Message);
      }

      if (stretch_ni > 0) {
        /*
         * Modify grid.dln to be consistent with stretch_ni,
         * before call to define_netcdf().
         */
        if (strcmp(grid.geometry,"globe") == 0) {
          grid.dln = (grid.globe_lontop-grid.globe_lonbot)/stretch_ni;
        }
        else if (strcmp(grid.geometry,"f-plane") == 0) {
          if (strcmp(grid.f_plane_map,"cartesian") == 0) {
            sprintf(Message,"-stretch_ni for case %s %s not yet implemented",grid.geometry,grid.f_plane_map);
            epic_error(dbmsname,Message);
          }
          else if (strcmp(grid.f_plane_map,"polar") == 0) {
            grid.dln = 360./(stretch_ni);
          }
        }
      }

      define_netcdf(planet,outfile,portion,stretch_ni,nc_id);
    }
    else {
      /*
       * Open existing netCDF file.
       */
      handle_file_compression(outfile);
      nc_err = nc_open(outfile,NC_WRITE,&nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"%s: \"%s\"",nc_strerror(nc_err),outfile);
        epic_error(dbmsname,Message);
      }

      /*
       * Put into define mode.
       */
      nc_err = nc_redef(nc_id);
      if (nc_err != NC_NOERR) {
        sprintf(Message,"nc_redef(), %s",nc_strerror(nc_err));
        epic_error(dbmsname,Message);
      }
    }

    switch (portion) {
      case SIZE_DATA:
        fprintf(stdout,"Writing SIZE_DATA to %s \n",outfile);
      break;
      case POST_SIZE_DATA:
        fprintf(stdout,"Writing POST_SIZE_DATA to %s \n",outfile);
      break;
      case HEADER_DATA:
        fprintf(stdout,"Writing HEADER_DATA to %s \n",outfile);
      break;
      case EXTRACT_HEADER_DATA:
        fprintf(stdout,"Writing EXTRACT_HEADER_DATA to %s \n",outfile);
      break;
      case POST_HEADER_DATA:
        fprintf(stdout,"Writing POST_HEADER_DATA to %s \n",outfile);
      break;
      case VAR_DATA:
        fprintf(stdout,"Writing VAR_DATA to %s for timestep %lu\n",outfile,grid.itime);
      break;
      case EXTRACT_DATA:
        fprintf(stdout,"Writing EXTRACT_DATA to %s for timestep %lu \n",outfile,grid.itime);
      break;
      case ALL_DATA:
        fprintf(stdout,"Writing ALL_DATA to %s for timestep %lu\n",outfile,grid.itime);
      break;
      default:
        sprintf(Message,"unrecognized portion = %d",portion);
        epic_error(dbmsname,Message);
      break;
    }
    fflush(stderr);
  }

  if (portion == SIZE_DATA           || 
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {

    /* 
     * Write global attributes.
     * We are following the Climate and Forecast (CF) Metadata Convention.
     */
    sprintf(Message,"CF-1.0");
    nc_put_att_text(nc_id,NC_GLOBAL,"Conventions",strlen(Message)+1,Message);
    if (portion == EXTRACT_HEADER_DATA) {
      sprintf(Message,"EPIC Model Version %4.2f: extract data file",grid.epic_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"title",strlen(Message)+1,Message);
      sprintf(Message,"EPIC Model Version %4.2f: specified fields are extracted from the running model at timestep interval -itextract",
                       grid.epic_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"source",strlen(Message)+1,Message);
    }
    else {
      sprintf(Message,"EPIC Model Version %4.2f: runable data file",grid.epic_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"title",strlen(Message)+1,Message);
      sprintf(Message,"EPIC Model Version %4.2f: model state (prognostic variables and tendencies)",
                       grid.epic_version);
      nc_put_att_text(nc_id,NC_GLOBAL,"source",strlen(Message)+1,Message);
    }
    sprintf(Message,"Developed at the University of Louisville, KY USA");
    nc_put_att_text(nc_id,NC_GLOBAL,"institution",strlen(Message)+1,Message);
   /* 
    * NOTE: Regarding time stamping the  global attribute "history," gcc 3.3.5 has a bug that generates a segmentation fault
    *       with time() and ctime(), so we are currently not including a timestamp.
    */   
    sprintf(Message,"This file produced by the EPIC Atmospheric model");
    nc_put_att_text(nc_id,NC_GLOBAL,"history",strlen(Message)+1,Message);
    sprintf(Message,"The EPIC Model is downloadable as open source from the NASA Planetary Data System (PDS) Atmospheres Node\n");
    nc_put_att_text(nc_id,NC_GLOBAL,"references",strlen(Message)+1,Message);
    
    WRITEF(&grid.epic_version,grid_epic_version,1);

    WRITETIME(&var.start_time,var_start_time);

    WRITEI(&planet->index,planet_index,1);
    WRITEC(planet->name,planet_name,32);
    WRITEC(planet->type,planet_type,16);
    WRITEC(planet->orbital_epoch,planet_orbital_epoch,8);
    WRITEF(&planet->re,planet_re,1);
    WRITEF(&planet->rp,planet_rp,1);
    WRITEF(&planet->obliquity,planet_obliquity,1);
    WRITEF(&planet->omega_sidereal,planet_omega_sidereal,1);
    WRITEF(&planet->omega_synodic,planet_omega_synodic,1);
    WRITEF(&planet->cp,planet_cp,1);
    WRITEF(&planet->rgas,planet_rgas,1);
    WRITEF(&planet->p0,planet_p0,1);
    WRITEF(&planet->GM,planet_GM,1);
    WRITEF(&planet->J2,planet_J2,1);
    WRITEF(&planet->x_he,planet_x_he,1);
    WRITEF(&planet->x_h2,planet_x_h2,1);
    WRITEF(&planet->x_3,planet_x_3,1);
    WRITEF(&planet->a,planet_a,1);
    WRITEF(&planet->e,planet_e,1);
    WRITEF(&planet->i,planet_i,1);
    WRITEF(&planet->lon_ascending_node,planet_lon_ascending_node,1);
    WRITEF(&planet->lon_perihelion,planet_lon_perihelion,1);
    WRITEF(&planet->mean_lon,planet_mean_lon,1);
    WRITEF(&planet->orbit_period,planet_orbit_period,1);
    WRITEF(&planet->vernal_equinox_anomaly,planet_vernal_equinox_anomaly,1);
    WRITEF(&planet->kinvisc,planet_kinvisc,1);
    WRITEF(&planet->dynvisc,planet_dynvisc,1);
    WRITEF(&planet->k_a,planet_k_a,1);

    WRITEI(&grid.nk,grid_nk,1);
    WRITEI(&grid.nj,grid_nj,1);

    if (stretch_ni > 0) {
      WRITEI(&stretch_ni,grid_ni,1);
    }
    else {
      WRITEI(&grid.ni,grid_ni,1);
    }

    WRITEI(&grid.jtp,grid_jtp,1);

    WRITEC(grid.geometry,grid_geometry,GEOM_STR);
    WRITEC(grid.vertical_coordinate,grid_vertical_coordinate,N_STR);
    WRITEI(&grid.coord_type,grid_coord_type,1);
    WRITEI(&grid.radiation_index,grid_radiation_index,1);

    WRITEC(grid.uv_timestep_scheme,grid_uv_timestep_scheme,N_STR);
    WRITEC(grid.radiation_scheme,grid_radiation_scheme,N_STR);
    WRITEC(grid.turbulence_scheme,grid_turbulence_scheme,N_STR);

    WRITEC(var.h.advection_scheme,var_h_advection_scheme,N_STR);
    WRITEC(var.theta.advection_scheme,var_theta_advection_scheme,N_STR);
    WRITEC(var.fpara.advection_scheme,var_fpara_advection_scheme,N_STR);
    WRITEC(var.nu_turb.advection_scheme,var_nu_turb_advection_scheme,N_STR);
    WRITEC(var.species[FIRST_SPECIES].advection_scheme,var_species_advection_scheme,N_STR);

    WRITEF(&grid.globe_lonbot,grid_globe_lonbot,1);
    WRITEF(&grid.globe_lontop,grid_globe_lontop,1);
    WRITEF(&grid.globe_latbot,grid_globe_latbot,1);
    WRITEF(&grid.globe_lattop,grid_globe_lattop,1);
    WRITEC(grid.f_plane_map,grid_f_plane_map,GEOM_STR);       
    WRITEF(&grid.f_plane_lat0,grid_f_plane_lat0,1);
    WRITEF(&grid.f_plane_half_width,grid_f_plane_half_width,1);

    WRITEI(grid.wrap,grid_wrap,TOPDIM);
    WRITEI(grid.pad,grid_pad,TOPDIM);
    WRITEI(&grid.jlo,grid_jlo,1);
    WRITEI(&grid.jfirst,grid_jfirst,1);
    WRITEI(&grid.jlast,grid_jlast,1);
    WRITEI(&grid.ilo,grid_ilo,1);
    WRITEI(&grid.k_sponge,grid_k_sponge,1);
    WRITEI(&grid.j_sponge,grid_j_sponge,1);
    WRITEI(&grid.extract_species_fraction_type,grid_extract_species_fraction_type,1);
    WRITEI(&grid.k_sigma,grid_k_sigma,1);
    WRITEF(&grid.du_vert,grid_du_vert,1);

    WRITEI(&grid.cloud_microphysics,grid_cloud_microphysics,1);
    WRITEI(&grid.include_nontrad_accel,grid_include_nontrad_accel,1);
    WRITEI(&grid.zonal_average_rt,grid_zonal_average_rt,1);
    WRITEI(var.on_list,var_on_list,LAST_INDEX-FIRST_INDEX+1);
    WRITEI(&var.extract_on,var_extract_on,1);
    WRITEI(var.extract_on_list,var_extract_on_list,LAST_INDEX-FIRST_INDEX+1);

    WRITEI(&var.ntp,var_ntp,1);
  }

  if (portion == SIZE_DATA) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /**********************************************
   *                                            *
   * Seam between SIZE_DATA and POST_SIZE_DATA. *
   *                                            *
   **********************************************/

  if (portion == POST_SIZE_DATA      ||
      portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA ||
      portion == ALL_DATA              ) {
    EPIC_FLOAT
      *buffer;

    nk = grid.nk;

    WRITEF(&grid.dlt,grid_dlt,1);
    WRITEF(&grid.dln,grid_dln,1);
    WRITEI(&grid.dt,grid_dt,1);

    WRITEF(&var.fpara_rate_scaling,var_fpara_rate_scaling,1);

    if (var.ntp > 0) {
      /*
       * Write reference temperature sounding profile.
       */
      buffer = fvector(0,var.ntp-1,dbmsname);
      for (i = 0; i < var.ntp; i++) {
        /*
         * Convert output pressure from Pa to hPa.
         */
        buffer[i] = var.pdat[i]/100.;
      }
      WRITEF(buffer,var_pdat,var.ntp);
      WRITEC("hPa",var_pdat_units,strlen("hPa")+1);
      WRITEF(var.tdat,var_tdat,var.ntp);
      WRITEC("K",var_tdat_units,strlen("K")+1);
      WRITEF(var.dtdat,var_dtdat,var.ntp);
      WRITEC("K",var_dtdat_units,strlen("K")+1);
      free_fvector(buffer,0,var.ntp-1,dbmsname);
    }

    WRITEF(grid.re,grid_re,nk+1);
    WRITEF(grid.rp,grid_rp,nk+1);

    WRITED(grid.sigmatheta,grid_sigmatheta,2*(nk+1)+2);

    WRITEF(grid.p_ref,grid_p_ref,2*(nk+1)+2);
    WRITEF(grid.t_ref,grid_t_ref,2*(nk+1)+2);
    WRITEF(grid.rho_ref,grid_rho_ref,2*(nk+1)+2);
    WRITEF(grid.theta_ref,grid_theta_ref,2*(nk+1)+2);
    WRITEF(grid.h_min,grid_h_min,nk+1);

    WRITEF(&grid.sgth_bot,grid_sgth_bot,1);
    WRITEF(&grid.sgth_top,grid_sgth_top,1);
    /*
     * Double precision to improve calculation of diagnostic theta.
     */
    WRITED(&grid.zeta0,grid_zeta0,1);
    WRITED(&grid.zeta1,grid_zeta1,1);
    WRITED(&grid.hybrid_alpha,grid_hybrid_alpha,1);
    WRITED(&grid.sigma_sigma,grid_sigma_sigma,1);

    WRITEI(&grid.newt_cool_adjust,grid_newt_cool_adjust,1);
    WRITEC(grid.eos,grid_eos,8);

    WRITEF(&grid.ptop,grid_ptop,1);
    WRITEF(&grid.pbot,grid_pbot,1);
    WRITEF(&grid.thetatop,grid_thetatop,1);
    WRITEF(&grid.thetabot,grid_thetabot,1);
    WRITEF(&grid.phi0,grid_phi0,1);

    WRITEI(&grid.nu_order,grid_nu_order,1);
    WRITEF(&grid.nudiv_nondim,grid_nudiv_nondim,1);
    WRITEF(&grid.nu_nondim,grid_nu_nondim,1);
    WRITEF(&grid.nu_hyper,grid_nu_hyper,1);
  }

  if (portion == HEADER_DATA         ||
      portion == EXTRACT_HEADER_DATA   ) {
    if (IAMNODE == NODE0) {
      nc_close(nc_id);
    }
    return;
  }

  /**********************************************************************
   *                                                                    *
   * Seam between HEADER_DATA/EXTRACT_HEADER_DATA and POST_HEADER_DATA. *
   *                                                                    *
   **********************************************************************/

  /*
   * Write non-array data that is time dependent.
   */
  WRITETIME(&var.model_time,var_model_time);

  grid.cfl_dt = cfl_dt(planet);
  WRITEI(&grid.cfl_dt,grid_cfl_dt,1);

  WRITEI(&grid.aux_a,grid_aux_a,1);
  WRITEI(&grid.aux_b,grid_aux_b,1);
  WRITEI(&grid.aux_c,grid_aux_c,1);
  WRITEF(&grid.aux_fa,grid_aux_fa,1);
  WRITEF(&grid.aux_fb,grid_aux_fb,1);
  WRITEF(&grid.aux_fc,grid_aux_fc,1);

  if (IAMNODE == NODE0) {
    /*
     * Leave define mode for netCDF file:
     */
    nc_err = nc_enddef(nc_id);
    if (nc_err != NC_NOERR) {
      sprintf(Message,"nc_enddef(), %s",nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Write time in days since var.start_time.
   */
  if (IAMNODE == NODE0) {
    t_index[ 0] = time_index;
    the_time[0] = TIME/86400.;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_err = nc_put_var1_double(nc_id,var.info[0].coorid[NETCDF_T_INDEX],t_index,the_time);
#else
    nc_err = nc_put_var1_float(nc_id,var.info[0].coorid[NETCDF_T_INDEX],t_index,the_time);
#endif

    if (nc_err != NC_NOERR) {
      sprintf(Message,"t_index=%lu, %s",(long unsigned)t_index[0],nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /*
   * Write solar longitude, L_s [deg].
   */
  if (IAMNODE == NODE0) {
    t_index[0] = time_index;
    the_L_s[0] = L_s;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_err = nc_put_var1_double(nc_id,var.l_s.info.id,t_index,the_L_s);
#else
    nc_err = nc_put_var1_float(nc_id,var.l_s.info.id,t_index,the_L_s);
#endif

    if (nc_err != NC_NOERR) {
      sprintf(Message,"Writing L_s, t_index=%lu, %s",(long unsigned)t_index[0],nc_strerror(nc_err));
      epic_error(dbmsname,Message);
    }
  }

  /* 
   * Loop over nodes to get all the information off one node before
   * moving on to the next. 
   */
  num_nodes = setup_write_array();

  /*
   * Load start and end vectors.
   */

  /*
   * Write layers 0 to grid.nk+1 for vertical dimension (grid.nk+1 is used for abyssal-layer
   * values for gas giants).
   */
  start[2] = 0;
  end[  2] = grid.nk+1;

  start[3] = time_index;
  end[  3] = time_index;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"  0%%");
    fflush(stdout);
  }

  for (node = 0; node < num_nodes; node++) {
    if (IAMNODE == NODE0) {
      fprintf(stdout,"\b\b\b\b%3d%%",(int)(100.*(EPIC_FLOAT)node/grid.we_num_nodes));
      fflush(stdout);
    }

    get_jlohi(node,num_nodes,start+1,end+1);
    get_ilohi(node,num_nodes,start+0,end+0);

    /*
     * Write parameter and key diagnostic arrays.
     */
    if (var.phi_surface.on) {
      if ((portion == EXTRACT_DATA && var.phi_surface.extract_on) ||
          (portion != EXTRACT_DATA)                                ) {
        write_array(node,TWODIM,start,end,stretch_ni,var.phi_surface.info[0].name,
                    var.phi_surface.info[0].index,var.phi_surface.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.pbot.on) {
      if ((portion == EXTRACT_DATA && var.pbot.extract_on) ||
          (portion != EXTRACT_DATA)                                ) {
        write_array(node,TWODIM,start,end,stretch_ni,var.pbot.info[0].name,
                    var.pbot.info[0].index,var.pbot.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.dzdt2.on) {
      if ((portion == EXTRACT_DATA && var.dzdt2.extract_on) ||
          (portion != EXTRACT_DATA)                           ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.dzdt2.info[0].name,
                    var.dzdt2.info[0].index,var.dzdt2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.gravity2.on) {
      /*
       * Suppress I dimension by passing -1 for stretch_ni.
       */
      write_array(node,THREEDIM,start,end,-1,var.gravity2.info[0].name,
                  var.gravity2.info[0].index,var.gravity2.value,EPIC_FLOAT_ARRAY,nc_id);
    }

    if (portion == EXTRACT_DATA) {
      /*
       * Write selected diagnostic variables to extract file.
       */

      /*
       * P3 is written as if it is the prognostic variable, rather than H,
       * so we handle H here instead of P3.
       */
      if (var.h.on && var.h.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.h.info[0].name,
                    var.h.info[0].index,var.h.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.h3.on && var.h3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.h3.info[0].name,
                    var.h3.info[0].index,var.h3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.hdry2.on && var.hdry2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.hdry2.info[0].name,
                    var.hdry2.info[0].index,var.hdry2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.hdry3.on && var.hdry3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.hdry3.info[0].name,
                    var.hdry3.info[0].index,var.hdry3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.pdry3.on && var.pdry3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.pdry3.info[0].name,
                    var.pdry3.info[0].index,var.pdry3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.p2.on && var.p2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.p2.info[0].name,
                    var.p2.info[0].index,var.p2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.theta2.on && var.theta2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.theta2.info[0].name,
                    var.theta2.info[0].index,var.theta2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.t2.on && var.t2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.t2.info[0].name,
                    var.t2.info[0].index,var.t2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.t3.on && var.t3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.t3.info[0].name,
                    var.t3.info[0].index,var.t3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.rho2.on && var.rho2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.rho2.info[0].name,
                    var.rho2.info[0].index,var.rho2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.rho3.on && var.rho3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.rho3.info[0].name,
                    var.rho3.info[0].index,var.rho3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.exner2.on && var.exner2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.exner2.info[0].name,
                    var.exner2.info[0].index,var.exner2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.exner3.on && var.exner3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.exner3.info[0].name,
                    var.exner3.info[0].index,var.exner3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.fgibb2.on && var.fgibb2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.fgibb2.info[0].name,
                    var.fgibb2.info[0].index,var.fgibb2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.z2.on && var.z2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.z2.info[0].name,
                    var.z2.info[0].index,var.z2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.phi2.on && var.phi2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.phi2.info[0].name,
                    var.phi2.info[0].index,var.phi2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.phi3.on && var.phi3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.phi3.info[0].name,
                    var.phi3.info[0].index,var.phi3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.mont2.on && var.mont2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.mont2.info[0].name,
                    var.mont2.info[0].index,var.mont2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.heat3.on && var.heat3.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.heat3.info[0].name,
                    var.heat3.info[0].index,var.heat3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.pv2.on && var.pv2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.pv2.info[0].name,
                    var.pv2.info[0].index,var.pv2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.ri2.on && var.ri2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.ri2.info[0].name,
                    var.ri2.info[0].index,var.ri2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_uv.on && var.diffusion_coef_uv.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_uv.info[0].name,
                    var.diffusion_coef_uv.info[0].index,var.diffusion_coef_uv.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_theta.on && var.diffusion_coef_theta.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_theta.info[0].name,
                    var.diffusion_coef_theta.info[0].index,var.diffusion_coef_theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.diffusion_coef_mass.on && var.diffusion_coef_mass.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.diffusion_coef_mass.info[0].name,
                    var.diffusion_coef_mass.info[0].index,var.diffusion_coef_mass.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.div_uv2.on && var.div_uv2.extract_on) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.div_uv2.info[0].name,
                    var.div_uv2.info[0].index,var.div_uv2.value,EPIC_FLOAT_ARRAY,nc_id);
      }
      if (var.w3.on && var.w3.extract_on) {
        switch(grid.coord_type) {
          case COORD_ISOBARIC:
            /*
             * Output pdot = D(p)/Dt [Pa/s], instead of the internal D(zeta)/Dt = W3 [K/s].
             */
            tmp = log(grid.ptop/grid.pbot)/(grid.zeta1-grid.zeta0);
            for (K = 0; K <= KHI; K++) {
              for (J = JLO; J <= JHI; J++) {
                for (I = ILO; I <= IHI; I++) {
                  /*
                   * NOTE: This formula for dp/dt is only valid in isobaric-coordinate case.
                   */
                  BUFF3D(K,J,I) = P3(K,J,I)*W3(K,J,I)*tmp;
                }
              }
              /* No need to apply bc_lateral() here. */
            }
            write_array(node,FOURDIM,start,end,stretch_ni,var.w3.info[0].name,
                        var.w3.info[0].index,buff3d,EPIC_FLOAT_ARRAY,nc_id);
          break;
          default:
            write_array(node,FOURDIM,start,end,stretch_ni,var.w3.info[0].name,
                        var.w3.info[0].index,var.w3.value,EPIC_FLOAT_ARRAY,nc_id);
          break;
        }
      }

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                         
       * Collect together here output options that do not have permanent memory.   *
       * These variables must be calculated here.                                  *
       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

      if (var.eddy_pv2.extract_on) {
        var.eddy_pv2.value = buff3d;
        for (K = 0; K <= KHI; K++) {
          for (J = JLOPAD; J <= JHIPADPV; J++) {
            avg = 0.;
            for (I = ILO; I <= IHI; I++) {
              avg += PV2(K,J,I);
            }

#if defined(EPIC_MPI)
            mpi_tmp = avg;
            MPI_Allreduce(&mpi_tmp,&avg,1,float_type,MPI_SUM,para.comm_JLO);
#endif

            avg /= grid.ni;

            for (I = ILOPAD; I <= IHIPAD; I++) {
              EDDY_PV2(K,J,I) = PV2(K,J,I)-avg;
            }
          }
        }
        /* No need to apply bc_lateral() here. */

        write_array(node,FOURDIM,start,end,stretch_ni,var.eddy_pv2.info[0].name,
                    var.eddy_pv2.info[0].index,var.eddy_pv2.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (var.rel_vort2.extract_on) {
        if (!var.rel_vort2.on) {
          /* Calculate variable here and store in BUFF3D memory. */
          var.rel_vort2.value = buff3d;
          for (K = 0; K <= KHI; K++) {
            vorticity(ON_SIGMATHETA,RELATIVE,2*K,
                      var.u.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      var.v.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      NULL,
                      var.rel_vort2.value+(K-Kshift)*Nelem2d);
          }
        }
        write_array(node,FOURDIM,start,end,stretch_ni,var.rel_vort2.info[0].name,
                    var.rel_vort2.info[0].index,var.rel_vort2.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (var.eddy_rel_vort2.extract_on) {
        if (!var.eddy_rel_vort2.on) {
          /* Calculate variable here and store in BUFF3D memory. */
          var.eddy_rel_vort2.value = buff3d;
          for (K = 0; K <= KHI; K++) {
            /* Temporarily store REL_VORT2(K,J,I) in EDDY_REL_VORT2(K,J,I). */
            vorticity(ON_SIGMATHETA,RELATIVE,2*K,
                      var.u.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      var.v.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      NULL,
                      var.eddy_rel_vort2.value+(K-Kshift)*Nelem2d);

            for (J = JLOPAD; J <= JHIPADPV; J++) {
              avg = 0.;
              for (I = ILO; I <= IHI; I++) {
                avg += EDDY_REL_VORT2(K,J,I);
              }

#if defined(EPIC_MPI)
              mpi_tmp = avg;
              MPI_Allreduce(&mpi_tmp,&avg,1,float_type,MPI_SUM,para.comm_JLO);
#endif

              avg /= grid.ni;

              for (I = ILOPAD; I <= IHIPAD; I++) {
                EDDY_REL_VORT2(K,J,I) -= avg;
              }
            }
          }
        }
        write_array(node,FOURDIM,start,end,stretch_ni,var.eddy_rel_vort2.info[0].name,
                    var.eddy_rel_vort2.info[0].index,var.eddy_rel_vort2.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (var.abs_vort2.extract_on) {
        if (!var.abs_vort2.on) {
          /* Calculate variable here and store in BUFF3D memory. */
          var.abs_vort2.value = buff3d;
          for (K = 0; K <= KHI; K++) {
            vorticity(ON_SIGMATHETA,ABSOLUTE,2*K,
                      var.u.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      var.v.value    +(K-Kshift)*Nelem2d+(grid.it_uv)*Nelem3d,
                      NULL,
                      var.abs_vort2.value+(K-Kshift)*Nelem2d);
          }
        }
        write_array(node,FOURDIM,start,end,stretch_ni,var.abs_vort2.info[0].name,
                    var.abs_vort2.info[0].index,var.abs_vort2.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (var.kinetic_energy2.extract_on) {
        if (!var.kinetic_energy2.on) {
          /* Calculate variable here and store in BUFF3D memory. */
          var.kinetic_energy2.value = buff3d;
          for (K = 0; K <= KHI; K++) {
            for (J = JLO; J <= JHI; J++) {
              for (I = ILO; I <= IHI; I++) {
                BUFF3D(K,J,I) = get_kin(planet,var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
                                               var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,2*K,J,I);
              }
            }
          }
          /* No need to apply bc_lateral() here */
        }
        write_array(node,FOURDIM,start,end,stretch_ni,var.kinetic_energy2.info[0].name,
                    var.kinetic_energy2.info[0].index,var.kinetic_energy2.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      if (var.molar_mass3.extract_on) {
        var.molar_mass3.value = buff3d;
        for (K = 0; K <= KHI; K++) {
          kk = 2*K+1;
          for (J = JLOPAD; J <= JHIPADPV; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              MOLAR_MASS3(K,J,I) = avg_molar_mass(planet,kk,J,I);
            }
          }
        }
        /* No need to apply bc_lateral() here. */
        write_array(node,FOURDIM,start,end,stretch_ni,var.molar_mass3.info[0].name,
                    var.molar_mass3.info[0].index,var.molar_mass3.value,EPIC_FLOAT_ARRAY,nc_id);
      }

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       * End of collection of output options that do not have permanent memory.      *
       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    }

    /* 
     * Write prognostic variables, and their associated tendency data as appropriate.
     */
    if (var.u.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if ((portion == EXTRACT_DATA && var.u.extract_on) ||
            (portion != EXTRACT_DATA                    )   ) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        if ((portion != VAR_DATA    ) &&
            (portion != EXTRACT_DATA)   ) {
          write_array(node,THREEDIM,start,end,stretch_ni,var.u.info_tend[0].name,
                      var.u.info_tend[0].index,var.u.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,THREEDIM,start,end,stretch_ni,var.u.info_tend[1].name,
                      var.u.info_tend[1].index,var.u.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion == EXTRACT_DATA && var.u.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else if (portion != EXTRACT_DATA) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[0].name,
                      var.u.info[0].index,var.u.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,FOURDIM,start,end,stretch_ni,var.u.info[1].name,
                      var.u.info[1].index,var.u.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    if (var.v.on) {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        if ((portion == EXTRACT_DATA && var.v.extract_on) ||
            (portion != EXTRACT_DATA                    )   ) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value,EPIC_FLOAT_ARRAY,nc_id);
        }
        if ((portion != VAR_DATA    ) &&
            (portion != EXTRACT_DATA)   ) {
          write_array(node,THREEDIM,start,end,stretch_ni,var.v.info_tend[0].name,
                      var.v.info_tend[0].index,var.v.tendency+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,THREEDIM,start,end,stretch_ni,var.v.info_tend[1].name,
                      var.v.info_tend[1].index,var.v.tendency+IT_MINUS2*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        if (portion == EXTRACT_DATA && var.v.extract_on) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
        else if (portion != EXTRACT_DATA) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[0].name,
                      var.v.info[0].index,var.v.value+IT_ZERO*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
          write_array(node,FOURDIM,start,end,stretch_ni,var.v.info[1].name,
                      var.v.info[1].index,var.v.value+IT_MINUS1*Nelem3d,EPIC_FLOAT_ARRAY,nc_id);
        }
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }

    /*
     * We write the closely related diagnostic variable P3 instead of the prognostic variable H, 
     * since pressure is more intuitive and H does not contain the pressure data at the top of the model.
     */
    if (var.p3.on) {
      if ((portion == EXTRACT_DATA && var.p3.extract_on) ||
          (portion != EXTRACT_DATA                     )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.p3.info[0].name,
                    var.p3.info[0].index,var.p3.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.theta.on) {
      if ((portion == EXTRACT_DATA && var.theta.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.theta.info[0].name,
                    var.theta.info[0].index,var.theta.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    if (var.fpara.on) {
      if ((portion == EXTRACT_DATA && var.fpara.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.fpara.info[0].name,
                    var.fpara.info[0].index,var.fpara.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }

    for (iq = 0; iq < grid.nq; iq++) {
      if ((portion == EXTRACT_DATA && var.species[grid.is[iq]].phase[grid.ip[iq]].extract_on) ||
          (portion != EXTRACT_DATA                                                          )   ) {
        /*
         * Map Q <= Q_MIN -> 0.0, so that extract file holds zero
         * instead of Q_MIN for voids, which makes plotting easier. 
         */
        for (K = KLOPAD; K <= KHIPAD; K++) {
          for (J = JLOPAD; J <= JHIPAD; J++) {
            for (I = ILOPAD; I <= IHIPAD; I++) {
              if (fcmp(Q(grid.is[iq],grid.ip[iq],K,J,I),Q_MIN) <= 0) {
                Q(grid.is[iq],grid.ip[iq],K,J,I) = 0.;
              }
            }
          }
        }
        /* No need to apply bc_lateral() here. */

        if ((portion == EXTRACT_DATA && grid.extract_species_fraction_type == MASS) ||
            (portion != EXTRACT_DATA                                              )   ) {
          write_array(node,FOURDIM,start,end,stretch_ni,var.species[grid.is[iq]].phase[grid.ip[iq]].info[MASS].name,
                      var.species[grid.is[iq]].phase[grid.ip[iq]].info[MASS].index,
                      var.species[grid.is[iq]].phase[grid.ip[iq]].q,
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        else if (grid.extract_species_fraction_type == MOLAR) {
          sync_x_to_q(planet);

          write_array(node,FOURDIM,start,end,stretch_ni,var.species[grid.is[iq]].phase[grid.ip[iq]].info[MOLAR].name,
                      var.species[grid.is[iq]].phase[grid.ip[iq]].info[MOLAR].index,
                      var.species[grid.is[iq]].phase[grid.ip[iq]].x,
                      EPIC_FLOAT_ARRAY,nc_id);
        }
        else {
          sprintf(Message,"grid.extract_species_fraction_type=%d unrecognized",grid.extract_species_fraction_type);
          epic_error(dbmsname,Message);
        }

        /*
         * Restore Q < Q_MIN to be Q_MIN.
         */
        restore_mass(planet,grid.is[iq],grid.ip[iq]);
        sync_x_to_q(planet);
      }
    }

    if (var.nu_turb.on) {
      if ((portion == EXTRACT_DATA && var.nu_turb.extract_on) ||
          (portion != EXTRACT_DATA                    )   ) {
        write_array(node,FOURDIM,start,end,stretch_ni,var.nu_turb.info[0].name,
                    var.nu_turb.info[0].index,var.nu_turb.value,EPIC_FLOAT_ARRAY,nc_id);
      }
    }
  }

  /* 
   * Close NODE0's data connection.
   */
  if (IAMNODE == NODE0) {
    nc_close(nc_id);
    fprintf(stdout,"\b\b\b\b%3d%%\n",100);
    fflush(stdout);
  }

  return;
}

/*======================= end of var_write() ==================================*/

/*======================= lookup_netcdf() =====================================*/

/*
 * Read numbers and names of global attributes and variables 
 * contained in infile, which must be in netCDF format.
 * Returns NC_NOERR if no read error, otherwise returns nc_err.
 *
 * NOTE: ngatts,**gattname,num_progs,**varname should be declared static in the calling
 *       function, with their input values equal to the last call, in order to
 *       properly reallocate memory.
 */

int lookup_netcdf(char   *infile,
                  int    *nc_id,
                  int    *ngatts,
                  char ***gattname,
                  int    *num_vars,
                  char ***varname)
{
  int
    i,
    ndims,unlimdimid,
    nc_err = NC_NOERR;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="lookup_netcdf";

  /*
   * Free previous memory:
   */
  for (i = 0; i < *ngatts; i++) {
    free((*gattname)[i]);
  }
  for (i = 0; i < *num_vars; i++) {
    free((*varname)[i]);
  }

  if (IAMNODE == NODE0) {
    /*
     * Decompress infile if necessary.
     */
    handle_file_compression(infile);
    nc_err = nc_open(infile,NC_NOWRITE,nc_id);
    if (nc_err == NC_NOERR) {
      nc_inq(*nc_id,&ndims,num_vars,ngatts,&unlimdimid);
    }
    else {
      return nc_err;
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(&nc_err,1,MPI_INT,NODE0,para.comm);
#endif

#if defined(EPIC_MPI)
  MPI_Bcast(ngatts,   1,MPI_INT,NODE0,para.comm);
  MPI_Bcast(num_vars, 1,MPI_INT,NODE0,para.comm);
#endif

  /*
   * Reallocate memory for character arrays. 
   * Look up and store names.
   */
  *gattname = (char **)realloc(*gattname,(*ngatts)*sizeof(char *));
  for (i = 0; i < *ngatts; i++) {
    (*gattname)[i] = (char *)calloc(NC_MAX_NAME,sizeof(char));
    if (IAMNODE == NODE0) {
      nc_inq_attname(*nc_id,NC_GLOBAL,i,(*gattname)[i]);
    }

#if defined (EPIC_MPI)
    MPI_Bcast((*gattname)[i],NC_MAX_NAME,MPI_CHAR,NODE0,para.comm);
#endif 

  }

  *varname = (char **)realloc(*varname,(*num_vars)*sizeof(char *));
  for (i = 0; i < *num_vars; i++) {
    (*varname)[i] = (char *)calloc(NC_MAX_NAME,sizeof(char));
    if (IAMNODE == NODE0) {
      nc_inq_varname(*nc_id,i,(*varname)[i]);
    }

#if defined (EPIC_MPI)
    MPI_Bcast((*varname)[i],NC_MAX_NAME,MPI_CHAR,NODE0,para.comm);
#endif

  }

  /* No errors. */
  return nc_err;
}

/*======================= end of lookup_netcdf() ==============================*/

/*======================= define_netcdf() =====================================*/
/*
 * Define variables and attributes for netCDF file.
 * We are following the CF conventions.
 */

#define DEFINE_NC_GRID(u) \
    var.u.info[0].dim = FOURDIM; \
    /* lon (I direction) */ \
    sprintf(coord_name,"lon_%s",var.u.info[0].name); \
    nc_err = nc_def_dim(nc_id,coord_name,idim,&var.u.info[0].dimid[NETCDF_I_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_def_var(nc_id,coord_name,float_type,ONEDIM, \
                        &var.u.info[0].dimid[NETCDF_I_INDEX],&var.u.info[0].coorid[NETCDF_I_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"standard_name", \
                             strlen("longitude")+1,"longitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"long_name", \
                             strlen("longitude")+1,"longitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"axis", \
                             strlen("X")+1,"X"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],"units", \
                             strlen("degrees_east")+1,"degrees_east"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    /* lat (J direction) */ \
    sprintf(coord_name,"lat_%s",var.u.info[0].name); \
    nc_err = nc_def_dim(nc_id,coord_name,jdim,&var.u.info[0].dimid[NETCDF_J_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_def_var(nc_id,coord_name,float_type,ONEDIM, \
                        &var.u.info[0].dimid[NETCDF_J_INDEX],&var.u.info[0].coorid[NETCDF_J_INDEX]); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"standard_name", \
                             strlen("latitude")+1,"latitude"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"long_name", \
                             strlen("latitude (planetographic)")+1,"latitude (planetographic)"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"axis", \
                             strlen("Y")+1,"Y"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],"units", \
                             strlen("degrees_north")+1,"degrees_north"); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
    } \
    /* vertical (K direction) */ \
    switch(grid.coord_type) { \
      case COORD_ISENTROPIC: \
        sprintf(coord_name,"theta_%s",var.u.info[0].name); \
        nc_err = nc_def_dim(nc_id,coord_name,kdim,&var.u.info[0].dimid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_def_var(nc_id,coord_name,NC_DOUBLE,ONEDIM, \
                            &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info[0].coorid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"standard_name", \
                                 strlen("air_potential_temperature")+1,"air_potential_temperature"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"long_name", \
                                 strlen("isentropic vertical coordinate")+1,"isentropic vertical coordinate"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"positive", \
                                 strlen("up")+1,"up"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"axis", \
                                 strlen("Z")+1,"Z"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"units", \
                                 strlen("K")+1,"K"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
      break; \
      case COORD_HYBRID: \
        sprintf(coord_name,"hybrid_sigma_theta%s",var.u.info[0].name); \
        nc_err = nc_def_dim(nc_id,coord_name,kdim,&var.u.info[0].dimid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_def_var(nc_id,coord_name,NC_DOUBLE,ONEDIM, \
                            &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info[0].coorid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"long_name", \
                                 strlen("hybrid sigma-isentropic vertical coordinate")+1,"hybrid sigma-isentropic vertical coordinate"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"positive", \
                                 strlen("up")+1,"up"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"axis", \
                                 strlen("Z")+1,"Z"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"units", \
                                 strlen("K")+1,"K"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
      break; \
      case COORD_ISOBARIC: \
        sprintf(coord_name,"p_%s",var.u.info[0].name); \
        nc_err = nc_def_dim(nc_id,coord_name,kdim,&var.u.info[0].dimid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_def_var(nc_id,coord_name,NC_DOUBLE,ONEDIM, \
                            &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info[0].coorid[NETCDF_K_INDEX]); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"standard_name", \
                                 strlen("air_pressure")+1,"air_pressure"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"long_name", \
                                 strlen("isobaric vertical coordinate")+1,"isobaric vertical coordinate"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"positive", \
                                 strlen("down")+1,"down"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"axis", \
                                 strlen("Z")+1,"Z"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],"units", \
                                 strlen("Pa")+1,"Pa"); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_GRID: %s\n",nc_strerror(nc_err)); \
        } \
      break; \
      default: \
        sprintf(Message,"unrecognized grid.coord_type=%d",grid.coord_type); \
        epic_error(dbmsname,Message); \
      break; \
    } \
    /* time */ \
    var.u.info[0].dimid[ NETCDF_T_INDEX] = var.info[0].dimid[NETCDF_T_INDEX]; \
    var.u.info[0].coorid[NETCDF_T_INDEX] = var.info[0].coorid[NETCDF_T_INDEX];

#define DEFINE_NC_VAR(p,h,num,on_array,has_standard_name) \
    for (in = 0; in < num; in++) { \
      if (on_array[in]) { \
        var.p.info[in].dim = var.h.info[0].dim;     /* rhs index should be 0 */ \
        for (i = 0; i < var.p.info[in].dim; i++) { \
          var.p.info[in].dimid[i]  = var.h.info[0].dimid[i]; \
          var.p.info[in].coorid[i] = var.h.info[0].coorid[i]; \
        } \
        nc_err = nc_def_var(nc_id,var.p.info[in].name,float_type,var.p.info[in].dim, \
                            var.p.info[in].dimid,&var.p.info[in].id); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
        } \
        if (has_standard_name) { \
          nc_err = nc_put_att_text(nc_id,var.p.info[in].id,"standard_name", \
                                   strlen(var.p.info[in].standard_name)+1,var.p.info[in].standard_name); \
          if (nc_err != NC_NOERR) { \
            fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
          } \
        } \
        nc_err = nc_put_att_text(nc_id,var.p.info[in].id,"long_name", \
                                 strlen(var.p.info[in].long_name)+1,var.p.info[in].long_name); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.p.info[in].id,"units", \
                                 strlen(var.p.info[in].units)+1,var.p.info[in].units); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_VAR: %s\n",nc_strerror(nc_err)); \
        } \
      } \
    }

#define DEFINE_NC_JI(phi_surface,h,has_standard_name) \
    nc_err = nc_def_var(nc_id,var.phi_surface.info[0].name,float_type,TWODIM, \
                        &var.h.info[0].dimid[NETCDF_J_INDEX],&var.phi_surface.info[0].id); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    } \
    if (has_standard_name) { \
      nc_err = nc_put_att_text(nc_id,var.phi_surface.info[0].id,"standard_name", \
                               strlen(var.phi_surface.info[0].standard_name)+1,var.phi_surface.info[0].standard_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
      } \
    } \
    nc_err = nc_put_att_text(nc_id,var.phi_surface.info[0].id,"long_name", \
                             strlen(var.phi_surface.info[0].long_name)+1,var.phi_surface.info[0].long_name); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.phi_surface.info[0].id,"units", \
                             strlen(var.phi_surface.info[0].units)+1,var.phi_surface.info[0].units); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_JI: %s\n",nc_strerror(nc_err)); \
    }

#define DEFINE_NC_KJ(gravity,h,has_standard_name) \
    nc_err = nc_def_var(nc_id,var.gravity.info[0].name,float_type,TWODIM, \
                        &var.h.info[0].dimid[NETCDF_K_INDEX],&var.gravity.info[0].id); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_KJ: %s\n",nc_strerror(nc_err)); \
    } \
    if (has_standard_name) { \
      nc_err = nc_put_att_text(nc_id,var.gravity.info[0].id,"standard_name", \
                               strlen(var.gravity.info[0].standard_name)+1,var.gravity.info[0].standard_name); \
      if (nc_err != NC_NOERR) { \
        fprintf(stderr,"DEFINE_NC_KJ: %s\n",nc_strerror(nc_err)); \
      } \
    } \
    nc_err = nc_put_att_text(nc_id,var.gravity.info[0].id,"long_name", \
                             strlen(var.gravity.info[0].long_name)+1,var.gravity.info[0].long_name); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_KJ: %s\n",nc_strerror(nc_err)); \
    } \
    nc_err = nc_put_att_text(nc_id,var.gravity.info[0].id,"units", \
                             strlen(var.gravity.info[0].units)+1,var.gravity.info[0].units); \
    if (nc_err != NC_NOERR) { \
      fprintf(stderr,"DEFINE_NC_KJ: %s\n",nc_strerror(nc_err)); \
    }

#define DEFINE_NC_TEND(u,num,on_array) \
    for (in = 0; in < num; in++) { \
      if (on_array[in]) { \
        nc_err = nc_def_var(nc_id,var.u.info_tend[in].name,float_type,THREEDIM, \
                            &var.u.info[0].dimid[NETCDF_K_INDEX],&var.u.info_tend[in].id); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"units", \
                                 strlen(var.u.info_tend[in].units)+1,var.u.info_tend[in].units); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"standard_name", \
                                 strlen(var.u.info_tend[in].standard_name)+1,var.u.info_tend[in].standard_name); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
        } \
        nc_err = nc_put_att_text(nc_id,var.u.info_tend[in].id,"long_name", \
                                 strlen(var.u.info_tend[in].long_name)+1,var.u.info_tend[in].long_name); \
        if (nc_err != NC_NOERR) { \
          fprintf(stderr,"DEFINE_NC_TEND: %s\n",nc_strerror(nc_err)); \
        } \
      } \
    }

/*
 * Specify maximum number of arrays associated with a given variable.
 * The array on_array[] is used to indicate which of these is to be defined for .nc files.
 */
#define MAX_NUM 2

void define_netcdf(planetspec   *planet,
                   char         *outfile,
                   int           portion,
                   int           stretch_ni,
                   int           nc_id)
{
  int
    nc_err,
    kdim,jdim,idim,
    K,J,I,
    is,ip,
    in,i,
    num,
    on_array[MAX_NUM];
  size_t
    index[1],
    message_len;
  char
    coord_name[VAR_NM_SZ+8];
  nc_type
    float_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="define_netcdf";

  if (IAMNODE != NODE0) {
    return;
  }

  /*
   * Layers 0 to grid.nk+1.
   */
  kdim = grid.nk+2;

  jdim = grid.nj-grid.jlo+1;

  if (stretch_ni > 0) {
    idim = stretch_ni;
  }
  else {
    idim = grid.ni;
  }

  /*
   * Specify floating-point type (precision).
   */
  if (EPIC_PRECISION == DOUBLE_PRECISION) {
    float_type = NC_DOUBLE;
  }
  else {
    float_type = NC_FLOAT;
  }

  /*
   * The variables all share the same time dimension, which is unlimited.
   *
   * time:
   */
  nc_def_dim(nc_id,"time",NC_UNLIMITED,&var.info[0].dimid[NETCDF_T_INDEX]); 
  nc_def_var(nc_id,"time",float_type,ONEDIM,
             &var.info[0].dimid[NETCDF_T_INDEX],&var.info[0].coorid[NETCDF_T_INDEX]);
  message_len = strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S 0",gmtime(&var.start_time));

  nc_err = nc_put_att_text(nc_id,var.u.info[0].coorid[NETCDF_T_INDEX],"axis",
                           strlen("T")+1,"T");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  } 

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"units",
                           message_len+1,Message);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"standard_name",
                           strlen("time")+1,"time");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"long_name",
                           strlen("elapsed time")+1,"elapsed time");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  nc_err = nc_put_att_text(nc_id,var.info[0].coorid[NETCDF_T_INDEX],"calendar",
                           strlen("julian")+1,"julian");
  if (nc_err != NC_NOERR) {
    sprintf(Message,"%s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  /*
   * Define solar longitude variable, L_s, which is a function of time.
   * L_s = 0 corresponds to the planet's vernal equinox.
   *
   * NOTE: EPIC uses dynamic North, not IAU North, consequently L_s is shifted by 180deg
   *       from the IAU values for planets with IAU obliquity greater than 90deg
   *       (notably Venus, Uranus and Pluto).
   */
  var.l_s.info.dim       = ONEDIM;
  var.l_s.info.dimid[0]  = var.info[0].dimid[NETCDF_T_INDEX];
  var.l_s.info.coorid[0] = var.info[0].coorid[NETCDF_T_INDEX];
  nc_err = nc_def_var(nc_id,var.l_s.info.name,float_type,var.l_s.info.dim,
                      var.l_s.info.dimid,&var.l_s.info.id);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }
  nc_err = nc_put_att_text(nc_id,var.l_s.info.id,"long_name",
                           strlen(var.l_s.info.long_name)+1,var.l_s.info.long_name);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }
  nc_err = nc_put_att_text(nc_id,var.l_s.info.id,"units",
                             strlen(var.l_s.info.units)+1,var.l_s.info.units);
  if (nc_err != NC_NOERR) {
    fprintf(stderr,"defining L_s: %s\n",nc_strerror(nc_err));
  }

  /*
   * The variables u, v, h, pv2, and p3 each reside on a different staggered
   * grid in our C-grid scheme. Thus, the grids for these variables are defined first.
   */
  DEFINE_NC_GRID(u);
  DEFINE_NC_GRID(v);
  DEFINE_NC_GRID(h);
  DEFINE_NC_GRID(pv2);
  DEFINE_NC_GRID(p3);

  /*
   * Define in netcdf the prognostic variables that are written to epic.nc files.
   */
  if (var.u.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.u.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(u,u,num,on_array,HAS_STANDARD_NAME);
      }
    } 
    else {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(u,u,num,on_array,HAS_STANDARD_NAME);
        if (portion == ALL_DATA) {
          num         = 2;
          on_array[0] = 1;
          on_array[1] = 1;
          DEFINE_NC_TEND(u,num,on_array);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        num         = 2;
        on_array[0] = 1;
        on_array[1] = 1;
        DEFINE_NC_VAR(u,u,num,on_array,HAS_STANDARD_NAME);
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  if (var.v.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.v.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(v,v,num,on_array,HAS_STANDARD_NAME);
      }
    } 
    else {
      if (strcmp(grid.uv_timestep_scheme,"3rd-order Adams-Bashforth") == 0) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(v,v,num,on_array,HAS_STANDARD_NAME);
        if (portion == ALL_DATA) {
          num         = 2;
          on_array[0] = 1;
          on_array[1] = 1;
          DEFINE_NC_TEND(v,num,on_array);
        }
      }
      else if (strcmp(grid.uv_timestep_scheme,"Leapfrog (Asselin filtered)") == 0) {
        num         = 2;
        on_array[0] = 1;
        on_array[1] = 1;
        DEFINE_NC_VAR(v,v,num,on_array,HAS_STANDARD_NAME);
      }
      else {
        sprintf(Message,"unrecognized grid.uv_timestep_scheme: %s",grid.uv_timestep_scheme);
        epic_error(dbmsname,Message);
      }
    }
  }

  /*
   * NOTE: P3 is treated as the prognostic variable for i/o rather than H
   */
  if (var.p3.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.p3.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(p3,p3,num,on_array,HAS_STANDARD_NAME);
      }
    }
    else {
      num         = 1;
      on_array[0] = 1;
      DEFINE_NC_VAR(p3,p3,num,on_array,HAS_STANDARD_NAME);
    }
  }

  if (var.theta.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.theta.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(theta,p3,num,on_array,HAS_STANDARD_NAME);
      }
    }
    else {
      num         = 1;
      on_array[0] = 1;
      DEFINE_NC_VAR(theta,p3,num,on_array,HAS_STANDARD_NAME);
    }
  }

  if (var.fpara.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.fpara.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(fpara,p3,num,on_array,NEEDS_STANDARD_NAME);
      }
    }
    else {
      num         = 1;
      on_array[0] = 1;
      DEFINE_NC_VAR(fpara,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
  }

  /*
   * NOTE: Not using grid.nq, grid.is[] and grid.ip[] here because of
   *       the extract_on conditional.
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
        if (var.species[is].extract_on) {
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            if (var.species[is].phase[ip].extract_on) {
              num = 2;
              if (grid.extract_species_fraction_type == MASS) {
                on_array[MASS ] = 1;
                on_array[MOLAR] = 0;
              }
              else if (grid.extract_species_fraction_type == MOLAR) {
                on_array[MASS ] = 0;
                on_array[MOLAR] = 1;
              }
              else {
                sprintf(Message,"unrecognized grid.extract_species_fraction_type=%d",grid.extract_species_fraction_type);
                epic_error(dbmsname,Message);
              }
              DEFINE_NC_VAR(species[is].phase[ip],p3,num,on_array,NEEDS_STANDARD_NAME);
            }
          }
        }
      }
      else {
        for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
          if (var.species[is].phase[ip].on) {
            num = 2;
            on_array[MASS ] = 1;
            on_array[MOLAR] = 0;
            DEFINE_NC_VAR(species[is].phase[ip],p3,num,on_array,NEEDS_STANDARD_NAME);
          }
        }
      }
    }
  }

  if (var.nu_turb.on) {
    if (portion == EXTRACT_HEADER_DATA || portion == EXTRACT_DATA) {
      if (var.nu_turb.extract_on) {
        num         = 1;
        on_array[0] = 1;
        DEFINE_NC_VAR(nu_turb,h,num,on_array,NEEDS_STANDARD_NAME);
      }
    }
    else {
      num         = 1;
      on_array[0] = 1;
      DEFINE_NC_VAR(nu_turb,h,num,on_array,NEEDS_STANDARD_NAME);
    }
  }

  /*
   * Define in netcdf the diagnostic variables that may optionally be written to extract.nc files.
   */
  num         = 1;
  on_array[0] = 1;

  if (portion == EXTRACT_HEADER_DATA) {
    if (var.h.extract_on) {
      DEFINE_NC_VAR(h,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.h3.extract_on) {
      DEFINE_NC_VAR(h3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.hdry2.extract_on) {
      DEFINE_NC_VAR(hdry2,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.hdry3.extract_on) {
      DEFINE_NC_VAR(hdry3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.pdry3.extract_on) {
      DEFINE_NC_VAR(pdry3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.p2.extract_on) {
      DEFINE_NC_VAR(p2,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.theta2.extract_on) {
      DEFINE_NC_VAR(theta2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.t2.extract_on) {
      DEFINE_NC_VAR(t2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.t3.extract_on) {
      DEFINE_NC_VAR(t3,p3,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.rho2.extract_on) {
      DEFINE_NC_VAR(rho2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.rho3.extract_on) {
      DEFINE_NC_VAR(rho3,p3,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.exner2.extract_on) {
      DEFINE_NC_VAR(exner2,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.exner3.extract_on) {
      DEFINE_NC_VAR(exner3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.fgibb2.extract_on) {
      DEFINE_NC_VAR(fgibb2,h,num,on_array,NEEDS_STANDARD_NAME); 
    }
    if (var.phi2.extract_on) {
      DEFINE_NC_VAR(phi2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.phi3.extract_on) {
      DEFINE_NC_VAR(phi3,p3,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.mont2.extract_on) {
      DEFINE_NC_VAR(mont2,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.ri2.extract_on) {
      DEFINE_NC_VAR(ri2,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.heat3.extract_on) {
      DEFINE_NC_VAR(heat3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.rel_vort2.extract_on) {
      DEFINE_NC_VAR(rel_vort2,pv2,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.eddy_rel_vort2.extract_on) {
      DEFINE_NC_VAR(eddy_rel_vort2,pv2,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.abs_vort2.extract_on) {
      DEFINE_NC_VAR(abs_vort2,pv2,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.kinetic_energy2.extract_on) {
      DEFINE_NC_VAR(kinetic_energy2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if(grid.coord_type == COORD_ISENTROPIC || 
       grid.coord_type == COORD_HYBRID) {
      /*
       * Ertel Potential vorticity is defined on potential temperature surfaces, but not
       * on pressure surfaces.  There is still a PV2 variable internally, but its denominator ("h")
       * is a constant on pressure surfaces, hence PV2 is proportional to ABS_VORT2 (absolute vorticity)
       * when using pressure coordinates.
       */
      if (var.pv2.extract_on) {
        DEFINE_NC_VAR(pv2,pv2,num,on_array,HAS_STANDARD_NAME);
      }
      if (var.eddy_pv2.extract_on) {
        DEFINE_NC_VAR(eddy_pv2,pv2,num,on_array,NEEDS_STANDARD_NAME);
      }
    }
    if (var.molar_mass3.extract_on) {
      DEFINE_NC_VAR(molar_mass3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.div_uv2.extract_on) {
      DEFINE_NC_VAR(div_uv2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.w3.extract_on) {
      DEFINE_NC_VAR(w3,p3,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.z2.extract_on) {
      DEFINE_NC_VAR(z2,h,num,on_array,HAS_STANDARD_NAME);
    }
    if (var.diffusion_coef_uv.extract_on) {
      DEFINE_NC_VAR(diffusion_coef_uv,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.diffusion_coef_theta.extract_on) {
      DEFINE_NC_VAR(diffusion_coef_theta,h,num,on_array,NEEDS_STANDARD_NAME);
    }
    if (var.diffusion_coef_mass.extract_on) {
      DEFINE_NC_VAR(diffusion_coef_mass,h,num,on_array,NEEDS_STANDARD_NAME);
    }
  }

  /*
   * Define for netCDF the parameter and key diagnostic arrays.
   */
  if (var.dzdt2.on) {
    if ((portion == EXTRACT_HEADER_DATA && var.dzdt2.extract_on) ||
        (portion != EXTRACT_HEADER_DATA)                           ) {
      DEFINE_NC_VAR(dzdt2,h,num,on_array,HAS_STANDARD_NAME);
    }
  }

  if (var.phi_surface.on) {
    if ((portion == EXTRACT_HEADER_DATA && var.phi_surface.extract_on) ||
        (portion != EXTRACT_HEADER_DATA)                                ) {
      DEFINE_NC_JI(phi_surface,h,HAS_STANDARD_NAME);
    }
  }

  if (var.pbot.on) {
    if ((portion == EXTRACT_HEADER_DATA && var.pbot.extract_on) ||
        (portion != EXTRACT_HEADER_DATA)                                ) {
      DEFINE_NC_JI(pbot,h,NEEDS_STANDARD_NAME);
    }
  }

  if (var.gravity2.on) {
    DEFINE_NC_KJ(gravity2,h,NEEDS_STANDARD_NAME);
  }

  /*
   * Leave define mode:
   */
  nc_enddef(nc_id);

  /*
   * Assign values to staggered C-grid coordinates.
   */

  /*
   * longitude:
   */
  if (stretch_ni > 0) {
    /*
     * Need to compute and output longitudes appropriate to stretch_ni.
     */
    EPIC_FLOAT
      *longitude;
    int
      ii;

    /*
     * Allocate memory.
     */
    longitude = fvector(0,2*(stretch_ni+1),dbmsname);

    /*
     * NOTE: The value of grid.dln should already be modified accordingly.
     */

    if (strcmp(grid.geometry,"globe") == 0) {
      longitude[0] = grid.globe_lonbot-grid.dln;
    }
    else {
      longitude[0] = -180.-grid.dln;
    }
    for (ii = 1; ii <= 2*(stretch_ni+1); ii++) {
      longitude[ii] = longitude[ii-1]+grid.dln*.5;
    }

    for (I = 1; I <= stretch_ni; I++) {
      index[0] = I-1;

#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I]));
      nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
      nc_put_var1_double(nc_id,var.h.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
      nc_put_var1_double(nc_id,var.pv2.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I]));
      nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(longitude[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I]));
      nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
      nc_put_var1_float(nc_id,var.h.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
      nc_put_var1_float(nc_id,var.pv2.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I]));
      nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(longitude[2*I+1]));
#endif
    }

    /*
     * Free allocated memory.
     */
    free_fvector(longitude,0,2*(stretch_ni+1),dbmsname);
  }
  else {
    for (I = 1; I <= grid.ni; I++) {
      index[0] = I-1;

#if EPIC_PRECISION == DOUBLE_PRECISION
      nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I]));
      nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
      nc_put_var1_double(nc_id,var.h.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
      nc_put_var1_double(nc_id,var.pv2.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I]));
      nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                         index,&(grid.lon[2*I+1]));
#else
      nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I]));
      nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
      nc_put_var1_float(nc_id,var.h.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
      nc_put_var1_float(nc_id,var.pv2.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I]));
      nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_I_INDEX],
                        index,&(grid.lon[2*I+1]));
#endif
    }
  }

  /*
   * latitude:
   */
  for (J = grid.jlo; J <= grid.nj; J++) {
    index[0] = J-grid.jlo;

#if EPIC_PRECISION == DOUBLE_PRECISION
    nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J]));
    nc_put_var1_double(nc_id,var.h.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
    nc_put_var1_double(nc_id,var.pv2.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J]));
    nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_J_INDEX],
                       index,&(grid.lat[2*J+1]));
#else
    nc_put_var1_float(nc_id,var.u.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
    nc_put_var1_float(nc_id,var.v.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J]));
    nc_put_var1_float(nc_id,var.h.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
    nc_put_var1_float(nc_id,var.pv2.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J]));
    nc_put_var1_float(nc_id,var.p3.info[0].coorid[NETCDF_J_INDEX],
                      index,&(grid.lat[2*J+1]));
#endif
  }

  /*
   * Vertical coordinate:
   */
  switch(grid.coord_type) {
    case COORD_ISENTROPIC:
    case COORD_HYBRID:
      for (K = 0; K <= grid.nk+1; K++) {
        index[0] = K;
        nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.sigmatheta[2*K  ]));
        nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.sigmatheta[2*K  ]));
        nc_put_var1_double(nc_id,var.h.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.sigmatheta[2*K  ]));
        nc_put_var1_double(nc_id,var.pv2.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.sigmatheta[2*K  ]));
        nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.sigmatheta[2*K+1]));
      }
    break;
    case COORD_ISOBARIC:
      for (K = 0; K <= grid.nk+1; K++) {
        index[0] = K;
        nc_put_var1_double(nc_id,var.u.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.p_ref[2*K  ]));
        nc_put_var1_double(nc_id,var.v.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.p_ref[2*K  ]));
        nc_put_var1_double(nc_id,var.h.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.p_ref[2*K  ]));
        nc_put_var1_double(nc_id,var.pv2.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.p_ref[2*K  ]));
        nc_put_var1_double(nc_id,var.p3.info[0].coorid[NETCDF_K_INDEX],
                           index,&(grid.p_ref[2*K+1]));
      }
    break;
    default:
      sprintf(Message,"unrecognized grid.coord_type=%d",grid.coord_type);
      epic_error(dbmsname,Message);
    break;
  }

  /*
   * NOTE: The values for the time dimension are written when
   *       the variables are written.
   */

  /*
   * Put back into define mode before returning:
   */
  nc_err = nc_redef(nc_id);
  if (nc_err != NC_NOERR) {
    sprintf(Message,"nc_redef(), %s",nc_strerror(nc_err));
    epic_error(dbmsname,Message);
  }

  return;
}

/*======================= end of define_netcdf() ==============================*/

/*======================= handle_file_compression() ===========================*/

/*
 * Removes any .gz suffix from file_name and  
 * decompresses the file if necessary.
 */
void handle_file_compression(char *file_name)
{
  char
    *ptr,
    alt_file_name[FILE_STR];
  FILE
    *fp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="handle_file_compression";

  ptr = strrchr(file_name,'.');
  if (ptr && strstr(ptr,".gz")) {
    /*
     * The file_name string has a .gz suffix.
     *
     * Check if the compressed file exists.
     */
    fp = fopen(file_name,"r");
    if (fp) {
      /*
       * The compressed file exists.
       * Decompress it.
       */
      fclose(fp);
      sprintf(Message,"gunzip %s",file_name);
      system(Message);
    }
    /*
     * Remove the .gz suffix from file_name and return.
     */
    *ptr = '\0';
    return;
  }
  else {
    /*
     * The file_name string does not have a .gz suffix.
     *
     * Check if the uncompressed file exists.
     */
    fp = fopen(file_name,"r");
    if (fp) {
      /*
       * The uncompressed file exists, so return.
       */
      fclose(fp);
      return;
    }
    else {
      /*
       * The uncompressed file does not exist.
       * Check if the compressed file exists.
       */
      sprintf(alt_file_name,"%s.gz",file_name);
      fp = fopen(alt_file_name,"r");
      if (fp) {
        /*
         * The compressed file exists.
         * Decompress the file and return.
         */
        sprintf(Message,"gunzip %s",alt_file_name);
        system(Message);
        fclose(fp);
        return;
      }
    }
  }

  return;
}

/*======================= end of handle_file_compression() ====================*/

/*======================= get_jlohi() =========================================*/

/*
 * Determines j range of specified node.
 * The value node = SETUP_GET_JLOHI signals that the data are being passed
 * in via jlo and jhi.
 */

void get_jlohi(int  node,
               int  num_nodes,
               int *jlo,
               int *jhi)
{
  int
    n;
  static int
    initialized = FALSE,
   *jlow,
   *jhigh;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_jlohi";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    jlow  = ivector(0,num_nodes-1,dbmsname);
    jhigh = ivector(0,num_nodes-1,dbmsname);

    initialized = TRUE;
  }

  if (node >= 0 && node < num_nodes) {
    *jlo = jlow[node];
    *jhi = jhigh[node];
  }
  else if (node == SETUP_GET_JLOHI) {
   /*
    * The information is passed in through the jlo, jhi arguments. 
    */
    for (n = 0; n < num_nodes; n++) {
      jlow[n]  = jlo[n];
      jhigh[n] = jhi[n];
    }
  }
  else {
    sprintf(Message,"node=%d,num_nodes=%d is invalid",node,num_nodes);
    epic_error(dbmsname,Message);
  } 

  return;
}

/*======================= end of get_jlohi() ==================================*/

/*======================= get_ilohi() =========================================*/

/*
 * Determines i range of specified node.
 * The value node = SETUP_GET_ILOHI signals that the data are being passed
 * in via ilo and ihi.
 */

void get_ilohi(int  node,
               int  num_nodes,
               int *ilo,
               int *ihi)
{
  int
    n;
  static int
    initialized = FALSE,
   *ilow,
   *ihigh;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_ilohi";

  if (!initialized) {
    /*
     * Allocate memory.
     */
    ilow  = ivector(0,num_nodes-1,dbmsname);
    ihigh = ivector(0,num_nodes-1,dbmsname);

    initialized = TRUE;
  }

  if (node >= 0 && node < num_nodes) {
    *ilo = ilow[node];
    *ihi = ihigh[node];
  }
  else if (node == SETUP_GET_ILOHI) {
   /*
    * The information is passed in through the ilo, ihi arguments. 
    */
    for (n = 0; n < num_nodes; n++) {
      ilow[n]  = ilo[n];
      ihigh[n] = ihi[n];
    }
  }
  else {
    sprintf(Message,"node=%d,num_nodes=%d is invalid",node,num_nodes);
    epic_error(dbmsname,Message);
  } 

  return;
}

/*======================= end of get_ilohi() ==================================*/

/*======================= prompt_extract_on() =================================*/

#define PRINT_EXTRACT_ON(u,U_INDEX) \
    if (var.extract_on_list[U_INDEX] != NOT_LISTED) { \
      if (var.extract_on_list[U_INDEX] == LISTED_AND_ON) { \
        fprintf(stdout,"     %2d >< %s \n",U_INDEX,var.u.info[0].name); \
      } \
      else if (var.extract_on_list[U_INDEX] == LISTED_AND_OFF) { \
        fprintf(stdout,"     %2d    %s \n",U_INDEX,var.u.info[0].name); \
      } \
      var.u.extract_on = var.extract_on_list[U_INDEX] = LISTED_AND_OFF; \
    }

void prompt_extract_on(char *def_extract_str,
                       int  *def_extract_species_fraction_type)
{
  int
    is,ip,index,
    inquire_species_fraction_type = FALSE;
  char
    extract_str[N_STR]="",
    *ptr;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="prompt_extract_on";

  /*
   * Use def_extract_str to turn on extract_on_list.
   */
  var.extract_on = FALSE;
  for (index = FIRST_INDEX; index <= LAST_INDEX; index++) {
    var.extract_on_list[index] = LISTED_AND_OFF;
  }

  for (ptr = def_extract_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      if (index >= FIRST_SPECIES && index <= LAST_SPECIES) {
        if (var.species[index].on != TRUE) {
          /*
           * Remove from def_extract_str.
           */
          sprintf(ptr," ");
        }
        else {
          var.extract_on             = TRUE;
          var.extract_on_list[index] = LISTED_AND_ON;
        }
      }
      else {
        var.extract_on             = TRUE;
        var.extract_on_list[index] = LISTED_AND_ON;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  /*
   * Print list of variable indices and names, and then reset extract_on to FALSE.
   */
  fprintf(stdout,"\n");
  fprintf(stdout,"   index   variable\n");
  PRINT_EXTRACT_ON(u,    U_INDEX);
  PRINT_EXTRACT_ON(v,    V_INDEX);
  PRINT_EXTRACT_ON(h,    H_INDEX);
  PRINT_EXTRACT_ON(theta,THETA_INDEX);
  if (var.nu_turb.on) PRINT_EXTRACT_ON(nu_turb,NU_TURB_INDEX);
  if (var.fpara.on  ) PRINT_EXTRACT_ON(fpara,    FPARA_INDEX);
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      PRINT_EXTRACT_ON(species[is],is);
    }
  }
  PRINT_EXTRACT_ON(p2,P2_INDEX);
  PRINT_EXTRACT_ON(p3,P3_INDEX);
  PRINT_EXTRACT_ON(pdry3,PDRY3_INDEX);
  PRINT_EXTRACT_ON(t2,T2_INDEX);
  PRINT_EXTRACT_ON(t3,T3_INDEX);
  PRINT_EXTRACT_ON(rho2,RHO2_INDEX);
  PRINT_EXTRACT_ON(rho3,RHO3_INDEX);
  PRINT_EXTRACT_ON(exner2,EXNER2_INDEX);
  PRINT_EXTRACT_ON(exner3,EXNER3_INDEX);
  if (var.fgibb2.on) {PRINT_EXTRACT_ON(fgibb2,FGIBB2_INDEX);}
  PRINT_EXTRACT_ON(phi2,PHI2_INDEX);
  PRINT_EXTRACT_ON(phi3,PHI3_INDEX);
  PRINT_EXTRACT_ON(mont2,MONT2_INDEX);
  PRINT_EXTRACT_ON(heat3,HEAT3_INDEX);
  if (grid.coord_type == COORD_ISENTROPIC ||
      grid.coord_type == COORD_HYBRID) {
    PRINT_EXTRACT_ON(pv2,PV2_INDEX);
    PRINT_EXTRACT_ON(eddy_pv2,EDDY_PV2_INDEX);
  }
  PRINT_EXTRACT_ON(molar_mass3,MOLAR_MASS3_INDEX);
  PRINT_EXTRACT_ON(ri2,RI2_INDEX);
  PRINT_EXTRACT_ON(rel_vort2,REL_VORT2_INDEX);
  PRINT_EXTRACT_ON(eddy_rel_vort2,EDDY_REL_VORT2_INDEX);
  PRINT_EXTRACT_ON(abs_vort2,ABS_VORT2_INDEX);
  PRINT_EXTRACT_ON(kinetic_energy2,KIN2_INDEX);
  PRINT_EXTRACT_ON(div_uv2,DIV_UV2_INDEX);
  PRINT_EXTRACT_ON(w3,W3_INDEX);
  PRINT_EXTRACT_ON(z2,Z2_INDEX);
  PRINT_EXTRACT_ON(dzdt2,DZDT2_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_uv,DIFFUSION_COEF_UV_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_theta,DIFFUSION_COEF_THETA_INDEX);
  PRINT_EXTRACT_ON(diffusion_coef_mass,DIFFUSION_COEF_MASS_INDEX);
  if (var.phi_surface.on) {PRINT_EXTRACT_ON(phi_surface,PHI_SURFACE_INDEX);}
  fprintf(stdout,"\n");
  sprintf(Message,"On one line, input the indices of variables to be included in extract.nc\n");
  input_string(Message,def_extract_str,extract_str);

  memset(def_extract_str,0,N_STR);
  strcpy(def_extract_str,extract_str);

  for (ptr = extract_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      var.extract_on             = TRUE;
      var.extract_on_list[index] = LISTED_AND_ON;

      /* Turn on extract_on. */
      switch(index) {
        case U_INDEX:                    var.u.extract_on                    = TRUE; break;
        case V_INDEX:                    var.v.extract_on                    = TRUE; break;
        case H_INDEX:                    var.h.extract_on                    = TRUE; break;
        case THETA_INDEX:                var.theta.extract_on                = TRUE; break;
        case FPARA_INDEX:                var.fpara.extract_on                = TRUE; break;
        case NU_TURB_INDEX:              var.nu_turb.extract_on              = TRUE; break;
        case H3_INDEX:                   var.h3.extract_on                   = TRUE; break;
        case HDRY2_INDEX:                var.hdry2.extract_on                = TRUE; break;
        case HDRY3_INDEX:                var.hdry3.extract_on                = TRUE; break;
        case PDRY3_INDEX:                var.pdry3.extract_on                = TRUE; break;
        case P2_INDEX:                   var.p2.extract_on                   = TRUE; break;
        case P3_INDEX:                   var.p3.extract_on                   = TRUE; break;
        case THETA2_INDEX:               var.theta2.extract_on               = TRUE; break;
        case T2_INDEX:                   var.t2.extract_on                   = TRUE; break;
        case T3_INDEX:                   var.t3.extract_on                   = TRUE; break;
        case RHO2_INDEX:                 var.rho2.extract_on                 = TRUE; break;
        case RHO3_INDEX:                 var.rho3.extract_on                 = TRUE; break;
        case EXNER2_INDEX:               var.exner2.extract_on               = TRUE; break;
        case EXNER3_INDEX:               var.exner3.extract_on               = TRUE; break;
        case FGIBB2_INDEX:               var.fgibb2.extract_on               = TRUE; break;
        case PHI2_INDEX:                 var.phi2.extract_on                 = TRUE; break;
        case PHI3_INDEX:                 var.phi3.extract_on                 = TRUE; break;
        case MONT2_INDEX:                var.mont2.extract_on                = TRUE; break;
        case HEAT3_INDEX:                var.heat3.extract_on                = TRUE; break;
        case PV2_INDEX:                  var.pv2.extract_on                  = TRUE; break;
        case EDDY_PV2_INDEX:             var.eddy_pv2.extract_on             = TRUE; break;
        case MOLAR_MASS3_INDEX:          var.molar_mass3.extract_on          = TRUE; break;
        case RI2_INDEX:                  var.ri2.extract_on                  = TRUE; break;
        case REL_VORT2_INDEX:            var.rel_vort2.extract_on            = TRUE; break;
        case EDDY_REL_VORT2_INDEX:       var.eddy_rel_vort2.extract_on       = TRUE; break;
        case ABS_VORT2_INDEX:            var.abs_vort2.extract_on            = TRUE; break;
        case KIN2_INDEX:                 var.kinetic_energy2.extract_on      = TRUE; break;
        case DIV_UV2_INDEX:              var.div_uv2.extract_on              = TRUE; break;
        case W3_INDEX:                   var.w3.extract_on                   = TRUE; break;
        case Z2_INDEX:                   var.z2.extract_on                   = TRUE; break;
        case DZDT2_INDEX:                var.dzdt2.extract_on                = TRUE; break;
        case DIFFUSION_COEF_UV_INDEX:    var.diffusion_coef_uv.extract_on    = TRUE; break;
        case DIFFUSION_COEF_THETA_INDEX: var.diffusion_coef_theta.extract_on = TRUE; break;
        case DIFFUSION_COEF_MASS_INDEX:  var.diffusion_coef_mass.extract_on  = TRUE; break;
        case PHI_SURFACE_INDEX:          var.phi_surface.extract_on          = TRUE; break;
        case PBOT_INDEX:                 var.pbot.extract_on                 = TRUE; break;
        default:
          if (index < FIRST_SPECIES || index > LAST_SPECIES) {
            sprintf(Message,"unrecognized index=%d",index);
            epic_error(dbmsname,Message);
          }
          var.species[index].extract_on = TRUE;
          for (ip = FIRST_PHASE; ip <= LAST_PHASE; ip++) {
            if (var.species[index].phase[ip].on == 1) {
              var.species[index].phase[ip].extract_on = TRUE;
              inquire_species_fraction_type           = TRUE;
            }
          }
        break;
      } /* end of switch */

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  if (inquire_species_fraction_type) {
    sprintf(Message,"Species fraction type to extract: [%d] => Mass mixing ratio, Q = mass_i/mass_dry_air ,\n"
                    "                                  [%d] => Mole fraction,     X = number_i/number_total\n",
                    MASS,MOLAR);
    *def_extract_species_fraction_type = input_int(Message,*def_extract_species_fraction_type);
  }

  return;
}

/*======================= end of prompt_extract_on() ==========================*/

/*======================= prompt_species_on() =================================*/

/*
 * Returns the number of optional species turned on.
 */

#define PRINT_SPECIES_ON(index) \
    if (var.on_list[index] != NOT_LISTED) { \
      if (var.on_list[index] == LISTED_AND_ON) { \
        fprintf(stdout,"  %2d  ><  %s \n",index,var.species[index].info[0].name); \
      } \
      else if (var.on_list[index] == LISTED_AND_OFF) { \
        fprintf(stdout,"  %2d      %s \n",index,var.species[index].info[0].name); \
      } \
    } \
    var.species[index].on = var.on_list[index] = LISTED_AND_OFF;

int prompt_species_on(planetspec *planet,
                      char       *def_species_str)
{
  int
    is,ip,index,
    count = 0;
  char
    species_str[N_STR],
    *ptr;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="prompt_species_on";

  /*
   * Prune optional-variable list as appropriate for each atmosphere.
   * Set on_list to NOT_LISTED to take off the printed list.
   */
  for (index = FPARA_INDEX; index <= LAST_SPECIES; index++) {
    var.on_list[index]    = NOT_LISTED;
    var.species[index].on = FALSE;
  }

  if (strcmp(planet->name,"Venus") == 0) {
    /*
     * Venus case.
     */
    ;
  }
  else if (strcmp(planet->name,"Venus_LLR05") == 0) {
    /*
     * Venus_LLR05 benchmark case.
     */
    ;
  }
  else if (strcmp(planet->name,"Earth") == 0) {
    /*
     * Earth case.
     */
    var.on_list[H_2O_INDEX] = LISTED_AND_OFF;
    var.on_list[CO_2_INDEX] = LISTED_AND_OFF;
    var.on_list[O_3_INDEX ] = LISTED_AND_OFF;
  }
  else if (strcmp(planet->name,"Held_Suarez") == 0) {
    /*
     * Held-Suarez benchmark case.
     */
    ;
  }
  else if (strcmp(planet->name,"Goldmine") == 0) {
    /*
     * Goldmine testbed case.
     */
    ;
  }
  else if (strcmp(planet->type,"gas-giant") == 0) {
    /*
     * Gas-giant case.
     */
    var.on_list[FPARA_INDEX ] = LISTED_AND_OFF;
    var.on_list[H_2O_INDEX  ] = LISTED_AND_OFF;
    var.on_list[NH_3_INDEX  ] = LISTED_AND_OFF;
    var.on_list[H_2S_INDEX  ] = LISTED_AND_OFF;
    var.on_list[CH_4_INDEX  ] = LISTED_AND_OFF;
    var.on_list[C_2H_2_INDEX] = LISTED_AND_OFF;
    var.on_list[C_2H_4_INDEX] = LISTED_AND_OFF;
    var.on_list[C_2H_6_INDEX] = LISTED_AND_OFF;
    var.on_list[NH_4SH_INDEX] = LISTED_AND_OFF;
    var.on_list[PH_3_INDEX  ] = LISTED_AND_OFF;
  }
  else if (strcmp(planet->name,"Titan") == 0) {
    /*
     * Titan case.
     */
    var.on_list[CH_4_INDEX  ] = LISTED_AND_OFF;
    var.on_list[C_2H_2_INDEX] = LISTED_AND_OFF;
    var.on_list[C_2H_6_INDEX] = LISTED_AND_OFF;
  }
  else {
    /*
     * Default to have them all appear on the list.
     */
    for (index = FPARA_INDEX; index <= LAST_SPECIES; index++) {
      var.on_list[index] = LISTED_AND_OFF;
    }
  }

  /*
   * Use def_species_str to turn on var.on_list[index].
   */
  for (ptr = def_species_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      if (index == FPARA_INDEX) {
        var.fpara.on       = TRUE;
        var.on_list[index] = LISTED_AND_ON;
      }
      else if (index >= FIRST_SPECIES && index <= LAST_SPECIES) {
        var.species[index].on = TRUE;
        var.on_list[index]    = LISTED_AND_ON;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  /*
   * Print list of variable indices and names, and then reset to FALSE.
   */
  fprintf(stdout," Optional prognostic variables:\n");

  if (var.fpara.on == 1) {
    fprintf(stdout,"  %2d  ><  %s \n",FPARA_INDEX,var.fpara.info[0].name);
  }
  else if (planet->x_h2 > 0.) {
    fprintf(stdout,"  %2d      %s \n",FPARA_INDEX,var.fpara.info[0].name);
  }
  var.fpara.on = var.on_list[FPARA_INDEX] = LISTED_AND_OFF;

  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    PRINT_SPECIES_ON(is);
  }

  fprintf(stdout,"\n");
  sprintf(Message,"On one line, input the indices of optional prognostic variables to be turned on [none => optional variables off]\n");
  input_string(Message,def_species_str,species_str);

  if (strcmp(species_str,"\n")   == 0 ||
      strcmp(species_str,"none") == 0   ) {
    return count;
  }

  for (ptr = species_str; strcmp(ptr,"\0") != 0; ptr++) {
    if (isdigit(*ptr)) {
      /* Read in integer. */
      sscanf(ptr,"%d",&index);

      /* Turn on species. */
      if (index == FPARA_INDEX) {
        var.fpara.on             = TRUE; 
        var.on_list[FPARA_INDEX] = LISTED_AND_ON;
        count++;
      }
      else if (index < FIRST_SPECIES || index > LAST_SPECIES) {
        sprintf(Message,"unrecognized index=%d",index);
        epic_error(dbmsname,Message);
      }
      else {
        var.species[index].on = TRUE;
        var.on_list[index]    = LISTED_AND_ON;
        count++;
      }

      /* Skip to end of integer just read. */
      while (isdigit(*ptr)) {
        ptr++;
      }
      ptr--;
    }
  }

  return count;
}

/*======================= end of prompt_species_on() ==========================*/

/*======================= bcast_char() ========================================*/

void bcast_char(int   node,
                char *str,
                int   num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(str,num,MPI_CHAR,node,para.comm);
#endif


  return;
}

/*======================= end of bcast_char() =================================*/

/*======================= bcast_int() =========================================*/

void bcast_int(int  node,
               int *val,
               int  num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,MPI_INT,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_int() ==================================*/

/*======================= bcast_float() =======================================*/

void bcast_float(int         node,
                 EPIC_FLOAT *val,
                 int         num)
{

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,float_type,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_float() ================================*/

/*======================= bcast_double() ======================================*/

void bcast_double(int     node,
                  double *val,
                  int     num)
{

#if defined(EPIC_MPI)
  MPI_Bcast(val,num,MPI_DOUBLE,node,para.comm);
#endif

  return;
}

/*======================= end of bcast_double() ===============================*/

/*====================== read_spacing_file() ==================================*/
/*
 * Read layer spacing data from the file def->layer_spacing_dat.
 * The following lines without asterisks illustrate the format of this file;
 * one can cut cut and paste these lines to make a template.
   nk         10
   K        p[hPa]
   0          1.0
   1          2.6
   2          6.4
   3         14.2
   4         35.1
   5        118.1
   6        375.6
   7       2416.9
   8       3624.2 
   9       4344.0
  10       6000.0
 */   
void read_spacing_file(init_defaultspec *def,
                       int               mode)
{
  register int
    K,kk;
  char
    token[32],
   *filename;
  EPIC_FLOAT
   *p_ref;
  FILE
    *spacing_file;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_spacing_file";

  filename = def->layer_spacing_dat;

  spacing_file = fopen(filename,"r");
  if (!spacing_file) {
    sprintf(Message,"cannot open %s",filename);
    epic_error(dbmsname,Message); 
  }

  fscanf(spacing_file,"%s %d",token,&def->nk);
  grid.nk = def->nk;

  /*
   * Allocate memory.
   */
  p_ref = fvector(0,2*(grid.nk+1),dbmsname);

  fscanf(spacing_file,"%s %s",token,token);

#if EPIC_PRECISION == DOUBLE_PRECISION
  for (K = 0; K <= grid.nk; K++) {
    fscanf(spacing_file,"%s %lf",token,p_ref+(2*K+1));
  }
#else
  for (K = 0; K <= grid.nk; K++) {
    fscanf(spacing_file,"%s %f",token,p_ref+(2*K+1));
  }
#endif

  /*
   * Fill in layer values.
   */
  for (K = KLO; K <= KHI; K++) {
    kk = 2*K;
    p_ref[kk] = onto_kk(planet,P2_INDEX,p_ref[kk-1],p_ref[kk+1],kk,JLO,ILO);
  }
  p_ref[      0] = p_ref[1]*p_ref[1]/p_ref[2];
  p_ref[2*KHI+2] = p_ref[2*KHI+1]*p_ref[2*KHI+1]/p_ref[2*KHI];

  /*
   * Convert from hPa to Pa.
   */
  for (kk = 0; kk <= 2*(grid.nk+1); kk++) {
    p_ref[kk] *= 100.;
  }

  def->ptop = grid.ptop = p_ref[1];
  def->pbot = grid.pbot = p_ref[2*grid.nk+1];

  if (mode == ALL_DATA) {
    for (kk = 0; kk <= 2*(grid.nk+1); kk++) {
      grid.p_ref[kk] = p_ref[kk];
    }
  }

  fclose(spacing_file);
  /*
   * Free allocated memory.
   */
  free_fvector(p_ref,0,2*(grid.nk+1),dbmsname);

  return;
}

/*====================== end of read_spacing_file() ===========================*/

/*====================== get_sounding() =======================================*/
/*
 * Smooth, monotonic interpolation of t_vs_p data table.
 *
 * Valid output_name values: "temperature" 
 *                           "theta"
 *
 * For example, to get a theta(p) via interpolation, use 
 *   get_sounding(planet,p,"theta",&theta);
 *
 * Interpolation uses -log(p).
 *
 * NOTE: If input_value is below the bottom of the data table, then the bottom
 *       value of the data table is used. 
 */
void get_sounding(planetspec *planet,
                  EPIC_FLOAT  pressure,
                  char       *output_name,
                  EPIC_FLOAT *pt_output_value)
{
  int
    ki;
  static int
    initialized = FALSE;
  EPIC_FLOAT
    fpara,fgibb,fpe,uoup,
    theta_ortho,theta_para,
    x,x_d,
    tmp;
  static float_triplet
    *tdat,
    *thetadat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="get_sounding";

  if (!initialized) {
    /* 
     * Allocate memory: 
     */
    tdat     = ftriplet(0,var.ntp-1,dbmsname);
    thetadat = ftriplet(0,var.ntp-1,dbmsname);

    /*
     * Set temperature data table.
     */
    for (ki = 0; ki < var.ntp; ki++) {
      tdat[ki].x = -log(var.pdat[ki]);
      tdat[ki].y = var.tdat[ki];
    }

    spline_pchip(var.ntp,tdat);

    /*
     * Set theta data table.
     */
    for (ki = 0; ki < var.ntp; ki++) {
      thetadat[ki].x = -log(var.pdat[ki]);
      fpara          = return_fpe(var.tdat[ki]);
      thetadat[ki].y = return_theta(planet,fpara,var.pdat[ki],var.tdat[ki],&theta_ortho,&theta_para);
    }

    spline_pchip(var.ntp,thetadat);

    /*
     * From the top down, iron out negative-slope regions in theta.
     * Adjust temperature accordingly.
     */
    if (var.ntp-2 >= 0) {
      for (ki = var.ntp-2; ki >= 0; ki--) {
        tmp = thetadat[ki+1].y;
        if (thetadat[ki].y > tmp) {
          thetadat[ki].y = tmp;
          fpara          = return_fpe(var.tdat[ki]);
          /*
           * Adjust tdat.
           */
          tdat[ki].y = return_temp(planet,fpara,var.pdat[ki],thetadat[ki].y);
        }
      }
    } 
    initialized = TRUE;
  } 
  /* End of initialization. */

  /*
   * Interpolate using -log(p).
   */
  x = -log(pressure);

  /*
   * Restrict y(x) to not fall below bottom of data table.
   */
  if (x < tdat[0].x) {
    if (strcmp(output_name,"temperature") == 0) {
      *pt_output_value = tdat[0].y;
    }
    else if (strcmp(output_name,"theta") == 0) {
      *pt_output_value = thetadat[0].y;
    }
    else {
      sprintf(Message,"unrecognized output_name=%s",output_name);
      epic_error(dbmsname,Message);
    }
    sprintf(Message,"-ln(p)=%e < %e; setting %s[0]=%e \n",
                     x,tdat[0].x,output_name,*pt_output_value);
    epic_warning(dbmsname,Message);
  }
  else {
    if (strcmp(output_name,"temperature") == 0) {
      ki               = find_place_in_table(var.ntp,tdat,x,&x_d);
      *pt_output_value = splint_pchip(x,tdat+ki,x_d);
    }
    else if (strcmp(output_name,"theta") == 0) {
      ki               = find_place_in_table(var.ntp,thetadat,x,&x_d);
      *pt_output_value = splint_pchip(x,thetadat+ki,x_d);
    }
  }

  return;
}

/*======================= end of get_sounding() ===============================*/

/*======================= read_t_vs_p() =======================================*/

/*
 * NOTE: This function should only be called during initialization, not when the 
 *       model is running, and hence does not need to be MPI ready.
 */

int read_t_vs_p(planetspec *planet,
                int         portion)
{
  char
    infile[FILE_STR],
    header[N_STR];
  int
    nn,ntp;
  EPIC_FLOAT
    p1,t1,dt1;
  FILE
    *t_vs_p;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_t_vs_p";

#if defined(EPIC_MPI)
  /* NOTE: not set up for MPI. */
  sprintf(Message,"not set up for MPI");
  epic_error(dbmsname,Message);
#endif

  /* Look in local directory first.*/
  sprintf(infile,"./t_vs_p.%s",planet->name);
  t_vs_p = fopen(infile,"r");
  if (!t_vs_p) {
    sprintf(infile,EPIC_PATH"/data/%s/t_vs_p.%s",planet->name,planet->name);
    t_vs_p = fopen(infile,"r");
    if (!t_vs_p) {
      sprintf(Message,"Cannot open file %s",infile);
      epic_error(dbmsname,Message);
    }
  }

  /* Skip over 6-line header: */
  for (nn = 0; nn < 6; nn++) {
    fgets(header,N_STR,t_vs_p);
  }
  /* input number of data points */
  fscanf(t_vs_p,"%d",&ntp); 
  if (portion == SIZE_DATA) {
    fclose(t_vs_p);
    return ntp;
  }

  /* 
   * Store in order of increasing sigmatheta.
   */
  for (nn = ntp-1; nn >= 0; nn--) { 

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(t_vs_p,"%lf %lf %lf",&p1,&t1,&dt1);
#else
    fscanf(t_vs_p,"%f %f %f",&p1,&t1,&dt1);
#endif

    /* convert from hPa to Pa */
    var.pdat[ nn] = 100.*p1; 
    var.tdat[ nn] = t1;
    var.dtdat[nn] = dt1;
  }

  fclose(t_vs_p);

  return ntp;
}

/*======================= end of read_t_vs_p() ================================*/

/*======================= read_meridional_plane() =============================*/

/*
 * Reads data from a text file arrayed in a pressure-latitude meridional plane.
 * An example of the format is the file epic/data/Jupiter/C_2H_2.Jupiter.CIRS.
 *
 * The input argument infile is the name of the file, including its path, 
 * and portion is SIZE_DATA or POST_SIZE_DATA.
 *
 * The output arguments np and nlat are the number of pressure and latitude levels
 * in the requested data set.  
 * The output arguments logp and lat are the respective 1D dimension values, where latitude
 * increases with index j and log pressure increases with index k.
 * The output argument value is the 2D array of data values associated with varname,
 * such that value[k][j] refers to pressure level k and latitude j.
 */

void read_meridional_plane(planetspec *planet,
                           char       *infile,
                           int         portion,
                           int        *np,
                           int        *nlat,
                           EPIC_FLOAT *logp,
                           EPIC_FLOAT *lat,
                           EPIC_FLOAT *value)
{
  int
    k,j;
  char
    header[N_STR];
  EPIC_FLOAT
    p;
  FILE
    *input;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="read_meridional_plane";

  input = fopen(infile,"r");
  if (!input) {
    sprintf(Message,"error opening %s",infile);
    epic_error(dbmsname,Message);
  }
  /* Skip over 6-line header: */
  for (k = 0; k < 6; k++) {
    fgets(header,N_STR,input);
  }
  /* input number of data points */
  fscanf(input,"%d %d",np,nlat);
  if (portion == SIZE_DATA) {
    fclose(input);
    return;
  }

  for (j = 0; j < *nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",lat+j);
#else
    fscanf(input,"%f", lat+j);
#endif

  }
  for (k = 0; k < *np; k++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",&p);
#else
    fscanf(input,"%f", &p);
#endif

    /* Convert hPa to Pa, then take log */
    logp[k] = log(p*100.);

    for (j = 0; j < *nlat; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
    fscanf(input,"%lf",value+(j+k*(*nlat)));
#else
    fscanf(input,"%f", value+(j+k*(*nlat)));
#endif

    }
  }

  fclose(input);
  return;
}

/*======================= end of read_meridional_plane() ======================*/
 
/*====================== inquire_radiation_scheme() ==========================*/

/*
 * Prompt the user to choose which radiation scheme to use.
 */

void inquire_radiation_scheme(planetspec *planet)
{
  int
    ii,
    num_radiation_schemes = 4;
  const char
    *radiation_scheme[]
      = {"off",
         "Correlated k",
         "Newtonian",
         "Heating from file"};
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="inquire_radiation_scheme";

  fprintf(stdout,"Choose radiation scheme: 0 => %s \n",radiation_scheme[0]);
  for (ii = 1; ii < num_radiation_schemes; ii++) {
    sprintf(Message,"                        %2d => %s \n",ii,radiation_scheme[ii]);
    if (ii < num_radiation_schemes-1) {
      fprintf(stdout,"%s",Message);
    }
  }

  grid.radiation_index = input_int(Message,grid.radiation_index);
  strcpy(grid.radiation_scheme,radiation_scheme[grid.radiation_index]);

  if (strcmp(grid.radiation_scheme,"Newtonian") == 0) {
    grid.newt_cool_adjust = input_int("Adjust layer average of Newtonian cooling to zero? [1=yes, 0=no]\n",
                                      grid.newt_cool_adjust);
  }
  else if (strcmp(grid.radiation_scheme,"Correlated k") == 0) {
    if (grid.zonal_average_rt != FALSE && 
        grid.zonal_average_rt != TRUE) {
      /* Avoid junk initial prompt. */
      grid.zonal_average_rt = TRUE;
    }
    grid.zonal_average_rt = input_int("Zonally average radiative transfer heating? [1=yes, 0=no]\n",
                                            grid.zonal_average_rt);
  }
}

/*====================== end of inquire_radiation_scheme() ====================*/

/*======================= input_float() =======================================*/

/* 
 * Read in floating-point datum, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

EPIC_FLOAT input_float(char       *prompt,
                       EPIC_FLOAT  def) 
{
  char  
    c,
    buffer[N_STR];
  int   
    len;
  EPIC_FLOAT 
    ans;

#if defined(EPIC_MPI)
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Datatype
       float_type = MPI_DOUBLE;
#  else
     MPI_Datatype
       float_type = MPI_FLOAT;
#  endif
#endif

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%g]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      ans = def;
    }
    else {

#if EPIC_PRECISION == DOUBLE_PRECISION
      sscanf(buffer,"%lf",&ans);
#else
      sscanf(buffer,"%f",&ans);
#endif

    }
  }

#if defined(EPIC_MPI)
   MPI_Bcast(&ans,1,float_type,NODE0,para.comm);
#endif

  return ans;
}

/*====================== end input_float() ====================================*/

/*====================== input_int() ========================================*/

/* 
 * Read in int, or set to default if input is a return ('\n').
 * C.Santori, T.Dowling 
 */

int input_int(char *prompt,
              int   def) 
{
  char  
    c,
    buffer[N_STR];
  int 
    ans,
    len;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%d]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      ans = def;
    }
    else {
      sscanf(buffer,"%d",&ans);
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(&ans,1,MPI_INT,NODE0,para.comm);
#endif

  return ans;
}

/*====================== end input_int() ====================================*/

/*====================== input_string() =====================================*/

/* 
 * Read in a string, or set to default if input is a return ('\n').
 *
 * C.Santori, T.Dowling 
 */

void input_string(char *prompt, 
                  char *def, 
                  char *ans) 
{
  char 
    c,
    buffer[N_STR];
  int  
    len;

  if (IAMNODE == NODE0) {
    fprintf(stdout,"%s[%s]: ",prompt,def);
    for (len = 0; (c = getchar()) != '\n' && len < N_STR; len++) {
      buffer[len]=c;
    }
    buffer[len] = '\0';
    if (len == 0) {
      strcpy(ans,def);
    }
    else {
      strcpy(ans,buffer);
      strcpy(def,buffer);
    }
  }

#if defined(EPIC_MPI)
  MPI_Bcast(ans,N_STR,MPI_CHAR,NODE0,para.comm);
#endif

}

/*====================== end input_string() ==================================*/

/*====================== print_model_description() ===========================*/
/*
 * Print to stdout a select listing of model parameters.
 */

void print_model_description(planetspec *planet)
{
  int
    is,index,
    dt_cfl;
  char
     header[N_STR],
    *ptr;
  double
    max_nu_horizontal[2+1],  /* NOTE: declared as double, not EPIC_FLOAT */
    tmp;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_model_description";

  /*
   * Compute full-strength viscosity coefficients.
   */
  set_max_nu(max_nu_horizontal);

  fprintf(stdout,"\n");
  fprintf(stdout,"                EPIC Version: %4.2f\n",grid.epic_version);

#if EPIC_PRECISION == DOUBLE_PRECISION
  fprintf(stdout,"              Floating point: double precision\n");
#else
  fprintf(stdout,"              Floating point: single precision\n");
#endif

  fprintf(stdout,"                      System: %s \n",planet->name);

  fprintf(stdout,"                    Geometry: %s \n",grid.geometry);
  fprintf(stdout,"         Vertical Coordinate: %s \n",grid.vertical_coordinate);
  if (strcmp(grid.vertical_coordinate,"hybrid") == 0) {
    fprintf(stdout,"           Hybrid parameters: k_sigma = %d, sigma_sigma = %.3f\n",grid.k_sigma,grid.sigma_sigma);

  }
  fprintf(stdout,"                   Lon range: %-.1f to %-.1f deg\n",
         grid.lon[2*1],grid.lon[2*(grid.ni+1)]);
  fprintf(stdout,"                   Lat range: %-.1f to %-.1f deg \n",
         grid.lat[2*(grid.jlo)],grid.lat[2*(grid.nj+1)]);

  fprintf(stdout,"           Momentum timestep: %s\n",
                 grid.uv_timestep_scheme);

  if (strcmp(grid.vertical_coordinate,"isobaric") != 0) {
    fprintf(stdout,"              Mass advection: %s \n",
                   var.h.advection_scheme);
  }

  if (strcmp(grid.vertical_coordinate,"isentropic") != 0) {
    fprintf(stdout,"             Theta advection: %s \n",
                   var.theta.advection_scheme);
  }

  if (strcmp(grid.radiation_scheme,"Newtonian") == 0) {
    if (grid.newt_cool_adjust) {
      fprintf(stdout,"            Radiation scheme: Newtonian, layer average adjusted to zero\n");
    }
    else {
      fprintf(stdout,"            Radiation scheme: Newtonian, no adjustment to layer average\n");
    }
  }
  else {
    fprintf(stdout,"            Radiation scheme: %s\n",grid.radiation_scheme);
  }

  if (grid.cloud_microphysics == OFF) {
    fprintf(stdout,"          Cloud microphysics: off\n");
  }
  else if (grid.cloud_microphysics == ACTIVE) {
    fprintf(stdout,"          Cloud microphysics: active\n");
  }
  else if (grid.cloud_microphysics == PASSIVE) {
    fprintf(stdout,"          Cloud microphysics: passive\n");
  }
  else if (grid.cloud_microphysics == STEADY) {
    fprintf(stdout,"          Cloud microphysics: steady\n");
  }
  else {
    sprintf(Message,"grid.cloud_microphysics=%d not recognized",grid.cloud_microphysics);
    epic_error(dbmsname,Message);
  }

  fprintf(stdout,"           Turbulence scheme: %s \n",grid.turbulence_scheme);

  strftime(Message,N_STR,"days since %Y-%m-%d %H:%M:%S (UTC)",gmtime(&var.start_time));
  fprintf(stdout,"                        Time: %g %s \n",TIME/86400.,Message);


  season_string(L_s,Message);
  fprintf(stdout,"                         L_s: %.1f deg %s\n",L_s,Message);

  /*
   * This CFL estimate is not currently a reliable guide for a 
   * numerically stable timestep, so we do not print it.
   */
  /****
  dt_cfl = cfl_dt(planet);
  fprintf(stdout,"                    Timestep: dt = %d s, CFL dt = %d s\n",
                  grid.dt,dt_cfl);
   ****/
  fprintf(stdout,"                    Timestep: dt = %d s\n",grid.dt);


  if (grid.ni == 1) {
    if (grid.nj == grid.jlo) {
      /* 1D model. */
      fprintf(stdout,"                        Size: nk = %d, 1D \n",
              grid.nk);
    }
    else {
      /* 2D model. */
      fprintf(stdout,"                        Size: nk = %d, nj = %d, 2D \n",
              grid.nk,grid.nj);
    }
  }
  else {
    if (grid.nj == grid.jlo) {
      /* 2D model. */
      fprintf(stdout,"                        Size: nk = %d, ni = %d, 2D \n",
              grid.nk,grid.ni);
    }
    else {
      /* 3D model. */
      fprintf(stdout,"                        Size: nk = %d, nj = %d, ni = %d \n",
              grid.nk,grid.nj,grid.ni);
    }
  }
  /* 
   * Print species list: 
   */
  for (is = FIRST_SPECIES; is <= LAST_SPECIES; is++) {
    if (var.species[is].on) {
      fprintf(stdout,"%28s: included \n",var.species[is].info[0].name);
    }
  }
  /*
   * Print chemical information:
   */
  if (var.fpara.on) {
    fprintf(stdout,"     ortho-para rate scaling: %4.2f\n",var.fpara_rate_scaling);
  }

  if (grid.nudiv_nondim > 0.) {
    fprintf(stdout,"          Divergence damping: %4.2f (%9.3e m^2/s)\n",
                    grid.nudiv_nondim,grid.nudiv_nondim*max_nu_horizontal[2]);
  }
  else {
    fprintf(stdout,"          Divergence damping: off \n");
  }

  sprintf(header,"              Hyperviscosity: ");
  if (grid.nu_order < 4) {
    strcat(header,"off");
  }
  else {
    index = grid.nu_order;
    sprintf(Message,"%dth order, %4.2f (%9.3e m^%d/s)",index,grid.nu_nondim,grid.nu_hyper,grid.nu_order);
    strcat(header,Message);
  }
  strcat(header,"\n");
  fprintf(stdout,"%s",header);
    
  if (grid.k_sponge < 0) {
    fprintf(stdout,"                    k_sponge: off\n");
  }
  else {
    fprintf(stdout,"                    k_sponge: %d\n",grid.k_sponge);
  }

  if (grid.j_sponge < 0) {
    fprintf(stdout,"                    j_sponge: off\n");
  }
  else {
    fprintf(stdout,"                    j_sponge: %d\n",grid.j_sponge);
  }

  return;
}

/*====================== end print_model_description() =======================*/


/*====================== print_zonal_info() ==================================*/

void print_zonal_info(planetspec *planet)
{
  int
    K,J,I;
  EPIC_FLOAT
    zeta,zf,zfy,fy,
    *buffji;
  FILE
    *u_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_zonal_info";

#if defined(EPIC_MPI)
  /* NOTE: not set up for MPI. */
  sprintf(Message,"not set up for MPI");
  epic_error(dbmsname,Message);
#endif

  /*
   * Allocate memory.
   */
  buffji = fvector(0,Nelem2d-1,dbmsname);

  K = grid.nk;
  I = grid.ilo;

  /*
   * Calculate absolute vorticity on pv-grid.
   */
  vorticity(ON_SIGMATHETA,ABSOLUTE,2*K+1,
            var.u.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            var.v.value+(K-Kshift)*Nelem2d+grid.it_uv*Nelem3d,
            NULL,
            buffji);

  u_dat = fopen("zonal_wind.dat","w");
  fprintf(u_dat," %s,  nk, nj, ni = %2d, %2d, %2d \n",
                 planet->name,grid.nk,grid.nj,grid.ni);
  fprintf(u_dat," zonal wind data at K,I = %2d, %2d \n",K,I);
  fprintf(u_dat,"      f,zeta,z+f units: 1.e-4 s^-1\n");
  fprintf(u_dat," df/dy,d(z+f)/dy units: 1.e-12 m^-1 s^-1 \n\n");
  fprintf(u_dat," lat(deg)  u(m/s)      f       zeta"
                "       zeta+f     df/dy     d(zeta+f)/dy \n");

  for (J = JHI; J >= JLO; J--) {
    /* Average vorticity onto u-grid. */
    zeta = .5*(BUFFJI(J,I)-grid.f[2*J]+BUFFJI(J+1,I)-grid.f[2*(J+1)]);
    zf   = .5*(BUFFJI(J,I)+BUFFJI(J+1,I));
    zfy  = grid.n[2*K+1][2*J+1]*(BUFFJI(J+1,I)-BUFFJI(J,I));
    fy   = grid.n[2*K+1][2*J+1]*(grid.f[2*J+2]-grid.f[2*J]);
    /* change units on zeta terms */
    /* change units on beta=df/dy terms */
    zeta *= 1.e+4;
    zf   *= 1.e+4;
    zfy  *= 1.e+12;
    fy   *= 1.e+12;
    fprintf(u_dat," %5.1f   %7.2f  %10.3e %10.3e %10.3e %10.3e %10.3e \n",
                  grid.lat[2*J+1],U(grid.it_uv,grid.nk,J,I),
                  grid.f[2*J+1]*1.e+4,zeta,zf,fy,zfy);
  }

  fclose(u_dat);

  /*
   * Free allocated memory.
   */
  free_fvector(buffji,0,Nelem2d-1,dbmsname);

  return;
}

/*====================== end of print_zonal_info() ===========================*/

/*====================== print_vertical_column() =============================*/

void print_vertical_column(planetspec *planet,
                           int         J,
                           int         I,
                           char       *filename)
{
  int
    K;
  EPIC_FLOAT
    pressure,theta,temperature,
    brunt2;
  FILE
    *vert_dat;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="print_vertical_column";

  if (IAMNODE != NODE0) {
    return;
  }

  /*
   * Check validity of J,I.
   */
  if (J < JLO || J > JHI) {
    sprintf(Message,"J=%d not in range [%d,%d]",J,JLO,JHI);
    epic_error(dbmsname,Message);
  }
  if (I < ILO || I > IHI) {
    sprintf(Message,"I=%d not in range [%d,%d]",I,ILO,IHI);
    epic_error(dbmsname,Message);
  }

  vert_dat = fopen(filename,"w");
  fprintf(stdout,"\n Column data at lat = %.1f, lon = %.1f: \n",
                 grid.lat[2*J+1],grid.lon[2*I+1]);
  fprintf(vert_dat,"  Vertical profile for lat=%.1f lon=%.1f \n",
                    grid.lat[2*J+1],grid.lon[2*I+1]);

  switch(grid.coord_type) {
    case COORD_HYBRID:
      fprintf(stdout,"\n     K  sigmatheta      press[hPa]     theta[K]  N2[1/s^2]  \n");
      fprintf(vert_dat,"   K  sigmatheta  press[hPa]     temp[K]  theta[K]  N2[1/s^2]  U[m/s]  re[km]   rp[km]\n");

      K = 0;
      fprintf(stdout,"-- %4.1f -- %7.2f -- %12.7g -- %7.2f ----------- \n",
             (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
      for (K = 1; K <= grid.nk; K++) {
        brunt2 = get_brunt2(planet,2*K,J,I);
        fprintf(stdout,"   %4.1f    %7.2f                    %7.2f  %9.6f \n", 
                       (EPIC_FLOAT)K,grid.sigmatheta[2*K],THETA2(K,J,I),brunt2);

        fprintf(vert_dat," %4.1f  %6.3f %13.6e  %9.1f %9.1f  %9.6f %7.2f %8.1f %8.1f\n",
                         (EPIC_FLOAT)K,grid.sigmatheta[2*K],P2(K,J,I)/100.,T2(K,J,I),THETA2(K,J,I),brunt2,U(grid.it_uv,K,J,I),
                          grid.re[K]/1000.,grid.rp[K]/1000.);
        if (K <= grid.k_sponge) {
          /*
           * The squiggles signify sponge layers.
           */
          fprintf(stdout,"~~ %4.1f ~~ %7.2f ~~ %12.7g ~~ %7.2f ~~~~~~~~~~~ \n",
                 (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
        }
        else {
          fprintf(stdout,"-- %4.1f -- %7.2f -- %12.7g -- %7.2f ----------- \n",
                 (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
        }
      }
    break;
    default:
      fprintf(stdout,"\n     K       sgth[K]      press[hPa]     theta[K]  N2[1/s^2]  \n");
      fprintf(vert_dat,"   K     sgth[K]    press[hPa]    temp[K]  theta[K]  N2[1/s^2]  U[m/s]   re[km]   rp[km]\n");

      K = 0;
      fprintf(stdout,"-- %4.1f -- %9.1f -- %12.7g -- %9.1f ----------- \n",
             (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
      for (K = 1; K <= grid.nk; K++) {
        brunt2 = get_brunt2(planet,2*K,J,I);
        fprintf(stdout,"   %4.1f    %9.1f                    %9.1f  %9.6f \n", 
                       (EPIC_FLOAT)K,grid.sigmatheta[2*K],THETA2(K,J,I),brunt2);

        fprintf(vert_dat," %4.1f  %9.1f %13.6e  %9.1f %9.1f  %9.6f %7.2f %8.1f %8.1f\n",
                         (EPIC_FLOAT)K,grid.sigmatheta[2*K],P2(K,J,I)/100.,T2(K,J,I),THETA2(K,J,I),brunt2,U(grid.it_uv,K,J,I),
                          grid.re[K]/1000.,grid.rp[K]/1000.);
        if (K <= grid.k_sponge) {
          /*
           * The squiggles signify sponge layers.
           */
          fprintf(stdout,"~~ %4.1f ~~ %9.1f ~~ %12.7g ~~ %9.1f ~~~~~~~~~~~ \n",
                 (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
        }
        else {
          fprintf(stdout,"-- %4.1f -- %9.1f -- %12.7g -- %9.1f ----------- \n",
                 (EPIC_FLOAT)K+0.5,grid.sigmatheta[2*K+1],P3(K,J,I)/100.,THETA(K,J,I));
        }
      }
    break;
  }

  fclose(vert_dat);

  return;
}

/*====================== end of print_vertical_column() ======================*/

/*====================== node0_barrier_print() ===============================*/

/*
 * Aaron Herrnstein.
 *
 * Routine prints a message using Processor 0.  MPI_Barrier() calls are placed
 * before and after the print statement.  This can be beneficial when tracking
 * down MPI bugs.
 */
void node0_barrier_print(char *calling_function,
                         char *Message, 
                         int  display_time)
{
  time_t rawtime;

#if defined(EPIC_MPI)
  MPI_Barrier(para.comm); 
#endif

  if (IAMNODE == NODE0) {
    /*
     * Print Message using processor 0
     */
    if (display_time) {
      time( &rawtime );
      fprintf(stdout, "%s(),  %s  at %s\n", calling_function, Message, ctime(&rawtime) );
    }
    else {
      fprintf(stdout, "%s(),  %s\n", calling_function, Message );
    }
    fflush(stdout);
  }

#if defined(EPIC_MPI)
  /* The barrier looks redundant here, but is necessary for proper io. */
  MPI_Barrier(para.comm);
#endif

}

/*====================== end of node0_barrier_print() ========================*/

/*====================== scdswap() ===========================================*/

/*
 *  Byte swapping for 8, 4 and 2 byte quantities.  
 */
#include <stdio.h>

void scdswap(char *arr, 
             int   nbytes, 
             int   cnt)
{
  char 
    buf[4];
  register int 
    nb=nbytes;
  register char 
   *parr, *pend;

  pend = arr+nb*cnt;

  switch(nb) {
    case 1:
    break;
    case 2:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        parr[0] = parr[1];
        parr[1] = buf[0];
      }
    break;
    case 4:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        buf[1]  = parr[1];
        parr[0] = parr[3];
        parr[1] = parr[2];
        parr[2] = buf[1];
        parr[3] = buf[0];
      }
    break;
    case 8:  
      for (parr=arr; parr<pend; parr+=nb) {
        buf[0]  = parr[0];
        buf[1]  = parr[1];
        buf[2]  = parr[2];
        buf[3]  = parr[3];
        parr[0] = parr[7];
        parr[1] = parr[6];
        parr[2] = parr[5];
        parr[3] = parr[4];
        parr[4] = buf[3];
        parr[5] = buf[2];
        parr[6] = buf[1];
        parr[7] = buf[0];
      }
    break;
    default: 
      fprintf(stderr," Bad length to scdswap()\n");
      exit(99);
    break;
  } 

  return;
} 

/*======================= end of scdswap() ===================================*/

/*======================= epic_error() =======================================*/
/*
 * Prints machine name, calling node rank, timestep, calling function name, 
 * and Message to stderr, then aborts.
 */
void epic_error(char *calling_function,
                char *Message)
{
  fprintf(stderr,"\n** EPIC Error: node=%d, timestep=%lu, %s(): %s\n",
                 IAMNODE,grid.itime,calling_function,Message);
  fflush(stderr);

#if defined(EPIC_MPI)
  MPI_Abort(MPI_COMM_WORLD,1);
#endif

  exit(1);
}

/*======================= end of epic_error() ================================*/

/*======================= epic_warning() =====================================*/
/*
 * Prints machine name, calling node rank, calling function name, 
 * and Message to stderr (does not exit).
 */
void epic_warning(char *calling_function,
                  char *Message)
{
  fprintf(stderr,"EPIC Warning: node=%d, timestep=%lu, %s(): %s\n",
                 IAMNODE,grid.itime,calling_function,Message);
  fflush(stderr);

  return;
}

/*======================= end of epic_warning() ==============================*/

/*======================= declare_copyright() ================================*/

void declare_copyright(void)
{
  
  static int
    oneshot=FALSE;

  if (!oneshot) {
    fprintf(stderr,"\n");
    fprintf(stderr," EPIC Model, Copyright (C) 1998-2023 Timothy E. Dowling \n");                                                                                         
    fprintf(stderr," This program is free software; you can redistribute it and/or \n");  
    fprintf(stderr," modify it under the terms of the GNU General Public License.  \n");    
    fprintf(stderr," This program is distributed WITHOUT ANY WARRANTY.             \n"); 
    fprintf(stderr,"\n");

    oneshot = TRUE;
  } 
                                                                 
  return;
}

/*======================= end of declare_copyright() =========================*/

/* * * * * * * * * * * * *  end of epic_funcs_io.c  * * * * * * * * * * * * * */




