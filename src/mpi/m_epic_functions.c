/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2019 Timothy E. Dowling                      *
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

/* * * * * * * * * * * * m_epic_functions() * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*====================== mpispec_init() ======================================*/ 

/*
 *  Sets up mpi bookkeeping.  
 *  The mpispec structure para is defined globally in epic.h.
 */

void mpispec_init(void)
{
  int   
    i,dim;
  MPI_Datatype
    float_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="mpispec_init";

  if (EPIC_PRECISION == DOUBLE_PRECISION) {
    float_type = MPI_DOUBLE;
  }
  else {
    float_type = MPI_FLOAT;
  }

 /* 
  * Define complex data type: 
  */
  MPI_Type_contiguous(2,float_type,&EPIC_MPI_COMPLEX);
  MPI_Type_commit(&EPIC_MPI_COMPLEX);

  /* 
   * Determine grid decompositions.
   *
   * Dimensions: 0 = zonal
   *             1 = meridional
   *             2 = vertical
   */
  para.nstart[0] = 1;
  /* set para.nstart[1] below */
  para.nstart[2] = 1;

  para.nend[0] = grid.ni;
  para.nend[1] = grid.nj;
  para.nend[2] = grid.nk;

  for (i = 0; i < TOPDIM; i++) {
    para.npad[i] = grid.pad[i];
    para.wrap[i] = grid.wrap[i];
  }

  /* 
   * A negative value flags automatic balancing.
   */
  para.nprocs[0] = -1;
  para.nprocs[1] = -1;
  para.nprocs[2] =  1;

  if (strcmp(grid.geometry,"globe") == 0) {
    /* 
     * J = 0 is an interior row (except for pv or v, 
     * but it is harmless to treat them the same). 
     */
    para.nstart[ 1] = 0;
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      para.nstart[1] = 1;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      para.nstart[1] = 0;
    }
    else {
      fprintf(stderr,"mpispec_init: Unrecognized f_plane_map \n");
      exit(1);
    }
  }
  else {
    fprintf(stderr,"m_epic_functions: Unrecognized geometry \n");
    exit(1);
  }

  for (dim = 0; dim < TOPDIM; dim++) {
    para.dimlen[dim] = para.nend[dim]-para.nstart[dim]+1;
  }

  /* Make 2D cartesian-decomposion communicator. */
  MPG_Cart_decomp(para.comm,2,para.dimlen,
                  para.nprocs,para.wrap,para.npad,
                  para.mylo,para.myhi,para.arraydims,&para.comm_ij);

  /* Make 3D cartesian-decomposition communicator. */
  MPG_Cart_decomp(para.comm,3,para.dimlen,
                  para.nprocs,para.wrap,para.npad,
                  para.mylo,para.myhi,para.arraydims,&para.comm_ijk);

  MPI_Comm_rank(para.comm_ijk,&para.iamnode);
  MPI_Comm_size(para.comm_ijk,&para.nproc);
  grid.we_num_nodes = para.nproc;
  para.ndim         = NINT(log((EPIC_FLOAT)para.nproc)/log(2.));

  for (dim = 0; dim < TOPDIM; dim++) {
    /* Shift array endpoints */
    para.mylo[dim] += para.nstart[dim];
    para.myhi[dim] += para.nstart[dim];
  }

  para.nelem2d = para.arraydims[0]*para.arraydims[1];
  para.nelem3d = para.nelem2d*para.arraydims[2];

  /*
   *  The parameters jfirst and jlast are used to handle the indexing of 
   *  boundaries in the latitude direction for C-grid staggering.
   */
  if (para.wrap[1]) {
    /* Periodic in J direction */
    para.jlast  = IMIN(grid.nj, JHI);
    para.jfirst = IMAX(grid.jlo,JLO);
  }
  else {
    /* Not periodic in J direction */
    para.jlast  = IMIN(grid.nj-1, JHI);
    para.jfirst = IMAX(grid.jlo+1,JLO);
  }

  if (strcmp(grid.geometry,"globe") == 0) {
    if (JLO == para.nstart[1]) {
      if (fcmp(grid.globe_latbot,-90.) == 0) {
        para.is_spole = TRUE;
      }
      else {
        para.is_spole = FALSE;
      }
    }
    else {
      para.is_spole = FALSE;
    }

    if (JHI == para.nend[1]) {
      if (fcmp(grid.globe_lattop,90.) == 0) {
        para.is_npole = TRUE;
      }
      else {
        para.is_npole = FALSE;
      }
    }
    else {
      para.is_npole = FALSE;
    }
  }
  else if (strcmp(grid.geometry,"f-plane") == 0) {
    para.is_spole = FALSE;
    if (strcmp(grid.f_plane_map,"cartesian") == 0) {
      para.is_npole = FALSE;
    }
    else if (strcmp(grid.f_plane_map,"polar") == 0) {
      if (JHI == para.nend[1]) {
        para.is_npole = TRUE;
      }
      else {
        para.is_npole = FALSE;
      }
    }
  }

  return;
}

/*====================== end of mpispec_init() ====================================*/

/*====================== bc_lateral() =============================================*/

void bc_lateral(EPIC_FLOAT *pt,
                int         dim)
{
  static int
    initialized = FALSE;
  static MPI_Datatype
    float_type;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="bc_lateral";

  if (!initialized) {
    if (EPIC_PRECISION == DOUBLE_PRECISION) {
      float_type = MPI_DOUBLE;
    }
    else if (EPIC_PRECISION == SINGLE_PRECISION) {
      float_type = MPI_FLOAT;
    }
    else {
      sprintf(Message,"unrecognized EPIC_PRECISION=%d",EPIC_PRECISION);
      epic_error(dbmsname,Message);
    }

    initialized = TRUE;

  } /* end initialization */

  if (dim == TWODIM) {
    MPG_Cart_edgeexch(para.comm_ij,dim,para.dimlen,para.npad,float_type,pt);
  }
  else if (dim == THREEDIM) {
    MPG_Cart_edgeexch(para.comm_ijk,dim,para.dimlen,para.npad,float_type,pt);
  }
  else {
    sprintf(Message,"unrecognized dim=%d",dim);
    epic_error(dbmsname,Message);
  }

  return;
}

/*====================== end of bc_lateral() ======================================*/

/* * * * * * * * * * * * end of m_epic_functions() * * * * * * * * * * * * * * * * */
