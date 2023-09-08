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

/* * * * * * * * * * * * * s_epic_functions.c  * * * * * * * * * * * * * * * * 
 *                                                                           *
 *       Single-processor specific functions.                                *
 *                                                                           *
 *       T. Dowling, R. LeBeau                                               *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*==========================  bc_lateral() ==================================*/

/*
 * Apply lateral boundary conditions.
 */

void bc_lateral(EPIC_FLOAT *pt,
                int         dim) {
  int  
    K,J,I;
  EPIC_FLOAT
    *buff3d,
    *buffji;  
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="bc_lateral";

  if (dim == TWODIM) {
    buffji = pt;
    /* Zonal periodicity */
    for (J = 0; J <= grid.nj+1; J++) {
      BUFFJI(J,        0) = BUFFJI(J,grid.ni);
      BUFFJI(J,grid.ni+1) = BUFFJI(J,      1);
    }

    if (grid.jlo == 1) {
      /* Doubly periodic. */
      for (I = 0; I <= grid.ni+1; I++) {
        BUFFJI(        0,I) = BUFFJI(grid.nj,I);
        BUFFJI(grid.nj+1,I) = BUFFJI(      1,I);
      }
    }
  }
  else if (dim == THREEDIM) {
    buff3d = pt;
    for (K = KLOPAD; K <= KHIPAD; K++) {
      /* Zonal periodicity */
      for (J = 0; J <= grid.nj+1; J++) {
        BUFF3D(K,J,        0) = BUFF3D(K,J,grid.ni);
        BUFF3D(K,J,grid.ni+1) = BUFF3D(K,J,      1);
      }

      if (grid.jlo == 1) {
        /* Doubly periodic. */
        for (I = 0; I <= grid.ni+1; I++) {
          BUFF3D(K,        0,I) = BUFF3D(K,grid.nj,I);
          BUFF3D(K,grid.nj+1,I) = BUFF3D(K,      1,I);
        }
      }
    }
  }
  else {
    sprintf(Message,"dim=%d unrecognized",dim);
    epic_error(dbmsname,Message);
  }

  return;
}

/*========================== end of bc_lateral() ============================*/

/* * * * * * * * * * * * * end of s_epic_functions.c * * * * * * * * * * * * */ 
