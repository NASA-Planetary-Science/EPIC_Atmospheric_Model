/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998 Joseph Matarese                              *
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

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpg.h"

void  index_row (int row, int *length, int *start, int dim, int *index) {
  int    i, hyperplane;

  for (hyperplane=1, i=0; i<dim; i++) hyperplane *= length[i];
  row *= length[0];
  index[0] = start[0];

  for (i=dim-1; i>0; i--) {
    hyperplane /= length[i];
    index[i] = row/hyperplane;
    row = row%hyperplane;
  }

  return;
}

void  offset_index (int *end, int *start, int *pad, int *index, int dim,
		    int *offset) {
  int    i, hyperplane;

  for (hyperplane=1, *offset=0, i=0; i<dim; i++) {
    *offset += (index[i]-start[i]+pad[i])*hyperplane;
    hyperplane *= (end[i]-start[i]+1+2*pad[i]);
  }

  return;
}

void  map_error (const char *s) {

  (void)fprintf (stderr, "%s: %s\n", __FILE__, s);
  MPI_Abort (MPI_COMM_WORLD, 0);
}
