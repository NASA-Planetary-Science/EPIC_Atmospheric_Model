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
#include "qa.h"

static char  
  **ptr_list=NULL;
static int     
  ptr_list_len=0;

void MPG_Cart_varcreate(int dim,int *start,int *end,int *pad,
                        MPI_Datatype etype,char **var,int *offset) 
{
  int       
    i, size;
  int  
    esize;
  char    
    **old_ptr_list;

  for (*offset=0, size=1, i=0; i<dim; i++) {
    *offset += pad[i]*size;
    size *= (end[i] - start[i] + 1 + 2*pad[i]);
  }
  MPI_Type_size (etype, &esize);
  size *= (int)esize;
  if ((*var = (char *)malloc(size)) == NULL) {
    fprintf (stderr, "Cannot allocate %d bytes for Cartesian variable\n",
	     size);
    MPI_Abort (MPI_COMM_WORLD, 0);
  }

  old_ptr_list = ptr_list;
  ptr_list_len++;
  if ((ptr_list = (char **)malloc(ptr_list_len*sizeof(char *))) == NULL) {
    fprintf (stderr, "Cannot allocate space\n");
    MPI_Abort (MPI_COMM_WORLD, 0);
  }
  (void)memcpy (ptr_list, old_ptr_list, (ptr_list_len-1)*sizeof(char *));
  ptr_list[ptr_list_len-1] = (char *)(*var);
  (void)free((char *)old_ptr_list);

  debug ('m', "Allocated %d bytes ", size);
  debug ('m', "at address %#x\n", *var);

  return;
}

void  MPG_Cart_varfree (int dim, int *start, int *end, int *pad,
			MPI_Datatype etype, char *var) {
  int       
    i, size;
  int  
    esize;

  for (size=1, i=0; i<dim; i++) {
    size *= (end[i] - start[i] + 1 + 2*pad[i]);
  }
  MPI_Type_size (etype, &esize);
  size *= (int)esize;

  debug ('m', "Freeing %d bytes ", size);
  debug ('m', "at address %#x...", var);

  for (i=0; i<ptr_list_len & (char *)var != ptr_list[i]; i++);
  if (i == ptr_list_len) {
    fprintf (stderr, "Cannot find record of allocation at address %p\n",
	     (void *)var);
    MPI_Abort (MPI_COMM_WORLD, 0);
  }
  
  (void)free (ptr_list[i]);
  ptr_list[i] = NULL;

  debug ('m', "%s\n", "done");

  return;
}
