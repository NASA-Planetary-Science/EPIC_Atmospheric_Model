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

#ifndef _mpg_h_
#define _mpg_h_

#include <stdio.h>

int MPE_Decomp1d(int n, 
                 int size, 
                 int rank, 
                 int *s, 
                 int *e);

extern void  MPG_Cart_decomp (
/* INPUT */   MPI_Comm comm, int dim, int *dimlen, int *dims, int *periods,
/* INPUT */   int *pad,
/* OUTPUT */  int *start, int *end, int *stride, MPI_Comm *comm_cart);

extern void  MPG_Cart_varcreate (
/* INPUT */   int dim, int *start, int *end, int *pad, MPI_Datatype etype,
/* OUTPUT */  char **var, int *offset);

extern void  MPG_Cart_varfree (
/* INPUT */   int dim, int *start, int *end, int *pad, MPI_Datatype etype,
/* INPUT */   char *var);

extern void  MPG_Cart_read (
/* INPUT */   FILE *stream, MPI_Comm comm, int dim, int *length, int *pad,
/* INPUT */   MPI_Datatype etype,
/* OUTPUT */  void *local);

extern void  MPG_Cart_write (
/* INPUT */   FILE *stream, MPI_Comm comm, int dim, int *length, int *pad,
/* INPUT */   MPI_Datatype etype,
/* INPUT */  void *local);

extern void  MPG_Cart_edgeexch (
/* INPUT */   MPI_Comm comm_cart, int dim, int *length, int *pad,
/* INPUT */   MPI_Datatype etype,
/* In/OUT */  void *input);

extern void  MPG_Cart_edgesend (
/* INPUT */   MPI_Comm comm_cart, int dim, int *length, int *pad,
/* INPUT */   MPI_Datatype etype, void *local,
/* OUTPUT */  MPI_Request *request);

extern void  MPG_Cart_edgerecv (
/* INPUT */   MPI_Comm comm_cart, int dim, int *length, int *pad,
/* INPUT */   MPI_Datatype etype,
/* OUTPUT */  void *local,
/* INPUT */   MPI_Request *request,
/* OUTPUT */  MPI_Status *stat);

#endif /* _mpg_h_ */
