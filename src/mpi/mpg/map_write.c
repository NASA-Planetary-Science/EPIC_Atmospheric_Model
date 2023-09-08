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

/*---------------------------------------------------------------------------*
 * Logic:
 *
 * 1. Node 0 reads row.
 * 2. Node 0 broadcasts row.
 * 3. Node i!=0 copy portion of row to local array, if applicable.
 * 
 *
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpg.h"
#include "qa.h"

#define M_NEED_ROW 240

extern void index_row(int row, 
                      int *length, 
                      int *start, 
                      int dim, 
                      int *index);
extern void offset_index(int *end, 
                         int *start, 
                         int *pad, 
                         int *index, 
                         int dim,
                         int *offset);
extern void map_error(const char *s);

void MPG_Cart_write(FILE         *stream,    /* stream from which to read data        */
                    MPI_Comm     comm_cart,  /* communicator for cartesian grid       */
                    int          ndim,       /* dimensions of Cartesian grid          */
                    int          *length,    /* number of elements along each axis    */
                    int          *pad,       /*                                       */
                    MPI_Datatype etype,      /* data type of grid elements to be read */
                    void         *local)     /* local array for storing data          */
{
  int                           
    inode,nnode,    /* ith node, total nodes          */                          
    idim,           /* ith dimension                  */                      
    tsize,          /* stores data type size          */                          
    esize,          /* stores data type size as int   */                          
    irow,nrow,      /* ith row, number of grid rows   */                        
    row_node,
    row_nodes,
    tag,                        
    ne;             /* number of elements in grid row */
  MPI_Status                    
    stat;
  char                         
    *buffer;        /* buffer to hold row             */
  int                          
    *periods, 
    *dims, 
    *coords,       /* dummy ptrs                      */                       
    *index,        /* index of starting row element   */                         
    *start, 
    *end,          /* start and end indices           */                       
    my_row,        /* row on my node                  */                         
    offset,i,
    count;
  double                        
    t1,t2;

  MPI_Comm_rank(comm_cart,&inode);
  MPI_Comm_size(comm_cart,&nnode);

  MPI_Type_size(etype,&tsize);
  esize = (int)tsize;

  /* Allocate memory: */
  if ((start   = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (end     = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (index   = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (dims    = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (periods = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (coords  = (int *)malloc(ndim*sizeof(int))) == NULL) {
    map_error ("Cannot allocate space");
  }

  ne = length[0];
  MPI_Cart_get(comm_cart,ndim,dims,periods,coords);
  nrow = 1;
  for (idim = 0; idim < ndim; idim++) {
    if (idim) {
      nrow *= length[idim];
    }
    MPE_Decomp1d(length[idim],dims[idim],coords[idim],
		 &start[idim],&end[idim]);
    start[idim]--;  
    end[  idim]--;
  }
  row_nodes = dims[0];  
  row_node  = coords[0];

/******
  fprintf(stderr, "node %d has row_node %d of %d\n",
	  inode,row_node,row_nodes);
*******/

/*- Allocate row buffer -----------------------------------------------------*/

  if ((buffer = malloc(ne*esize)) == NULL) {
    map_error ("Cannot allocate space");
  }
  
  MPI_Barrier(comm_cart);
  t1 = MPI_Wtime();

  for (irow = 0; irow < nrow; irow++) {

/*- Node i determines whether it has row ------------------------------------*/

    index_row(irow,length,start,ndim,index);

    my_row = 1;
    for (idim = 1; idim < ndim; idim++) {
      if (index[idim] < start[idim] || index[idim] > end[idim]) {
	my_row = 0;  
        break;
      }
    }

/*- Node i identifies portion of local memory belonging to row. -------------*/

    count = offset = 0;
    if (my_row) {
      offset_index (end,start,pad,index,ndim,&offset);
      count = end[0]-start[0]+1;
    }

/*- On node !=0 if my_row ---------------------------------------------------*/

    if (inode && my_row) {
      /*- Send portion size to node 0 with message tag = row position -------*/
      tag = row_node;
/******
      fprintf(stderr,"node %d sending tag %d to node 0\n", inode, tag);
*******/
      MPI_Send(&count,1,MPI_INT,0,tag,comm_cart);
      /*- Node !=0 sends portion to node 0 with tag = nnodes + row position -*/
      tag += row_nodes;
      MPI_Send((char *)local+offset*esize,count,etype,0,tag,comm_cart);
    }

/*- On node 0, for each of nnodes messages containing portion sizes ---------*/

    if (inode == 0) {
      if (my_row && fwrite ((char *)local + offset*esize, esize, count, stream) != count) {
	map_error ("Cannot write expected number of elements");
      }
      for (i = my_row; i < row_nodes; i++) {
	/*- Read portion size message ---------------------------------------*/
	tag = i;
	MPI_Recv(&count,1,MPI_INT,MPI_ANY_SOURCE,tag,comm_cart,&stat);
	/*- Read portion message --------------------------------------------*/
	tag += row_nodes;
	MPI_Recv(buffer,count,etype,MPI_ANY_SOURCE,tag,comm_cart,&stat);
	/*- Write portion to file -------------------------------------------*/
	if (fwrite(buffer,esize,count,stream) != count) {
	  map_error ("Cannot write expected number of elements");
        }
      }
    }

    MPI_Barrier(comm_cart);

    if (inode == 0) {
      debug('i',"%d of ",irow+1);
      debug('i',"%d rows written\r",nrow);
    }
  }

  if (inode == 0) {
    debug('i',"%s","\n");
  }

  MPI_Barrier(comm_cart);
  t2 = MPI_Wtime();

  if (inode == 0) {
    debug('t',"MPI_Cart_write:  elapsed time = %g\n",t2-t1);
  }

  (void)free((char *)buffer);
  (void)free((char *)start);
  (void)free((char *)end);
  (void)free((char *)index);
  (void)free((char *)dims);
  (void)free((char *)periods);
  (void)free((char *)coords);

  return;
}
