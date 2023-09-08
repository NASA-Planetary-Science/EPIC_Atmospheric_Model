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
extern void  map_error(const char *s);

void MPG_Cart_rslab(FILE         *stream,    /* stream from which to read data        */
                    MPI_Comm     comm_cart,  /* communicator for cartesian grid       */
                    int          ndim,       /* dimensions of Cartesian grid          */
                    int          *slab_hi,   /* index of maximum value to write       */
                    int          *slab_lo,   /* index of minimum value to write       */
                    int          *length,    /* number of elements along each axis    */
                    int          *pad,       /*                                       */
                    MPI_Datatype etype,      /* data type of grid elements to be read */
                    void         *local)     /* local array for storing data */
{
  int                           
    inode,nnode,   /* ith node, total nodes          */                           
    idim,          /* ith dimension                  */                      
    tsize,         /* stores data type size          */                        
    esize,         /* stores data type size as int   */                      
    irow,nrow,     /* ith row, number of grid rows   */                                            
    ne;            /* number of elements in grid row */
  char                         
    *buffer;       /* buffer to hold row             */
  int                          
    *periods, 
    *dims, 
    *coords,       /* dummy ptrs                     */                         
    *index,        /* index of starting row element  */                          
    *start, 
    *end,          /* start and end indices          */                         
    my_row,        /* row on my node                 */                         
    my_slice, 
    slice_offset, 
    slice_row, 
    count,
    offset;
  double                        
    t1,t2;
  int
    orow=0,                    
    bytes_read=0;

  MPI_Comm_rank(comm_cart,&inode);
  MPI_Comm_size(comm_cart,&nnode);

  MPI_Type_size(etype,&tsize);
  esize = (int)tsize;

  /* Allocate memory */
  if ((start   = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (end     = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (index   = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (dims    = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (periods = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (coords  = (int *)malloc(ndim*sizeof(int))) == NULL) {
    map_error ("Cannot allocate space");
  }

  ne = slab_hi[0]-slab_lo[0]+1;
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

  MPI_Barrier(comm_cart);
  t1 = MPI_Wtime();

/*- Allocate row buffer -----------------------------------------------------*/

  if ((buffer = malloc(ne*esize)) == NULL) {
    map_error("Cannot allocate space");
  }

  for (irow = 0; irow < nrow; irow++) {

/*- Node i determines whether it needs row ----------------------------------*/
/*- Node i determines whether row intersects slice --------------------------*/
/*- Node 0 determines whether node i has slice ------------------------------*/

    index_row(irow,length,start,ndim,index);
    my_row = 1;
    for (idim = 1; idim < ndim; idim++) {
      if (index[idim] < start[idim] || index[idim] > end[idim]) {
	my_row = 0; 
        break;
      }
    }

    slice_row = 1;
    for (idim = 1; idim < ndim; idim++) {
      if (index[idim] < slab_lo[idim] || index[idim] > slab_hi[idim]) {
	slice_row = 0;  
        break;
      }
    }
    if (slice_row && (slab_hi[0] >= start[0] && slab_lo[0] <= end[0])) {
      my_slice = my_row;
    }
    else {
      my_slice = 0;
    }

    if (slice_row) {

/*- Node 0 reads and broadcasts row -----------------------------------------*/

      if (inode == 0) {
	if (fread (buffer, esize, ne, stream) != ne) {
	  map_error("Cannot read expected number of elements");
        }
      }
      MPI_Bcast(buffer,ne,etype,0,comm_cart);

/*- Node i identifies portion of local memory belonging to row. -------------*/

      if (my_slice) {
	count = offset = slice_offset = 0;
	/* debug ('i', "node %d has slice\n", inode); */
	offset_index(end,start,pad,index,ndim,&offset);
	slice_offset = slab_lo[0]-start[0];
	if (slice_offset < 0) {
          slice_offset = 0;
        }

	if (slab_hi[0] > end[0]) {
          count = end[0]-(start[0]+slice_offset)+1;
        }
	else {
          count = slab_hi[0]-(start[0]+slice_offset)+1;
        }
      
	if(start[0]-slab_lo[0] > 0) {
	  memcpy((char *)local+(offset+slice_offset)*esize,
		  &buffer[(start[0]-slab_lo[0])*esize],count*esize);
        }
	else {
	  memcpy((char *)local+(offset+slice_offset)*esize,
		  &buffer[0],count*esize);
        }
      
	debug('f',"%d bytes read ",count*esize);
	debug('f',"at offset %d bytes ",(offset+slice_offset)*esize);
	debug('f',"from address %#x\n",local);
      
	bytes_read += count*esize;
      }
      
      if (inode == 0) {
	debug('i',"%d of ",++orow);
	debug('i',"%d rows read\r",nrow);
      }
    }
  }

  if (inode == 0) {
    debug('i',"%s","\n");
  }

  MPI_Barrier(comm_cart);
  t2 = MPI_Wtime();

  if (inode == 0) {
    debug('t',"MPI_Cart_read:  elapsed time = %g\n",t2-t1);
  }

  (void)free((char *)buffer);

  debug ('f',"%d bytes read\n",bytes_read);

  (void)free((char *)start);
  (void)free((char *)end);
  (void)free((char *)index);
  (void)free((char *)dims);
  (void)free((char *)periods);
  (void)free((char *)coords);

  return;
}
