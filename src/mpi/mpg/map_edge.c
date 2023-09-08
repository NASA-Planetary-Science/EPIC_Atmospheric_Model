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

/*
 * See README.txt for instuctions on how to turn on the debugging.
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpg.h"
#include "qa.h"

#define M_EDGE_SEND 250

extern void  map_error (const char *s);

struct s_commlist 
{
  MPI_Status     
    *stat;
  MPI_Datatype   
    *stride;
  int            
    *hi_offset_src,
    *hi_offset_rcv,
    *hi,           
    *lo_offset_src, 
    *lo_offset_rcv, 
    *lo,            
    ref_offset;
  MPI_Comm        
    comm;
};

typedef struct s_commlist 
  Commlist;

void MPG_Cart_edge(MPI_Comm      comm_cart, 
                   int          *length, 
                   int          *pad, 
                   int           ndim,
		   MPI_Datatype  etype,
		   int          *ref_offset, 
                   int          *lo_offset_src, 
                   int          *lo_offset_rcv,
		   int          *hi_offset_src, 
                   int          *hi_offset_rcv, 
                   int          *lo, 
                   int          *hi,
		   MPI_Datatype *stride) 
{
  int       
    inode,nnode,idim,   
    *start,*end,
    *dims,*coords,
    *periods,*padded_length,
    padded_size,bstride,bcount,blength,
    esize;

  for (idim = 0; idim < ndim && !pad[idim]; idim++) {
    ;
  }
  if (idim == ndim) {
    map_error("edge exchange routine called with all pads equal to zero");
  }
  MPI_Comm_rank(comm_cart, &inode);
  MPI_Comm_size(comm_cart, &nnode);

  MPI_Type_size(etype,&esize);

  /* Allocate memory: */
  if ((start         = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (end           = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (dims          = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (padded_length = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (periods       = (int *)malloc(ndim*sizeof(int))) == NULL ||
      (coords        = (int *)malloc(ndim*sizeof(int))) == NULL) {
    map_error ("Cannot allocate space");
  }

  MPI_Cart_get(comm_cart,ndim,dims,periods,coords);
  for (idim = 0; idim < ndim; idim++) {
    MPI_Cart_shift(comm_cart,idim,1,&lo[idim],&hi[idim]);

    debug('p',"I am node %d; ",inode);
    debug('p',"above me in direction %d ",idim);
    debug('p',"is node %d, ",hi[idim]);
    debug('p',"below me is node %d\n",lo[idim]);
    MPE_Decomp1d(length[idim],dims[idim],coords[idim],
		 &start[idim],&end[idim]);
    start[idim]--;  
    end[  idim]--;
    padded_length[idim] = end[idim]-start[idim]+1+2*pad[idim];
  }

/*---------------------------------------------------------------------------*
 * Get the offset of the first element in the grid - (pad,pad,pad,...) w.r.t.
 * the padded dimension lengths.  Also get the offset of elements, (nx+pad-1,
 * pad,pad,...), (pad,ny+pad-1,pad,...) ...
 *---------------------------------------------------------------------------*/

  *ref_offset = 0;
  padded_size = 1;
  for (idim = 0; idim < ndim; idim++) {
    *ref_offset += pad[idim]*padded_size;
    padded_size *= padded_length[idim];
  }
  *ref_offset *= (int)esize;

  padded_size = (int)esize;
  for (idim = 0; idim < ndim; idim++) {
    lo_offset_src[idim] = pad[idim]*padded_size;
    lo_offset_rcv[idim] = 0;
    hi_offset_rcv[idim] = (end[idim]-start[idim]+1+pad[idim])*padded_size;
    hi_offset_src[idim] = (end[idim]-start[idim]+1)*padded_size;
    padded_size        *= padded_length[idim];
  }

/*---------------------------------------------------------------------------*
 * Create a new, "strided" datatype for the exchange in the "non-contiguous"
 * direction
 *---------------------------------------------------------------------------*/

  padded_size /= (int)esize;

  bstride = 1;  
  bcount  = padded_size;
  for (idim = 0; idim < ndim; idim++) {
    bcount  /= padded_length[idim];
    blength  = pad[idim]*bstride;  
    bstride *= padded_length[idim];
    MPI_Type_vector(bcount,blength,bstride,etype,&stride[idim]);
    MPI_Type_commit(&stride[idim]);

    debug('v',"For dimension %d, ",idim);
    debug('v',"type block count = %d, ",bcount);
    debug('v',"length = %d, ",blength);
    debug('v',"and stride = %d\n",bstride);
  }

  /* Free allocated memory */
  free(coords), free(periods), free(padded_length), free(dims), free(end), free(start);

  return;
}

void MPG_Cart_edgeexch(MPI_Comm      comm_cart,  /* communicator for cartesian grid       */
                       int           ndim,       /* dimensions of Cartesian grid          */
                       int          *length,     /* number of elements along each axis    */
                       int          *pad,        /* pad depth for each axis               */
                       MPI_Datatype  etype,      /* data type of grid elements to be read */
                       void         *input)      /* local array for storing data          */
{
  static int             
    ncomms=0;
  static Commlist       
    *commlist;
  Commlist              
    *cp;
  int                    
    inode,idim,icomm;

  for (icomm = 0; icomm < ncomms; icomm++) {
    if (commlist[icomm].comm == comm_cart) {
      break;
    }
  }

  if (icomm >= ncomms) {     /* create new communicator info */
    cp = commlist;
    ncomms++;
    debug('v',"Creating new comm. attributes - no. comm's = %d\n",ncomms);
    /* Allocate memory: */
    if ((commlist = (Commlist *)malloc(ncomms*sizeof(Commlist))) == NULL) {
      map_error ("Cannot allocate space");
    }
    for (icomm = 0; icomm < ncomms-1; icomm++) {
      commlist[icomm] = cp[icomm];
    }
    (void)free((char *)cp);

    commlist[icomm].comm = comm_cart;
    
    /* Allocate memory: */
    if ((commlist[icomm].stride =
        (MPI_Datatype *)malloc(ndim*sizeof(MPI_Datatype))) == NULL ||
        (commlist[icomm].stat =
        (MPI_Status *)malloc(ndim*2*sizeof(MPI_Status)))   == NULL ||
        (commlist[icomm].lo_offset_src =
        (int *)malloc(ndim*sizeof(int)))                   == NULL ||
        (commlist[icomm].lo_offset_rcv =
        (int *)malloc(ndim*sizeof(int)))                   == NULL ||
        (commlist[icomm].hi_offset_src =
        (int *)malloc(ndim*sizeof(int)))                   == NULL ||
        (commlist[icomm].hi_offset_rcv =
        (int *)malloc(ndim*sizeof(int)))                   == NULL ||
        (commlist[icomm].lo = 
        (int *)malloc(ndim*sizeof(int)))                   == NULL ||
        (commlist[icomm].hi = 
        (int *)malloc(ndim*sizeof(int)))                   == NULL) {
      map_error ("map_edge: cannot allocate space");
    }

    MPG_Cart_edge(comm_cart,length,pad,ndim,etype,
                &(commlist[icomm].ref_offset),
                  commlist[icomm].lo_offset_src,
                  commlist[icomm].lo_offset_rcv,
                  commlist[icomm].hi_offset_src,
                  commlist[icomm].hi_offset_rcv,
                  commlist[icomm].lo, 
                  commlist[icomm].hi,
                  commlist[icomm].stride);
  } 
  else {
    debug('v',"Using old comm. attributes - no. comm. = %d\n",icomm);
  }

  /* 
   * At this point (commlist+icomm) points to an existing comm in commlist
   * or points to the comm (comm_cart) just appended to the commlist.
   */

  cp = commlist+icomm;

  MPI_Comm_rank(MPI_COMM_WORLD,&inode);

  debug('z',"On node %d ",inode);
  for (idim = 0; idim < ndim; idim++) {
    debug('z',"lo_offset_src[%d] = ",idim);
    debug('z',"%d; ",cp->lo_offset_src[idim]);
  }
  debug('z',"%s","\n");

  for (idim = ndim-1; idim >= 0; idim--) {
    MPI_Sendrecv((char *)input+cp->lo_offset_src[idim],1,
                  cp->stride[idim],cp->lo[idim],M_EDGE_SEND+idim*2,
                  (char *)input+cp->hi_offset_rcv[idim],1,
                  cp->stride[idim],cp->hi[idim],M_EDGE_SEND+idim*2,
                  cp->comm,&(cp->stat[idim*2]));

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",cp->lo_offset_src[idim]);
    debug('p',"to node %d ",cp->lo[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2);

    debug('p',"Node %d ",inode);
    debug('p',"received edge from address %#x, ",input);
    debug('p',"offset %d ",cp->hi_offset_rcv[idim]);
    debug('p',"from node %d ",cp->hi[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2);

    MPI_Sendrecv((char *)input+cp->hi_offset_src[idim],1,
                  cp->stride[idim],cp->hi[idim],M_EDGE_SEND+idim*2+1,
                  (char *)input+cp->lo_offset_rcv[idim],1,
                  cp->stride[idim],cp->lo[idim],M_EDGE_SEND+idim*2+1,
                  cp->comm,&(cp->stat[idim*2+1]));

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",cp->hi_offset_src[idim]);
    debug('p',"to node %d ",cp->hi[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2+1);

    debug('p',"Node %d ",inode);
    debug('p',"received edge from address %#x, ",input);
    debug('p',"offset %d ",cp->lo_offset_rcv[idim]);
    debug('p',"from node %d ",cp->lo[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2+1);
  }
  
  return;
}

void MPG_Cart_edgesend(MPI_Comm      comm_cart,   /* communicator for cartesian grid       */
                       int           ndim,        /* dimensions of Cartesian grid          */
                       int          *length,      /* number of elements along each axis    */
                       int          *pad,         /* number of elements along each axis    */
                       MPI_Datatype  etype,       /* data type of grid elements to be read */
                       void         *input,       /* local array for storing data          */
                       MPI_Request  *request)
{
  static int            
    initialized=0,           
    *hi_offset_src, 
    *hi_offset_rcv, 
    *hi,            
    *lo_offset_src,
    *lo_offset_rcv, 
    *lo,             
    ref_offset;
  int                    
    idim,inode;
  static MPI_Datatype   
    *stride;

  if (!initialized) {
    /* Allocate memory: */
    if ((stride        = (MPI_Datatype *)malloc(ndim*sizeof(MPI_Datatype))) == NULL ||
	(lo_offset_src = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(lo_offset_rcv = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi_offset_src = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi_offset_rcv = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(lo            = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi            = (int          *)malloc(ndim*sizeof(int         ))) == NULL) {
      map_error ("Cannot allocate space");
    }

    MPG_Cart_edge(comm_cart,length,pad,ndim,etype,&ref_offset,
		  lo_offset_src,lo_offset_rcv,hi_offset_src,hi_offset_rcv,
		  lo,hi,stride);

    initialized = 1;
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&inode);

  for (idim = ndim-1; idim >= 0; idim--) {
    MPI_Isend((char *)input+lo_offset_src[idim],1,stride[idim],
	       lo[idim],M_EDGE_SEND+idim*2+1,comm_cart,&request[idim*2]);

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",lo_offset_src[idim]);
    debug('p',"to node %d ",lo[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2+1);

    MPI_Isend((char *)input+hi_offset_src[idim],1,stride[idim],
	      hi[idim],M_EDGE_SEND+idim*2,comm_cart,&request[idim*2+1]);

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",hi_offset_src[idim]);
    debug('p',"to node %d ",hi[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2+1);
  }
  
  return;
}

void MPG_Cart_edgerecv(MPI_Comm      comm_cart,  /* communicator for cartesian grid       */
                       int           ndim,       /* dimensions of Cartesian grid          */
                       int          *length,     /* number of elements along each axis    */
                       int          *pad,        /* number of elements along each axis    */
                       MPI_Datatype  etype,      /* data type of grid elements to be read */
                       void         *input,      /* local array for storing data          */
                       MPI_Request  *request,
                       MPI_Status   *stat)
{
  static int             
    initialized=0,          
    *hi_offset_src, 
    *hi_offset_rcv, 
    *hi,           
    *lo_offset_src, 
    *lo_offset_rcv, 
    *lo,            
    ref_offset;
  int                    
    idim,inode;
  static MPI_Datatype   
    *stride;

  if (!initialized) {
    /* Allocate memory */
    if ((stride        = (MPI_Datatype *)malloc(ndim*sizeof(MPI_Datatype))) == NULL ||
	(lo_offset_src = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(lo_offset_rcv = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi_offset_src = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi_offset_rcv = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(lo            = (int          *)malloc(ndim*sizeof(int         ))) == NULL ||
	(hi            = (int          *)malloc(ndim*sizeof(int         ))) == NULL) {
      map_error ("Cannot allocate space");
    }
    MPG_Cart_edge(comm_cart,length,pad,ndim,etype,&ref_offset,
		  lo_offset_src,lo_offset_rcv,hi_offset_src,hi_offset_rcv,
		  lo,hi,stride);

    initialized = 1;
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&inode);
  MPI_Waitall(ndim*2,request,stat);

  for (idim = ndim-1; idim >= 0; idim--) {
    MPI_Recv((char *)input+lo_offset_rcv[idim],1,stride[idim],
	     lo[idim],M_EDGE_SEND+idim*2,comm_cart,&stat[idim*2]);

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",lo_offset_rcv[idim]);
    debug('p',"to node %d ",lo[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2);

    MPI_Recv((char *)input+hi_offset_rcv[idim],1,stride[idim],
	     hi[idim],M_EDGE_SEND+idim*2+1,comm_cart,&stat[idim*2+1]);

    debug('p',"Node %d ",inode);
    debug('p',"sent edge from address %#x, ",input);
    debug('p',"offset %d ",hi_offset_rcv[idim]);
    debug('p',"to node %d ",hi[idim]);
    debug('p',"via tag %d\n",M_EDGE_SEND+idim*2+1);
  }

  return;
}
