/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2009 Timothy E. Dowling                      *
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

#ifndef EPIC_IO_MACROS_H
#define EPIC_IO_MACROS_H

/*======================= input/output macros =====================*/

/*
 * NOTE: With gcc 3.3, preprocessor constructions like 
 *          planet ## -> ## g
 *       generate a fatal error, because pasting is checked in a 
 *       pairwise fashion, and  neither planet-> nor ->g
 *       is a valid preprocessing token.  In previous versions of
 *       gcc, a warning was issued when such constructions were encountered, 
 *       but gcc 3.3 has elevated this warning to a compilation error.
 *       As a workaround, we have rewritten these input/output macros to
 *       not use the ## preprocessor operator.
 *
 * NOTE: The precision of floating-point data is given by EPIC_FLOAT
 */

/*
 * Because netCDF assumes that the last dimension of a variable varies the 
 * fastest, the number associated with each dimension is backwards from 
 * the usual EPIC convention:
 */
#define NETCDF_T_INDEX  0
#define NETCDF_K_INDEX  1
#define NETCDF_J_INDEX  2
#define NETCDF_I_INDEX  3

/*
 * Input macros:
 */

/*
 * Macro for reading character data.
 */
#define READC(iname,ename,num) \
  if (IAMNODE == NODE0) { \
    nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
    if (nc_err == NC_NOERR) { \
      if (the_nc_type == NC_CHAR) { \
        nc_get_att_text(nc_id,NC_GLOBAL,#ename,iname); \
      } \
      else { \
        sprintf(Message,"READC, %s, the_nc_type != NC_CHAR",#iname); \
        epic_error(dbmsname,Message); \
      } \
    } \
    else { \
      sprintf(Message,"READC, %s, %s",#iname,nc_strerror(nc_err)); \
      epic_warning(dbmsname,Message); \
    } \
  } \
  bcast_char(NODE0,iname,(num));

/*
 * Macro for reading integer data.
 */
#define READI(iname,ename,num) \
  if (IAMNODE == NODE0) { \
    nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
    if (nc_err == NC_NOERR) { \
      if (the_nc_type == NC_INT) { \
        nc_get_att_int(nc_id,NC_GLOBAL,#ename,iname); \
      } \
      else { \
        sprintf(Message,"READI, %s, the_nc_type != NC_INT",#iname); \
        epic_error(dbmsname,Message); \
      } \
    } \
    else { \
      sprintf(Message,"READI, %s, %s",#iname,nc_strerror(nc_err)); \
      epic_error(dbmsname,Message); \
    } \
  } \
  bcast_int(NODE0,iname,(num));

/*
 * Macro for reading time_t datum.
 *
 * NOTE: When 64-bit integers are easy to handle in MPI and netCDF,
 *       change the buffer to 64 bit.
 */
#define READTIME(iname,ename) \
  if (1) { \
    int ibuff; \
    if (IAMNODE == NODE0) { \
      nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
      if (nc_err == NC_NOERR) { \
        if (the_nc_type == NC_INT) { \
          nc_get_att_int(nc_id,NC_GLOBAL,#ename,&ibuff); \
        } \
        else { \
          sprintf(Message,"READTIME, %s, the_nc_type != NC_INT",#iname); \
          epic_error(dbmsname,Message); \
        } \
      } \
      else { \
        sprintf(Message,"READTIME, %s, %s",#iname,nc_strerror(nc_err)); \
        epic_error(dbmsname,Message); \
      } \
    } \
    bcast_int(NODE0,&ibuff,1); \
    *iname = (time_t)ibuff; \
  }

/*
 * Macro for reading floating-point data (precision controlled by EPIC_PRECISION).
 */
#if EPIC_PRECISION == DOUBLE_PRECISION
#  define READF(iname,ename,num) \
     if (IAMNODE == NODE0) { \
       nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
       if (nc_err == NC_NOERR) { \
         if (the_nc_type == NC_DOUBLE) { \
           nc_get_att_double(nc_id,NC_GLOBAL,#ename,iname); \
         } \
         else { \
           sprintf(Message,"READF, %s, the_nc_type != NC_DOUBLE",#iname); \
           epic_error(dbmsname,Message); \
         } \
       } \
       else { \
         sprintf(Message,"READF, %s, %s",#iname,nc_strerror(nc_err)); \
         epic_error(dbmsname,Message); \
       } \
     } \
     bcast_float(NODE0,iname,(num));
#else
#  define READF(iname,ename,num) \
     if (IAMNODE == NODE0) { \
       nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
       if (nc_err == NC_NOERR) { \
         if (the_nc_type == NC_FLOAT) { \
           nc_get_att_float(nc_id,NC_GLOBAL,#ename,iname); \
         } \
         else { \
           sprintf(Message,"READF, %s, the_nc_type != NC_FLOAT",#iname); \
           epic_error(dbmsname,Message); \
         } \
       } \
       else { \
         sprintf(Message,"READF, %s, %s",#iname,nc_strerror(nc_err)); \
         epic_error(dbmsname,Message); \
       } \
     } \
     bcast_float(NODE0,iname,(num));
#endif

/*
 * Macro for reading double precision data.
 */
#define READD(iname,ename,num) \
   if (IAMNODE == NODE0) { \
     nc_err = nc_inq_atttype(nc_id,NC_GLOBAL,#ename,&the_nc_type); \
     if (nc_err == NC_NOERR) { \
       if (the_nc_type == NC_DOUBLE) { \
         nc_get_att_double(nc_id,NC_GLOBAL,#ename,iname); \
       } \
       else { \
         sprintf(Message,"READD, %s, the_nc_type != NC_DOUBLE",#iname); \
         epic_error(dbmsname,Message); \
       } \
     } \
     else { \
       sprintf(Message,"READD, %s, %s",#iname,nc_strerror(nc_err)); \
       epic_error(dbmsname,Message); \
     } \
   } \
   bcast_double(NODE0,iname,(num));

/* 
 * Output macros:
 */

/*
 * Macro for writing character data.
 */
#define WRITEC(iname,ename,num) \
  if (IAMNODE == NODE0) { \
    nc_err = nc_put_att_text(nc_id,NC_GLOBAL,#ename,(num),iname); \
    if (nc_err != NC_NOERR) { \
      sprintf(Message,"WRITEC, %s, %s",#iname,nc_strerror(nc_err)); \
      epic_error(dbmsname,Message); \
    } \
  }

/*
 * Macro for writing integer data.
 */
#define WRITEI(iname,ename,num) \
  if (IAMNODE == NODE0) { \
    nc_err = nc_put_att_int(nc_id,NC_GLOBAL,#ename,NC_INT,(num),iname); \
    if (nc_err != NC_NOERR) { \
      sprintf(Message,"WRITEI, %s, %s",#iname,nc_strerror(nc_err)); \
      epic_error(dbmsname,Message); \
    } \
  }

/*
 * Macro for writing time_t datum.
 *
 * NOTE: When 64-bit integers are easy to handle in MPI and netCDF,
 *       change the buffer to 64 bit.
 */
#define WRITETIME(iname,ename) \
  if (1) { \
    int ibuff; \
    if (IAMNODE == NODE0) { \
      ibuff = (int)*iname; \
      nc_err = nc_put_att_int(nc_id,NC_GLOBAL,#ename,NC_INT,1,&ibuff); \
      if (nc_err != NC_NOERR) { \
        sprintf(Message,"WRITETIME, %s, %s",#iname,nc_strerror(nc_err)); \
        epic_error(dbmsname,Message); \
      } \
    } \
  }

/*
 * Macro for writing floating-point data (precision controlled by EPIC_PRECISION).
 */
#if EPIC_PRECISION == DOUBLE_PRECISION
#  define WRITEF(iname,ename,num) \
     if (IAMNODE == NODE0) { \
         nc_err = nc_put_att_double(nc_id,NC_GLOBAL,#ename,NC_DOUBLE,(num),iname); \
       if (nc_err != NC_NOERR) { \
         sprintf(Message,"WRITEF, %s, %s",#iname,nc_strerror(nc_err)); \
         epic_error(dbmsname,Message); \
       } \
     }
#else
#  define WRITEF(iname,ename,num) \
     if (IAMNODE == NODE0) { \
         nc_err = nc_put_att_float(nc_id,NC_GLOBAL,#ename,NC_FLOAT,(num),iname); \
       if (nc_err != NC_NOERR) { \
         sprintf(Message,"WRITEF, %s, %s",#iname,nc_strerror(nc_err)); \
         epic_error(dbmsname,Message); \
       } \
     }
#endif

/*
 * Macro for writing double precision data.
 */
#define WRITED(iname,ename,num) \
   if (IAMNODE == NODE0) { \
       nc_err = nc_put_att_double(nc_id,NC_GLOBAL,#ename,NC_DOUBLE,(num),iname); \
     if (nc_err != NC_NOERR) { \
       sprintf(Message,"WRITEF, %s, %s",#iname,nc_strerror(nc_err)); \
       epic_error(dbmsname,Message); \
     } \
   }

/*============================ end of input/output macros =================================*/

#endif
