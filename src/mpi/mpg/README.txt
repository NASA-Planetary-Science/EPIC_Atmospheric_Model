/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2010 Timothy E. Dowling                      *
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

The mpg library provides functions for domain decomposition of 
gridded models.

To turn on debugging, edit the file $EPIC_PATH/include/qa.h
and uncomment the debug() macro. Set the QA_DEBUG environment
variable in your .cshrc file (or equivalent) as described.

The mpg library was written by Joe Matarese in the late 1990s, 
in his capacity as a system supporter for the MIT Earth Resources 
Laboratory's nCUBE 2. Matarese is no longer supporting the
mpg code, but the good news is the code has held up over time.

The original mpg code calls one function from the MPICH
MPE library, MPE_Decomp1d(). To allow the MPG library to compile
on implementations of MPI other than MPICH, the source code for 
MPE_Decomp1d() is included here with permission as mpich_decomp.c.

The Fortran interface source code has been separated into
a subdirectory because it caused problems that are
described in ./fortran_interface/README.txt. The original
C interface code is  preserved in the subdirectory ./c_original.
