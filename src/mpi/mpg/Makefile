#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# *                                                                 *
# * Copyright (C) 1998-2010 Timothy E. Dowling                      *
# *                                                                 *
# * This program is free software; you can redistribute it and/or   *
# * modify it under the terms of the GNU General Public License     *
# * as published by the Free Software Foundation; either version 2  *
# * of the License, or (at your option) any later version.          *
# * A copy of this License is in the file:                          *
# *   $EPIC_PATH/License.txt                                        *
# *                                                                 *
# * This program is distributed in the hope that it will be useful, *
# * but WITHOUT ANY WARRANTY; without even the implied warranty of  *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            *
# *                                                                 *
# * You should have received a copy of the GNU General Public       *
# * License along with this program; if not, write to the Free      *
# * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     *
# * Boston, MA 02110-1301, USA.                                     *
# *                                                                 *
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

# Supported MPI implementations: openmpi, LAM, mpich


all:
	@if [ $(MPI_TYPE) = openmpi ] ; then \
	  echo "Making mpg..." ; \
	  make -f Makefile.openmpi ; \
	elif [ $(MPI_TYPE) = LAM ] ; then \
	  echo "Making mpg..." ; \
	  make -f Makefile.lam ; \
	else \
	  mpireconfig Makefile.mpich ; \
	  echo "Making mpg..." ; \
	  make -f Makefile.mpich ; \
	fi ; \

install:
	@if [ $(MPI_TYPE) = openmpi ] ; then \
	  echo "Install mpg..." ; \
	  make install -f Makefile.openmpi ; \
	elif [ $(MPI_TYPE) = LAM ] ; then \
	  echo "Installing mpg..." ; \
	  make install -f Makefile.lam ; \
	else \
	  mpireconfig Makefile.mpich ; \
	  echo "Installing mpg..." ; \
	  make install -f Makefile.mpich ; \
	fi

clean:
	@rm -f core *o *BAK *bak *~ *% *.log ._* ./c_original/._* ./fortran_interface/._*

clear:
	@rm -f core *o *BAK *bak *~ *% *.log *.a ._* ./c_original/._* ./fortran_interface/._*
	
depend: 
	@if [ $(MPI_TYPE) = openmpi ] ; then \
	  make depend -f Makefile.openmpi ; \
	elif [ $(MPI_TYPE) = LAM ] ; then \
	  make depend -f Makefile.lam ; \
	else \
	  mpireconfig Makefile.mpich ; \
	  make depend -f Makefile.mpich ; \
	fi

