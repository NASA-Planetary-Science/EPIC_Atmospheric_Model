/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2011 Timothy E. Dowling                      *
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

#ifndef EPIC_PV_SCHEMES_H
#define EPIC_PV_SCHEMES_H
/* * * * * * * * * * * * * epic_pv_schemes.h * * * * * * * * * * * * * * * 
 *                                                                       *
 *  Timothy E. Dowling                                                   *
 *                                                                       *
 *  How potential vorticity is averaged onto the u and v grids controls  *
 *  important conservation properties.                                   *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#if PV_SCHEME == ARAKAWA_LAMB_1981
#  define AL_U   (2.*PV2(K,J+1,I+1)+PV2(K,J+1,I)+2.*PV2(K,J,I)+PV2(K,J,I+1))
#  define BE_U   (PV2(K,J+1,I)+2.*PV2(K,J+1,I-1)+PV2(K,J,I-1)+2.*PV2(K,J,I))
#  define GA_U   (2.*PV2(K,J+1,I)+PV2(K,J+1,I-1)+2.*PV2(K,J,I-1)+PV2(K,J,I))
#  define DE_U   (PV2(K,J+1,I+1)+2.*PV2(K,J+1,I)+PV2(K,J,I)+2.*PV2(K,J,I+1))
#  define EP1    (PV2(K,J+1,I+1)-PV2(K,J,I+1)+PV2(K,J+1,I)-PV2(K,J,I))
#  define EP2    (PV2(K,J+1,I)-PV2(K,J,I)+PV2(K,J+1,I-1)-PV2(K,J,I-1))
#  define AL_V   (2.*PV2(K,J,I+1)+PV2(K,J,I)+2.*PV2(K,J-1,I)+PV2(K,J-1,I+1))
#  define BE_V   (PV2(K,J,I+1)+2.*PV2(K,J,I)+PV2(K,J-1,I)+2.*PV2(K,J-1,I+1))
#  define GA_V   (2.*PV2(K,J+1,I+1)+PV2(K,J+1,I)+2.*PV2(K,J,I)+PV2(K,J,I+1))
#  define DE_V   (PV2(K,J+1,I+1)+2.*PV2(K,J+1,I)+PV2(K,J,I)+2.*PV2(K,J,I+1))
#  define PH1    (-PV2(K,J+1,I+1)+PV2(K,J+1,I)+PV2(K,J,I)-PV2(K,J,I+1))
#  define PH2    (-PV2(K,J,I+1)+PV2(K,J,I)+PV2(K,J-1,I)-PV2(K,J-1,I+1))
#  define PV_COEF (1./24.)
#elif PV_SCHEME == SADOURNY_1975
#  define AL_U   (PV2(K,J,I)+PV2(K,J+1,I)+PV2(K,J+1,I+1))
#  define BE_U   (PV2(K,J+1,I-1)+PV2(K,J+1,I)+PV2(K,J,I))
#  define GA_U   (PV2(K,J,I-1)+PV2(K,J,I)+PV2(K,J+1,I))
#  define DE_U   (PV2(K,J,I)+PV2(K,J+1,I)+PV2(K,J,I+1))
#  define EP1    0.
#  define EP2    0.
#  define AL_V   (PV2(K,J-1,I)+PV2(K,J,I)+PV2(K,J,I+1))
#  define BE_V   (PV2(K,J,I)+PV2(K,J,I+1)+PV2(K,J-1,I+1))
#  define GA_V   (PV2(K,J,I)+PV2(K,J,I+1)+PV2(K,J+1,I+1))
#  define DE_V   (PV2(K,J,I)+PV2(K,J+1,I)+PV2(K,J,I+1))
#  define PH1    0.
#  define PH2    0.
#  define PV_COEF (1./12.)
#endif

/* * * * * * * * * * end of epic_pv_schemes.h * * * * * * * * * * * * * * * */ 
#endif
