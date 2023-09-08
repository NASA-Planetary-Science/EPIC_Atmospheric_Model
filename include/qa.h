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

#include <stdlib.h>
#include <string.h>

/*---------------------------------------------------------------------------*
 * qa.h - Quality Assurance macros
 *
 * QA macros "debug" and "assert" are controlled by the environment variables
 * "QA_DEBUG" and "QA_ASSERT", respectively.
 *
 * If QA_ASSERT is set, then the assert macro is enabled.  Otherwise, the
 * assert macro is disabled.
 *
 * QA_DEBUG may be set to a string containing one or more keywords to
 * determine the behaviour of the debug macro.  The keywords are:
 *
 * f - for file i/o related diagnostics
 * i - for purely informational diagnostics
 * m - for memory related diagnostics
 * p - for parallel code diagnostics
 * v - for variable value diagnostics
 * z - for miscellaneous diagnostics not covered by the above
 * t - for timing info
 *
 *---------------------------------------------------------------------------*/

/********
#define debug(key,format,var) if (strchr(getenv("QA_DEBUG"), key)) { \
                                fprintf(stderr,format,var); \
                              }
*********/
/**********/
#define debug(key,format,var) ;
/***********/

/********
#define assert(code) if (getenv("QA_ASSERT") && !(code)) { \
                       fprintf (stderr, "Assertion \"%s\" failed\n", #code); \
                     }
*********/
#define assert(code) ;
