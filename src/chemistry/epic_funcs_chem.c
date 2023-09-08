/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 * Copyright (C) 1998-2018 Timothy E. Dowling                      *
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

/* * * * * * * * *  epic_funcs_chem.c  * * * * * * * * * * * * * * * * *
 *                                                                     *
 *  Timothy E. Dowling                                                 *
 *                                                                     *
 *  Chemistry functions that do not reference EPIC model variables     *
 *  by name or index (may reference planetspec).                       *
 *                                                                     *
 *  This file includes the following:                                  *
 *                                                                     *
 *      parse_species_name()                                           *
 *      b_vir(), b1_vir(), b2_vir()                                    *
 *      sum_xx()                                                       *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <epic.h>

/*====================== parse_species_name() ================================*/

/*
 * Takes a species name string (e.g., "NH_4SH") and returns  
 * the number of distinct elements, their symbols, and how many atoms 
 * of each are present.  Reallocates the necessary memory to hold this
 * information.
 *
 * NOTE: num_elements and **symbols shoud be declared static in the calling 
 *       function, with their input values equal to the last call, 
 *       in order properly reallocate memory.
 *
 */
void parse_species_name(char   *name,
                        int    *num_elements,
                        char ***symbols,
                        int   **counts)
{
  int
    num_caps,
    i,ii;
  char
    *ptr,
    subscript[4],
    format[4];
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="parse_species_name";

  /*
   * Free previous memory:
   */
  for (i = 0; i < *num_elements; i++) {
    free((*symbols)[i]);
  }

  /*
   * Count number of capital letters to determine working array sizes:
   */
  num_caps = 0;
  if (name) {
    ptr = name;
  }
  else {
    sprintf(Message,"input name is null");
    epic_error(dbmsname,Message);
  }
  while (*ptr != '\0') {
    if (isupper(*ptr)) num_caps++;
    ptr++;
  }

  /*
   * Reallocate memory:
   */
  *counts  = (int   *)realloc(*counts, num_caps*sizeof(int   ));
  *symbols = (char **)realloc(*symbols,num_caps*sizeof(char *));
  for (i = 0; i < num_caps; i++) {
    (*symbols)[i] = (char *)calloc(4,sizeof(char));
  }

  /*
   * Determine symbols and counts:
   */
  i   = 0;
  ptr = name;
  while (*ptr != '\0') {
    if (isupper(*ptr)) {
      if (islower(*(ptr+1))) {
        /* Element symbol has two letters: */
        strncpy((*symbols)[i],ptr,2);
        ptr+=2;
      }
      else {
        /* Element symbol has one letter: */
        strncpy((*symbols)[i],ptr,1);
        ptr+=1;
      }
      if (*ptr == '_') {
        subscript[0] = '\0';
        /* Determine subscript's number of places: */
        ii = 0;
        while(isdigit(*(++ptr))) {
          ii++;
        }
        if (ii == 0) {
          sprintf(Message," \"_\" not followed by digit: %s",name);
          epic_error(dbmsname,Message);
        }
        else {
          sprintf(format,"%%%dd",ii);
          sscanf(ptr-ii,format,*counts+i);
        }
      }
      else {
        (*counts)[i] = 1;
      }
      i++;
    }
    else {
      ptr++;
    }
  }
  /*
   * Trim arrays to refer to distinct elements:
   */
  *num_elements = num_caps;
  for (i = 0; i < num_caps; i++) {
    for (ii = i+1; ii < num_caps; ii++) {
      if ((*counts)[ii] > 0 && strcmp((*symbols)[i],(*symbols)[ii]) == 0) {
        (*counts)[i ] += (*counts)[ii];
        (*counts)[ii]  = 0;
        (*num_elements)--;
      }
    }
  }
  /* Remove zero entries by shifting: */
  for (i = 0; i < num_caps; i++) {
    if ((*counts)[i] == 0) {
      for (ii = i; ii < num_caps-1; ii++) {
        (*counts)[ii] = (*counts)[ii+1];
        strcpy((*symbols)[ii],(*symbols)[ii+1]);
      }
      (*counts)[num_caps-1] = 0;
    }
  }
  /* Trim allocated memory */
  for (i = (*num_elements); i < num_caps; i++) {
    free((*symbols)[i]);
  }
  *counts  = (int   *)realloc(*counts, (*num_elements)*sizeof(int   ));
  *symbols = (char **)realloc(*symbols,(*num_elements)*sizeof(char *));

  return;
}

/*====================== end of parse_species_name() =========================*/

/*====================== b_vir() =============================================*/

/*
 * Returns 2nd virial coefficient B_{ab} as a function of temperature.
 * For example, b_vir("H_2","H_2",temperature) returns B for pure molecular
 * hydrogen, H_2, whereas b_vir("H_2","He",temperature) returns the cross-term  
 * B for H_2+He. Data are taken from "The virial coefficients of pure gases
 * and mixtures," by J.H. Dymond and E.B. Smith (1980, Oxford), and converted
 * for a pressure expansion in mks units.
 */

#define MAX_CHEM_PAIRS 8

EPIC_FLOAT b_vir(char       *chem_a,
                 char       *chem_b,
                 EPIC_FLOAT  temperature) 
{
  char
    header[N_STR],
    infile[N_STR],
    chem_ab[16];
  /* Memory for up to MAX_CHEM_PAIRS different chemical pairs: */
  static char
    list[MAX_CHEM_PAIRS][16];
  int
    i,j;
  static int
    ndat[MAX_CHEM_PAIRS],
    count=0;
  EPIC_FLOAT 
    b,t_d,
    *buffer;
  static float_triplet
    *b_table[MAX_CHEM_PAIRS];
  static EPIC_FLOAT
    first_dbdt[MAX_CHEM_PAIRS],
    last_dbdt[ MAX_CHEM_PAIRS];
  FILE
    *input;
  /* 
   * The following are part of DEBUG_MILESTONE(.) statements: 
   */
  int
    idbms=0;
  static char
    dbmsname[]="b_vir";

  /* 
   * Determine list index for chemical pair: 
   */
  sprintf(chem_ab,"%s%s",chem_a,chem_b);
  i = 0;
  while (strcmp(chem_ab,list[i]) != 0 && i < count && i < MAX_CHEM_PAIRS) {
    i++;
  }

  if (i >= MAX_CHEM_PAIRS) {
    sprintf(Message,"b_vir() exceeded MAX_CHEM_PAIRS = %d",MAX_CHEM_PAIRS);
    epic_error(dbmsname,Message);
  }
  else if (i == count) {
    /* 
     * Chemical pair not on list. 
     * Add to list and input data. 
     */
    if (IAMNODE == 0) {
      sprintf(list[count],"%s",chem_ab);
      if (strcmp(chem_a,chem_b) == 0) {
        sprintf(infile,EPIC_PATH"/data/chemistry/virial/b_vs_t.%s",chem_a);
        input = fopen(infile,"r");
        if (!input) {
          fprintf(stderr,"Warning: b_vir() cannot find %s \n",infile);
          fprintf(stderr,"         Defaulting to ideal equation of state for %s. \n",chem_a);
        }
      }
      else {
        sprintf(infile,EPIC_PATH"/data/chemistry/virial/b_vs_t.%sx%s",chem_a,chem_b);
        input = fopen(infile,"r");
        if (!input) {
          /* Try the names reversed: */
          sprintf(infile,EPIC_PATH"/data/chemistry/virial/b_vs_t.%sx%s",chem_b,chem_a);
          input = fopen(infile,"r");
          if (!input) {
            fprintf(stderr,"Warning: b_vir() cannot find %s \n",infile);
            fprintf(stderr,"         Defaulting to ideal equation of state for %s.\n",chem_ab);
          }
        }
      }
      if (input) {
        /* Skip over header: */
        for (j = 0; j < 6; j++) {
          fgets(header,100,input);  
        }
        /* Input number of data points: */
        fscanf(input,"%d",ndat+count); 
        /* 
         * Allocate memory: 
         */
        b_table[count] = ftriplet(0,ndat[count]-1,dbmsname);

        /* Input B(T): */
        for (j = 0; j < ndat[count]; j++) {

#if EPIC_PRECISION == DOUBLE_PRECISION
          fscanf(input,"%lf %lf %*f",&b_table[count][j].x,&b_table[count][j].y);
#else
          fscanf(input,"%f %f %*f",&b_table[count][j].x,&b_table[count][j].y);
#endif

          /* Convert for pressure expansion in mks units: */
          b_table[count][j].y /= 1.e+3*R_GAS*b_table[count][j].x;
        }
        fclose(input);
      }
      else {
        /*
         * When a data file is not available, default to ideal equation of state
         * by setting bdat to zero:
         */
        ndat[count] = 3;
        /*
         * Allocate memory.
         */
        b_table[count] = ftriplet(0,ndat[count]-1,dbmsname);
        for (j = 0; j < ndat[count]; j++) {
          b_table[count][j].x = (EPIC_FLOAT)(j+1)*100.;
          b_table[count][j].y = 0.;
        }
      }
      count++;
    }

#if defined(EPIC_MPI)
    MPI_Bcast(&count,                   1,MPI_INT,   NODE0,para.comm);
    MPI_Bcast(list[count-1],            8,MPI_CHAR,  NODE0,para.comm);
    MPI_Bcast(ndat+(count-1),           1,MPI_INT,   NODE0,para.comm);
    /* pack buffer */
    buffer = fvector(0,2*ndat[count-1]-1,dbmsname);
    for (j = 0; j < ndat[count-1]; j++) {
      buffer[j              ] = b_table[count-1][j].x;
      buffer[j+ndat[count-1]] = b_table[count-1][j].y;
    }
#  if EPIC_PRECISION == DOUBLE_PRECISION
     MPI_Bcast(buffer,2*ndat[count-1],MPI_DOUBLE,NODE0,para.comm);
#  else
     MPI_Bcast(buffer,2*ndat[count-1],MPI_FLOAT,NODE0,para.comm);
#  endif
    /* unpack buffer */
    for (j = 0; j < ndat[count-1]; j++) {
      b_table[count-1][j].x = buffer[j              ];
      b_table[count-1][j].y = buffer[j+ndat[count-1]];
    }
    free_fvector(buffer,0,2*ndat[count-1]-1,dbmsname);
#endif

    /* Set endpoint slopes: */
    first_dbdt[count-1] = 
          (b_table[count-1][1].y-b_table[count-1][0].y)/
          (b_table[count-1][1].x-b_table[count-1][0].x);
    last_dbdt[count-1]  = 
          (b_table[count-1][ndat[count-1]-1].y-b_table[count-1][ndat[count-1]-2].y)/
          (b_table[count-1][ndat[count-1]-1].x-b_table[count-1][ndat[count-1]-2].x);
    /* Prepare for monotonic spline interpolation: */
    spline_pchip(ndat[count-1],b_table[count-1]);
  }

  /* 
   * Main function evaluation:
   */

  /* Use cubic-spline interpolation: */
  if (temperature <= b_table[i][0].x) {
    /* At or before start of table. */
    b = b_table[i][0].y+first_dbdt[i]*(temperature-b_table[i][0].x);
  }
  else if (temperature >= b_table[i][ndat[i]-1].x) {
    /* At or after end of table. */
    b = b_table[i][ndat[i]-1].y+last_dbdt[i]*(temperature-b_table[i][ndat[i]-1].x);
  }
  else {
    j = find_place_in_table(ndat[i],b_table[i],temperature,&t_d);
    b = splint_pchip(temperature,b_table[i]+j,t_d);
  }

  return b;
}

/*====================== end of b_vir() ======================================*/

/*====================== b1_vir() ============================================*/

/*
 * Returns B1 = T dB/dT.
 */
EPIC_FLOAT b1_vir(char       *chem_a,
                  char       *chem_b,
                  EPIC_FLOAT  temperature)
{
  EPIC_FLOAT
    b1,tt;
  static EPIC_FLOAT
    dt=1.;

  tt = temperature/dt;
  b1 = (b_vir(chem_a,chem_b,temperature+.5*dt)-
        b_vir(chem_a,chem_b,temperature-.5*dt))*tt;

  return b1;
}

/*====================== end of b1_vir() =====================================*/

/*====================== b2_vir() ============================================*/

/*
 * Returns B2 = T^2 (d/dT)^2 B.
 */
EPIC_FLOAT b2_vir(char       *chem_a,
                  char       *chem_b,
                  EPIC_FLOAT  temperature)
{
  EPIC_FLOAT
    b2,tt;
  static EPIC_FLOAT
    dt=1.;

  tt = temperature/dt;
  b2 = (    b_vir(chem_a,chem_b,temperature+dt)
        -2.*b_vir(chem_a,chem_b,temperature   )
           +b_vir(chem_a,chem_b,temperature-dt))*tt*tt;

  return b2;
}

/*====================== end of b2_vir() =====================================*/

/*====================== sum_xx() ============================================*/
/*
 * Returns sum of 2nd virial coefficient, or related function,
 * with quadratic mole-fraction weighting appropriate to specified planet.
 */
EPIC_FLOAT sum_xx(planetspec *planet,
                  EPIC_FLOAT (*b_func)(char *,
                               char *,
                               EPIC_FLOAT),
                  EPIC_FLOAT   temperature)
{
  static int
    initialized = FALSE;
  EPIC_FLOAT
    b_sum,x_sum;
  static EPIC_FLOAT
    x_H_2,x_He;

  if (strcmp(grid.eos,"ideal") == 0) {
    /*
     * B = 0 for ideal case:
     */
    b_sum = 0.;
  }
  else if (strcmp(grid.eos,"virial") == 0) {
    if (strcmp(planet->type,"gas-giant") == 0) {
      if (!initialized) {
        /* NOTE: currently including only H_2+He: */
        x_sum  = planet->x_h2+planet->x_he;
        x_H_2  = planet->x_h2/x_sum;
        x_He   = planet->x_he/x_sum;

        initialized = TRUE;
      }
      b_sum =    (*b_func)("H_2","H_2",temperature)*x_H_2*x_H_2
             +2.*(*b_func)("H_2","He", temperature)*x_H_2*x_He
                +(*b_func)("He", "He", temperature)*x_He *x_He;
    }
    else {
      switch(planet->index) {
        case VENUS_INDEX:
        case VENUS_LLR05_INDEX:
        case MARS_INDEX:
          b_sum = (*b_func)("CO_2","CO_2",temperature);
        break;
        case EARTH_INDEX:
        case TITAN_INDEX:
          b_sum = (*b_func)("N_2","N_2",temperature);
        break;
        default:
          /* Default to ideal case. */
          if (!initialized) {
            if (IAMNODE == 0) {
              fprintf(stderr,"Warning: sum_xx(): equation of state = %s not defined for %s, "
                             "using ideal equation of state.\n",grid.eos,planet->name);
            }
            initialized = TRUE;
          }
          b_sum = 0.;
        break;
      }
    }
  }
  else {
    fprintf(stderr,"Unrecognized equation of state = %s in sum_xx()\n",grid.eos);
    exit(1);
  }

  return b_sum;
}

/*====================== end of sum_xx() =====================================*/

/* * * * * * * * * * * end of epic_funcs_chem.c * * * * * * * * * * * * * * * */
