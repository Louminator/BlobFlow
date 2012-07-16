/* ESTAT.C */
/* Copyright (c) 2000 Louis F. Rossi                                    *
 * Elliptical Corrected Core Spreading Vortex Method (ECCSVM) for the   *
 * simulation/approximation of incompressible flows.                    *

 * This program is free software; you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation; either version 2 of the License, or    *
 * (at your option) any later version.					*

 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
 * GNU General Public License for more details.                         *

 * You should have received a copy of the GNU General Public License    *
 * along with this program; if not, write to the Free Software          *
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 *
 * USA                                                                  *
 *                                                                      *
 * or until 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Massachusetts Lowell                                   *
 * One University Ave                                                   *
 * Lowell, MA 01854                                                     *
 * 
 * or after 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19715-2553                                                */

#include <math.h>
#include <stdio.h>
#include <string.h>

#define  SQR(X) ((X) * (X))

#define  Title  160
#define  NMax   100000

int  grid,i,j,N,k;
char filename[Title],gridname[Title];
double posx,posy,amp,var,a2,theta,maxa2,avga2,wtavga2,wt;

void init(void)
{
   int  dummy;
   FILE *datafile,*fopen();
   
   N = 0;
   strcpy(gridname,filename);
   if ( (filename[strlen(filename)-3] != 'v') ||
       (filename[strlen(filename)-2] != 't') ||
       (filename[strlen(filename)-1] != 'x') )
     strcat(filename,".vtx");
   
   datafile = fopen(filename,"r");
   while ((dummy = getc(datafile)) != EOF)
     {
	ungetc( (char) dummy,datafile);
	fscanf(datafile,"%lf%lf%lf%lf%lf%lf\n",
	       &posx,&posy,&amp,&var,&a2,&theta);
	if (a2>1.0)
	  {
	     if (a2 > maxa2) maxa2=a2;
	     avga2 += a2;
	     wtavga2 += fabs(amp)*a2;
	  }
	else
	  {
	     if (1.0/a2 > maxa2) maxa2=1.0/a2;
	     avga2 += 1.0/a2;
	     wtavga2 += fabs(amp)/a2;
	  }
	wt += fabs(amp);
	++N;
     }
   avga2 /= N;
   wtavga2 /= wt;
}

int main(int argc, char *argv[])
{

   if (argc == 2)
     {
	strcpy(filename,argv[1]);
	init();
	printf("Max a2: %12.4e  Avg a2: %12.4e  Wt'd avg a2: %12.4e\n",
	       maxa2,avga2,wtavga2);
     }
   return(1);
}
