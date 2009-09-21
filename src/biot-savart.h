/* BIOT-SAVART.H */
/* Copyright (c) 2009 Louis F. Rossi                                    *
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
 * or                                                                   *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

/* This header describes computation of the Biot-Savart integral for
   elliptical Gaussian basis functions. */


/* This requires knowing what a particle is.*/

/* Direct evaluation */

extern void vort_vort_interaction(int,int);

extern void induced_v(Blob_external *,Blob_internal *,Blob_parms *, 
		      double, double, double[]);

/* Rodrigo's spectral velocity code */
extern void read_data(void);
extern void eval_biot(double, double, double, double[]);

/* Asymptotic approximation headers. */

#define psi2coeffA 1.0/(1.0+a2)
#define psi2coeffB a2/(1.0+a2)

/* PSI_ORDER sets the order of the streamfunction 
   in epsilon in the near field. */

#define PSI_ORDER 6

#if PSI_ORDER == 4
#define MAXEXP 6
#define MAXPOLYN 5
#endif

#if PSI_ORDER == 5
#define MAXEXP 7
#define MAXPOLYN 6
#endif

#if PSI_ORDER == 6
#define MAXEXP 8
#define MAXPOLYN 7
#endif

extern double expint(double[], double, double);

extern void build_rt(double,double,double,double[],double[]);
extern void build_psi(double,double,double,
		      double[],double[],double[],double[][MAXEXP]);

extern double build_psi_x(double,double,double,double,
			  double,double,double,
			  double[],double[],double[],
			  double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_y(double,double,double,double,
			  double,double,double,
			  double[],double[],double[],
			  double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_xx(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_yy(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_xy(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_xxx(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_xxy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_xyy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][MAXEXP],double[][MAXEXP]);
extern double build_psi_yyy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][MAXEXP],double[][MAXEXP]);

/* Fast multipole summation */

#define PMAX 30          /* Maximum number of multipole coefficients. */
#define LMAX 7           /* Maximum number of levels of refinement in 
                            multipole summation */

typedef struct EltLink
{
   int element;
   struct EltLink *next;
}
FineGridLink;

extern double       minX,maxX,minY,maxY,distX,distY;
extern Complex      *Level_Ptr[LMAX];
extern int          gridx[LMAX][NMAX], gridy[LMAX][NMAX], mplevels;
extern FineGridLink **FineGridLinks;
extern int          numk2;


extern void Create_Hierarchy();
extern void Release_Links(int);

extern int  Set_Level();
extern void partition(int);

extern void MP_Sum(int, int);
extern void MP_Direct(int, int);
extern void MP_Direct2(int, int);
extern void MP_Direct3(int, int);

extern int  C(int,int);
extern Complex cmult(Complex, Complex);
extern Complex cpowi(Complex, int);
extern void Calc_Coeffs(int, Complex, Complex[]);


/* Complete velocity information. */

extern void   dpos_vel(int);
extern void   dpos_vel_fast(int);
extern Vector dpos_vel_gen(Vector,Blob_internal,Blob_parms);
extern Vector vel_gen(double,double);

extern void correct_vel_4();
