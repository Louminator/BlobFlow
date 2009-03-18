/* GLOBAL.H */
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
 * or after 15 January 2001                                             *
 * Louis Rossi                                                          *
 * Department of Mathematical Sciences                                  *
 * University of Delaware                                               *
 * Newark, DE 19716                                                     */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define SQR(X) ((X)*(X))
#define a2diff(X) (1.0/(X).a2-(X).a2)
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define SGN(X) ((X) < (0.0) ? (-1) : (1))

#define CORRECTVEL4
#undef  NOFASTMP
#define NOBOUNDARY

/* Set LINEAR to impose your own velocity field. */
#undef  LINEAR
#undef  cashkarp

#define NMAX 10000       /* Maximum number of computational elements. */
#define PMAX 30          /* Maximum number of multipole coefficients. */
#define LMAX 7           /* Maximum number of levels of refinement in 
                            multipole summation */
#define BMAX 250
#define MAX_SPLIT_CONF 8 /* Maximum number of elements that a single element */
                         /* can split into. */
#define FILENAME_LEN 500 /* Maximum file name length */

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

#define psi2coeffA 1.0/(1.0+a2)
#define psi2coeffB a2/(1.0+a2)

typedef struct
{
    double x,y;
}
Vector;

typedef struct
{
    double x,y;
    double strength;
    double dx,dy;
}
Blob_external;

typedef struct
{
    double a2,s2,th;
}
Blob_internal;

typedef struct
{
  double sin2,cos2,sincos,costh,sinth,du11,du12,du21;
  double u_xx,u_xy,u_yy,v_xx;
}
Blob_parms;

typedef struct
{
   double du11,du12,du21;  
}
Tensor;

typedef struct
{
    Blob_external blob0,blob1;
}
Metablob;

typedef struct
{
  double m,x,y,nx,ny,l;
}
Panel;

typedef struct
{
   double re,im;
}
Complex;

typedef struct EltLink
{
   int element;
   struct EltLink *next;
}
FineGridLink;

extern FILE          *comp_log,*diag_log,*mpi_log,*cpu_log;

extern int           N,oldN;
extern Metablob      mblob[NMAX];
extern Blob_internal blobguts[NMAX];
extern Blob_parms    tmpparms[NMAX];
extern double        visc;                       /* Physical constants */
extern double        alpha,l2tol,dtth_delta;     /* Numerical parameters */

/* Boundaries */
extern int    B,Bpiv[BMAX];
extern Panel  walls[BMAX];
extern double BdyMat[BMAX][BMAX];

/* Dynamic memory allocation might be better here. Nah! */
extern Vector       refinevels[NMAX][3][MAX_SPLIT_CONF];
extern int          refineindex[NMAX],refinestack;

/* Time integration parameters */
extern double TimeStep,PrefStep,FrameStep,EndTime,SimTime;
extern int    Frame,MaxOrder;

/*Fusion parameters */
extern double merge_budget,merge_p,merge_a2tol,merge_thtol;
extern double merge_mom3wt,merge_mom4wt,clusterR,aM2;
extern double MergeStep,merge_c,merge_growth_rate;
extern int    merge_estimator,nsplit,nmerge,totsplit,totmerge,MergeFrame;

extern double BoundaryStep;
extern int    BoundaryFrame;

/*Splitting parameters */
extern int    split_method;
extern double *split_parm_ptr;

/*Fast multipole variables*/
extern double       minX,maxX,minY,maxY,distX,distY;
extern Complex      *Level_Ptr[LMAX];
extern int          gridx[LMAX][NMAX], gridy[LMAX][NMAX], mplevels;
extern FineGridLink **FineGridLinks;
extern int          numk2;

/*Data export variables */
extern char filename[];

extern clock_t tot_cputime_ref,tot_cputime,
  vel_cputime_ref,vel_cputime,velsum_cputime_ref,velsum_cputime,
  veldirect_cputime_ref,veldirect_cputime,
  mp_cputime,mp_cputime_ref;

/* Function headers */

extern void init(int, char**);
extern int  main(int, char**);
extern void rectify();

extern void flip_blob();
extern void set_blob();

extern double dta2(Blob_internal*,Blob_parms*);
extern double dts2(Blob_internal*);
extern double dtth(Blob_internal*,Blob_parms*,double,double);

extern Vector induced_vel();
extern Tensor induced_veldev();

extern void rk4(double);
extern void split14(Blob_external *, Blob_internal *, Blob_parms *,
		    Blob_external *, Blob_internal *, Blob_parms *,
		    Blob_external *, Blob_internal *, Blob_parms *,
		    Blob_external *, Blob_internal *, Blob_parms *);
extern void splitvels(Vector[4],Blob_external,Blob_internal,Blob_parms);

/* Streamfunction approximation routines. */

extern double expint(double[], double, double);

extern void build_rt(double,double,double,double[],double[]);
extern void build_psiRT(double,double,double,double[],double[],double[][MAXEXP]);
extern void build_psi(double,double,double,
		      double[],double[],double[],double[][MAXEXP]);

extern void dx_coeff(double[], double[], int, double);
extern void dy_coeff(double[], double[], int, double);
extern void dxx_coeff(double[], double[], int, double);
extern void dyy_coeff(double[], double[], int, double);
extern void dxy_coeff(double[], double[], int, double);

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

extern void eval_biot(double, double, double, double[]);
extern void induced_v(Blob_external *,Blob_internal *,Blob_parms *,
		      double, double, double[]);
extern void induced_v_asympt(Blob_external *,Blob_internal *,Blob_parms *,
			     double, double, double[]);
extern void induced_v_platte(Blob_external *,Blob_internal *,Blob_parms *,
			     double, double, double[]);

extern double searchlist_uniform(int *, int , double,int *);
extern double searchlist_distribution(int *, int, double,int *);

extern int  Set_Level();
extern void partition(int);
extern void InitFineGrid(int);
extern int  C(int,int);
extern Complex cmult(Complex, Complex);
extern Complex cpowi(Complex, int);
extern void Calc_Coeffs(int, Complex, Complex[]);
extern void correct_vel_4();
extern void dpos_vel_fast(int);
extern void dpos_vel_exact(int);
extern void dpos_vel_exact_pi(int);
extern void dpos_vel_mix(int);
extern void Create_Hierarchy();
extern void vort_vort_interaction(int,int);
extern void MP_Sum(int, int);
extern void MP_Direct(int, int);
extern void MP_Direct2(int, int);
extern void MP_Direct3(int, int);
extern void resort();
extern void Release_Links(int);
extern void write_vorts(int);
extern void write_pmvorts(int);
extern void write_partition(int);
extern void split12(Blob_external*, Blob_internal*, Blob_parms*, Blob_external*, Blob_internal*, Blob_parms*);
extern void dpos_vel(int);
extern Vector dpos_vel_gen(Vector,Blob_internal,Blob_parms);
extern Vector dpos_vel_gen_linear(Vector,Blob_internal,Blob_parms);
extern Vector vel_gen(double,double);
extern void dpos_vel_linear(int);
extern void potential_flow(int);
extern Vector potential_vel(double,double);
extern void reflect_X();
extern void cache_resort();
extern void vel_field();

extern void BoundaryConstrain();

extern void stop(int);

extern void chk_nan(int);

extern void biot(double, int, double[], double[], int, int, double[]);

extern void solve_bdy_matrix(Panel[],int[],double[][BMAX]);
extern void factor_bdy_matrix(Panel[],int[],double[][BMAX]);
extern void bdy_vel(double, double, double*, double*);
extern void clip(double);

/* Rodrigo's code */
extern void read_data(void);

/* Lapack headers */
#ifndef NOBOUNDARY
extern int dgesv_(int*,int*,double[][BMAX],int*,int[],double[],int*,int*);
extern int dgetrf_( int*, int*, double[][BMAX], int*, int[], int*);
void dgetrs_( char*, int*, int*, double[][BMAX], int*, int[], double[], int*, int*);
#endif
