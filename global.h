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

#define CORRECTVEL4
#undef  NOFASTMP

/* Set LINEAR to impose your own velocity field. */
#undef  LINEAR
#undef  cashkarp

#define NMax 50000   /* Maximum number of computational elements. */
#define PMax 30      /* Maximum number of multipole coefficients. */
#define LMax 7     /* Maximum number of levels of refinement in 
                        multipole summation */
#define MaxSplitConf 8 /*Maximum number of elements that a single element */
                       /* can split into. */
#define Title 500     /* Maximum file name length */

/* PSI_ORDER sets the order of the streamfunction 
   in epsilon in the near field. */

#define PSI_ORDER 6

#if PSI_ORDER == 4
#define maxexp 6
#define maxpolyn 5
#endif

#if PSI_ORDER == 5
#define maxexp 7
#define maxpolyn 6
#endif

#if PSI_ORDER == 6
#define maxexp 8
#define maxpolyn 7
#endif

#define psi2coeffA 1.0/(1.0+a2)
#define psi2coeffB a2/(1.0+a2)

typedef struct
{
    double x,y;
}
vector;

typedef struct
{
    double x,y;
    double strength;
    double dx,dy;
    int    order;
}
blob_external;

typedef struct
{
    double a2,s2,th;
}
blob_internal;

typedef struct
{
  double sin2,cos2,sincos,costh,sinth,du11,du12,du21;
  double u_xx,u_xy,u_yy,v_xx;
  int    refinecnt,nint;
}
blobparms;

typedef struct
{
   double du11,du12,du21;  
}
tensor;

typedef struct
{
    blob_external blob0,blob1,blob2,blob3,blob4;
}
metablob;

typedef struct
{
   double re,im;
}
complex;

typedef struct EltLink
{
   int element;
   struct EltLink *next;
}
FineGridLink;

extern FILE          *comp_log,*diag_log,*mpi_log,*cpu_log;

extern int           N,oldN;
extern metablob      mblob[NMax];
extern blob_internal blobguts[NMax];
extern blobparms     tmpparms[NMax];
extern double        visc;                       /* Physical constants */
extern double        alpha,l2tol,dtth_delta;     /* Numerical parameters */

/* Dynamic memory allocation might be better here. Nah! */
extern vector       refinevels[NMax][3][MaxSplitConf];
extern int          refineindex[NMax],refinestack;

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
extern complex      *Level_Ptr[LMax];
extern int          gridx[LMax][NMax], gridy[LMax][NMax], mplevels;
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

extern double dta2(blob_internal*,blobparms*);
extern double dts2(blob_internal*);
extern double dtth(blob_internal*,blobparms*);

extern vector induced_vel();
extern tensor induced_veldev();

extern void push();
extern void ab2();
extern void ab3();
extern void ab4();
extern void am4();
extern void ab2half();
extern void ab3half();
extern void ab4half();
extern void am4half();
extern void split14(blob_external *, blob_internal *, blobparms *,
		    blob_external *, blob_internal *, blobparms *,
		    blob_external *, blob_internal *, blobparms *,
		    blob_external *, blob_internal *, blobparms *);
extern void splitvels(vector[4],blob_external,blob_internal,blobparms);

/* Streamfunction approximation routines. */

extern double expint(double[], double, double);

extern void build_rt(double,double,double,double[],double[]);
extern void build_psiRT(double,double,double,double[],double[],double[][]);
extern void build_psi(double,double,double,
		      double[],double[],double[],double[][]);

extern void dx_coeff(double[], double[], int, double);
extern void dy_coeff(double[], double[], int, double);
extern void dxx_coeff(double[], double[], int, double);
extern void dyy_coeff(double[], double[], int, double);
extern void dxy_coeff(double[], double[], int, double);

extern double build_psi_x(double,double,double,double,
			  double,double,double,
			  double[],double[],double[],
			  double[],double[][],double[][]);
extern double build_psi_y(double,double,double,double,
			  double,double,double,
			  double[],double[],double[],
			  double[],double[][],double[][]);
extern double build_psi_xx(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][],double[][]);
extern double build_psi_yy(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][],double[][]);
extern double build_psi_xy(double,double,double,double,
			   double,double,double,
			   double[],double[],double[],
			   double[],double[][],double[][]);
extern double build_psi_xxx(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][],double[][]);
extern double build_psi_xxy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][],double[][]);
extern double build_psi_xyy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][],double[][]);
extern double build_psi_yyy(double,double,double,double,
			    double,double,double,
			    double[],double[],double[],
			    double[],double[][],double[][]);

extern void induced_v(blob_external *,blob_internal *,blobparms *,
		      double, double, double[]);

extern double searchlist_uniform(int *, int , double,int *);
extern double searchlist_distribution(int *, int, double,int *);

extern int  Set_Level();
extern void partition(int);
extern void InitFineGrid(int);
extern int  C(int,int);
extern complex cmult(complex, complex);
extern complex cpowi(complex, int);
extern void Calc_Coeffs(int, complex, complex[]);
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
extern void rk4internal(blob_internal*, blobparms*, double);
extern void internal_march(blob_internal*, blobparms*, double);
extern void rkckmarch(blob_internal*, blobparms*, double,double);
extern void chksplit();
extern void split14vels(vector*, blob_external, blob_internal, blobparms);
extern void split14(blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*);
extern void split15vels(vector*, blob_external, blob_internal, blobparms);
extern void split15(blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*,
		    blob_external*, blob_internal*, blobparms*);
extern void split15asymvels(vector*, blob_external, blob_internal, blobparms, double*);
extern void split15asym(blob_external*, blob_internal*, blobparms*,
			blob_external*, blob_internal*, blobparms*,
			blob_external*, blob_internal*, blobparms*,
			blob_external*, blob_internal*, blobparms*,
			blob_external*, blob_internal*, blobparms*, double*);
extern void merge();
extern void resort();
extern void Release_Links(int);
extern void write_vorts(int);
extern void write_pmvorts(int);
extern void write_partition(int);
extern void split12(blob_external*, blob_internal*, blobparms*, blob_external*, blob_internal*, blobparms*);
extern void dpos_vel(int);
extern vector dpos_vel_gen(vector,blob_internal,blobparms);
extern vector dpos_vel_gen_linear(vector,blob_internal,blobparms);
extern vector vel_gen(double,double);
extern void dpos_vel_linear(int);
extern void potential_flow(int);
extern vector potential_vel(double,double);
extern void reflect_X();
extern void cache_resort();
extern void vel_field();

extern void BoundaryConstrain();

extern void stop(int);

extern void chk_nan(int);


