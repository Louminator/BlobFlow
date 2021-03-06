#include <math.h>
#include <stdio.h>
#include <string.h>

#define  SQR(X) ((X) * (X))

#define  Title  160
#define  NMax   300000

int  grid,i,j,N,k,quiet=0,xanti=0;
char filename[Title],gridname[Title];
double posx[NMax],posy[NMax],amp[NMax],var[NMax],a2[NMax],theta[NMax],field,
       xbegin,xend,ybegin,yend,xstep,ystep;

init()
{
   int  dummy;
   FILE *datafile,*deffile,*fopen();
   
   N = 0;
   strcpy(gridname,filename);
   if ( (filename[strlen(filename)-3] == 'v') &&
       (filename[strlen(filename)-2] == 't') &&
       (filename[strlen(filename)-1] == 'x') )
     {
	gridname[strlen(gridname)-3] = 'g';
	gridname[strlen(gridname)-2] = 'r';
	gridname[strlen(gridname)-1] = 'd';
     }
   else
     {    
	strcat(gridname,".grd");
	strcat(filename,".vtx");
     }
   
   if ((deffile = fopen("egrid.default","r")) == NULL)
     {
	printf("No defaults file found.\n");
	printf("Lower left: ");
	scanf("%lf %lf",&xbegin,&ybegin);
	printf("Upper right: ");
	scanf("%lf %lf",&xend,&yend);
	printf("Grid size: ");
	scanf("%d",&grid);
     }
   else
     {
	fscanf(deffile,"%lf %lf",&xbegin,&ybegin);
	fscanf(deffile,"%lf %lf",&xend,&yend);
	fscanf(deffile,"%d",&grid);
	fclose(deffile);
	if (!quiet)
	  {
	     printf("Using egrid defaults file for parameters...\n");
	     printf("Lower left: (%lf,%lf)\n",xbegin,ybegin);
	     printf("Upper right: (%lf,%lf)\n",xend,yend);
	     printf("Grid: %d\n",grid);
	  }
     }

  xstep=(xend-xbegin)/(grid-1.0);
  ystep=(yend-ybegin)/(grid-1.0);

  datafile = fopen(filename,"r");
  while ((dummy = getc(datafile)) != EOF)
    {
      ungetc( (char) dummy,datafile);
      fscanf(datafile,"%lf%lf%lf%lf%lf%lf\n",
	     &posx[N],&posy[N],&amp[N],&var[N],&a2[N],&theta[N]);
      ++N;
    }
  fclose(datafile);

  if (xanti==1)
    {
      for (i=0; i<N; ++i)
	{
	  posx[i+N]  = posx[i];
	  posy[i+N]  = -posy[i];
	  amp[i+N]   = -amp[i];
	  var[i+N]   = var[i];
	  a2[i+N]    = a2[i];
	  theta[i+N] = -theta[i];
	}
      N *= 2;
    }
}

main(int argc, char *argv[])
{
   double dx,dy,c2,s2,c,s,temp;
   FILE *outfile,*fopen();
   
   quiet = 0;
   
   if (argc > 1)
     {
	for (i=1; i<argc; ++i)
	  {
	     if (strcmp(argv[i],"-q") == 0)
	       quiet = 1;
	     if (strcmp(argv[i],"-x") == 0)
	       xanti = 1;
	     if (memcmp("-",argv[i],1) != 0)
	       strcpy(filename,argv[i]);
	  }
	
	init();
	
	outfile = fopen(gridname,"w");
	
	for (j=0; j<grid; ++j)
	  {
	     
	     if (!quiet) printf("Row %d...\n",j);
	     
	     for (i=0; i<grid; ++i)
	       {
		  field = 0.0;
		  for (k=0; k<N; ++k)
		    {
		       c = cos(theta[k]);
		       s = sin(theta[k]);
		       c2 = SQR(c);
		       s2 = SQR(s);
		       dx = (xbegin + i*xstep) - posx[k];
		       dy = (ybegin + j*ystep) - posy[k];
		       temp = -(SQR(dx)*(c2/a2[k]+s2*a2[k])+
				SQR(dy)*(s2/a2[k]+c2*a2[k])+
				2.0*dx*dy*c*s*(1.0/a2[k]-a2[k]))/
			 (4.0*var[k]);
		       /* A little code to avoid underflow handling. */
		       if (-temp<16)
			 field += 
			   amp[k]*exp(temp)/(2.0*var[k]);
		    }
		  fprintf(outfile,"%14.8e\n",field);
	       }
	  }
     }
}
