#include <math.h>
#define SQR(X) ((X)*(X))

void find_coeffs (double a, double coeffs[289], double Md[289][700])
{
 
	int j,k;
	double c1,c2,a1,a2;
	
	if (a>10)
	{
		a=10;
	}
	
	j=(int) floor(699*sqrt((a-1.+1e-15)/9));
	a1=SQR(j*sqrt(9.0)/699)+1;
	a2=SQR((j+1)*sqrt(9.0)/699)+1;
	c2=(a-a1)/(a2-a1);
	c1=(a2-a)/(a2-a1);
	
	for (k=0; k<289; k++)
		coeffs[k]=c1*Md[k][j]+c2*Md[k][j+1];
	
}


