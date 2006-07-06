#include <math.h>

// -----------------------------------------------------------------------
//  Chebyshev polynomial evaluation
// -----------------------------------------------------------------------
// Returns the value of the polynomial (or its derivative)
// double z = z value (i.e. evaluation point)
// int d = degree of the polynomial
// int nder = order of the deirvative

// by Rodrigo Platte 06/14/2006

double Cheb(double z, double acosz, double sqrt1z, int d, int nder)

{
  switch (nder)
    {
    case 0: // dth Chebyshev polynomial evaluated at z
      return cos(d*acosz);
      break;

    case 1: // 1st deriv of dth Chebyshev polynomial evaluated at z
      return sin(d*acosz)*d/(sqrt1z);
      break;

    case 2: // 2nd deriv of dth Chebyshev polynomial evaluated at z
      return (sin(d*acosz)*z-cos(d*acosz)*d*sqrt1z)*d/((1-z*z)*sqrt1z);
      break;

    case 3: // 3rd deriv of dth Chebyshev polynomial evaluated at z
      return ( ((d*d+2)*z*z-d*d+1)*d*sin(d*acosz)
	       -3*cos(d*acosz)*d*d*sqrt1z*z ) / ((1-z*z)*(1-z*z)*sqrt1z);
      break;

    }
}

//------------------------------------------------------------------------
