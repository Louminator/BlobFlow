function strfn()

xmin = -1.25;
ymin = -1.25;
xmax = 1.25;
ymax = 1.25;
n    = 100;

xi = 0.1;
xi1 = 3*xi;
xi2 = -3.5*xi;
xi3 = xi;

xi4 =  6*xi;
xi5 = -4*xi;

dx = (xmax-xmin)/n;
dy = (ymax-ymin)/n;

[xa,ya] = meshgrid(xmin:(xmax-xmin)/(n-1):xmax,ymin:(ymax-ymin)/(n-1):ymax);

R = xa.^2+ya.^2;

psi = pi*R/4 + xi1*R + xi2*R.^2 + xi3*R.^3 + xi4*xa.*ya;

contour(xa,ya,psi,20,'r');
hold on;


psi = pi*R/4 + xi1*R + xi2*R.^2 + xi3*R.^3;

contour(xa,ya,psi,20,'b');
hold off;


