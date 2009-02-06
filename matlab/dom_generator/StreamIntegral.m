function [phi]=StreamIntegral(xs,ys,a)

Rs=sqrt(xs^2/a^2+ys^2*a^2);

if Rs<1e-14
    xs=1e-14;  ys=1e-14;  Rs=1e-14;
end

phi2=0.5*( (xs^2/a+ys^2*a)/(a+1/a) );

sstep=1e-10;

if Rs<.2
   phi = - quadl(@(R)integrand(R,xs,ys,a), 1e-15, Rs, 1e-12);
elseif Rs<=5.5
   phi = - quadl(@(R)integrand(R,xs,ys,a), 1e-15, .1, 1e-12) ...
         - quadl(@(R)integrand(R,xs,ys,a), .1, Rs, 1e-12);
elseif Rs<=11
   phi = - quadl(@(R)integrand(R,xs,ys,a), 1e-15, .1, 1e-12) ...
         - quadl(@(R)integrand(R,xs,ys,a), .1, 5 , 1e-12) ...
         - quadl(@(R)integrand(R,xs,ys,a), 5 , Rs, 1e-12);
else
phi = - quadl(@(R)integrand(R,xs,ys,a), 1e-15, .1, 1e-12) ...
      - quadl(@(R)integrand(R,xs,ys,a), .1, 5 , 1e-12) ...
      - quadl(@(R)integrand(R,xs,ys,a), 5 , 10, 1e-12) ...
      - quadl(@(R)integrand(R,xs,ys,a), 10 , Rs, 1e-12);  
end

phi=phi+phi2*exp(-Rs^2/4)/(4*pi); 
    
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=integrand(R,xs,ys,a)

rhos=sqrt(xs^2+ys^2);

Rs=sqrt(xs^2/a^2+ys^2*a^2);

xi=0.5*( rhos^2-R.^2.*(a^2+1/a^2)+sqrt( (R.^2*(a^2+1/a^2)-rhos^2).^2+4*R.^2.*(Rs^2-R.^2) ) );

beta=sqrt(R.^2/a^2+xi);

alpha=sqrt(R.^2*a^2+xi);

phi1=0.5.*( (xs^2./alpha + ys^2./beta)./(alpha + beta) + log( (alpha + beta)./ (R.*a+R./a) ) );

y=phi1.*(-R.^3.*exp(-R.^2/4)/(8*pi));