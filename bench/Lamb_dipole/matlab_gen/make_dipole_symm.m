clear

% Domain (=[-1,1]^2)
N= 50;
x0 = -1;
x1 =  1;

y0 = 0;
y1 = 2;

h = (x1-x0)/N
var = h*h

x = x0+((1:N)-0.5)*h;
y = y0+((1:N)-0.5)*h;

% Dipole radius (=1).
R = 1;
k = fminsearch(@(s) besselj(1,s*R).^2,3.8/R);

% Dipole velocity (=1).
U = 1;
C = 2*U/k/besselj(0,k*R);

[xa,ya] = meshgrid(x,y);

ra = sqrt(xa.^2+ya.^2);

wa = - C*k^2*(besselj(1,k*ra).*ya./ra).*(ra<=R);

% Clean up pesky problems at the origin.
for k=1:N
    for l=1:N
        if (isnan(wa(k,l)))
                wa(k,l) = 0.0;
        end
    end
end

surf(xa,ya,wa);

w = wa';
w = w(:);
save w.grd w -ASCII