function [xa,ya,y] = MZM(x,y,Ri,Ro,a,b)

[xa,ya] = meshgrid(x,y);

dx = x(2)-x(1);
dy = y(2)-y(1);

y = w(xa/a,ya/b,Ri,Ro);
%y = y(:)/dx/dy/sum(sum(y));

function out = w(x,y,Ri,Ro);

out = 0*x;
r = sqrt(x.^2+y.^2);

r = r(:);
out = out(:);

idx = find(r<=Ri);
out(idx) = 1;

idx = find(r>=Ro);
out(idx) = 0;

idx = find(r>Ri & r<Ro);
out(idx) = 1-fk((r(idx)-Ri)/(Ro-Ri));

out = reshape(out,size(x));

function y = fk(r)

kappa = 1/2*log(2)*exp(2);
y = exp(-kappa./r.*exp(1./(r-1)));

