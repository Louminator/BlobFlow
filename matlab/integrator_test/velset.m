function [v,dv] = velset(blob)

x=blob.x;
y=blob.y;

r2 = x^2 + y^2;

r = sqrt(r2);

% Define f() an f'() where
%
% v = f(r^2) [-y/r ; x/y] = f(r^2) theta-hat

f = r2;

fp = 1;

v = f*[-y/r ; x/r];

dv = zeros(2,2);

dv(1,1) = -fp*2*x*y/r + f*x*y/(r^(3/2));

dv(1,2) = -fp*2*y^2/r - f/r + f*y^2/(r^(3/2));

dv(2,1) = fp*2*x^2/r + f/r - f*x^2/(r^(3/2));

dv(2,2) = fp*2*x*y/r - f*x*y/(r^(3/2));