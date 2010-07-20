clear

% Domain (=[-1,1]^2)
N= 150;
x0 = -1;
x1 =  1;

h = (x1-x0)/N
var = h*h

x = x0+((1:N)-0.5)*h;
y = x;

[xa,ya] = meshgrid(x,y);

ra2 = xa.^2+ya.^2;
ra2A = (xa-0.5).^2+ya.^2;
ra2B = (xa+0.5).^2+ya.^2;

sig = 0.06;

wa = -pi*(1/2)^2*(20)*exp(-ra2/4/sig^2)/4/pi/sig^2 + ...
    pi*(1/2)^2*(40)*exp(-ra2A/4/sig^2)/4/pi/sig^2 + ...
    pi*(1/2)^2*(40)*exp(-ra2B/4/sig^2)/4/pi/sig^2;

% Clean up pesky problems at the origin.
for k=1:N
    for l=1:N
        if (isnan(wa(k,l)))
                wa(k,l) = 0.0;
        end
    end
end

contourf(xa,ya,wa,30);
xlabel('x'); ylabel('y'); title('\omega');
colorbar;

w = wa';
w = w(:);
save w.grd w -ASCII
