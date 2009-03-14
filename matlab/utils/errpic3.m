function [w,ex,dx,dy] = errpic3(name,fileno,time)

figno  = fileno;

[fid,message] = fopen('egrid.default','r');
if fid == -1
disp(message)
end

grdname = sprintf('%s%04d.grd',name,fileno);

xmin = fscanf(fid,'%lf',1);
ymin = fscanf(fid,'%lf',1);
xmax = fscanf(fid,'%lf',1);
ymax = fscanf(fid,'%lf',1);
n    = fscanf(fid,'%d',1);

dx = (xmax-xmin)/n;
dy = (ymax-ymin)/n;

fclose(fid);

[fid,message] = fopen(grdname,'r');
if fid == -1
disp(message)
end

for i = 1:n
tmp(:,i) = fscanf(fid,'%lf',n);
end
fclose(fid);
w = tmp.';

[xa,ya] = meshgrid(xmin:(xmax-xmin)/(n-1):xmax,ymin:(ymax-ymin)/(n-1):ymax);

R = xa.^2+ya.^2;

ex = (1/4)*(1/((0.25)^2+time))*exp(-R/(4.0*((0.25)^2+time)));

err = (w-ex).*(R<1);



figure;

%%subplot(1,3,1);

mesh(xa,ya,err);

%%subplot(1,3,2);

%%surf(xa,ya,((err)./ex).*(R<1));

%%subplot(1,3,3);

%%mesh(xa,ya,w);

disp('sup norm (rel):');
disp(max(max(abs(err)))*4*((0.25)^2+time));

disp('sup norm (abs):');
disp(max(max(abs(err))));

disp('l2 norm:');
disp(sqrt(dx*dy*sum(sum((err).^2))));
