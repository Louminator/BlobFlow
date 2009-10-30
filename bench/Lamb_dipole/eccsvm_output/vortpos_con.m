function [w,levels] = vortpos(name,fileno,l)

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

minw = min(min(w));
maxw = max(max(w));

if (nargin == 2)
    levels = minw:(maxw-minw)/19:maxw;
else
    levels = l;
end

contourf(xa,ya,w,levels);

