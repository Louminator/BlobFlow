function blob = th_slave(blob,dv)

if (dv(1,1)*(1/blob.a2-blob.a2)<0)
    th = atan2(-dv(1,1)+sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
else
    th = atan2(-dv(1,1)-sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
end
    
% Add an epsilon correction.
th = th - (1/blob.a2-blob.a2)/(1/blob.a2+blob.a2)*(dv(2,1)-dv(1,2))/2/...
    (4*sin(th)*cos(th)*(dv(1,2)+dv(2,1))/2*(1+(dv(1,1)*2/(dv(1,2)+dv(2,1)))^2));

blob = set_blob(blob);