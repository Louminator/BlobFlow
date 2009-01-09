function blob = th_slave(blob,dv)

if ((1/blob.a2-blob.a2)<0)
    blob.th = atan2(-dv(1,1)+sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
else
    blob.th = atan2(-dv(1,1)-sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
end

% Add an epsilon correction.
blob.th = blob.th - (1/blob.a2-blob.a2)/(1/blob.a2+blob.a2)*(dv(2,1)-dv(1,2))/2/...
     (4*sin(blob.th)*cos(blob.th)*(dv(1,2)+dv(2,1))/2*(1+(dv(1,1)*2/(dv(1,2)+dv(2,1)))^2));
 
blob = set_blob(blob);
