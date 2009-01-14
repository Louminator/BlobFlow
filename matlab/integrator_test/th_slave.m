function blob = th_slave(blob,dv)

if ((1/blob.a2-blob.a2)<0)
    blob.th = atan2(-dv(1,1)+sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
else
    blob.th = atan2(-dv(1,1)-sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);
end

% Add an epsilon correction.
epsilon = (1/blob.a2-blob.a2)/(1/blob.a2+blob.a2);

th1 = -(dv(2,1)-dv(1,2))/2/...
     (4*sin(blob.th)*cos(blob.th)*(dv(1,2)+dv(2,1))/2*(1+(dv(1,1)*2/(dv(1,2)+dv(2,1)))^2));

blob.th = blob.th + epsilon*th1 + 2/3*epsilon^3*th1^3;
 
blob = set_blob(blob);
