function blob = th_slave(blob,dv)

temp = dv(1,1)*(blob.a2+1/blob.a2);

arg = max([0 ...
    (dv(2,1)/blob.a2+dv(1,2)*blob.a2)*...
    (dv(1,2)/blob.a2+dv(2,1)*blob.a2)]);

th = atan2((dv(1,2)*blob.a2+dv(2,1)/blob.a2),-(temp+sqrt(arg)));

blob = set_blob(blob);

if ((blob.a2-1/blob.a2)*...
        ((blob.sin2-blob.cos2)*dv(1,1)-(dv(1,2)+dv(2,1))*...
        blob.sincos)>0)
    th = atan2((dv(1,2)*blob.a2+dv(2,1)/blob.a2),-(temp-sqrt(arg)));
    
end

th = atan2(-dv(1,1)+sqrt(dv(1,1)^2+(dv(1,2)+dv(2,1))^2/4),(dv(1,2)+dv(2,1))/2);

blob = set_blob(blob);