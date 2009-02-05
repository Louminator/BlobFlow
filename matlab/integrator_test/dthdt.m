function dthdt = dthdt(blob,dv,axisymmtol,step)

dthdt = (dv(2,1)-dv(1,2))/2 + ...
    ( (dv(2,1)+dv(1,2))/2*(blob.sin2-blob.cos2) + ...
    2*dv(1,1)*blob.sincos )*...
    (1/blob.a2+blob.a2)/(1/blob.a2-blob.a2);

if (abs(blob.a2-1/blob.a2)<axisymmtol*step)
    dthdt = (dv(2,1)-dv(1,2))/2;
end