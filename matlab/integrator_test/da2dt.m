function da2dt = da2dt(blob,dv,visc)

da2dt = 2*(dv(1,1)*(blob.cos2-blob.sin2) + ...
    (dv(1,2)+dv(2,1))*blob.sincos)*blob.a2 + ...
    (visc/(2*blob.s2))*(1-blob.a2^2);