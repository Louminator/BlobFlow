function blob4 = rk4step_slave(blob,step,visc)

blob = set_blob(blob);

[v,dv] = velset(blob);

% FCE

blob1 = blob;

blob = th_slave(blob,dv);
blob1.a2 = blob.a2 + 0.5*step*da2dt(blob,dv,visc);
blob1.s2 = blob.s2 + 0.5*step*ds2dt(blob,visc);
blob1.x = blob.x + 0.5*step*v(1);
blob1.y = blob.y + 0.5*step*v(2);

blob1 = set_blob(blob1);

% BCE

[v1,dv1] = velset(blob1);

blob2 = blob1;

blob1 = th_slave(blob1,dv1);
blob2.a2 = blob.a2 + 0.5*step*da2dt(blob1,dv1,visc);
blob2.s2 = blob.s2 + 0.5*step*ds2dt(blob1,visc);
blob2.x = blob.x + 0.5*step*v1(1);
blob2.y = blob.y + 0.5*step*v1(2);

blob2 = set_blob(blob2);

% Midpoint rule

[v2,dv2] = velset(blob2);

blob3 = blob2;

blob2 = th_slave(blob2,dv2);
blob3.a2 = blob.a2 + step*da2dt(blob2,dv2,visc);
blob3.s2 = blob.s2 + step*ds2dt(blob2,visc);
blob3.x = blob.x + step*v2(1);
blob3.y = blob.y + step*v2(2);

blob3 = set_blob(blob3);

% Simpson's Rule

[v3,dv3] = velset(blob3);

blob4 = blob3;

blob3 = th_slave(blob3,dv3);

blob4.a2 = blob.a2 + step*...
    (da2dt(blob,dv,visc) + 2*da2dt(blob1,dv1,visc)+ ...
    2*da2dt(blob2,dv2,visc) + da2dt(blob3,dv3,visc))/6;

blob4.s2 = blob.s2 + step*...
    (ds2dt(blob,visc) + 2*ds2dt(blob1,visc)+ ...
    2*ds2dt(blob2,visc) + ds2dt(blob3,visc))/6;

blob4.x = blob.x + step*(v(1) + 2*v1(1) + 2*v2(1) + v3(1))/6;
blob4.y = blob.y + step*(v(2) + 2*v1(2) + 2*v2(2) + v3(2))/6;

blob4 = set_blob(blob4);
[v4,dv4] = velset(blob4);
blob4 = th_slave(blob4,dv4);


