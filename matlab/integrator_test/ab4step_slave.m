function blob = ab4step_slave(blob1,blob2,blob3,blob4,step,visc)

[v1,dv1] = velset(blob1);
[v2,dv2] = velset(blob2);
[v3,dv3] = velset(blob3);
[v4,dv4] = velset(blob4);

blob = blob1;

blob.x = blob1.x + step/24*...
    (55*v1(1)-59*v2(1)+37*v3(1)-9*v4(1));
blob.y = blob1.y + step/24*...
    (55*v1(2)-59*v2(2)+37*v3(2)-9*v4(2));
blob.s2 = blob1.s2 + step/24*...
    (55*ds2dt(blob1,visc)-59*ds2dt(blob2,visc)+...
    37*ds2dt(blob3,visc)-9*ds2dt(blob4,visc));
blob.a2 = blob1.a2 + step/24*...
    (55*da2dt(blob1,dv1,visc)-59*da2dt(blob2,dv2,visc)+...
    37*da2dt(blob3,dv3,visc)-9*da2dt(blob4,dv4,visc));

[v,dv] = velset(blob);
blob = th_slave(blob,dv);

