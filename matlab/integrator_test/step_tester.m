%%

visc = 1.0e-4;

axisymmtol = 1.0e-6;

step = 1.0e-1/4;

blob.x = 1;
blob.y = 0;
blob.s2 = 0.01;
blob.a2 = 1;
blob.th = 0;

blob=set_blob(blob);

%% RK4

for k=1:40
    
    if (abs(blob.a2-1/blob.a2)<axisymmtol)
        blob = rk4step_slave(blob,step,visc);
    else
        blob = rk4step(blob,step,visc);
    end
end

disp(blob.x^2+blob.y^2-1);


%% AB 4
    
blob4 = blob;

if (abs(blob4.a2-1/blob4.a2)<axisymmtol)
    blob3 = rk4step_slave(blob4,step,visc);
else
    blob3 = rk4step(blob4,step,visc);
end

if (abs(blob3.a2-1/blob3.a2)<axisymmtol)
    blob2 = rk4step_slave(blob3,step,visc);
else
    blob2 = rk4step(blob3,step,visc);
end

if (abs(blob2.a2-1/blob2.a2)<axisymmtol)
    blob1 = rk4step_slave(blob2,step,visc);
else
    blob1 = rk4step(blob2,step,visc);
end

for k=1:100-3
    if ( (abs(blob1.a2-1/blob1.a2)<axisymmtol) | ...
            (abs(blob2.a2-1/blob2.a2)<axisymmtol) | ...
            (abs(blob3.a2-1/blob3.a2)<axisymmtol) | ...
            (abs(blob4.a2-1/blob4.a2)<axisymmtol) )
        blob = ab4step_slave(blob1,blob2,blob3,blob4,step,visc);
    else
        blob = ab4step(blob1,blob2,blob3,blob4,step,visc);
    end
    
    blob4 = blob3;
    blob3 = blob2;
    blob2 = blob1;
    blob1 = blob;
end

disp(blob.x^2+blob.y^2-1);
