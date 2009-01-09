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

%% Build a reference solution.

visc = 1.0e-4;

axisymmtol = 1.0e-6;

% Set j to a big number like 8.

j = 10;

blob.x = 1;
blob.y = 0;
blob.s2 = 0.01;
blob.a2 = 1.0;
blob.th = 0;

blob=set_blob(blob);

N = 10*2^(j-1);
step = 1.0/N;

for k=1:N
    
    if (abs(blob.a2-1/blob.a2)<axisymmtol)
        blob = rk4step_slave(blob,step,visc);
    else
        blob = rk4step(blob,step,visc);
    end
end

refblob=blob;    

%% RK4

visc = 1.0e-4;

axisymmtol = 1.0e-6;

err_rk4 = [];
err_rk4_th = [];
err_rk4_a2 = [];
err_rk4_s2 = [];
err_rk4_pos = [];
err_rk4_w = [];

numpts = 8;

    a2 = max(refblob.a2,1/refblob.a2);
    x = linspace(refblob.x-4*sqrt(a2*refblob.s2),...
        refblob.x+4*sqrt(a2*refblob.s2),20);
    y = linspace(refblob.y-4*sqrt(a2*refblob.s2),...
        refblob.y+4*sqrt(a2*refblob.s2),20);
    dA = (x(2)-x(1))*(y(2)-y(1));
    [xa,ya] = meshgrid(x,y);
    tmpx = xa-refblob.x;
    tmpy = ya-refblob.y;
    tmp = -((tmpx*refblob.costh+tmpy*refblob.sinth).^2/refblob.a2 + ...
        (tmpy*refblob.costh-tmpx*refblob.sinth).^2*refblob.a2);
    refw = 1./(4*pi*refblob.s2).*exp(tmp./4./refblob.s2);


for j = 1:numpts
    
    blob.x = 1;
    blob.y = 0;
    blob.s2 = 0.01;
    blob.a2 = 1.0;
    blob.th = 0;
    
    blob=set_blob(blob);
    
    N = 10*2^(j-1);
    step = 1.0/N;
    
    for k=1:N
        
        if (abs(blob.a2-1/blob.a2)<axisymmtol)
            blob = rk4step_slave(blob,step,visc);
        else
            blob = rk4step(blob,step,visc);
        end
    end
    
    disp(j);
%     err_rk4(end+1) = 1-blob.x^2-blob.y^2;
%     
%     err_rk4_th(end+1) = sqrt((blob.th-refblob.th)^2);
%     err_rk4_s2(end+1) = sqrt((blob.s2-refblob.s2)^2);
%     if ( (refblob.a2-1)*(blob.a2-1) > 0)
%         err_rk4_a2(end+1) = ...
%             sqrt((sqrt(blob.cos2*blob.a2)-sqrt(refblob.cos2*refblob.a2))^2+...
%             (sqrt(blob.sin2*blob.a2)-sqrt(refblob.sin2*refblob.a2))^2);
%     else
%         err_rk4_a2(end+1) = ...
%             sqrt((sqrt(blob.sin2*1/blob.a2)-sqrt(refblob.cos2*refblob.a2))^2+...
%             (sqrt(blob.cos2*1/blob.a2)-sqrt(refblob.sin2*refblob.a2))^2);
%     end
    err_rk4_pos(end+1) = sqrt((blob.x-refblob.x)^2+(blob.y-refblob.y)^2);
    tmpx = xa-blob.x;
    tmpy = ya-blob.y;
    tmp = -((tmpx*blob.costh+tmpy*blob.sinth).^2/blob.a2 + ...
        (tmpy*blob.costh-tmpx*blob.sinth).^2*blob.a2);
    w = 1./(4*pi*blob.s2).*exp(tmp./4./blob.s2);
    err_rk4_w(end+1) = sqrt(sum(sum((w-refw).^2))*dA);
end

loglog(1./2.^(1:numpts),err_rk4_w,'o-');
grid on;

%% AB 4
    
visc = 1.0e-4;

axisymmtol = 1.0e-6;

err_ab4 = [];
err_ab4_th = [];
err_ab4_a2 = [];
err_ab4_s2 = [];
err_ab4_pos = [];

numpts = 8;

for j = 1:numpts
    
    blob.x = 1;
    blob.y = 0;
    blob.s2 = 0.01;
    blob.a2 = 1.0;
    blob.th = 0;
    
    blob=set_blob(blob);
    
    N = 10*2^(j-1);
    step = 1.0/N;
    
    
    blob4 = blob;
    [v4,dv4] = velset(blob4);

    if (abs(blob4.a2-1/blob4.a2)<axisymmtol)
        blob4 = th_slave(blob4,dv4);
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
    
    for k=1:N-3
        if (    (abs(blob1.a2-1/blob1.a2)<axisymmtol) | ...
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
    
    disp(j);
    err_ab4(end+1) = 1-blob.x^2-blob.y^2;
    
    err_ab4_th(end+1) = sqrt((blob.th-refblob.th)^2);
    err_ab4_s2(end+1) = sqrt((blob.s2-refblob.s2)^2);
    if ( (refblob.a2-1)*(blob.a2-1) > 0)
        err_ab4_a2(end+1) = ...
            sqrt((sqrt(blob.cos2*blob.a2)-sqrt(refblob.cos2*refblob.a2))^2+...
            (sqrt(blob.sin2*blob.a2)-sqrt(refblob.sin2*refblob.a2))^2);
    else
        err_ab4_a2(end+1) = ...
            sqrt((sqrt(blob.sin2*1/blob.a2)-sqrt(refblob.cos2*refblob.a2))^2+...
            (sqrt(blob.cos2*1/blob.a2)-sqrt(refblob.sin2*refblob.a2))^2);
    end
    err_ab4_pos(end+1) = sqrt((blob.x-refblob.x)^2+(blob.y-refblob.y)^2);
    
end

loglog(1./2.^(1:numpts),err_ab4,'o-');
grid on;
