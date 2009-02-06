%%% Data generator for Biot-Savart evaluation code
%%% Polynomial coefficiets are computed through interpolation
%%% and saved for future used. 
%%%
%%% More details can  be found in section 3 of
%%% Platte, Rossi, Mitchell, Using global interpolation to evaluate the
%%% Biot-Savart integral for deformable elliptical Gaussian vortex
%%% elements.

clear all

%%% Core aspect ratio samples %%%
aa=1+linspace(0,sqrt(9),700).^2;

Dom1=[]; Dom2=[];Dom3=[];Dom4=[];Dom5=[];Dom6=[];Dom7=[];Dom8=[];Dom9=[];
Dom10=[]; Dom11=[]; Dom12=[]; Dom13=[];

%%% Generate grid points on [-1,1]x[-1,1]%%%
%%% Chebyshev nodes are used %%%%
Nx=16; Ny=16;
x1=-cos((0:Nx)*pi/(Nx)); 
y1=-cos((0:Ny)*pi/(Ny)); 
[xxp,yyp]=meshgrid(x1,y1);
xx=xxp(:); yy=yyp(:);

%%% Generate interpolation matrix %%%%
%%% Monomials are used -- restriceted to low degree %%
A=[];
for nx=1:length(x1)
    for ny=1:length(y1)
        A(:,end+1) =xx.^(nx-1).*yy.^(ny-1);
    end
end

%%% Interpolate on Domains 1 - 13 %%%
%%% Domains are scaled and shifted to [-1,1]x[-1,1]%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=1, a=', num2str(a)])
    ymin=0; ymax=0.5;
    xmin=0; xmax=5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    %pause
    Dom1(:,end+1)=A\b;    
end

save 'Dom1.dat' Dom1 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=2, a=', num2str(a)])
    ymin=0.5; ymax=1;
    xmin=0; xmax=5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom2(:,end+1)=A\b;    
end

save 'Dom2.dat' Dom2 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=3, a=', num2str(a)])
    ymin=1; ymax=3;
    xmin=0; xmax=5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom3(:,end+1)=A\b;    
end

save 'Dom3.dat' Dom3 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=4, a=', num2str(a)])
    ymin=3; ymax=7*sqrt(a);
    xmin=0; xmax=5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom4(:,end+1)=A\b;    
end

save 'Dom4.dat' Dom4 -ascii -double 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=5, a=', num2str(a)])
    ymin=7*sqrt(a); ymax=11.5*a;
    xmin=0; xmax=5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom5(:,end+1)=A\b;    
end

 save 'Dom5.dat' Dom5 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=6, a=', num2str(a)])
    ymin=0; ymax=0.5;
    xmin=5*a; xmax=11.5*a;
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom6(:,end+1)=A\b;    
end

save 'Dom6.dat' Dom6 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=7, a=', num2str(a)])
    ymin=0.5; ymax=1;
    xmin=5*a; xmax=11.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom7(:,end+1)=A\b;    
end

save 'Dom7.dat' Dom7 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=8, a=', num2str(a)])
    ymin=1; ymax=3;
    xmin=5*a; xmax=11.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom8(:,end+1)=A\b;    
end

save 'Dom8.dat' Dom8 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=9, a=', num2str(a)])
    ymin=3; ymax=7*sqrt(a);
    xmin=5*a; xmax=11.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom9(:,end+1)=A\b;    
end

save 'Dom9.dat' Dom9 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for na=1:length(aa)
    a=aa(na);
    disp(['Dom=10, a=', num2str(a)])
    ymin=7*sqrt(a); ymax=11.5*a;
    xmin=5*a; xmax=11.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom10(:,end+1)=A\b;    
end

save 'Dom10.dat' Dom10 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for na=1:length(aa)
    a=aa(na);
    disp(['Dom=11, a=', num2str(a)])
    ymin=0; ymax=11.5*a;
    xmin=11.5*a; xmax=40.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom11(:,end+1)=A\b;    
end

save 'Dom11.dat' Dom11 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for na=1:length(aa)
    a=aa(na);
    disp(['Dom=12, a=', num2str(a)])
    ymin=11.5*a; ymax=40.5*a;
    xmin=0; xmax=11.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom12(:,end+1)=A\b;    
end

save 'Dom12.dat' Dom12 -ascii -double

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for na=1:length(aa)
    a=aa(na);
    disp(['Dom=13, a=', num2str(a)])
    ymin=11.5*a; ymax=40.5*a;
    xmin=11.5*a; xmax=40.5*a; 
    b=EvaluatePhi((xmax-xmin)*xx/2+(xmax+xmin)/2, (ymax-ymin)*yy/2+(ymax+ymin)/2, a);
    Dom13(:,end+1)=A\b;    
end

save 'Dom13.dat' Dom13 -ascii -double
