%% Creation cell.

% User supplied basis function core size for computations.
sigma = 0.05/2/sqrt(2);

% Base symmetric MZM disk.
Ri = 0.0;
Ro = 1.6/sqrt(2);

% Aspect ratio of the deformed disk.

aspect = 2.0;

a = sqrt(aspect);
b = 1/sqrt(aspect);
X = sigma*(ceil(Ro*a/sigma)+2)
Y = sigma*(ceil(Ro*b/sigma)+2);

if (X>Y)
    Y=X;
else
    X=Y;
end

x = -X:sigma:X;
size(x)
y = -Y:sigma:Y;

% Create the desired function.
[xa,ya,w] = MZM(x,y,Ri,Ro,a,b);
w = 20*w(:);
%subplot(1,2,1);
surf(xa,ya,reshape(w,size(xa)));
% Make sure the entire function is within the specified domain.

save w.grd w -ASCII
