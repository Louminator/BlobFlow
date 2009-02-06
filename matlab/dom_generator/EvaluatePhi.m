function [Phi]=EvaluatePhi(x,y,a)
% [Phi]=evaluatePhi(x,y,a)
% Evalates Phi at points given in the vectors x and y

xx=x(:);
yy=y(:);
Phi=zeros(length(xx),1);
for i=1:length(xx)
    [Phi(i)]=StreamIntegral(xx(i),yy(i),a);
end
Phi=reshape(Phi,size(x));