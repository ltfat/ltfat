function [coef]=comp_dwiltiii(f,g,M)
%COMP_DWILTIII  Compute Discrete Wilson transform type III.
%   

L=size(f,1);
a=M;
N=L/a;
W=size(f,2);

coef=zeros(M,N,W,assert_classname(f,g));

fwasreal=isreal(f);

if 0
  coef2 = comp_gdgt(f,g,a,2*M,L,0,.5,0,0);
else
  halfmod=exp(-pi*i*(0:L-1).'/(2*M));
  f=f.*repmat(halfmod,1,W);
  coef2=comp_dgt(f,g,a,2*M,[0 1],0,0,0);
  

end;
  
% --- m is even ---------
coef(1:2:M,1:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(1:2:M,1:2:N,:)+exp(i*pi/4)*coef2(2*M:-2:M+1,1:2:N,:));
coef(1:2:M,2:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(1:2:M,2:2:N,:)+exp(-i*pi/4)*coef2(2*M:-2:M+1,2:2:N,:));

% --- m is odd ----------
coef(2:2:M,1:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(2:2:M,1:2:N,:)+exp(-i*pi/4)*coef2(2*M-1:-2:M+1,1:2:N,:));
coef(2:2:M,2:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(2:2:M,2:2:N,:)+exp(i*pi/4)*coef2(2*M-1:-2:M+1,2:2:N,:));

coef=reshape(coef,M,N,W);

% Clean coefficients if they are known to be real
% If input is real, and window is real, output must be real as well.
if (isreal(g) && fwasreal);
  coef=real(coef);
end;


