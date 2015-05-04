function [coef]=comp_dwiltiii(f,g,M)
%COMP_DWILTIII  Compute Discrete Wilson transform type III.
%   

L = size(f,1);
a=M;
N=L/a;
W=size(f,2);

fwasreal = isreal(f);
coef=zeros(M,N,W,assert_classname(f,g));

halfmod=exp(-pi*i*(0:L-1).'/(2*M));
f=f.*repmat(halfmod,1,W);

coef2=comp_sepdgt(f,g,a,2*M,0);
  
if (isreal(g) && fwasreal)
   % --- m is even ---------
   coef(1:2:M,1:2:N,:)= real(coef2(1:2:M,1:2:N,:)) + imag(coef2(1:2:M,1:2:N,:));
   coef(1:2:M,2:2:N,:)= real(coef2(1:2:M,2:2:N,:)) - imag(coef2(1:2:M,2:2:N,:));

   % --- m is odd ----------
   coef(2:2:M,1:2:N,:)= real(coef2(2:2:M,1:2:N,:)) - imag(coef2(2:2:M,1:2:N,:));
   coef(2:2:M,2:2:N,:)= real(coef2(2:2:M,2:2:N,:)) + imag(coef2(2:2:M,2:2:N,:));    
else
   % --- m is even ---------
   coef(1:2:M,1:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(1:2:M,1:2:N,:)+exp(i*pi/4)*coef2(2*M:-2:M+1,1:2:N,:));
   coef(1:2:M,2:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(1:2:M,2:2:N,:)+exp(-i*pi/4)*coef2(2*M:-2:M+1,2:2:N,:));

   % --- m is odd ----------
   coef(2:2:M,1:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(2:2:M,1:2:N,:)+exp(-i*pi/4)*coef2(2*M-1:-2:M+1,1:2:N,:));
   coef(2:2:M,2:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(2:2:M,2:2:N,:)+exp(i*pi/4)*coef2(2*M-1:-2:M+1,2:2:N,:));
end
