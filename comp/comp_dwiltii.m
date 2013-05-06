function [coef]=comp_dwiltii(coef2,a)
%COMP_DWILT  Compute Discrete Wilson transform.
%   

M=size(coef2,1)/2;
N=size(coef2,2);
W=size(coef2,3);
L=N*a;

coef=zeros(2*M,N/2,W,assert_classname(coef2));


% ---- m is zero ---------
coef(1,:,:)=coef2(1,1:2:N,:);

% --- m is odd ----------
coef(2:2:M,:,:)    = i/sqrt(2)*(coef2(2:2:M,1:2:N,:)+coef2(2*M:-2:M+2,1:2:N,:));
coef(M+2:2:2*M,:,:)= 1/sqrt(2)*(coef2(2:2:M,2:2:N,:)-coef2(2*M:-2:M+2,2:2:N,:));

% --- m is even ---------
coef(3:2:M,:,:)=     1/sqrt(2)*(coef2(3:2:M,1:2:N,:)-coef2(2*M-1:-2:M+2,1:2:N,:));
coef(M+3:2:2*M,:,:)= i/sqrt(2)*(coef2(3:2:M,2:2:N,:)+coef2(2*M-1:-2:M+2,2:2:N,:));

% --- m is nyquest ------
if mod(M,2)==0
  coef(M+1,:,:) = i*coef2(M+1,2:2:N,:);
else
  coef(M+1,:,:) = i*coef2(M+1,1:2:N,:);
end;

coef=reshape(coef,M*N,W);




