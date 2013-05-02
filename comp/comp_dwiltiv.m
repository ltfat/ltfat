function [coef]=comp_dwiltiv(coef2,a)
%COMP_DWILTIV  Compute Discrete Wilson transform type IV.
%   

M=size(coef2,1)/2;
N=size(coef2,2);
W=size(coef2,3);
L=N*a;

coef=zeros(M,N,W,assert_classname(coef2));

% --- m is even ---------
coef(1:2:M,1:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(1:2:M,1:2:N,:)+exp(-i*pi*3/4)*coef2(2*M:-2:M+1,1:2:N,:));
coef(1:2:M,2:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(1:2:M,2:2:N,:)+exp(i*pi*3/4)*coef2(2*M:-2:M+1,2:2:N,:));

% --- m is odd ----------
coef(2:2:M,1:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(2:2:M,1:2:N,:)+exp(i*pi*3/4)*coef2(2*M-1:-2:M+1,1:2:N,:));
coef(2:2:M,2:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(2:2:M,2:2:N,:)+exp(-i*pi*3/4)*coef2(2*M-1:-2:M+1,2:2:N,:));

coef=reshape(coef,M*N,W);

