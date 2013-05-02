function [coef2]=comp_idwiltiv(coef,a,M)
%COMP_IDWILTIV  Compute Inverse discrete Wilson transform type IV.
% 
%   This is a computational routine. Do not call it
%   directly.

%   AUTHOR : Peter L. Søndergaard

N=size(coef,1)/M;
W=size(coef,2);
L=N*a;

coef=reshape(coef,M,N,W);

coef2=zeros(2*M,N,W,assert_classname(coef));

coef2(1:2:M,1:2:N,:)        = exp(i*pi/4)/sqrt(2)*coef(1:2:M,1:2:N,:);
coef2(2*M:-2:M+1,1:2:N,:)   = exp(i*pi*3/4)/sqrt(2)*coef(1:2:M,1:2:N,:);

coef2(1:2:M,2:2:N,:)        = exp(-i*pi/4)/sqrt(2)*coef(1:2:M,2:2:N,:);
coef2(2*M:-2:M+1,2:2:N,:)   = exp(-i*pi*3/4)/sqrt(2)*coef(1:2:M,2:2:N,:);

coef2(2:2:M,1:2:N,:)        = exp(-i*pi/4)/sqrt(2)*coef(2:2:M,1:2:N,:);
coef2(2*M-1:-2:M+1,1:2:N,:) = exp(-i*pi*3/4)/sqrt(2)*coef(2:2:M,1:2:N,:);

coef2(2:2:M,2:2:N,:)        = exp(i*pi/4)/sqrt(2)*coef(2:2:M,2:2:N,:);
coef2(2*M-1:-2:M+1,2:2:N,:) = exp(i*pi*3/4)/sqrt(2)*coef(2:2:M,2:2:N,:);


