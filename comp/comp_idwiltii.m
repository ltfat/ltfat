function [coef2]=comp_idwiltii(coef,a,M)
%COMP_IDWILTII  Compute Inverse discrete Wilson transform type II
% 
%   This is a computational routine. Do not call it
%   directly.

%   AUTHOR : Peter L. Søndergaard

N=size(coef,1)/M;
W=size(coef,2);
L=N*a;

coef=reshape(coef,M*2,N/2,W);

coef2=zeros(2*M,N,W,assert_classname(coef));

% First and middle modulation are transferred unchanged.
coef2(1,1:2:N,:) = coef(1,:,:);

coef2(2:2:M,1:2:N,:)        = -i/sqrt(2)*coef(2:2:M,:,:);
coef2(2*M:-2:M+2,1:2:N,:)   = -i/sqrt(2)*coef(2:2:M,:,:);

coef2(2:2:M,2:2:N,:)        =  1/sqrt(2)*coef(M+2:2:2*M,:,:);
coef2(2*M:-2:M+2,2:2:N,:)   = -1/sqrt(2)*coef(M+2:2:2*M,:,:);

if M>2
  coef2(3:2:M,1:2:N,:)        = 1/sqrt(2)*coef(3:2:M,:,:);
  coef2(2*M-1:-2:M+2,1:2:N,:) = -1/sqrt(2)*coef(3:2:M,:,:);
  
  coef2(3:2:M,2:2:N,:)        = -i/sqrt(2)*coef(M+3:2:2*M,:,:);
  coef2(2*M-1:-2:M+2,2:2:N,:) = -i/sqrt(2)*coef(M+3:2:2*M,:,:);
end;

if mod(M,2)==0
  coef2(M+1,2:2:N,:) = -i*coef(M+1,:,:);
else
  coef2(M+1,1:2:N,:) = -i*coef(M+1,:,:);
end;







