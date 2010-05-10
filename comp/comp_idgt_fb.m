function [f]=comp_idgt_fb(coef,g,L,a,M)
%COMP_IDGT_FB  Filter bank IDGT.
%   Usage:  f=comp_idgt_fb(c,g,L,a,M);
%       
%   This is a computational routine. Do not call it directly.
%
%   See also: idgt

%   AUTHOR : Peter Soendergaard.
%   TESTING: OK
%   REFERENCE: OK

% Calculate the parameters that was not specified.
N=L/a;
b=L/M;

R=size(g,2);

W=prod(size(coef))/(M*N*R);

N=L/a;
b=L/M;

gl=length(g);
glh=floor(gl/2);  % gl-half

% Apply ifft to the coefficients.
coef=ifft(reshape(coef,M,N*W))*sqrt(M);

coef=reshape(coef,M,N,W);

% The fftshift actually makes some things easier.
g=fftshift(g);

f=zeros(L,W);

% Make multicolumn g by replication.
gw=repmat(g,1,W);

% Shift the coefficients.
for n=0:N-1
  coef(:,n+1,:)=circshift(coef(:,n+1,:),-n*a+glh);
end;

% ----- Handle the first boundary using periodic boundary conditions. ---
for n=0:ceil(glh/a)-1
  
  % Rotate the coefficients, duplicate them until they have same
  % length as g, and multiply by g.
  %ff=repmat(circshift(squeeze(coef(:,n+1,:)),-n*a+glh),gl/M,1).*gw;
  ff=repmat(squeeze(coef(:,n+1,:)),gl/M,1).*gw;
  
  sp=mod(n*a-glh,L);
  ep=mod(n*a-glh+gl-1,L);

  % Add the ff vector to f at position sp.
  f(sp+1:L,:)=f(sp+1:L,:)+ff(1:L-sp,:);
  f(1:ep+1,:)=f(1:ep+1,:)+ff(L-sp+1:gl,:);
    
end;

% ----- Handle the middle case. ---------------------
for n=ceil(glh/a):floor((L-ceil(gl/2))/a)
    
  % Rotate the coefficients, duplicate them until they have same
  % length as g, and multiply by g.
  %ff=repmat(circshift(squeeze(coef(:,n+1,:)),-n*a+glh),gl/M,1).*gw;
  ff=repmat(squeeze(coef(:,n+1,:)),gl/M,1).*gw;

  sp=mod(n*a-glh,L);
  ep=mod(n*a-glh+gl-1,L);
  
  % Add the ff vector to f at position sp.
  f(sp+1:ep+1,:)=f(sp+1:ep+1,:)+ff;
  
end;

% ----- Handle the last boundary using periodic boundary conditions. ---
% This n is one-indexed, to avoid to many +1
for n=floor((L-ceil(gl/2))/a)+1:N-1
  
  % Rotate the coefficients, duplicate them until they have same
  % length as g, and multiply by g.
  %ff=repmat(circshift(squeeze(coef(:,n+1,:)),-n*a+glh),gl/M,1).*gw;
  ff=repmat(squeeze(coef(:,n+1,:)),gl/M,1).*gw;
 
  sp=mod(n*a-glh,L);
  ep=mod(n*a-glh+gl-1,L);
 
  % Add the ff vector to f at position sp.
  f(sp+1:L,:)=f(sp+1:L,:)+ff(1:L-sp,:);
  f(1:ep+1,:)=f(1:ep+1,:)+ff(L-sp+1:gl,:);
    
end;

% Scale correctly.
f=sqrt(M)*f;




