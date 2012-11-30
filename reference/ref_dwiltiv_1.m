function [coef]=ref_dwiltiv_1(f,g,a,M)
%COMP_DWILT  Compute Discrete Wilson transform.
%   
%   Do not call this function directly, use DWILT instead.

%   Author : Peter L. Søndergaard.

L=size(g,1);
N=L/a;
W=size(f,2);

c=ref_gdgt(f,g,a,2*M,0.5,.5,0);

coef2=reshape(c,2*M,N,W); % ----- Type IV ------

coef=zeros(M,N,W);

if 0
  % --- Loop version ---
  for n=0:N-1
    for m=0:M-1
      if rem(m+n,2)==0
	coef(m+1,n+1,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(m+1,n+1,:)+exp(-i*pi*3/4)*coef2(2*M-m,n+1,:));
      else
	coef(m+1,n+1,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(m+1,n+1,:)+exp(i*pi*3/4)*coef2(2*M-m,n+1,:));
      end;
    end;
  end;

else
  % --- Vector version---

  % --- m is even ---------
  coef(1:2:M,1:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(1:2:M,1:2:N,:)+exp(-i*pi*3/4)*coef2(2*M:-2:M+1,1:2:N,:));
  coef(1:2:M,2:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(1:2:M,2:2:N,:)+exp(i*pi*3/4)*coef2(2*M:-2:M+1,2:2:N,:));
  
  % --- m is odd ----------
  coef(2:2:M,1:2:N,:)= 1/sqrt(2)*(exp(i*pi/4)*coef2(2:2:M,1:2:N,:)+exp(i*pi*3/4)*coef2(2*M-1:-2:M+1,1:2:N,:));
  coef(2:2:M,2:2:N,:)= 1/sqrt(2)*(exp(-i*pi/4)*coef2(2:2:M,2:2:N,:)+exp(-i*pi*3/4)*coef2(2*M-1:-2:M+1,2:2:N,:));

end;


coef=reshape(coef,M*N,W);





