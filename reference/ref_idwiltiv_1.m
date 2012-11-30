function [f]=ref_idwiltiv_1(coef,g,a,M)
%REF_IDWILTIV_1  Reference IDWILTIV by IDGT type IV
% 

%   Author : Peter L. Søndergaard

L=size(g,1);
N=L/a;
W=size(coef,2);

coef=reshape(coef,M,N,W);

coef2=zeros(2*M,N,W);

if 0

  % --- loop version ---
  for n=0:N-1
    for m=0:M-1
      if rem(m+n,2)==0
	coef2(m+1,n+1,:)   =  exp(i*pi/4)/sqrt(2)*coef(m+1,n+1,:);
	coef2(2*M-m,n+1,:) =  exp(i*pi*3/4)/sqrt(2)*coef(m+1,n+1,:);
      else
	coef2(m+1,n+1,:)   =  exp(-i*pi/4)/sqrt(2)*coef(m+1,n+1,:);
	coef2(2*M-m,n+1,:) =  exp(-i*pi*3/4)/sqrt(2)*coef(m+1,n+1,:);
      end;
    end;
  end;
    
else

  % --- Vector version ---

  coef2(1:2:M,1:2:N,:)        = exp(i*pi/4)/sqrt(2)*coef(1:2:M,1:2:N,:);
  coef2(2*M:-2:M+1,1:2:N,:)   = exp(i*pi*3/4)/sqrt(2)*coef(1:2:M,1:2:N,:);
  
  coef2(1:2:M,2:2:N,:)        = exp(-i*pi/4)/sqrt(2)*coef(1:2:M,2:2:N,:);
  coef2(2*M:-2:M+1,2:2:N,:)   = exp(-i*pi*3/4)/sqrt(2)*coef(1:2:M,2:2:N,:);

  coef2(2:2:M,1:2:N,:)        = exp(-i*pi/4)/sqrt(2)*coef(2:2:M,1:2:N,:);
  coef2(2*M-1:-2:M+1,1:2:N,:) = exp(-i*pi*3/4)/sqrt(2)*coef(2:2:M,1:2:N,:);
  
  coef2(2:2:M,2:2:N,:)        = exp(i*pi/4)/sqrt(2)*coef(2:2:M,2:2:N,:);
  coef2(2*M-1:-2:M+1,2:2:N,:) = exp(i*pi*3/4)/sqrt(2)*coef(2:2:M,2:2:N,:);
  
end;

f=ref_igdgt(reshape(coef2,2*M*N,W),g,a,2*M,0.5,0.5,0);


