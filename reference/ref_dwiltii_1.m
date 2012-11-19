function [coef]=ref_dwiltii_1(f,g,a,M)
%COMP_DWILT  Compute Discrete Wilson transform.
%   
%   Do not call this function directly, use DWILT instead.

%   Author : Peter L. Søndergaard.

L=size(g,1);
N=L/a;
W=size(f,2);

c=ref_gdgt(f,g,a,2*M,.5,0,0);

coef2=reshape(c,2*M,N,W); % ----- Type II ------
%coef2=dgt(f,g,a,2*M);

coef=zeros(2*M,N/2,W);

if 0
  % --- Loop version ---
  for n=0:N/2-1

    % ---- m is zero ---------
    coef(1,n+1,:)=coef2(1,2*n+1,:);
    
    for m=1:2:M-1
      % --- m is odd ----------
      coef(m+1,n+1,:)=   i/sqrt(2)*(coef2(m+1,2*n+1,:)+coef2(2*M-m+1,2*n+1,:));
      coef(M+m+1,n+1,:)= 1/sqrt(2)*(coef2(m+1,2*n+2,:)-coef2(2*M-m+1,2*n+2,:));
    end;
    for m=2:2:M-1
      % --- m is even ---------
      coef(m+1,n+1,:)=   1/sqrt(2)*(coef2(m+1,2*n+1,:)-coef2(2*M-m+1,2*n+1,:));
      coef(M+m+1,n+1,:)= i/sqrt(2)*(coef2(m+1,2*n+2,:)+coef2(2*M-m+1,2*n+2,:)); 
    end;

    % --- m is nyquest ------
    if mod(M,2)==0
      coef(M+1,n+1,:) = i*coef2(M+1,2*n+2,:);
    else
      coef(M+1,n+1,:) = i*coef2(M+1,2*n+1,:);
    end;
    
  end;


else
  % --- Vector version---

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

end;


coef=reshape(coef,M*N,W);





