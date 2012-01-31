function [coef]=ref_dwilt_1(f,g,a,M)
%COMP_DWILT  Compute Discrete Wilson transform by DGT
%   

L=size(g,1);
N=L/a;
W=size(f,2);

coef2=dgt(f,g,a,2*M);

coef=zeros(2*M,N/2,W);

if 1

  % Loop version
  
  for n=0:N/2-1

    % ---- m is zero ---------
    coef(1,n+1,:)=coef2(1,2*n+1,:);
    
    for m=1:2:M-1
      % --- m is odd ----------
      coef(m+1,n+1,:)=  i/sqrt(2)*(coef2(m+1,2*n+1,:)-coef2(2*M-m+1,2*n+1,:));
      coef(M+m+1,n+1,:)=1/sqrt(2)*(coef2(m+1,2*n+2,:)+coef2(2*M-m+1,2*n+2,:));
    end;
    for m=2:2:M-1
      % --- m is even ---------
      coef(m+1,n+1,:)=  1/sqrt(2)*(coef2(m+1,2*n+1,:)+coef2(2*M-m+1,2*n+1,:));
      coef(M+m+1,n+1,:)=i/sqrt(2)*(coef2(m+1,2*n+2,:)-coef2(2*M-m+1,2*n+2,:)); 
    end;

    % --- m is nyquest ------
    if mod(M,2)==0
      coef(M+1,n+1,:) = coef2(M+1,2*n+1,:);
    else
      coef(M+1,n+1,:) = coef2(M+1,2*n+2,:);
    end;
    
  end;


else

  % Vector version
  % ---- m is zero ---------
  
  coef(1,:,:)=coef2(1,2*n+1,:);
  
  
  % --- m is odd ----------
  % sine, first column.
  coef(2:2:M,:,:)=1/sqrt(2)*i*(coef2(2:2:M,1:2:N,:)-coef2(2*M:-2:M+2,1:2:N,:));
  
  % cosine, second column
  coef(M+2:2:2*M,:,:)=1/sqrt(2)*(coef2(2:2:M,2:2:N,:)+coef2(2*M:-2:M+2,2:2:N,:));
  
  % --- m is even ---------
  
  % cosine, first column.
  coef(3:2:M,:,:)=1/sqrt(2)*(coef2(3:2:M,1:2:N,:)+coef2(2*M-1:-2:M+2,1:2:N,:));
  
  % sine, second column
  coef(M+3:2:2*M,:,:)=1/sqrt(2)*i*(coef2(3:2:M,2:2:N,:)-coef2(2*M-1:-2:M+2,2:2:N,:));
  
  % --- m is nyquest ------
  if mod(M,2)==0
    coef(M+1,:,:) = coef2(M+1,1:2:N,:);
  else
    coef(M+1,:,:) = coef2(M+1,2:2:N,:);
  end;

end;

coef=reshape(coef,M*N,W);




