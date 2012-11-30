function [f]=ref_idwilt_1(coef,g,a,M)
%REF_IDWILT_1  Reference IDWILT by IDGT
% 

%   Author : Peter L. Søndergaard

L=size(g,1);
N=L/a;
W=size(coef,2);

coef=reshape(coef,M*2,N/2,W);

coef2=zeros(2*M,N,W);

if 1

  % ----- Loop version ---

  for n=0:N/2-1

    % m=0
    coef2(1,2*n+1,:) = coef(1,n+1,:);

    for m=0:M-1
  
      % m odd
      for m=1:2:M-1
	coef2(m+1,2*n+1,:)     = -i/sqrt(2)*coef(m+1,n+1,:);
	coef2(m+1,2*n+2,:)     =  1/sqrt(2)*coef(M+m+1,n+1,:);
      
	coef2(2*M-m+1,2*n+1,:) =  i/sqrt(2)*coef(m+1,n+1,:);
	coef2(2*M-m+1,2*n+2,:) =  1/sqrt(2)*coef(M+m+1,n+1,:);
      end;      

      % m even
      for m=2:2:M-1
	coef2(m+1,2*n+1,:)     =  1/sqrt(2)*coef(m+1,n+1,:);
	coef2(m+1,2*n+2,:)     = -i/sqrt(2)*coef(M+m+1,n+1,:);
	
	coef2(2*M-m+1,2*n+1,:) =  1/sqrt(2)*coef(m+1,n+1,:);
	coef2(2*M-m+1,2*n+2,:) =  i/sqrt(2)*coef(M+m+1,n+1,:);
      end;    

    end;

    % m=nyquest
    if mod(M,2)==0
      coef2(M+1,2*n+1,:) = coef(M+1,n+1,:);
    else
      coef2(M+1,2*n+2,:) = coef(M+1,n+1,:);
    end;

  end;

else
  % ----- Vector version ---

  % First and middle modulation are transferred unchanged.
  coef2(1,1:2:N,:) = coef(1,:,:);
  if mod(M,2)==0
    coef2(M+1,1:2:N,:) = coef(M+1,:,:);
  else
    coef2(M+1,2:2:N,:) = coef(M+1,:,:);
  end;
  
  if M>2
    coef2(3:2:M,1:2:N,:)        = 1/sqrt(2)*coef(3:2:M,:,:);
    coef2(3:2:M,2:2:N,:)        = -1/sqrt(2)*i*coef(M+3:2:2*M,:,:);
    
    coef2(2*M-1:-2:M+2,1:2:N,:) = 1/sqrt(2)*coef(3:2:M,:,:);
    coef2(2*M-1:-2:M+2,2:2:N,:) =  1/sqrt(2)*i*coef(M+3:2:2*M,:,:);
  end;
  
  
  % sine, first column.
  coef2(2:2:M,1:2:N,:)        = -1/sqrt(2)*i*coef(2:2:M,:,:);
  coef2(2:2:M,2:2:N,:)        = 1/sqrt(2)*coef(M+2:2:2*M,:,:);
  
  coef2(2*M:-2:M+2,1:2:N,:)   =  1/sqrt(2)*i*coef(2:2:M,:,:);
  coef2(2*M:-2:M+2,2:2:N,:)   = 1/sqrt(2)*coef(M+2:2:2*M,:,:);
  

end;


f=idgt(reshape(coef2,2*M,N,W),g,a);



