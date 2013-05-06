function [coef]=comp_dwilt_fb(f,g,M,L)
%COMP_DWILT_FB  Compute Discrete Wilson transform.
%   

a=M;
N=L/a;
W=size(f,2);


coef=zeros(2*M,N/2,W,assert_classname(f,g));

if (isreal(f) && isreal(g))

  coef2=comp_dgt_fb(f,g,a,2*M);

  % If the input coefficients are real, the calculations can be
  % be simplified. The complex case code also works for the real case.

  % Unmodulated case.
  coef(1,:,:)=coef2(1,1:2:N,:);

  % cosine, first column.
  coef(3:2:M,:,:)=sqrt(2)*real(coef2(3:2:M,1:2:N,:));
  
  % sine, second column
  coef(M+3:2:2*M,:,:)=-sqrt(2)*imag(coef2(3:2:M,2:2:N,:));
  
  % sine, first column.
  coef(2:2:M,:,:)=-sqrt(2)*imag(coef2(2:2:M,1:2:N,:));
  
  % cosine, second column
  coef(M+2:2:2*M,:,:)=sqrt(2)*real(coef2(2:2:M,2:2:N,:));

  % Nyquest case
  if mod(M,2)==0
    coef(M+1,:,:) = coef2(M+1,1:2:N,:);
  else
    coef(M+1,:,:) = coef2(M+1,2:2:N,:);
  end;


else
  % Complex valued case
  
  coef2=comp_dgt_fb(f,g,a,2*M);
  
  % Unmodulated case.
  coef(1,:,:)=coef2(1,1:2:N,:);

  % odd value of m
  coef(2:2:M,:,:)=1/sqrt(2)*i*(coef2(2:2:M,1:2:N,:)-coef2(2*M:-2:M+2,1:2:N,:));
  coef(M+2:2:2*M,:,:)=1/sqrt(2)*(coef2(2:2:M,2:2:N,:)+coef2(2*M:-2:M+2,2:2:N,:));

  % even value of m
  coef(3:2:M,:,:)=1/sqrt(2)*(coef2(3:2:M,1:2:N,:)+coef2(2*M-1:-2:M+2,1:2:N,:));
  coef(M+3:2:2*M,:,:)=1/sqrt(2)*i*(coef2(3:2:M,2:2:N,:)-coef2(2*M-1:-2:M+2,2:2:N,:));

  % Nyquest case
  if mod(M,2)==0
    coef(M+1,:,:) = coef2(M+1,1:2:N,:);
  else
    coef(M+1,:,:) = coef2(M+1,2:2:N,:);
  end;

end;



