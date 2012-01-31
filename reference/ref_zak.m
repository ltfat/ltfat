function coef=ref_zak(f,K)
%REF_ZAK   Reference Zak-transform.
%   Usage:  c=ref_zak(f,K);
%
%   This function computes a reference Zak-transform by
%   an explicit summation.
%
%   It returns the coefficients in the rectangular layout.
%

L=size(f,1);


% Slow, explicit method

% Workspace
coef=zeros(K,L/K);

for jj=0:K-1
  for kk=0:L/K-1
    for ll=0:L/K-1
      coef(jj+1,kk+1)=coef(jj+1,kk+1)+f(mod(jj-ll*K,L)+1)*exp(2*pi*i*kk*ll*K/L);
    end;
  end;
end;

coef=coef*sqrt(K/L);



