function C=ref_dsft(F)
%REF_DSFT Reference DSFT

%   AUTHOR : Peter L. Søndergaard

K=size(F,1);
L=size(F,2);

C=zeros(L,K);

for m=0:L-1
  for n=0:K-1
    for l=0:L-1
      for k=0:K-1
	C(m+1,n+1)=C(m+1,n+1)+F(k+1,l+1)*exp(2*pi*i*(k*n/K-l*m/L));
      end;
    end;
  end;
end;

C=C/sqrt(K*L);

