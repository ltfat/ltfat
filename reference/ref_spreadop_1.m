function fout=ref_spreadop_1(f,c,a);
%REF_SPREADOP_1  Ref. Spreading function
%   Usage: h=ref_spreadfun_1(f,c,a);
%

%XXX Wrong sign convention
  
W=size(f,2);
M=size(c,1);
N=size(c,2);
L=N*a; 

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

fout=zeros(L,W);

cf1=fft(coef);

for j=0:M-1
  for m=0:b-1    
    for n=0:N-1
      fout(j+m*M+1,:)=fout(j+m*M+1,:)+cf1(j+1,n+1)*f(mod(j+m*M-n*a,L)+1,:);
    end;
  end;
end;



