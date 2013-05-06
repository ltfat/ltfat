function fout=ref_spreadop(f,c,a);
%REF_SPREADOP  Reference Spreading operator.
%


W=size(f,2);
M=size(c,1);
N=size(c,2);
L=N*a;

fout=zeros(L,W,assert_classname(f,c));

for l=0:L-1
  for m=0:M-1
    for n=0:N-1
      fout(l+1,:)=fout(l+1,:)+c(m+1,n+1)*exp(2*pi*i*l*m/M)*f(mod(l-n*a,L)+1,:);
    end;
  end;
end;


