function c=ref_dgt_3(f,g,a,M)
%REF_DGT_5  DGT algorithm 5
%
%  This algorithm uses explicit convolution inside the loops.

L=size(g,1);
N=L/a;
b=L/M;

[c,h_a,h_m]=gcd(-a,M);
p=a/c;
q=M/c;
d=N/q;

w=zeros(M,N);

% This version uses matrix-vector products and ffts

F=zeros(c,d,p,q);
G=zeros(c,d,p,q);
C=zeros(c,d,q,q);

% Set up the matrices
for r=0:c-1    
  for s=0:d-1
    for k=0:p-1
      for l=0:q-1
        F(r+1,s+1,k+1,l+1)=f(mod(r+k*M+s*p*M-l*h_a*a,L)+1);
        G(r+1,s+1,k+1,l+1)=sqrt(M*d)*g(mod(r+k*M-l*a+s*p*M,L)+1);
      end;
    end;
  end;
end;

% fft them
F=dft(F,[],2);
G=dft(G,[],2);

% Multiply them
for r=0:c-1    
  for s=0:d-1
    GM=reshape(G(r+1,s+1,:,:),p,q);
    FM=reshape(F(r+1,s+1,:,:),p,q);
    C(r+1,s+1,:,:)=GM'*FM;
  end;
end;

% Inverse fft
C=idft(C,[],2);

% Place the result
for r=0:c-1    
  for l=0:q-1
    for u=0:q-1
      for s=0:d-1
        w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)=C(r+1,s+1,u+1,l+1);
      end;
    end;
  end;
end; 

c=dft(w);
