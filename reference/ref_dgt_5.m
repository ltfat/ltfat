function c=ref_dgt_5(f,g,a,M)
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

if 0
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
end;

fv=zeros(d,1);
gv=zeros(d,1);
fvres=zeros(d,1);

% Place the result
for r=0:c-1    
  for l=0:q-1
    for u=0:q-1
      
      fvres=zeros(d,1);
      
      for k=0:p-1
        for s=0:d-1
          fv(s+1)=f(mod(r+k*M+s*p*M-l*h_a*a,L)+1);
          gv(s+1)=sqrt(M)*g(mod(r+k*M-l*a+s*p*M,L)+1);
        end;
        
        fv=pconv(fv,gv,'rr');
        
        fvres=fv+fvres;
      end;
      
      for s=0:d-1
        w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)=fvres(s+1);
      end;
    end;
  end;
end; 

c=dft(w);
