function coef=ref_dgt_6(f,g,a,M)
%REF_DGT_6  DGT algorithm 6
%
%  This algorithm assumes an FIR window.
%
%  This algorithm uses OLA to compute the small convolution

L=size(f,1);
Lg=size(g,1);
N=L/a;
b=L/M;

[c,h_a,h_m]=gcd(-a,M);
p=a/c;
q=M/c;
d=N/q;

w=zeros(M,N);

gwork=fir2long(g,L);

% As above, but the inner convolution is now explicitly expressed
for r=0:c-1    
  for l=0:q-1
    for u=0:q-1
      for k=0:p-1     
        psi=zeros(d,1);
        phi=zeros(d,1);
        for s=0:d-1
          psi(s+1)=f(mod(r+k*M+s*p*M-l*h_a*a,L)+1);   
        end;

        for s=0:d-1
          offset=r+k*q*c-u*p*c;
          phi(s+1)=gwork(mod(offset+s*p*q*c,L)+1);
        end;

        innerconv = pconv(psi,phi,'r');
        
        for s=0:d-1
          w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)=w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)+innerconv(s+1);
        end;
      end;
    end;
  end;    
end;

coef=fft(w);


  
  
  
  
  








  
  

