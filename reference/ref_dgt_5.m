function coef=ref_dgt_5(f,g,a,M)
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

if 0
  % As ref_dgt_2, but substitution mt=k+st*p
  % and                            n =u+s*q
  % and apply reductions as described in the paper.
  for r=0:c-1    
    for l=0:q-1
      for u=0:q-1
        for k=0:p-1            
          for s=0:d-1
            for st=0:d-1
              w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)=w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)+f(mod(r+k*M+st*p*M-l*h_a*a,L)+1)*...
                  conj(g(mod(r+k*M-u*a+(st-s)*p*M,L)+1));
            end;
          end;
        end;
      end;
    end;    
  end;
end;

if 1

  % As above, but the inner convolution is now explicitly expressed
  for r=0:c-1    
    for l=0:q-1
      for u=0:q-1
        for k=0:p-1     
          psi=zeros(d,1);
          phi=zeros(d,1);
          for s=0:d-1
            psi(s+1)=f(mod(r+k*M+s*p*M-l*h_a*a,L)+1);              
            phi(s+1)=g(mod(r+k*M-u*a+s*p*M,L)+1);
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
end;

  
  
  
  
  
  

