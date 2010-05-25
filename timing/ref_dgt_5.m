function c=ref_dgt_5(f,g,a,M)
%REF_DGT_4  DGT algorithm 4
%
%  This algorithm makes the r-loop be the outermost loop, in order to
%  reduce the size of the intermediate buffers.
%
%  It integrates the DFTs into the factorizations and the matrix
%  computation.
%
L=size(g,1);
R=1;
W=size(f,2);
N=L/a;
b=L/M;

[c,h_a,h_m]=gcd(-a,M);
p=a/c;
q=M/c;
d=N/q;

w=zeros(M,N);

% This version uses matrix-vector products and ffts

F=zeros(p,q,d);
G=zeros(p,q,c,d);
C=zeros(q,q,d);
sbuf=zeros(d,1);

% Setup the window factorization.
% This is done for both r and s at once to simulate precomputing the
% window factorization.
for s=0:d-1
  for r=0:c-1    
    for k=0:p-1
      for l=0:q-1
        G(k+1,l+1,r+1,s+1)=sqrt(M*d)*g(mod(r+k*M-l*a+s*p*M,L)+1);
      end;
    end;
  end;
end;

if d>1
  G=dft(G,[],4);
end;

% -------- Main loop -------------
for r=0:c-1
  
  % Set up the signal factorization
  for k=0:p-1
    for l=0:q-1
      for s=0:d-1
	sbuf(s+1)=f(mod(r+k*M+s*p*M-l*h_a*a,L)+1);
      end;      

      sbuf=dft(sbuf);

      F(k+1,l+1,:)=sbuf;
    end;
  end;

  % Multiply them
  for s=0:d-1    
    GM=reshape(G(:,:,r+1,s+1),p,q);
    FM=reshape(F(:,:,s+1),p,q);
    C(:,:,s+1)=GM'*FM;
  end;

  for l=0:q-1
    for u=0:q-1        	
      sbuf(:)=C(u+1,l+1,:);
      sbuf=idft(sbuf);
      
      for s=0:d-1
	w(r+l*c+1,mod(u+s*q-l*h_a,N)+1)=sbuf(s+1);
      end;
    end;
  end;
  
end;  % Main loop ends here.

c=dft(w);
