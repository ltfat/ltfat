function c=ref_dgt_1(f,g,a,M)
%REF_DGT_1  DGT by Poisson summation.

L=size(g,1);
N=L/a;
b=L/M;

w=zeros(M,N);

if 0
  
  % This version uses the definition
  
  for jj=0:M-1  
    for nn=0:N-1
      for kk=0:b-1
	w(jj+1,nn+1)=w(jj+1,nn+1)+f(jj+kk*M+1)*conj(g(mod(jj+kk*M-nn*a,L)+1));
      end;
    end;
  end;

else

  % This version uses matrix-vector products.
  
  W=zeros(b,N);
  v=zeros(b,1);
  for jj=0:M-1

    % Setup the matrix and vector
    for kk=0:b-1
      for nn=0:N-1
	W(kk+1,nn+1)=g(mod(jj+kk*M-nn*a,L)+1);
      end;

      v(kk+1)=f(mod(jj+kk*M,L)+1);
    end;

    % do the product.
    s1=W'*v;

    % Arrange in w
    for nn=0:N-1
      w(jj+1,nn+1)=s1(nn+1);
    end;
   
  end;


end;

c=fft(w);


