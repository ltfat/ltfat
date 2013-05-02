function cout=comp_dgt_walnut(f,gf,a,M)
%COMP_DGT_WALNUT  First step of full-window factorization of a Gabor matrix.
%   Usage:  c=comp_dgt_walnut(f,gf,a,M);
%
%   Input parameters:
%         f      : Factored input data
%         gf     : Factorization of window (from facgabm).
%         a      : Length of time shift.
%         M      : Number of channels.
%   Output parameters:
%         c      : M x N*W*R array of coefficients, where N=L/a
%
%   Do not call this function directly, use DGT instead.
%   This function does not check input parameters!
%
%   The length of f and gamma must match.
%
%   If input is a matrix, the transformation is applied to
%   each column.
%
%   This function does not handle the multidim case. Take care before
%   calling this.
%
%   References: so07-2 st98-8

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

L=size(f,1);
W=size(f,2);
LR=numel(gf);
R=LR/L;

N=L/a;

[c,h_a,h_m]=gcd(a,M);
h_a=-h_a;
p=a/c;
q=M/c;
d=N/q;

ff=zeros(p,q*W,c,d,assert_classname(f,gf));

if p==1
  % --- integer oversampling ---

  if (c==1) && (d==1) && (W==1) && (R==1)
    % --- Short time Fourier transform of single signal ---
    % This is used for spectrograms of short signals.      
      ff(1,:,1,1)=f(:);
  else
    for s=0:d-1
      for r=0:c-1    
	for l=0:q-1
	  ff(1,l+1:q:W*q,r+1,s+1)=f(r+s*M+l*c+1,:);
	end;
      end;
    end;    
  end;

else
  % --- rational oversampling ---
  % Set up the small matrices
  % The r-loop (runs up to c) has been vectorized
  for w=0:W-1
    for s=0:d-1
      for l=0:q-1
	for k=0:p-1	    	  
	  ff(k+1,l+1+w*q,:,s+1)=f((1:c)+mod(k*M+s*p*M-l*h_a*a,L),w+1);
	end;
      end;
    end;
  end;
end;

% This version uses matrix-vector products and ffts

% fft them
if d>1
  ff=fft(ff,[],4);
end;

C=zeros(q*R,q*W,c,d,assert_classname(f,gf));

for r=0:c-1    
  for s=0:d-1
    GM=reshape(gf(:,r+s*c+1),p,q*R);
    FM=reshape(ff(:,:,r+1,s+1),p,q*W);
    
    C(:,:,r+1,s+1)=GM'*FM;
  end;
end;

% Inverse fft
if d>1
  C=ifft(C,[],4);
end;

% Place the result

cout=zeros(M,N,R,W,assert_classname(f,gf));

if p==1
  % --- integer oversampling ---

  if (c==1) && (d==1) && (W==1) && (R==1)
    
    % --- Short time Fourier transform of single signal ---
    % This is used for spectrograms of short signals.      
    for l=0:q-1
      cout(l+1,mod((0:q-1)+l,N)+1,1,1)=C(:,l+1,1,1);
    end;
    
  else

    % The r-loop (runs up to c) has been vectorized
    for rw=0:R-1
      for w=0:W-1    
	for s=0:d-1
	  for l=0:q-1
	    for u=0:q-1
	      cout((1:c)+l*c,mod(u+s*q+l,N)+1,rw+1,w+1)=C(u+1+rw*q,l+1+w*q,:,s+1);
	    end;
	  end;
	end;
      end; 
    end;
  end;

else

  % Rational oversampling
  % The r-loop (runs up to c) has been vectorized
  for rw=0:R-1
    for w=0:W-1    
      for s=0:d-1
	for l=0:q-1
	  for u=0:q-1
	    cout((1:c)+l*c,mod(u+s*q-l*h_a,N)+1,rw+1,w+1)=C(u+1+rw*q,l+1+w*q,:,s+1);
	  end;
	end;
      end;
    end; 
  end;

end;

cout=reshape(cout,M,N*W*R);




