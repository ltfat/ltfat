function gf=comp_wfac(g,a,M)
%COMP_WFAC  Compute window factorization
%  Usage: gf=comp_wfac(g,a,M);
%
%   References: so07-2 st98-8

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK
  
L=size(g,1);
R=size(g,2);

N=L/a;
b=L/M;

c=gcd(a,M);
p=a/c;
q=M/c;
d=N/q;

gf=zeros(p,q*R,c,d,assert_classname(g));

% Set up the small matrices
% The w loop is only used for multiwindows, which should be a rare occurence.
% Therefore, we make it the outermost
if p==1  
  % Integer oversampling
  if (c==1) && (d==1) && (R==1)
    % --- Short time Fourier transform of single signal ---
    % This is used for spectrograms of short signals.            
    for l=0:q-1	  
      gf(1,l+1,1,1)=g(mod(-l,L)+1);
    end;

  else

    for w=0:R-1
      for s=0:d-1
	for l=0:q-1	  
	  gf(1,l+1+q*w,:,s+1)=g((1:c)+mod(-l*a+s*p*M,L),w+1);
	end;
      end;
    end;

  end;
else
  % Rational oversampling

  for w=0:R-1
    for s=0:d-1
      for l=0:q-1
	for k=0:p-1	    
	  gf(k+1,l+1+q*w,:,s+1)=g((1:c)+c*mod(k*q-l*p+s*p*q,d*p*q),w+1);
	end;
      end;
    end;
  end;

end;

% dft them
if d>1
  gf=fft(gf,[],4);
end;

% Scale by the sqrt(M) comming from Walnuts representation
gf=gf*sqrt(M);

gf=reshape(gf,p*q*R,c*d);

