function f=comp_idgtreal(coef,g,a,M,L,phasetype)
%COMP_IDGTREAL  Compute IDGTREAL
%   Usage:  f=comp_idgtreal(c,g,a,M,L);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : length of transform.
%   Output parameters:
%         f     : Signal.
%

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

  
b=L/M;
N=L/a;

Lwindow=size(g,1);
W=size(coef,3);
M2=floor(M/2)+1;
M2short=ceil(M/2);

if phasetype==1
    TimeInd = (0:(N-1))/N;
    FreqInd = (0:(M2-1))*b;
    
    phase = FreqInd'*TimeInd;
    phase = exp(-2*i*pi*phase);
    
    % Handle multisignals
    for w=1:W
        coef(:,:,w) = coef(:,:,w).*phase;
    end;
    
end;

if L==Lwindow
  % Do full-window algorithm.

  % Get the factorization of the window.
  gf = comp_wfac(g,a,M);      

  % Call the computational subroutine.
  f = comp_idgtreal_fac(reshape(coef,M2,N*W),gf,L,a,M);
  
else
  % Do filter bank algorithm.
  % Call the computational subroutine.
  coef2=zeros(M,N,W);
  coef2(1:M2,:,:)=coef;
  if rem(M,2)==0
    coef2(M2+1:M,:,:)=conj(coef(M2-1:-1:2,:,:));
  else
    coef2(M2+1:M,:,:)=conj(coef(M2:-1:2,:,:));
  end;
  f = comp_idgt_fb(coef2,g,L,a,M);

end;