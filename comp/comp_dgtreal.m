function c=comp_dgtreal(f,g,a,M,L,phasetype)
%COMP_DGTREAL  Compute a DGTREAL
%   Usage:  c=comp_dgt_real(f,g,a,M,L);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : M*N array of coefficients.
%

%   AUTHOR : Peter Soendergaard.

Lwindow=size(g,1);

W=size(f,2);
N=L/a;

M2=floor(M/2)+1;

if Lwindow<L
  % Do the filter bank algorithm
  % Periodic boundary conditions
  c=comp_dgtreal_fb(f,g,a,M,0);
  
else
  % Do the factorization algorithm 
  c=comp_dgtreal_long(f,g,a,M);
  
end;

c=reshape(c,M2,N,W);

if phasetype==1
    
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M2-1))/M;
    
    phase = FreqInd'*TimeInd;
    phase = exp(2*i*pi*phase);
    
    for w=1:W
        c(:,:,w) = c(:,:,w).*phase;
    end;

end;




