function c=comp_dgtreal(f,g,a,M,L,wilson,callfun)
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
  c=c(1:M2,:);
  
else
  % Do the factorization algorithm 
  c=comp_dgtreal_long(f,g,a,M);
  
end;

c=reshape(c,M2,N,W);




