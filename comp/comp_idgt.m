function f=comp_idgt(coef,g,a,M,L,phasetype)
%COMP_IDGT  Compute IDGT
%   Usage:  f=comp_idgt(c,g,a,M,L,phasetype);
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

% AUTHOR : Peter Soendergaard.

b=L/M;
N=L/a;

Lwindow=size(g,1);
W=size(coef,3);

if phasetype==1
    coef=phaseunlock(coef,a);
end;

% FIXME: This line is necessary because the mex and oct interfaces expect
% a matrix as input.
coef=reshape(coef,M,prod(size(coef))/M);

if L==Lwindow
  % Do full-window algorithm.

  % Get the factorization of the window.
  gf = comp_wfac(g,a,M);      

  % Call the computational subroutine.
  f  = comp_idgt_fac(coef,gf,L,a,M);
  
 else
   
  % Do filter bank algorithm.
  % Call the computational subroutine.
  f=comp_idgt_fb(coef,g,L,a,M);

end;
%OLDFORMAT
