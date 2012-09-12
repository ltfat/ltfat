function G=frsynmat(F,L);
%FRSYNMAT  Frame synthesis operator matrix
%   Usage: G=frsynmat(F,L);
%
%   `G=frsynmat(F,L)` returns the matrix representation *G* of the frame
%   synthesis operator for a frame *F* of length *L*. The frame object *F*
%   must have been created using |newframe|_.
%
%   The frame matrix contains all the frame atoms as column vectors. It has
%   dimensions $L \times Ncoef$, where $Ncoef$ is the number of
%   coefficients. The number of coefficients can be found as
%   `Ncoef=framered(F)*L`. This means than the frame matrix is usually
%   **very** large, and this routine should only be used for small values
%   of *L*.
%
%   The action of the frame synthesis operator |frsyn|_ is equal to
%   multiplication with the frame synthesis operator
%   matrix. Consider the following simple example:::
%
%     L=200;
%     f=randn(L,1);
%     F=newframe('dgt','gauss','dual',10,20);
%     G=frsynmat(F,L);
%     testsig  = randn(L,1);
%     testcoef = frana(F,f); 
%     res = frsyn(F,testcoef)-G*testcoef;
%     norm(res)
%
%   See also: newframe, frana, frsyn, franaadj, franamat

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

switch(F.type)
 case {'dgtreal','fftreal'}
  error('%s: The synthesis operator of this frame does not have a matrix representation.',upper(mfilename));
 otherwise
  
  Lcheck=framelength(F,L);
  if Lcheck~=L
    error('%s: Incompatible frame length.',upper(mfilename));
  end;
  
  % Generic code handles all frames where there are no extra coefficients
  % in the representation
  Ncoef=framered(F)*L;
  coef=eye(Ncoef);
  G = frsyn(F,coef);  
end;

