function h=gabmul(f,c,p3,p4,p5)
%GABMUL  Gabor multiplier.
%   Usage:  h=gabmul(f,c,a);
%           h=gabmul(f,c,g,a);
%           h=gabmul(f,c,ga,gs,a);
%
%   Input parameters:
%         f     : Input signal
%         c     : symbol of Gabor multiplier
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%   Output parameters:
%         h     : Output signal
%
%   GABMUL(f,c,g,a) will filter f by a Gabor multiplier determined by
%   the symbol c over the rectangular time-frequency lattice determined by 
%   _a and _M, where _M is deduced from the size of c. The rows of c
%   correspond to frequency, the columns to temporal sampling points.
%   The window g will be used for both analysis and synthesis.
%
%   GABMUL(f,c,a) will do the same using an optimally concentrated, tight
%   Gaussian as window function.
%
%   GABMUL(f,c,ga,gs,a) will do the same using the window ga for
%   analysis and gs for synthesis.
%
%   The adjoint operator of GABMUL(f,c,ga,gs,a) is given by
%   GABMUL(f,conj(c),gs,ga,a).
%
%   See also: dgt, idgt, gabdual, gabtight
%
%   Demos: demo_gabmul


error(nargchk(3,5,nargin));

M=size(c,1);
N=size(c,2);

if nargin==3
  a=p3;
  L=a*N;
  ga=gabtight(a,M,L);
  gs=ga;
end;

if nargin==4;
  ga=p3;
  gs=p3;
  a=p4;
 end;

if nargin==5;  
  ga=p3;
  gs=p4;
  a=p5;
end;

if numel(c)==1
  error('Size of symbol is too small. You probably forgot to supply it.');
end;


assert_squarelat(a,M,'GABMUL',0);

% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'DGT',0);

[coef,Ls]=dgt(f,ga,a,M);

for ii=1:W
  coef(:,:,ii)=coef(:,:,ii).*c;
end;


h=idgt(coef,gs,a,Ls);

% Change h to have same shape as f originally had.
h=comp_sigreshape_post(h,Ls,wasrow,remembershape);


%OLDFORMAT
