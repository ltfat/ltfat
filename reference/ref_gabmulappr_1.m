function sym=ref_gabmulappr_1(T,p2,p3,p4,p5);
%GABMULAPPR_1  Best Approximation by a Gabor multiplier.
%   Usage:  sym=gabmulappr(T,a,M);
%           sym=gabmulappr(T,g,a,M);
%           sym=gabmulappr(T,ga,gs,a,M);
%
%   Input parameters:
%         T     : matrix to be approximated
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%         M     : Number of channels.
%
%   Output parameters:
%         sym   : symbol
%
%   GABMULAPPR(T,g,a,M) will calculate the best approximation of the given
%   matrix T in the frobenius norm by a Gabor multiplier determined by the
%   symbol sym over the rectangular time-frequency lattice determined by a
%   and M.  The window g will be used for both analysis and synthesis.
%   IMPORTANT: The chosen Gabor system has to be a frame sequence!
%
%   GABMULAPPR(T,a,M) will do the same using an optimally concentrated,
%   tight Gaussian as window function.
%
%   GABMULAPPR(T,gs,ga,a) will do the same using the window ga for analysis
%   and gs for synthesis.
%
%   See also: gabmul
%
%   Demos: demo_gabmulappr
%
%   References: ba07 ba06 fehakr04
  
%   AUTHOR: P. Balazs (XXL)

% ---------- Verify the input -----------------
  
complainif_argnonotinrange(nargin,3,5,mfilename);

L=size(T,1);

if size(T,2)~=L
  error('T must be square.');
end;

if nargin==3
  % Usage: sym=gabmulappr(T,a,M);
  a=p2;
  M=p3;
  ga=gabtight(a,M,L);
  gs=ga;
end;

if nargin==4
  % Usage: sym=gabmulappr(T,g,a,M);
  ga=p2;
  gs=p2;
  a=p3;
  M=p4;
end;
  
if nargin==5
  % Usage: sym=gabmulappr(T,ga,gm,a,M);
  ga=p2;
  gs=p3;
  a=p4;
  M=p5;
end;

N=L/a;
b=L/M;

if size(ga,2)>1
  if size(ga,1)>1
    error('Input g/ga must be a vector');
  else
    % ga was a row vector.
    ga=ga(:);
  end;
end;

if size(gs,2)>1
  if size(gs,1)>1
    error('Input g/gs must be a vector');
  else
    % gs was a row vector.
    gs=gs(:);
  end;
end;

% -------- Algorithm starts here --------------------------

% Calculate the lower symbol. This is basically linear algebra

part1=reshape(dgt(T',ga,a,M),M*N,L);
part2=reshape(dgt(part1',gs,a,M),M*N,M*N).';
lowsym = reshape(diag(part2),M,N);

% Change from lower symbol to upper symbol. This is a quick calculation
% of a 2D convolution.

Gramfirst = conj(dgt(ga,ga,a,M)).*dgt(gs,gs,a,M);

Gramfirst = fft2(Gramfirst);
lowsym    = fft2(lowsym);
lowsym    = lowsym./Gramfirst;
sym       = ifft2(lowsym);
 




