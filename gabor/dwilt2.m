function [c,Ls]=dwilt2(f,g1,p3,p4,p5)
%DWILT2  2-D Discrete Wilson transform
%   Usage: c=dwilt2(f,g,M); 
%          c=dwilt2(f,g1,g2,[M1,M2]);
%          c=dwilt2(f,g1,g2,[M1,M2],[L1,L2]);
%          [c,L]=dwilt2(f,g1,g2,[M1,M2],[L1,L2]);
%
%   Input parameters:
%         f        : Input data, matrix.
%         g,g1,g2  : Window functions.
%         M,M1,M2  : Number of bands.
%         L1,L2    : Length of transform to do.
%   Output parameters:
%         c        : array of coefficients.
%         Ls       : Original size of input matrix.
%
%   `dwilt2(f,g,M)` calculates a two dimensional discrete Wilson transform
%   of the input signal *f* using the window *g* and parameter *M* along each
%   dimension.
%
%   For each dimension, the length of the transform will be the smallest
%   possible that is larger than the length of the signal along that dimension.
%   f will be appropriately zero-extended.
%
%   All windows must be whole-point even.
% 
%   `dwilt2(f,g,M,L)` computes a Wilson transform as above, but does
%   a transform of length *L* along each dimension. *f* will be cut or
%   zero-extended to length *L* before the transform is done.
%
%   `[c,Ls]=dwilt(f,g,M)` or `[c,Ls]=dwilt(f,g,M,L)` additionally returns the
%   length of the input signal *f*. This is handy for reconstruction.
%
%   `c=dwilt2(f,g1,g2,M)` makes it possible to use a different window along the
%   two dimensions. 
%
%   The parameters *L*, *M* and *Ls* can also be vectors of length 2. In
%   this case the first element will be used for the first dimension and the
%   second element will be used for the second dimension.
%
%   The output *c* has 4 or 5 dimensions. The dimensions index the
%   following properties:
%
%     1. Number of translation along 1st dimension of input.
%
%     2. Number of channel along 1st dimension  of input
%
%     3. Number of translation along 2nd dimension of input.
%
%     4. Number of channel along 2nd dimension  of input
%
%     5. Plane number, corresponds to 3rd dimension of input. 
% 
%   See also:  dwilt, idwilt2, dgt2, wildual

%   AUTHOR : Peter L. Søndergaard.

error(nargchk(3,5,nargin));

L=[];

if prod(size(p3))>2
  % Two windows was specified.
  g2=p3;
  M=p4;
  if nargin==5
    L=p5;
  end;
else
  g2=g1;
  M=p3;
  if nargin==4
    L=p4;
  end;
end;
  
if size(g1,2)>1
  if size(g1,1)>1
    error('g1 must be a vector');
  else
    % g1 was a row vector.
    g1=g1(:);
  end;
end;

if size(g2,2)>1
  if size(g2,1)>1
    error('g2 must be a vector');
  else
    % g2 was a row vector.
    g2=g2(:);
  end;
end;

if prod(size(M))>2 || prod(size(M))==0
  error('M must be a scalar or 1x2 vector');
end;

% If input is real, and window is real, output must be real as well.
inputwasreal = (isreal(g1) && isreal(g2) && isreal(f));

if length(M)==2
  M1=M(1);
  M2=M(2);
else
  M1=M;
  M2=M;
end;

assert_squarelat(M1,M1,1,'DWILT2');
assert_squarelat(M2,M2,1,'DWILT2');


Lwindow1=size(g1,1);
Lwindow2=size(g2,1);

Ls(1)=size(f,1);
Ls(2)=size(f,2);
W=size(f,3);

if ndims(f)<2 || ndims(f)>3
  error('f must be 2 or 3 dimensional.');
end;

if ~isempty(L) && (numel(L)>2 || numel(L)==0)
  error('L must be a scalar or 1x2 vector');
end;

if length(L)==2
  L1=L(1);
  L2=L(2);
else
  L1=L;
  L2=L;
end;

[b1,N1,L1]=assert_L(Ls(1),Lwindow1,L1,M1,M1,'DWILT2');
[b2,N2,L2]=assert_L(Ls(2),Lwindow2,L2,M2,M2,'DWILT2');

a1=M1;
a2=M2;

% --- Do transform along first dimension ---
% Change f to correct shape.
[f,fl,Wtemp,wasrow,remembershape]=comp_sigreshape_pre(f,'DWILT2',3);

% Zero-extend if neccesary.
f=postpad(f,L1);

c=comp_dwilt(f,g1,M1,L1);

% Combine dimensions again.
c=reshape(c,L1,Ls(2),W);

% Exchange first and second dimension.
c=permute(c,[2,1,3]);

% --- do transform along second dimension ---

% Change f to correct shape.
[c,cl,Wtemp,wasrow,remembershape]=comp_sigreshape_pre(c,'DWILT2',3);

% Zero-extend if neccesary.
c=postpad(c,L2);

c=comp_dwilt(c,g2,M2,L2);
c=comp_sigreshape_post(c,L2,wasrow,remembershape);

% Combine dimensions again.
c=reshape(c,L2,L1,W);

% Exchange first and second dimension.
c=permute(c,[2,1,3]);

% Reshape to final layout.
c=reshape(c,M1*2,N1/2,M2*2,N2/2,W);

% Clean coefficients if they are known to be real
if inputwasreal
  c=real(c);
end;

