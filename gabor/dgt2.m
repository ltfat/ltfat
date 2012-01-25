function [c,Ls]=dgt2(f,g1,p3,p4,p5,p6)
%DGT2  2-D Discrete Gabor transform.
%   Usage: c=dgt2(f,g,a,M);
%          c=dgt2(f,g1,g2,[a1,a2],[M1,M2]);
%          c=dgt2(f,g1,g2,[a1,a2],[M1,M2],[L1,L2]);
%          [c,Ls]=dgt2(f,g1,g2,[a1,a2],[M1,M2]);
%          [c,Ls]=dgt2(f,g1,g2,[a1,a2],[M1,M2],[L1,L2]);
%
%   Input parameters:
%         f       : Input data, matrix.
%         g,g1,g2 : Window functions.
%         a,a1,a2 : Length of time shifts.
%         M,M1,M2 : Number of modulations.
%         L1,L2   : Length of transform to do 
%
%   Output parameters:
%         c       : array of coefficients.
%         Ls      : Original size of input matrix.
%
%   `dgt2(f,g,a,M)` will calculate a separable two-dimensional discrete
%   Gabor transformation of the input signal *f* with respect to the window
%   *g* and parameters *a* and *M*.
%
%   For each dimension, the length of the transform will be the smallest
%   possible that is larger than the length of the signal along that dimension.
%   f will be appropriately zero-extended.
%
%   `dgt2(f,g,a,M,L)` computes a Gabor transform as above, but does
%   a transform of length *L* along each dimension. *f* will be cut or
%   zero-extended to length *L* before the transform is done.
%
%   `[c,Ls]=dgt2(f,g,a,M)` or `[c,Ls]=dgt2(f,g,a,M,L)` additionally returns
%   the length of the input signal *f*. This is handy for reconstruction::
%
%                [c,Ls]=dgt2(f,g,a,M);
%                fr=idgt2(c,gd,a,Ls);
%
%   will reconstruct the signal *f* no matter what the size of *f* is, provided
%   that *gd* is a dual window of *g*. 
%
%   `dgt2(f,g1,g2,a,M)` makes it possible to use a different window along the
%   two dimensions. 
%
%   The parameters *a*, *M*, *L* and *Ls* can also be vectors of length 2.
%   In this case the first element will be used for the first dimension
%   and the second element will be used for the second dimension.
%
%   The output *c* will be have *4* or *5* dimensions. The dimensions index the
%   following properties:
%
%      1. Number of translation along 1st dimension of input.
%
%      2. Number of channel along 1st dimension  of input
%
%      3. Number of translation along 2nd dimension of input.
%
%      4. Number of channel along 2nd dimension  of input
%
%      5. Plane number, corresponds to 3rd dimension of input. 
% 
%   See also:  dgt, idgt2, gabdual

error(nargchk(4,6,nargin));

L=[];

% Refuse to do nd-arrays.
%if ndims(f)>2
%  error('DGT2 does not work for nd-arrays.');
%end;

if prod(size(p3))>2
  % Two windows was specified.
  g2=p3;
  a=p4;
  M=p5;
  if nargin==6
    L=p6;
  end;
else
  g2=g1;
  a=p3;
  M=p4;
  if nargin==5
    L=p5;
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

if prod(size(a))>2 || prod(size(a))==0
  error('a must be a scalar or 1x2 vector');
end;

if prod(size(M))>2 || prod(size(M))==0
  error('M must be a scalar or 1x2 vector');
end;

if length(a)==2
  a1=a(1);
  a2=a(2);
else
  a1=a;
  a2=a;
end;

if length(M)==2
  M1=M(1);
  M2=M(2);
else
  M1=M;
  M2=M;
end;

assert_squarelat(a1,M1,1,'DGT2');
assert_squarelat(a2,M2,1,'DGT2');


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

[b1,N1,L1]=assert_L(Ls(1),Lwindow1,L1,a1,M1,'DGT2');
[b2,N2,L2]=assert_L(Ls(2),Lwindow2,L2,a2,M2,'DGT2');

% --- Do transform along first dimension ---

% Change f to correct shape.
[f,fl,Wtemp,wasrow,remembershape]=comp_sigreshape_pre(f,'DWILT2',3);

% Zero-extend if neccesary.
f=postpad(f,L1);

c=comp_dgt(f,g1,a1,M1,L1,0);

% Combine dimensions again.
c=reshape(c,M1*N1,Ls(2),W);

% Exchange first and second dimension.
c=permute(c,[2,1,3]);

% --- do transform along second dimension ---

% Change f to correct shape.
[c,cl,Wtemp,wasrow,remembershape]=comp_sigreshape_pre(c,'DWILT2',3);

% Zero-extend if neccesary.
c=postpad(c,L2);

c=comp_dgt(c,g2,a2,M2,L2,0);

c=comp_sigreshape_post(c,M2*N2,wasrow,remembershape);

% Combine dimensions again.
c=reshape(c,M2*N2,M1*N1,W);

% Exchange first and second dimension.
c=permute(c,[2,1,3]);

% Reshape to final layout.
c=reshape(c,M1,N1,M2,N2,W);
