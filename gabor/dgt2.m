function [c,Ls]=dgt2(f,g1,p3,p4,p5,p6)
%DGT2  2-D Discrete Gabor transform
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
%   The output *c* has *4* or *5* dimensions. The dimensions index the
%   following properties:
%
%   1. Number of translation along 1st dimension of input.
%
%   2. Number of channel along 1st dimension  of input
%
%   3. Number of translation along 2nd dimension of input.
%
%   4. Number of channel along 2nd dimension  of input
%
%   5. Plane number, corresponds to 3rd dimension of input. 
% 
%   See also:  dgt, idgt2, gabdual

complainif_argnonotinrange(nargin,4,6,mfilename);

L=[];

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
  
if isempty(L)
  L1=[];
  L2=[];
else
  L1=L(1);
  L2=L(2);
end;

% Expand 'a' and M if necessary to two elements
a=bsxfun(@times,a,[1 1]);
M=bsxfun(@times,M,[1 1]);

Ls=size(f);
Ls=Ls(1:2);

c=dgt(f,g1,a(1),M(1),L1);
c=dgt(c,g2,a(2),M(2),L2,'dim',3);
