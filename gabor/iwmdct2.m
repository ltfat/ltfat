function [f]=iwmdct2(c,g1,p3,p4)
%IWMDCT2  2D Inverse windowed MDCT transform
%   Usage: f=iwmdct2(c,g);
%          f=iwmdct2(c,g1,g2);
%          f=iwmdct2(c,g1,g2,Ls);
%
%   Input parameters:
%         c       : Array of coefficients.
%         g,g1,g2 : Window functions.
%         Ls      : Size of reconstructed signal.
%   Output parameters:
%         f       : Output data, matrix.
%
%   `iwmdct2(c,g)` calculates a separable two dimensional inverse |wmdct|
%   transform of the input coefficients *c* using the window *g*. The number of
%   channels is deduced from the size of the coefficients *c*.
%
%   `iwmdct2(c,g1,g2)` does the same using the window *g1* along the first
%   dimension, and window *g2* along the second dimension.
%
%   `iwmdct2(c,g1,g2,Ls)` cuts the signal to size *Ls* after the transform
%   is done.
%
%   See also:  wmdct2, dgt2, wildual

%   AUTHOR : Peter L. Søndergaard

complainif_argnonotinrange(nargin,2,4,mfilename);

Ls=[];

switch nargin
  case 2
    g2=g1;
  case 3
    if prod(size(p3))>2
      % Two windows was specified.
      g2=p3;
    else
      g2=g1;
      Ls=p3;
    end;
  case 4
    g2=p3;
    Ls=p4;
end;
  

if ndims(c)<4 || ndims(c)>5
  error('c must be 4 or 5 dimensional.');
end;

M1=size(c,1);
N1=size(c,2);
M2=size(c,3);
N2=size(c,4);
W=size(c,5);

L1=M1*N1;
L2=M2*N2;

[g1,info]=wilwin(g1,M1,L1,'IWMDCT2');
[g2,info]=wilwin(g2,M2,L2,'IWMDCT2');

% If input is real, and window is real, output must be real as well.
inputwasreal = (isreal(g1) && isreal(g2) && isreal(c));


if isempty(Ls)
  Ls(1)=L1;
  Ls(2)=L2;
else
  Ls=bsxfun(@times,Ls,[1 1]);
end;

% --- first dimension

% Change c to correct shape.
c=reshape(c,M1,N1,L2*W);

c=comp_idwiltiii(c,g1);

c=postpad(c,Ls(1));

c=reshape(c,Ls(1),L2,W);

c=permute(c,[2,1,3]);

% --- second dimension

% Change c to correct shape.
c=reshape(c,M2,N2,Ls(1)*W);

c=comp_idwiltiii(c,g2);

c=postpad(c,Ls(2));

c=reshape(c,Ls(2),Ls(1),W);

f=permute(c,[2,1,3]);

% Clean signal if it is known to be real
if inputwasreal
  f=real(f);
end;

