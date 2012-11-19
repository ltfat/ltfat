function [f]=idwilt2(c,g1,p3,p4)
%IDWILT2  2D Inverse Discrete Wilson transform
%   Usage: f=idwilt2(c,g);
%          f=idwilt2(c,g1,g2);
%          f=idwilt2(c,g1,g2,Ls);
%
%   Input parameters:
%         c       : Array of coefficients.
%         g,g1,g2 : Window functions.
%         Ls      : Size of reconstructed signal.
%   Output parameters:
%         f       : Output data, matrix.
%
%   `idwilt2(c,g)` calculates a separable two dimensional inverse
%   discrete Wilson transformation of the input coefficients *c* using the
%   window *g*. The number of channels is deduced from the size of the
%   coefficients *c*.
%
%   `idwilt2(c,g1,g2)` does the same using the window *g1* along the first
%   dimension, and window *g2* along the second dimension.
%
%   `idwilt2(c,g1,g2,Ls)` cuts the signal to size *Ls* after the transformation
%   is done.
%
%   See also:  dwilt2, dgt2, wildual

%   AUTHOR : Peter L. Søndergaard

error(nargchk(2,4,nargin));

doLs=0;

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
      doLs=1;
    end;
  case 4
    g2=p3;
    Ls=p4;
    doLs=1;
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

Lwindow1=size(g1,1);
Lwindow2=size(g2,1);

if ndims(c)<4 || ndims(c)>5
  error('c must be 4 or 5 dimensional.');
end;

M1=size(c,1)/2;
N1=size(c,2)*2;
M2=size(c,3)/2;
N2=size(c,4)*2;
W=size(c,5);

a1=M1;
a2=M2;

L1=a1*N1;
L2=a2*N2;

% Length of window must be dividable by M.
% We cannot automically zero-extend the window, as it can
% possible break some symmetry properties of the window, and we don't
% know which symmetries to preserve.
if rem(Lwindow1,M1)~=0
  error('Length of window no. 1 must be dividable by M1.')
end;

if rem(Lwindow2,M2)~=0
  error('Length of window no. 2 must be dividable by M2.')
end;

% If input is real, and window is real, output must be real as well.
inputwasreal = (isreal(g1) && isreal(g2) && isreal(c));


% --- first dimension

% Change c to correct shape.
c=reshape(c,2*M1,N1/2,L2*W);

c=comp_idwilt(c,g1);

% Check if Ls was specified.
if doLs
  c=postpad(c,Ls(1));
else
  Ls(1)=L1;
end;

% Change to correct size
c=reshape(c,Ls(1),L2,W);

% Exchange first and second dimension.
c=permute(c,[2,1,3]);

% --- second dimension

% Change c to correct shape.
c=reshape(c,2*M2,N2/2,Ls(1)*W);

c=comp_idwilt(c,g2);

% Check if Ls was specified.
if doLs
  c=postpad(c,Ls(2));
else
  Ls(2)=L2;
end;

% Change to correct size
c=reshape(c,Ls(2),Ls(1),W);

% Exchange first and second dimension.
f=permute(c,[2,1,3]);

% Clean signal if it is known to be real
if inputwasreal
  f=real(f);
end;

