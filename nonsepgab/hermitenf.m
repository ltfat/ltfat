function [B,err]=hermitenf(V)
%HERMITENF  Hermite normal form of 2x2 matrix
%   Usage:  A=hermitenf(V);
%
%   `hermitenf(V)` coverts the $2\times 2$ matrix `V` to Hermite normal form.
%   The output is lower triangular.  
%
%   For more information, see
%   `<http://en.wikipedia.org/wiki/Hermite_normal_form>`_.
%
%   See also: smithnf

% Code by Arno J. van Leest, 1999.
% and Peter Soendergaard, 2004.

if nargin~=1
  error('Wrong number of input arguments.');
end;

% Check if matrix has correct size.
if size(V,1)~=2 || size(V,2)~=2
  error('V must be a 2x2 matrix.');
end;

% Integer values
if norm(mod(V,1))~=0
  error('V must have all integer values.');
end;


% Convert to Arnos normal form.
gcd1=gcd(V(1,1),V(1,2));
gcd2=gcd(V(2,1),V(2,2));

A=zeros(2);
A(1,:)=V(1,:)/gcd1;
A(2,:)=V(2,:)/gcd2;

D = det(A);

% Exchange vectors if determinant is negative.
if D<0
  D=-D;
  A=fliplr(A);
end;

[g,h0,h1] = gcd(A(1,1),A(1,2));

x = A(2,:)*[h0;h1];

x = mod(x,D);
B = [gcd1 0;x*gcd2 D*gcd2];

% The following is sometimes needed for Octave.
% Octave mistakenly converts the numbers to floats,
% and later complains, so we round to integers.
B=round(B);

