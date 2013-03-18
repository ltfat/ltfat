function f=shah(L,a);
%SHAH  Discrete Shah-distribution
%   Usage: f=shah(L,a);
% 
%   `shah(L,a)` computes the discrete, normalized Shah-distribution of
%   length *L* with a distance of *a* between the spikes.
%
%   The Shah distribution is defined by 
%
%   .. f(n*a+1)=1/sqrt(L/a) 
%
%   .. math:: f(n\cdot a+1)=\frac{1}{\sqrt(L/a)} 
%
%   for integer *n*, otherwise *f* is zero.
% 
%   This is also known as an impulse train or as the comb function, because
%   the shape of the function resembles a comb. It is the sum of unit
%   impulses ('diracs') with the distance *a*.
% 
%   If *a* divides *L*, then the |dft| of `shah(L,a)` is `shah(L,L/a)`.
% 
%   The Shah function has an extremely bad time-frequency localization.
%   It does not generate a Gabor frame for any *L* and *a*.
%
%   Examples:
%   ---------
%
%   A simple spectrogram of the Shah function (includes the negative
%   frequencies to display the whole TF-plane):::
%
%     sgram(shah(256,16),'dynrange',80,'nf')
%
  
%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

if nargin~=2
  error('Wrong number of input parameters.');
end;

%if mod(L,a)~=0
%  error('a must divide L.');
%end;

f=zeros(L,1);

f(1:a:L)=1/sqrt(L/a);

