function c = symphase(c,a,varargin)
%SYMPHASE  Change Gabor coefficients to symmetric phase
%   Usage:  c=symphase(c,a);
%
%   `symphase(c,a)` alters the phase of the Gabor coefficients *c* so as if
%   they were obtained from a Gabor transform based on symmetric
%   time/frequency shifts. The coefficient must have been obtained from a
%   |dgt| with parameter *a*.
%
%   See also: dgt, phaselock

%   AUTHORS : Peter Balazs, Peter L. Søndergaard.

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.lt=[0 1];
[flags,kv]=ltfatarghelper({},definput,varargin);

if  (prod(size(a))~=1 || ~isnumeric(a))
  error('a must be a scalar');
end;

if rem(a,1)~=0
  error('a must be an integer');
end;

M=size(c,1);
N=size(c,2);
L=N*a;
b=L/M;

if rem(b,1)~=0
  error('Lattice error. The a parameter is probably incorrect.');
end;

TimeInd = (0:(N-1))*a;
FreqInd = (0:(M-1));

phase = FreqInd'*TimeInd;
phase = mod(phase,M);
phase = exp(1i*pi*phase/M);

if kv.lt(1)>0 
    % truly non-separable case
    
    for n=0:(N-1)
        w = mod(n*kv.lt(1)/kv.lt(2),1);
        phase(:,n+1) = phase(:,n+1)*exp(pi*1i*a*w*n/M);
    end
end

% Handle multisignals
c=bsxfun(@times,c,phase);

