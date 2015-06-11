function c = symphase(c,a,varargin)
%SYMPHASE  Change Gabor coefficients to symmetric phase
%   Usage:  c=symphase(c,a);
%
%   `symphase(c,a)` alters the phase of the Gabor coefficients *c* so as if
%   they were obtained from a Gabor transform based on symmetric
%   time/frequency shifts. The coefficient must have been obtained from a
%   |dgt| with parameter *a*.
%
%   Gabor coefficients with symmetric phase correspond to the following 
%   transform:
%   Consider a signal *f* of length *L* and define $N=L/a$.
%   The output from `c=symphase(dgt(f,g,a,M),a)` is given by
%
%   ..             L-1 
%     c(m+1,n+1) = sum f(l+1)*exp(-2*pi*i*m*(l-n*a/2)/M)*conj(g(l-a*n+1)), 
%                  l=0  
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)e^{-2\pi im(l-na/2)/M}\overline{g(l-an+1)}
%
%   where $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$ and $l-an$ is computed modulo *L*.
%
%   `symphase(c,a,'lt',lt)` does the same for a non-separable lattice
%   specified by *lt*. Please see the help of |matrix2latticetype| for a
%   precise description of the parameter *lt*.
%
%   See also: dgt, phaselock, phaseunlock
%
%   References: cmdaaufl97

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

