function c = phaseunlockreal(c,a,M,varargin)
%PHASEUNLOCKREAL  Undo phase lock of Gabor coefficients
%   Usage:  c=phaseunlockreal(c,a,M);
%
%   `phaseunlockreal(c,a,M)` removes phase locking from the Gabor coefficients *c*.
%   The coefficient must have been obtained from a |dgtreal| with parameter *a*
%   and using the 'timeinv' flag.
%
%   Phase locking the coefficients modifies them so as if they were obtained
%   from a frequency-invariant Gabor system. A filter bank produces phase locked
%   coefficients. 
%
%   See also: dgt, phaselockreal
%
%   References: puc95

%   AUTHOR:    Peter Balazs, Peter L. Søndergaard, Zdenek Prusa
%   TESTING:   OK
%   REFERENCE: OK

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if  ~isscalar(a) || ~isnumeric(a) || rem(a,1)~=0
  error('a must be integer');
end;

if  ~isscalar(M) || ~isnumeric(M) || rem(M,1)~=0
  error('M must be integer');
end;

definput.flags.accuracy={'normal', 'precise'};
flags=ltfatarghelper({},definput,varargin);

M2=size(c,1);
M2user = floor(M/2) + 1;

if M2~=M2user
    error('%s: Size of s does not comply with M.',upper(mfilename));
end

N=size(c,2);

if flags.do_normal
    TimeInd = (0:(N-1))*a;
    FreqInd = (0:(M2-1));

    phase = FreqInd'*TimeInd;
    phase = mod(phase,M);
    phase = exp(-2*1i*pi*phase/M);

    % Handle multisignals
    c=bsxfun(@times,c,phase);
elseif flags.do_precise
    
    frames = ifftreal(c,M);
    for n=0:N-1
        frames(:,n+1,:) = circshift(frames(:,n+1,:),n*a);
    end
    c = fftreal(frames);
    
end
