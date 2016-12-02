function c = phaselockreal(c,a,M,varargin)
%PHASELOCKREAL  Phaselock Gabor coefficients
%   Usage:  c=phaselockreal(c,a,M);
%
%   `phaselockreal(c,a,M)` phaselocks the Gabor coefficients *c*. 
%   The coefficients must have been obtained from a |dgtreal| with 
%   parameter *a*.
%
%   Phaselocking the coefficients modifies them so as if they were obtained
%   from a time-invariant Gabor system. A filter bank produces phase locked
%   coefficients.
%
%   Please see help of |phaselock| for more details.
%
%   See also: dgtreal, phaseunlockreal
%
%   References: puc95

%   AUTHOR:    Christoph Wiesmeyr, Peter L. Søndergaard, Zdenek Prusa
%   TESTING:   test_phaselock
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
    phase = exp(2*1i*pi*phase/M);

    % Handle multisignals
    c=bsxfun(@times,c,phase);
elseif flags.do_precise
    
    frames = ifftreal(c,M);
    for n=0:N-1
        frames(:,n+1,:) = circshift(frames(:,n+1,:),-n*a);
    end
    c = fftreal(frames);
    
end


