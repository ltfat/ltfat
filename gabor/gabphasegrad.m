function [tgrad,fgrad,c]=gabphasegrad(method,varargin)
%GABPHASEGRAD   Phase gradient of the DGT
%   Usage:  [tgrad,fgrad,c] = gabphasegrad('dgt',f,g,a,M);
%           [tgrad,fgrad]   = gabphasegrad('phase',cphase,a);
%           [tgrad,fgrad]   = gabphasegrad('abs',s,g,a);
%
%   `[tgrad,fgrad]=gabphasegrad(method,...)` computes the relative 
%   time-frequency gradient of the phase of the |dgt| of a signal. 
%   The derivative in time *tgrad* is the relative instantaneous 
%   frequency while the frequency derivative *fgrad* is the negative
%   of the local group delay.
%
%   *tgrad* is a measure the deviation from the current channel frequency,
%   so a value of zero means that the instantaneous frequency is equal to 
%   the center frequency of the considered channel, a positive value means
%   the true absolute intantaneous frequency is higher than the current 
%   channel frequency and vice versa. 
%   Similarly, *fgrad* is a measure of deviation from the current time 
%   positions.
%
%   *fgrad* is scaled such that distances are measured in samples. Similarly,
%   *tgrad* is scaled such that the Nyquist frequency (the highest possible
%   frequency) corresponds to a value of L/2. The absolute time and 
%   frequency positions can be obtained as
%
%      tgradabs = bsxfun(@plus,tgrad,fftindex(M)*L/M);
%      fgradabs = bsxfun(@plus,fgrad,(0:L/a-1)*a);
%
%   Please note that neither *tgrad* and *fgrad* nor *tgradabs* and 
%   *fgradabs* are true derivatives of the |dgt| phase. To obtain the true
%   phase derivatives, one has to explicitly pass either 'freqinv' or 
%   'timeinv' flags and scale both *tgrad* and *fgrad* by 2*pi/L.
%
%   The computation of *tgrad* and *fgrad* is inaccurate when the absolute
%   value of the Gabor coefficients is low. This is due to the fact the the
%   phase of complex numbers close to the machine precision is almost
%   random. Therefore, *tgrad* and *fgrad* may attain very large random values
%   when `abs(c)` is close to zero.
%
%   The computation can be done using four different methods.
%
%     'dgt'    Directly from the signal using algorithm by Auger and
%              Flandrin.
%
%     'phase'  From the phase of a DGT of the signal. This is the
%              classic method used in the phase vocoder.
%
%     'abs'    From the absolute value of the DGT. Currently this
%              method works only for Gaussian windows.
%
%     'cross'  Directly from the signal using algorithm by Nelson.
%
%   `[tgrad,fgrad]=gabphasegrad('dgt',f,g,a,M)` computes the time-frequency
%   gradient using a DGT of the signal *f*. The DGT is computed using the
%   window *g* on the lattice specified by the time shift *a* and the number
%   of channels *M*. The algorithm used to perform this calculation computes
%   several DGTs, and therefore this routine takes the exact same input
%   parameters as |dgt|.
%
%   The window *g* may be specified as in |dgt|. If the window used is
%   'gauss', the computation will be done by a faster algorithm.
%
%   `[tgrad,fgrad,c]=gabphasegrad('dgt',f,g,a,M)` additionally returns the
%   Gabor coefficients *c*, as they are always computed as a byproduct of the
%   algorithm.
%
%   `[tgrad,fgrad]=gabphasegrad('cross',f,g,a,M)` does the same as above
%   but this time using algorithm by Nelson which is based on computing 
%   several DGTs.
%
%   `[tgrad,fgrad]=gabphasegrad('phase',cphase,a)` computes the phase
%   gradient from the phase *cphase* of a DGT of the signal. The original DGT
%   from which the phase is obtained must have been computed using a
%   time-shift of *a* using the default phase convention (`'freqinv'`) e.g.::
%
%        [tgrad,fgrad]=gabphasegrad('phase',angle(dgt(f,g,a,M)),a)
%
%   `[tgrad,fgrad]=gabphasegrad('abs',s,g,a)` computes the phase gradient
%   from the spectrogram *s*. The spectrogram must have been computed using
%   the window *g* and time-shift *a* e.g.::
%
%        [tgrad,fgrad]=gabphasegrad('abs',abs(dgt(f,g,a,M)),g,a)
%
%   `[tgrad,fgrad]=gabphasegrad('abs',s,g,a,difforder)` uses a centered finite
%   diffence scheme of order *difforder* to perform the needed numerical
%   differentiation. Default is to use a 4th order scheme.
%
%   Currently the `'abs'` method only works if the window *g* is a Gaussian
%   window specified as a string or cell array.
%
%   See also: resgram, gabreassign, dgt
%
%   References: aufl95 cmdaaufl97 fl65 ltfatnote042


% AUTHOR: Peter L. Søndergaard, 2008; Zdenek Průša 2015

%narginchk(4,6);

% If no phaseconv flag was passed, add 'relative'
definput = arg_gabphasederivconv;
if ~any(cellfun(@(el) any(strcmpi(el,varargin)),definput.flags.phaseconv))
    varargin{end+1} = 'relative';
end

if nargout<3
    phased = gabphasederiv({'t','f'},method,varargin{:});
else
    [phased,c] = gabphasederiv({'t','f'},method,varargin{:});
end

[tgrad,fgrad] = deal(phased{:});
