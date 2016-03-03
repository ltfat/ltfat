function [c,newphase,usedmask,tgrad,fgrad]=constructphase(s,g,a,varargin)
%CONSTRUCTPHASE  Construct phase for DGT
%   Usage:  c=constructphase(s,g,a);
%           c=constructphase(s,g,a,tol);
%           c=constructphase(c,g,a,tol,mask);
%           c=constructphase(c,g,a,tol,mask,usephase);
%           [c,newphase,usedmask,tgrad,fgrad] = constructphase(...);
%
%   Input parameters:
%         s        : Initial coefficients.
%         g        : Analysis Gabor window.
%         a        : Hop factor.
%         tol      : Relative tolerance.
%         mask     : Mask for selecting known phase.
%         usephase : Explicit known phase.
%   Output parameters:
%         c        : Coefficients with the constructed phase.
%         newphase : Just the (unwrapped) phase.
%         usedmask : Mask for selecting coefficients with the new phase.
%         tgrad    : Relative time phase derivative.
%         fgrad    : Relative frequency phase derivative.
%
%   `constructphase(s,g,a)` will construct a suitable phase for the postive
%   valued coefficients *s*. 
%  
%   If *s* is the absolute values of the Gabor coefficients of a signal
%   obtained using the window *g* and time-shift *a*, i.e.:
%
%     c=dgt(f,g,a,M);
%     s=abs(c);
%
%   then `constuctphase(s,g,a)` will attempt to reconstruct *c*.
%
%   The window *g* must be Gaussian, i.e. *g* must have the value `'gauss'`
%   or be a cell array `{'gauss',tfr}`. 
%
%   `constructphase(s,g,a,tol)` does as above, but sets the phase of
%   coefficients less than *tol* to random values. 
%   By default, *tol* has the value 1e-10.
%
%   `constructphase(c,g,a,M,tol,mask)` accepts real or complex valued
%   *c* and real valued *mask* of the same size. Values in *mask* which can
%   be converted to logical true (anything other than 0) determine
%   coefficients with known phase which is used in the output. Only the
%   phase of remaining coefficients (for which mask==0) is computed.
%
%   `constructphasereal(c,g,a,M,tol,mask,usephase)` does the same as before
%   but uses the known phase values from *usephase* rather than from *c*.
%
%   In addition, *tol* can be a vector containing decreasing values. In 
%   that case, the algorithm is run `numel(tol)` times, initialized with
%   the result from the previous step in the 2nd and the further steps.
%
%   Further, the function accepts the following flags:
%
%      'freqinv'  The constructed phase complies with the frequency
%                 invariant phase convention such that it can be directly
%                 used in |idgtreal|.
%                 This is the default.
%
%      'timeinv'  The constructed phase complies with the time-invariant
%                 phase convention. The same flag must be used in the other
%                 functions e.g. |idgtreal|
%
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%   See also:  dgt, gabphasegrad, ltfatmex
%
%   References: ltfatnote040
%

% AUTHOR: Peter L. Søndergaard, Zdenek Prusa

thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);

definput.keyvals.tol=[1e-1 1e-10];
definput.keyvals.mask=[];
definput.keyvals.usephase=[];
definput.flags.phase={'freqinv','timeinv'};
[flags,~,tol,mask,usephase]=ltfatarghelper({'tol','mask','usephase'},definput,varargin);

if ~isnumeric(s) 
    error('%s: *s* must be numeric.',thismfilename);
end

if ~isempty(usephase) && isempty(mask)
    error('%s: Both mask and usephase must be used at the same time.',...
          thismfilename);
end

if isempty(mask) 
    if ~isreal(s) || any(s(:)<0)
        error('%s: *s* must be real and positive when no mask is used.',...
              thismfilename);
    end
else 
    if any(size(mask) ~= size(s)) || ~isreal(mask)
        error(['%s: s and mask must have the same size and mask must',...
               ' be real.'],thismfilename)
    end
    % Sanitize mask (anything other than 0 is true)
    mask = cast(mask,'double');
    mask(mask~=0) = 1;
end

if ~isempty(usephase)
    if any(size(mask) ~= size(s)) || ~isreal(usephase)
        error(['%s: s and usephase must have the same size and usephase must',...
               ' be real.'],thismfilename)        
    end
else
    usephase = angle(s);
end

if ~isnumeric(tol) || ~isequal(tol,sort(tol,'descend'))
    error(['%s: *tol* must be a scalar or a vector sorted in a ',...
           'descending manner.'],thismfilename);
end

abss = abs(s);
% Compute phase gradient, check parameteres
[tgrad,fgrad] = gabphasegrad('abs',abss,g,a,2);

absthr = max(abss(:))*tol;
if isempty(mask)
    usedmask = zeros(size(s));
else
    usedmask = mask;
end

if isempty(mask)
    % Build the phase (calling a MEX file)
    newphase=comp_heapint(abss,tgrad,fgrad,a,tol(1),flags.do_timeinv);
    % Set phase of the coefficients below tol to random values
    bigenoughidx = abss>absthr(1);
    usedmask(bigenoughidx) = 1;
else
    newphase=comp_maskedheapint(abss,tgrad,fgrad,mask,a,tol(1),...
                                flags.do_timeinv,usephase);
    % Set phase of small coefficient to random values
    % but just in the missing part
    % Find all small coefficients in the unknown phase area
    missingidx = find(usedmask==0);
    bigenoughidx = abss(missingidx)>absthr(1);
    usedmask(missingidx(bigenoughidx)) = 1;
end

% Do further tol
for ii=2:numel(tol)
    newphase=comp_maskedheapint(abss,tgrad,fgrad,usedmask,a,tol(ii),...
                                flags.do_timeinv,newphase);
    missingidx = find(usedmask==0);
    bigenoughidx = abss(missingidx)>absthr(ii);
    usedmask(missingidx(bigenoughidx)) = 1;                  
end

% Convert the mask so it can be used directly for indexing
usedmask = logical(usedmask);
% Assign random values to coefficients below tolerance
zerono = numel(find(~usedmask));
newphase(~usedmask) = rand(zerono,1)*2*pi;

% Combine the magnitude and phase
c=abss.*exp(1i*newphase);

