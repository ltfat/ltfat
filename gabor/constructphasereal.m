function [c,newphase,tgrad,fgrad]=constructphasereal(s,g,a,M,varargin)
%CONSTRUCTPHASEREAL  Construct the phase of a DGTREAL
%   Usage:  c=constructphasereal(s,g,a,M);
%           c=constructphasereal(s,g,a,M,tol);
%           c=constructphasereal(c,g,a,M,tol,mask);
%           [c,newphase,tgrad,fgrad] = constructphasereal(...);
%
%   `constructphasereal(s,g,a,M)` will construct a suitable phase for the positive
%   valued coefficients *s*.
%
%   If *s* is the absolute values of the Gabor coefficients of a signal
%   obtained using the window *g*, time-shift *a* and number of channels *M*, i.e.:
%
%     c=dgtreal(f,g,a,M);
%     s=abs(c);
%
%   then `constuctphasereal(s,g,a,M)` will attempt to reconstruct *c*.
%
%   The window *g* must be Gaussian, i.e. *g* must have the value `'gauss'`
%   or be a cell array `{'gauss',tfr}`.
%
%   `constructphasereal(s,g,a,M,tol)` does as above, but sets the phase of
%   coefficients less than *tol* to random value.
%   By default, *tol* has the value 1e-10.
%
%   `constructphasereal(c,g,a,M,tol,mask)` accepts real or complex valued
%   *c* and real valued *mask* of the same size. Values in *mask* which can
%   be converted to logical true (anything other than 0) determine
%   coefficients with known phase which is used in the output. Only the
%   phase of remaining coefficients (for which mask==0) is computed.
%
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%
%   See also:  dgt, gabphasegrad, ltfatmex
%

% AUTHOR: Peter L. Søndergaard, Zdenek Prusa
thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
complainif_notposint(M,'M',thismfilename);

definput.keyvals.tol=1e-10;
definput.keyvals.mask=[];
definput.flags.phase={'freqinv','timeinv'};
[flags,~,tol,mask]=ltfatarghelper({'tol','mask'},definput,varargin);

if ~isnumeric(s) 
    error('%s: *s* must be numeric.',thismfilename);
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

if ~isscalar(tol)
    error('%s: *tol* must be scalar.',thismfilename);
end

[M2,N,W] = size(s);

if W>1
    error('%s: *s* must not be 3 dimensional.',thismfilename);
end

M2true = floor(M/2) + 1;

if M2true ~= M2
    error('%s: Mismatch between *M* and the size of *s*.',thismfilename);
end

L=N*a;

[~,info]=gabwin(g,a,M,L,'callfun',upper(mfilename));

if ~info.gauss
    error(['%s: The window must be a Gaussian window (specified ',...
           'as a string or as a cell array)'],upper(mfilename));
end

% Here we try to avoid calling gabphasegrad as it only works with full
% dgts.
abss = abs(s);
logs=log(abss+realmin);
tt=-11;
logs(logs<max(logs(:))+tt)=tt;

fgrad = info.tfr*pderiv(logs,2,2)/(2*pi);
% Undo the scaling done by pderiv and scale properly
tgrad = pderiv(logs,1,2)/(2*pi*info.tfr)*(M/M2);

% Fix the first and last rows .. the
% borders are symmetric so the centered difference is 0
tgrad(1,:) = 0;
tgrad(end,:) = 0;

absthr = max(abss(:))*tol;

% Build the phase
if isempty(mask)
    % s is real and positive
    newphase=comp_heapintreal(abss,tgrad,fgrad,a,M,tol,flags.do_timeinv);

    % Set phase of small coefficient to random values
    toosmallidx = abss<absthr;
    zerono = numel(find(toosmallidx));
    newphase(toosmallidx) = rand(zerono,1)*2*pi;
else
    newphase=comp_maskedheapintreal(s,tgrad,fgrad,mask,a,M,tol,flags.do_timeinv);
    % Set phase of small coefficient to random values
    % but just in the missing part
    missingidx = find(mask==0);
    toosmallidx = abss(missingidx)<absthr;
    zerono = numel(find(toosmallidx));
    newphase(missingidx(toosmallidx)) = rand(zerono,1)*2*pi;
end

c=abss.*exp(1i*newphase);
