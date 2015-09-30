function [c,newphase,tgrad,fgrad]=constructphasereal(s,g,a,M,tol)
%CONSTRUCTPHASEREAL  Construct the phase of a DGTREAL
%   Usage:  c=constructphasereal(s,g,a,M);
%           c=constructphasereal(s,g,a,M,tol);
%
%   `constructphasereal(s,g,a,M)` will construct a suitable phase for the postive
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
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%   Examples:
%   ---------
%
%   
%
%   See also:  dgt, gabphasegrad, ltfatmex
%

% AUTHOR: Peter L. Søndergaard, Zdenek Prusa
thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
complainif_notposint(M,'M',thismfilename);

if ~isnumeric(s) || ~isreal(s)
    error('%s: *s* must be a real matrix.',thismfilename);
end

if nargin<5
    tol=1e-10;
else
    if ~isscalar(tol)
        error('%s: *tol* must be scalar.',thismfilename);
    end
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
b=L/M;

[~,info]=gabwin(g,a,M,L,'callfun',upper(mfilename));

if ~info.gauss
    error(['%s: The window must be a Gaussian window (specified ',...
           'as a string or as a cell array)'],upper(mfilename));
end

% Here we try to avoid calling gabphasegrad as it only works with full
% dgts.

logs=log(s+realmin);
tt=-11;
logs(logs<max(logs(:))+tt)=tt;

fgrad = info.tfr*pderiv(logs,2,2)/(2*pi);
% Undo the scaling done by pderiv and scale properly
tgrad = pderiv(logs,1,2)/(2*pi*info.tfr)*(M/M2);

% Fix the first and last rows .. the
% borders are symmetric so the centered difference is 0
tgrad(1,:) = 0;
tgrad(end,:) = 0;

% Build the phase
newphase=comp_heapintreal(s,tgrad,fgrad,a,M,tol);

% Set phase of small coefficient to random values
absthr = max(s(:))*tol;
toosmallidx = s<absthr;
zerono = numel(find(toosmallidx));
newphase(toosmallidx) = rand(zerono,1)*2*pi;

c=s.*exp(1i*newphase);
