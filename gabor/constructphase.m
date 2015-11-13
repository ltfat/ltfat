function [c,newphase,tgrad,fgrad]=constructphase(s,g,a,varargin)
%CONSTRUCTPHASE  Construct the phase of a DGT
%   Usage:  c=constructphase(s,g,a);
%           c=constructphase(s,g,a,tol);
%           c=constructphase(c,g,a,M,tol,mask);
%           [c,newphase,tgrad,fgrad] = constructphase(...);
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
%   This function requires a computational subroutine that is only
%   available in C. Use |ltfatmex| to compile it.
%
%   See also:  dgt, gabphasegrad, ltfatmex
%

% AUTHOR: Peter L. Søndergaard, Zdenek Prusa

thismfilename = upper(mfilename);
complainif_notposint(a,'a',thismfilename);
complainif_notposint(M,'M',thismfilename);

definput.keyvals.tol=1e-10;
definput.keyvals.mask=[];
[~,~,tol,mask]=ltfatarghelper({'tol','mask'},definput,varargin);

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

if ~iscalar(tol)
    error('%s: tol must be scalar.',thismfilename);
end

abss = abs(s);
% Compute phase gradients, check parameteres
[tgrad,fgrad] = gabphasegrad('abs',abss,g,a,2);

absthr = max(abss(:))*tol;

if isempty(mask)
    % Build the phase (calling a MEX file)
    newphase=comp_heapint(abss,tgrad,fgrad,a,tol);
    % Set phase of the coefficients below tol to random values
    toosmallidx = abss<absthr;
    zerono = numel(find(toosmallidx));
    newphase(toosmallidx) = rand(zerono,1)*2*pi;
else
    newphase=comp_maskedheapint(s,tgrad,fgrad,mask,a,tol);
    % Set phase of small coefficient to random values
    % but just in the missing part
    missingidx = find(mask==0);
    toosmallidx = abss(missingidx)<absthr;
    zerono = numel(find(toosmallidx));
    newphase(missingidx(toosmallidx)) = rand(zerono,1)*2*pi;    
end

% Combine the magnitude and phase
c=abss.*exp(1i*newphase);

