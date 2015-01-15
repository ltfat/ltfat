function sr=filterbankreassign(s,tgrad,fgrad,a,var)
%FILTERBANKREASSIGN  Reassign filterbank spectrogram
%   Usage:  sr = filterbankreassign(s,a,tgrad,fgrad,cfreq);
%           sr = filterbankreassign(s,a,tgrad,fgrad,g);
%
%   Input parameters:
%      s     : Spectrogram to be reassigned.
%      tgrad : Group delay relative to original position.
%      fgrad : Instantaneous frequency relative to original position.
%      a     : Vector of time steps.
%      cfreq : Vector of relative center frequencies in ]-1,1].
%      g     : Set of filters.
%   Output parameters:
%      sr    : Reassigned filterbank spectrogram.
%
%   `filterbankreassign(s,a,fgrad,tgrad,cfreq)` will reassign the values of 
%   the filterbank spectrogram *s* using the group *tgrad* delay and 
%   instantaneous frequency  and *fgrad*. The time-frequency sampling 
%   pattern is determined from the time steps *a* and the center 
%   frequencies *cfreq*. 
%
%   `filterbankreassign(s,a,fgrad,tgrad,g)` will do the same thing except
%   the center frequencies are estimated from a set of filters *g*.
%
%   The arguments *s*, *tgrad* and *fgrad* must be cell-arrays of vectors
%   of the same lengths. Arguments *a* and *cfreq* or *g* must have the
%   same number of elements as the cell arrays with coefficients.
%   

%   AUTHOR : Nicki Holighaus.

% Sanity checks
complainif_notenoughargs(nargin,5,'FILTERBANKREASSIGN');

if isempty(s) || ~iscell(s) 
    error('%s: s should be a nonempty cell array.',upper(mfilename));
end

if isempty(tgrad) || ~iscell(tgrad) 
    error('%s: tgrad should be a nonempty cell array.',upper(mfilename));
end

if isempty(fgrad) || ~iscell(fgrad) 
    error('%s: fgrad should be a nonempty cell array.',upper(mfilename));
end

if any(cellfun(@(sEl,tEl,fEl) ~isvector(sEl) || ~isvector(tEl) || ...
                              ~isvector(fEl), s,tgrad,fgrad))
   error('%s: s, tgrad, frad must be cell arrays of numeric vectors.',...
         upper(mfilename)); 
end

if any(size(s)~=size(tgrad)) || any(size(s)~=size(fgrad)) || ...
   any(cellfun(@(sEl,tEl,fEl) any(size(sEl)~=size(tEl)) ...
                           || any(size(sEl)~=size(fEl)),...
                           s,tgrad,fgrad))
   error('%s: s, tgrad, frad does not have the same format.',upper(mfilename));   
end


W = cellfun(@(sEl)size(sEl,2),s);
if any(W>1)
   error('%s: Only one-channel signals are supported.',upper(mfilename)); 
end

% Number of channels
M = numel(s);

% Number of elements in channels
N = cellfun(@(sEl)size(sEl,1),s);

% Sanitize a
a=comp_filterbank_a(a,M);
a = a(:,1)./a(:,2);

% Check if a comply with subband lengths
L = N.*a;
if any(abs(L-L(1))>1e-6)
   error(['%s: Subsampling factors and subband lengths do not ',...
          'comply.'],upper(mfilename));
end    

% Determine center frequencies
if isempty(var) || numel(var)~=M || ~isvector(var) && ~iscell(var)
   error(['%s: cfreq must be length-M numeric vector or a cell-array ',...
          'containg M filters.'],upper(mfilename)); 
else
    if iscell(var)
       cfreq = cent_freqs(g,L);
    else
       cfreq = var;
    end
end

% Do the computations
sr = comp_filterbankreassign(s,tgrad,fgrad,a,cfreq);