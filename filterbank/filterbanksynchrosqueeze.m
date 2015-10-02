function [cr,repos,Lc]=filterbanksynchrosqueeze(c,tgrad,var)
%FILTERBANKREASSIGN  Reassign filterbank spectrogram
%   Usage:  cr = filterbankreassign(c,tgrad,cfreq);
%           cr = filterbankreassign(c,tgrad,g);
%           [cr,repos,Lc] = filterbankreassign(...);
%
%   Input parameters:
%      c     : Coefficients to be synchrosqueezed.
%      tgrad : Instantaneous frequency relative to original position.
%      cfreq : Vector of relative center frequencies in ]-1,1].
%      g     : Set of filters.
%   Output parameters:
%      cr    : Synchrosqueezed filterbank coefficients.
%      repos : Reassigned positions.
%      Lc    : Subband lengths.
%
%   `filterbanksynchrosqueeze(c,tgrad,cfreq)` will reassign the values of 
%   the filterbank coefficients *c* instantaneous frequency *tgrad*. 
%   The frequency center frequencies of filters are given by *cfreq*.
%   The filterbank coefficients *c* are assumed to be obtained from a
%   non-subsampled filterbank (a=1).
%
%   `filterbanksynchrosqueeze(s,tgrad,g)` will do the same thing except
%   the center frequencies are estimated from a set of filters *g*.
%
%   `[sr,repos,Lc]=filterbanksynchrosqueeze(...)` does the same thing, but 
%   in addition returns a vector of subband lengths *Lc* (`Lc = cellfun(@numel,s)`)
%   and cell array *repos* with `sum(Lc)` elements. Each element corresponds 
%   to a single coefficient obtained by `cell2mat(sr)` and it is a vector 
%   of indices identifying coefficients from `cell2mat(s)` assigned to 
%   the particular time-frequency position.
%
%   The arguments *s*, *tgrad* must be cell-arrays of vectors
%   of the same lengths. Arguments *cfreq* or *g* must have the
%   same number of elements as the cell arrays with coefficients.
%
%   Examples:
%   ---------
%
%   This example shows how to reassign a ERB filterbank spectrogram:::
%
%     % Genrate 3 chirps 1 second long
%     L = 44100; fs = 44100; l = 0:L-1;
% 
%     f = sin(2*pi*(l/35+(l/300).^2)) + ...
%         sin(2*pi*(l/10+(l/300).^2)) + ...
%         sin(2*pi*(l/5-(l/450).^2));
%     f = 0.7*f';
%     
%     % Create ERB filterbank
%     [g,a,fc]=erbfilters(fs,L,'fractional','spacing',1/12,'warped');
%     
%     % Compute phase gradient
%     [tgrad,fgrad,cs,c]=filterbankphasegrad(f,g,a);
%     % Do the reassignment
%     sr=filterbankreassign(cs,tgrad,fgrad,a,cent_freqs(fs,fc));
%     figure(1); subplot(211);
%     plotfilterbank(cs,a,fc,fs,60);
%     title('ERBlet spectrogram of 3 chirps');
%     subplot(212);  
%     plotfilterbank(sr,a,fc,fs,60);
%     title('Reassigned ERBlet spectrogram of 3 chirps');
%
%   See also: filterbankphasegrad, gabreassign
%
%   References: ltfatnote041

%   AUTHOR : Nicki Holighaus.

% Sanity checks
complainif_notenoughargs(nargin,3,'FILTERBANKREASSIGN');

if isempty(c) || ~iscell(c) 
    error('%s: s should be a nonempty cell array.',upper(mfilename));
end

if isempty(tgrad) || ~iscell(tgrad) || any(~cellfun(@isreal,tgrad))
    error('%s: tgrad should be a nonempty cell array.',upper(mfilename));
end

if any(cellfun(@(sEl,tEl) ~isvector(sEl) || ~isvector(tEl) ...
                              , c,tgrad))
   error('%s: s, tgrad, must be cell arrays of numeric vectors.',...
         upper(mfilename)); 
end

if ~isequal(size(c),size(tgrad)) || ...
   any(cellfun(@(sEl,tEl) ~isequal(size(sEl),size(tEl)), ...
               c,tgrad))
   error('%s: s, tgrad does not have the same format.',upper(mfilename));   
end


W = cellfun(@(sEl)size(sEl,2),c);
if any(W>1)
   error('%s: Only one-channel signals are supported.',upper(mfilename)); 
end

% Number of channels
M = numel(c);

% Number of elements in channels
Lc = cellfun(@(sEl)size(sEl,1),c);

% Check if a comply with subband lengths
L = Lc;
if any(abs(L-L(1))>1e-6)
   error(['%s: Subsampling factors and subband lengths do not ',...
          'comply.'],upper(mfilename));
end
L = L(1);

% Determine center frequencies
if isempty(var) || numel(var)~=M || ~isvector(var) && ~iscell(var)
   error(['%s: cfreq must be length-M numeric vector or a cell-array ',...
          'containg M filters.'],upper(mfilename)); 
else
    if iscell(var)
       cfreq = cent_freqs(var,L);
    else
       cfreq = var;
    end
end

% Dummy fgrad
fgrad = tgrad;
for m=1:numel(fgrad);
    fgrad{m} = zeros(L,1);
end

a = ones(M,1);

% Do the computations
if nargout>1
   [cr,repos] = comp_filterbankreassign(c,tgrad,fgrad,a,cfreq);
else
   cr = comp_filterbankreassign(c,tgrad,fgrad,a,cfreq);
end

