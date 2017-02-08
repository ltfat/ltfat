function cfreq = cent_freqs(g,L)
%CENT_FREQS   Determine relative center frequencies
%   Usage:  cfreq = cent_freqs(g);
%           cfreq = cent_freqs(g,L);
%           cfreq = cent_freqs(fs,fc);
%           cfreq = cent_freqs(g,fc);
%
%   Input parameters:
%      g     : Set of filters.
%      L     : Signal length.
%      fs    : Sampling rate (in Hz).
%      fc    : Vector of center frequencies (in Hz).
%   Output parameters:
%      cfreq : Vector of relative center frequencies in ]-1,1].
%
%   `cent_freqs(g)` will compute the center frequencies of the filters 
%   contained in *g* by determining their circular center of gravity. To
%   that purpose, the transfer function of each filter will be computed for
%   a default signal length on 10000 samples. For improved accuracy, the 
%   factual signal length *L* can be supplied as an optional parameter.
%   Alternatively, the center frequencies can be obtained from a set of
%   center frequencies *fc* (in Hz) and the sampling rate *fs*. The
%   sampling rate can also be determined from the field *fs* of the filter
%   set *g*.
%
%   Note: If g.H contains full-length, numeric transfer functions, *L*
%   must be specified for correct results.

complainif_notenoughargs(nargin,1,'FILTERBANKCENTFREQS');

if nargin > 1 % Try to determine cfreq from fc and fs in Hz
    if numel(g) == 1 && isnumeric(g) && numel(L) > 1
        cfreq = modcent(2*L/g,2);
        return
    elseif numel(g) == numel(L) && ( numel(g) > 1 || isfield(g,'fs') ) 
        cfreq = cellfun(@(fcEl,gEl) modcent(2*fcEl/gEl.fs,2),num2cell(L),g.');
        return
    end
else
    L = 10000; % Default value
end

g = filterbankwin(g,1,L,'normal');

% Compute l1-normalized absolute value of the transfer functions 
gH = cellfun(@(gEl) comp_transferfunction(gEl,L),g,'UniformOutput',false);
gH = cellfun(@(gHEl) abs(gHEl)./norm(gHEl,1),gH,'UniformOutput',false);

% Compute circular center of gravity
circInd = exp(2*pi*1i*(0:L-1)/L).';
cfreq = cellfun(@(gHEl) sum(circInd.*gHEl),gH.');
cfreq = real((pi*1i)\log(cfreq));   

cfreq(isnan(cfreq)) = 0;
 
