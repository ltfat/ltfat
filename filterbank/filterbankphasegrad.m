function [tgrad,fgrad,c_s,c]=filterbankphasegrad(f,g,a,L,minlvl)
%FILTERBANKPHASEGRAD   Phase gradient of a filterbank representation
%   Usage:  [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,L,minlvl);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,L);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,minlvl);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a);
%           [tgrad,fgrad,c_s] = ...
%           [tgrad,fgrad]  = ...
% 
%   Input parameters:
%      f     : Signal to be analyzed.
%      g     : Cell array of filters
%      a     : Vector of time steps.
%      L     : Signal length (optional).
%      minlvl: Regularization parameter (optional, required < 1).
%   Output parameters:
%      tgrad : Group delay relative to original position.
%      fgrad : Instantaneous frequency relative to original position.
%      c_s   : Filterbank spectrogram.
%      c     : Filterbank coefficients.
%
%   `[fgrad,tgrad,c_s,c] = filterbankphasegrad(f,g,a,L)` computes the group
%   delay *tgrad* and instantaneous frequency *fgrad* of the filterbank
%   spectrogram *c_s* obtained from the signal *f* and filterbank
%   parameters *g* and *a*. Both quantities are specified relative to the
%   original coefficient position. *tgrad* is given in samples, while
%   *fgrad* is given as values on the unit circle, easily converted into 
%   relative frequencies by $\log(tgrad)/(pi*i)$.
%   This routine uses the equivalence of the filterbank coefficients in 
%   each channel with coefficients obtained from an STFT obtained with a
%   certain window (possibly different for every channel). As a consequence
%   of this equivalence, the formulas derived in the reference apply. 
%
%   References: aufl95

%   AUTHOR : Nicki Holighaus.

if nargin < 5 
    if isempty(L) 
        Ls = size(f,1);
        L=filterbanklength(Ls,a);
        minlvl = eps;
    elseif L >= 1
        minlvl = eps;
    else
        minlvl = L;
    end;
end;

% Sanity checks
complainif_notposint(L,'L','FILTERBANKPHASEGRAD');

% Reshape input signal
[f,~,W]=comp_sigreshape_pre(f,'FILTERBANKPHASEGRAD',0);

if W>1
    error('%s: Only one-channel signals supported.',upper(mfilename));
end

% Unify format of coefficients
[g,info]=filterbankwin(g,a,L,'normal');

% Sanitize format of a
info.a = comp_filterbank_a(info.a, info.M);

% Precompute filters
[hg, dg, g] = comp_phasegradfilters(g, info.a, L);

f=postpad(f,L);

% Run the computation
[tgrad,fgrad,c_s,c] = comp_filterbankphasegrad(f,g,hg,dg,info.a,minlvl);
