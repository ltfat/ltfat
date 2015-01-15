function [fgrad,tgrad,c_s,c]=filterbankphasegrad(f,g,a,L,minlvl)
%FILTERBANKPHASEGRAD   Phase gradient of a filterbank representation
%   Usage:  [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,L,minlvl);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,L);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a,minlvl);
%           [tgrad,fgrad,c_s,c] = filterbankphasegrad(f,g,a);
%           [tgrad,fgrad,c_s] = filterbankphasegrad(...)
%           [tgrad,fgrad]  = filterbankphasegrad(...)
% 
%   Input parameters:
%      f     : Signal to be analyzed.
%      g     : Cell array of filters
%      a     : Vector of time steps.
%      L     : Signal length (optional).
%      minlvl: Regularization parameter (optional, required < 1).
%   Output parameters:
%      fgrad : Instantaneous frequency relative to original position.
%      tgrad : Group delay relative to original position.
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

complainif_notenoughargs(nargin,3,'FILTERBANKPHASEGRAD');

% Reshape input signal
[f,~,W]=comp_sigreshape_pre(f,'FILTERBANKPHASEGRAD',0);
Ls = size(f,1);

if W>1
    error('%s: Only one-channel signals supported.',upper(mfilename));
end

if nargin < 5 
    if nargin < 4
        L = filterbanklength(Ls,a);
        minlvl = eps;
    else
        if ~(isscalar(L) && isnumeric(L) ) && L>0
            error('%s: Fourth argument shoud be a positive number.',...
                  upper(mfilename));
        end
        if L >= 1
            minlvl = eps;
        else
            minlvl = L;
        end;
    end
end;

complainif_notposint(L,'L','FILTERBANKPHASEGRAD');

Luser = filterbanklength(L,a);
if Luser~=L
    error(['%s: Incorrect transform length L=%i specified. ', ...
           'Next valid length is L=%i. See the help of ',...
           'FILTERBANKLENGTH for the requirements.'],...
           upper(mfilename),L,Luser);
end


% Unify format of coefficients
[g,asan]=filterbankwin(g,a,L,'normal');

% Precompute filters
[hg, dg, g] = comp_phasegradfilters(g, asan, L);

f=postpad(f,L);

% Run the computation
[fgrad,tgrad,c_s,c] = comp_filterbankphasegrad(f,g,hg,dg,asan,minlvl);
