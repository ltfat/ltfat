function [tgrad,fgrad,s,c]=filterbankphasegrad(f,g,a,L,minlvl)
%FILTERBANKPHASEGRAD   Phase gradient of a filterbank representation
%   Usage:  [tgrad,fgrad,s,c] = filterbankphasegrad(f,g,a,L,minlvl);
%           [tgrad,fgrad,s,c] = filterbankphasegrad(f,g,a,L);
%           [tgrad,fgrad,s,c] = filterbankphasegrad(f,g,a,minlvl);
%           [tgrad,fgrad,s,c] = filterbankphasegrad(f,g,a);
%           [tgrad,fgrad,s] = filterbankphasegrad(...)
%           [tgrad,fgrad]  = filterbankphasegrad(...)
% 
%   Input parameters:
%      f     : Signal to be analyzed.
%      g     : Cell array of filters
%      a     : Vector of time steps.
%      L     : Signal length (optional).
%      minlvl: Regularization parameter (optional, required < 1).
%   Output parameters:
%      tgrad : Instantaneous frequency relative to original position.
%      fgrad : The negative of the local group delay. 
%      cs    : Filterbank spectrogram.
%      c     : Filterbank coefficients.
%
%   `[tgrad,fgrad,s,c] = filterbankphasegrad(f,g,a,L)` computes the 
%   relative instantaneous frequency *tgrad* and the negative of the group
%   delay *fgrad* of the filterbank spectrogram *s* obtained from the 
%   signal *f* and filterbank parameters *g* and *a*. 
%   Both *tgrad* and *fgrad* are specified relative to the original 
%   coefficient position entirely similar to |gabphasegrad|.
%   *fgrad* is given in samples, while *tgrad* is given in normalised
%   frequencies such that the absolute frequencies are in the range of ]-1,1]. 
%
%   This routine uses the equivalence of the filterbank coefficients in 
%   each channel with coefficients obtained from an STFT obtained with a
%   certain window (possibly different for every channel). As a consequence
%   of this equivalence, the formulas derived in the reference apply. 
%
%   See also: gabphasegrad
%
%   References: aufl95 ltfatnote041

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
[gh, gd, g] = comp_phasegradfilters(g, asan, L);

f=postpad(f,L);

c=comp_filterbank(f,g,asan); 
% Compute filterbank coefficients with frequency weighted window
ch=comp_filterbank(f,gh,asan);
% Compute filterbank coefficients with time weighted window
cd=comp_filterbank(f,gd,asan);

% Run the computation
[tgrad,fgrad,s] = comp_filterbankphasegrad(c,ch,cd,L,minlvl);
