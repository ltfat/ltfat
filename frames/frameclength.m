function [Ncoef, L]=frameclength(F,Ls)
%FRAMECLENGTH  Number of coefficients from length of signal
%   Usage: Ncoef=frameclength(F,Ls);
%          [Ncoef,L]=frameclength(...);
%
%   `Ncoef=frameclength(F,Ls)` returns the total number of coefficients 
%   obtained by applying the analysis operator of frame *F* to a signal
%   of length *Ls* i.e. `size(frana(F,f),1)` for `Ls=length(f)`. 
%
%   `[Ncoef,L]=frameclength(F,Ls)` additionally returns *L*, which is the 
%   same as returned by |framelength|.
%
%   If the frame length *L* is longer than the signal length *Ls*, the 
%   signal will be zero-padded to *L* by |frana|.
%
%   See also: frame, framelengthcoef

callfun = upper(mfilename);
complainif_notposint(Ls,'Ls',callfun);
complainif_notvalidframeobj(F,callfun);

L = F.length(Ls);

% Some frames need special function
if isfield(F,'clength')
    Ncoef = F.clength(L);
else
    % Generic, works for any non-realonly frame and for
    % all representaions not having any extra coefficients

    Ncoef = L*F.red;
    
    if F.realinput
        Ncoef=Ncoef/2;
    end

    assert(abs(Ncoef-round(Ncoef))<1e-3,...
           sprintf('%s: There is a bug. L=%d should be an integer.',...
           upper(mfilename),Ncoef));

    Ncoef=round(Ncoef);
end
