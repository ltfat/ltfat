function sr=gabreassignreal(s,tgrad,fgrad,a,M)
%GABREASSIGNREAL  Reassign time-frequency distribution for real signals
%   Usage:  sr = gabreassignreal(s,tgrad,fgrad,a,M);
%
%   `gabreassignreal(s,tgrad,fgrad,a,M)` reassigns the values of the positive
%   time-frequency distribution *s* using the phase gradient given by *fgrad*
%   and *tgrad*. The lattice is determined by the time shift *a* and the 
%   number of channels *M*.
%
%   *fgrad* and *tgrad* can be obtained by the routine |gabphasegrad|.
%
%   Examples:
%   ---------
%
%   The following example demonstrates how to manually create a
%   reassigned spectrogram. An easier way is to just call |resgram|:::
%
%     % Create reassigned vector field of the bat signal.
%     a=4; M=100;
%     [phased,c] = gabphasederivreal({'t','f'},'dgt',bat,'gauss',a,M,'relative');
%     [tgrad, fgrad] = deal(phased{:});
%
%     % Perform the actual reassignment
%     sr = gabreassignreal(abs(c).^2,tgrad,fgrad,a,M);
%
%     % Display it using plotdgt
%     plotdgt(sr,a,143000,50);
%  
%   See also: resgram, gabphasederivreal, gabreassign
%
%   References: aufl95

% AUTHOR: Peter L. Søndergaard, 2008.
%         Nicki Holighaus, 2023.


thisname = upper(mfilename);
complainif_notenoughargs(nargin,5,thisname);
complainif_notposint(a,'a',thisname);
complainif_notposint(M,'M',thisname);


% Basic checks
if any(cellfun(@(el) isempty(el) || ~isnumeric(el),{s,tgrad,fgrad}))
    error('%s: s, tgrad, fgrad must be non-empty and numeric.',...
          upper(mfilename));
end

% Check if argument sizes are consistent
if ~isequal(size(s),size(tgrad),size(fgrad))
   error('%s: s, tgrad, fgrad must all have the same size.',...
          upper(mfilename));
end

% Check if any argument is not real
if any(cellfun(@(el) ~isreal(el),{tgrad,fgrad}))
   error('%s: tgrad, fgrad must be real.',...
          upper(mfilename));
end

% if any(s<0)
%     error('%s: s must contain positive numbers only.',...
%         upper(mfilename));
% end

sr=comp_gabreassignreal(s,tgrad,fgrad,a,M);