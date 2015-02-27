function [tgrad,fgrad,cs] = comp_filterbankphasegrad(c,ch,cd,L,minlvl)
%
%

% Compute spectrogram and
% remove small values because we need to divide by cs
temp = cell2mat(c);
minlvl = minlvl*max(abs(temp(:)).^2);
cs = cellfun(@(x) max(abs(x).^2,minlvl),c,'UniformOutput',false);

% Compute instantaneous frequency
fgrad=cellfun(@(x,y,z) real(x.*conj(y)./z)/L*2,ch,c,cs,'UniformOutput',false);
%fgrad=cellfun(@(x) x/L*2,fgrad,'UniformOutput',false);
% Limit 
fgrad = cellfun(@(fEl) fEl.*(abs(fEl)<=2) ,fgrad,'UniformOutput',0);

% Compute group delay
tgrad=cellfun(@(x,y,z) imag(x.*conj(y)./z),cd,c,cs,'UniformOutput',false);
