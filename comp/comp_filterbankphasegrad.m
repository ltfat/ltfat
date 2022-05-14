function [tgrad,fgrad,s] = comp_filterbankphasegrad(c,ch,cd,L,minlvl)
%
%this function is called by filterbankconstphase

% Compute spectrogram and
% remove small values because we need to divide by cs
temp = cell2mat(c);
minlvl = minlvl*max(abs(temp(:)).^2);
s = cellfun(@(x) max(abs(x).^2,minlvl),c,'UniformOutput',false);

% Compute instantaneous frequency
tgrad=cellfun(@(x,y,z) real(x.*conj(y)./z)/L*2,cd,c,s,'UniformOutput',false);

% Limit 
tgrad = cellfun(@(fEl) fEl.*(abs(fEl)<=2) ,tgrad,'UniformOutput',0);

% Compute group delay
fgrad=cellfun(@(x,y,z) imag(x.*conj(y)./z),ch,c,s,'UniformOutput',false);
