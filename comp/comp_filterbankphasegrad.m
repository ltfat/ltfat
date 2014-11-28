function [fgrad,tgrad,c_s,c] = comp_filterbankphasegrad(f,g,hg,dg,a,minlvl)
%
%


L = size(f,1);

c=comp_filterbank(f,g,a); 
% Compute filterbank coefficients with frequency weighted window
c_h=comp_filterbank(f,hg,a);
% Compute filterbank coefficients with time weighted window
c_d=comp_filterbank(f,dg,a);

% Compute spectrogram and
% remove small values because we need to divide by c_s
temp = cell2mat(c);
minlvl = minlvl*max(abs(temp(:)).^2);
c_s = cellfun(@(x) max(abs(x).^2,minlvl),c,'UniformOutput',false);

% Compute instantaneous frequency
tgrad=cellfun(@(x,y,z) real(x.*conj(y)./z)/L*2,c_h,c,c_s,'UniformOutput',false);
%fgrad=cellfun(@(x) x/L*2,fgrad,'UniformOutput',false);
% Limit 
tgrad = cellfun(@(fEl) fEl.*(abs(fEl)<=2) ,tgrad,'UniformOutput',0);

% Compute group delay
fgrad=cellfun(@(x,y,z) imag(x.*conj(y)./z),c_d,c,c_s,'UniformOutput',false);