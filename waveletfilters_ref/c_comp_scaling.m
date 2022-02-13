function [scal, a] = c_comp_scaling(afull, delayvec, lowpass_at_zero, flags)
%calculate the scaling of the filterbank

%if size(afull, 2) == 2
    scal=sqrt(afull(:,1)./afull(:,2));
%else
%    scal = afull.^(0.5);
%end

if flags.do_real && lowpass_at_zero
    % Scale the lowpass filters
    scal(1)=scal(1)/sqrt(2);
    a = afull;
    %delayvec = delayvec;
elseif flags.do_complex
    % Replicate the scaling, sampling rates and delays, except at zero
    % frequency
    if lowpass_at_zero
       a=[afull;flipud(afull(2:end,:))];
       scal=[scal;flipud(scal(2:end))];
       %delayvec=[delayvec;flipud(delayvec(2:end))];
    else
        a=[afull;flipud(afull)];
        scal=[scal;flipud(scal)];
        %delayvec=[delayvec;flipud(delayvec)];
    end
else
    a = afull;
    %delayvec = delayvec;
end