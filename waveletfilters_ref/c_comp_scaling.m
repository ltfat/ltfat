function [a, scal, delayvec] = c_comp_scaling(afull, delayvec, lowpass_at_zero, flags)

%scal=sqrt(afull(:,1)./afull(:,2));
scal = afull.^(0.5);
if flags.do_real && lowpass_at_zero
    % Scale the lowpass filters
    scal(1)=scal(1)/sqrt(2);
    a = afull;
    delayvec = delayvec;
elseif flags.do_complex
    % Replicate the scaling, sampling rates and delays, except at zero
    % frequency
    if lowpass_at_zero
       a=[afull;flipud(afull(2:end,:))];
       scal=[scal;flipud(scal(2:end))];
       delayvec=[delayvec;flipud(delayvec(2:end))];
    else
        a=[afull;flipud(afull)];
        scal=[scal;flipud(scal)];
        delayvec=[delayvec;flipud(delayvec)];
    end
else
    a = afull;
    delayvec = delayvec;
end