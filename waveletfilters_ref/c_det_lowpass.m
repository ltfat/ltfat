function [aprecise, M2, lowpass_number, lowpass_at_zero] = c_det_lowpass(Ls, scales, basea, flags, kv)
%get the number of lowpasses
%aprecise = [scales + #lowpass, 1]
%M2 = M + #lowpass
%#lowpass
%lowpass_at_zero...flag


if numel(scales) < 4 && flags.do_single
    error('%s: Lowpass generation requires at least 4 scales.',upper(mfilename));
elseif numel(scales) < 2 && flags.do_repeat
    error('%s: Lowpass generation requires at least two scales.',upper(mfilename));
end

% Get number of scales and sort them
M = numel(scales);
scales_sorted = sort(scales,'descend');
%% Determine total number of filters and natural subsampling factor for lowpass
if flags.do_repeat 
    lowpass_number = scales_sorted(2)/(scales_sorted(1)-scales_sorted(2)); % Maybe adjust this to not guarantee some distance between first filter and zero frequency.
    if abs(lowpass_number - round(lowpass_number)) < eps*10^3
        lowpass_number = round(lowpass_number);
        lowpass_at_zero = 1;
    else
        lowpass_at_zero = 0;
    end
    lowpass_number = floor(lowpass_number);
    if lowpass_number == 0
        lowpass_number = 1;
    end
    M2 = M + lowpass_number;
    aprecise = (basea.*scales_sorted(1))*ones(lowpass_number,1);
elseif flags.do_single
    lowpass_number = 1;
    lowpass_at_zero = 1;
    M2 = M+1;
    aprecise = (0.2./scales_sorted(4))*Ls; % This depends on how single lowpass is called (l.195ff). Maybe automate. Adapt if necessary!!!
else
    lowpass_number = 0;
    lowpass_at_zero = 0;
    M2 = M;
    aprecise = [];
end

%% Get subsampling factors
aprecise = [aprecise;basea.*scales];

if any(aprecise<1)
    error(['%s: Bandwidth of at least one of the filters is bigger than fs. '],upper(mfilename));
end

aprecise=aprecise/kv.redmul;
if any(aprecise<1)
    error('%s: The maximum redundancy mult. for this setting is %5.2f',...
         upper(mfilename), min(basea./scales));
end