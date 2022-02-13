function delayvec = c_comp_delay(kv, M2, afull, flags, lowpass_at_zero)
%if delay is a function handle, assign delay for each filter
%if delay is a scalar, make a vector from it

if isa(kv.delay,'function_handle')
    delayvec = zeros(M2,1);
    for kk = 1:M2
        delayvec(kk) = kv.delay(kk-1,afull(kk,1)./afull(kk,2));
    end
elseif numel(kv.delay) == 1
    delayvec = repmat(kv.delay,M2,1);
%elseif any(size(kv.delay,2)) && numel(kv.delay) >= M2
%    delayvec = kv.delay(:);
elseif ~isempty(kv.delay) && size(kv.delay,2) > 1
    delayvec = kv.delay(:);
else
    error('%s: delay must be scaler or have enough elements to cover all channels.',upper(mfilename));
end

% if flags.do_complex
%     if lowpass_at_zero
%         delayvec=[delayvec;flipud(delayvec(2:end))];
%     else
%         delayvec=[delayvec;flipud(delayvec)];
%     end
% end
