function delayvec = c_comp_delay(kv, M2, afull)

if isa(kv.delay,'function_handle')
    delayvec = zeros(M2,1);
    for kk = 1:M2
        delayvec(kk) = kv.delay(kk-1,afull(kk,1)./afull(kk,2));
    end
elseif numel(kv.delay) == 1
    delayvec = repmat(kv.delay,M2,1);
elseif any(size(kv.delay,2)) == 1 && numel(kv.delay) >= M2
    delayvec = kv.delay(:);
else
    error('%s: delay must be scaler or have enough elements to cover all channels.',upper(mfilename));
end