function [a, L] = c_comp_downsampling(Ls, M2, scales, aprecise, lowpass_at_zero, flags, kv)
%account for downsampling type and target redundancy

if flags.do_regsampling % This should only be used for lowpass = single!!!
    a = ones(M2,1);
    
    [lower_scale,~] = max(scales);
    [upper_scale,~] = min(scales);
    lower_scale = floor(log2(1/lower_scale));
    upper_scale = floor(log2(1/upper_scale));
    
    % Find minimum a in each octave and floor23 it.
    ct=1;
    for kk = lower_scale:upper_scale
        tempidx = find( floor(log2(1./scales)) == kk );
        [~,tempminidx] = min(1/scales(tempidx));
        idx = tempidx(tempminidx);
        
        % Deal the integer subsampling factors
        a(tempidx) = floor23(aprecise(idx));
        ct=ct+1;
    end   
    
    % Determine the minimal transform length lcm(a)
    L = filterbanklength(Ls,a);
    
    % Heuristic trying to reduce lcm(a)
    while L>2*Ls && ~(all(a==a(1)))
        maxa = max(a);
        a(a==maxa) = 0;
        a(a==0) = max(a);
        L = filterbanklength(Ls,a);
    end

elseif flags.do_fractional
    L = Ls;
    N=ceil(Ls./aprecise);
    a=[repmat(Ls,M2,1),N];
elseif flags.do_fractionaluniform
    L = Ls;
    if lowpass_at_zero
        aprecise(2:end)= min(aprecise(2:end));
    else 
        aprecise= repmat(min(aprecise),numel(aprecise),1);
    end
    N=ceil(Ls./aprecise);
    a=[repmat(Ls,M2,1),N];
elseif flags.do_uniform
    a=floor(min(aprecise));
    L=filterbanklength(Ls,a);
    a = repmat(a,M2,1);
end
afull=comp_filterbank_a(a,M2,struct());


%==========================================================================
%% Adjust the downsampling rates in order to achieve 'redtar'    

if ~isempty(kv.redtar)
   % if size(afull,2) == 2
        a = afull(:,1)./afull(:,2);
   % else
   %     a = afull;
   % end

    if ~flags.do_real
        org_red = sum(1./a);
    elseif lowpass_at_zero
        org_red = 1./a(1) + sum(2./a(2:end));
    else
        org_red = sum(2./a);
    end
    
    a = floor(a*org_red/kv.redtar);
    a(a==0) = 1;
    
    if ~flags.do_uniform
        N_new=ceil(L./a);
        if flags.do_complex
            N_new = [N_new;N_new(end:-1:2)];
        end
        a=[repmat(L,numel(N_new),1),N_new];
    else 
        L = filterbanklength(L,a);
        a=[repmat(L,numel(a),1),ones(numel(a), 1)];
    end
else
    a = afull;
end

%if size(a,2) == 2
%    a = a(:,1)./a(:,2);
%end