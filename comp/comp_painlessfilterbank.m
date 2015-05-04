function gout = comp_painlessfilterbank(g,a,L,type,do_real)
%COMP_PAINLESSFILTERBANK
% 
%   Function computes filterbank dual or tight frame for the painless case.
%

M = numel(g);

F = comp_filterbankresponse(g,a,L,do_real);
% inf here ensures FIR filters will stay FIR
g = comp_filterbank_pre(g,a,L,inf);

if strcmpi(type,'tight')
    F = sqrt(F);  
elseif strcmpi(type,'dual')
    % Do nothing
else
    %Fail
    error('%s: Internal error. Unrecognized frame type.',upper(mfilename));
end

gout=cell(1,M);
for m=1:M
    thisgd=struct();
    if isfield(g{m},'H')
       H=circshift(comp_transferfunction(g{m},L)./F,-g{m}.foff);
       thisgd.H=H(1:numel(g{m}.H));
       thisgd.foff=g{m}.foff;
       thisgd.realonly=0;
       thisgd.delay=0;
       thisgd.L = L;
    elseif isfield(g{m},'h')
       H=comp_transferfunction(g{m},L)./F; 
       thisgd = ifft(H);
    end

    gout{m}=thisgd;
end;

% Convert appropriate filters to structs with .h fields.
gId = cellfun(@(gEl) isfield(gEl,'h'),g);
if any(gId)
    gout(gId) = filterbankwin(gout(gId),a(gId));
end
