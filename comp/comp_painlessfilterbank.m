function gout = comp_painlessfilterbank(g,a,L,type,do_real)
%COMP_PAINLESSFILTERBANK
% 
%   Function computes filterbank dual or tight frame for the painless case
%

M = numel(g);

F=comp_filterbankresponse(g,a,L,do_real);

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
    elseif isfield(g{m},'h')
       H=comp_transferfunction(g{m},L)./F; 
       thisgd.h = ifft(H);
       thisgd.offset = 0;
    end

    gout{m}=thisgd;
end;