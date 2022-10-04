function H=comp_transferfunction(g,L)
%COMP_TRANSFERFUNCTION  Compute the transfer function
%
%  `comp_transferfunction(g,L)` computes length *L* transfer function 
%  (frequency response) of a single filter *g*. This function can only
%  handle filters in a proper internal format i.e. already processed by
%  |filterbankwin|.

% Setting crossover to 0 ensures FIR filters to be transformed to 
% full-length Frequency-domain defined filters with g.H and g.foff fields.
g = comp_filterbank_pre({g},1,L,0);
% Band-limited filters have to be made full-length
H = circshift(postpad(g{1}.H(:),L),g{1}.foff);

% Realonly has to be treated separately for band-limited filters
if isfield(g,'realonly') && g.realonly
     H=(H+involute(H))/2;
end;

