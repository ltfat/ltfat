function [proj,F] = comp_fuana(F,f,varargin)
%FUANA fusion frame analysis operator
%
%   Input parameters:
%     F        fusion frame
%     f        input signal
%
%   Output parameters:
%     proj   projection
%     F      loop back of fusion frames
%
%
%   `fuana` does the fusion frame analysis, i.e 
%   $ p_j = (w_j \pi_{W_j} f )_j.$, where $(f)$ 
%   is the signal/ function, $w_j$ teh weights, 
%   $W_j$ the subspaces (spanned by the local frames).
%
%
%   Optional parameters:
%     'range',r   vector or scalar, only calculate a specific range of projections, i.e. the projection on $W_j$
%
% Started: 23.12.2021
% Last: 8.11.2022
%
% Author: P. Balazs
 
complainif_notvalidframeobj(F,'Fuana');

if ~strcmp(F.type, 'fusion')
    error('fuana only works for fusion frames');
end


definput.keyvals.range = [];
[~,kv]=ltfatarghelper({},definput,varargin);

if isfield(F, 'localdual')
    FD = F.localdual;
else    
    FD = framedual(F);
end
L = length(f);
F = comp_checkfudim(F, L);

if isnan(F.cdim)
    error('Local frames have to have the same dimension (for now).');
elseif F.cdim ~= length(f)
    error('Local frames and signal do not fit.')
end

if ~isempty(kv.range)
    proj = cell(1, 1);
else
    kv.range = 1:F.Nframes;
    proj = cell(F.Nframes, 1);
end

cnt = 1;
for ii=1:numel(kv.range)
    jj = kv.range(ii);
    c = frana(F.frames{jj},f);
    c= flip(F.w(jj))*c;
    proj{cnt} = frsyn(FD.frames{jj}, c);    
    cnt = ii+1;
end
