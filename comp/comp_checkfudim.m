function F = comp_checkfudim(F, L)
%CHECKFUDIM   check Fusion Frame Dimension
%
%   Usage: checkfudim(G);
%
%   Input parameters:
%     F     frame
%     L     desired framelength
%
%   Output parameters:
%     F     frame with additional parameter .cdim
%
%
%     `checkfudim` checks if all local frames in a 
%     fusion frame share the same dimension.
%
%     Example:
%     F1 = [1 0; 0 1];
%     F2 = [1 0 0; 0 1 0; 0 0 1];
%     FF1 = frame('gen',F1);
%     FF2 = frame('gen',F2);
%     G = frame('fusion',1,FF1,FF2);
%     G = checkfudim(G);
%
%
%   Started: 23.12.2021
%   Current: 04.11.2022
%
% Author: P. Balazs


complainif_notvalidframeobj(F,'MatrixRep');

%if nargin==1
%L = 160;
%end

if isfield(F, 'cdim')
    return
elseif ~isfield(F, 'type') || ~strcmp(F.type, 'fusion')
    error('Checkfudim only works for fusion frames');
end

for ii=1:F.Nframes
    localframe = frameaccel(F.frames{ii}, framelength(F.frames{ii}, L));
    FF=frsynmatrix(localframe,localframe.L);
    M = size(FF, 1);
    if ii > 2 && M ~= Mold
        F.cdim = nan;
        break
    else
        Mold = M;
    end
end

F.cdim = Mold;