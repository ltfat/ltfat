function FfD = comp_fudual(F, L)
%FUDUAL dual fusion frame
%
%   Input parameters:
%     F     frame of which to calculate the dual
%
%   Output parameters:
%     FfD   the dual frame
%
%   `fudual` calculates the dual of the fusion frame operator.
%   it differs from framedual, which calculates the dual of each single
%   frame and inverts the weights to create a fusion frame (and stores
%   them in a separate variable
%
% Started 15.01.2022
% Last 22.02.2022

complainif_notvalidframeobj(F,'Fudual');

F = comp_checkfudim(F, L);

if ~strcmp(F.type, 'fusion')
    error('fudual only works for fusion frames');
end

if isfield(F, 'frameoperator')
    Sf = F.frameoperator;
else
   Id = eye(F.cdim);
   c = comp_fuana(F,Id);
   Sf = comp_fusyn(F,c);
   F.frameoperator = Sf;
end


for ii=1:F.Nframes
    G = frsynmatrix(F.frames{ii}, F.cdim);
    G = cell2mat(Sf)\G; %supposedly more efficient than G = inv(Sf)*G;
    eval(sprintf("G_%i = frame('gen',G);",ii));
    if ii == 1 
         fusionstring = "G_1";
    else
        fusionstring = sprintf("%s, G_%i",fusionstring, ii);
    end
end
%weights could probably be passed directly
FfD = eval(sprintf('frame(''fusion'',1,%s)',fusionstring));
FfD.w = F.w; %weights are just the same
