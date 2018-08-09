function [varargincell,winCell] = arghelper_filterswinparser(windows,varargincell)

% Search for window given as cell array
candCellId = cellfun(@(vEl) iscell(vEl) && any(strcmpi(vEl{1},windows)),varargincell);

winCell = {};
% If there is such window, replace cell with function name so that 
% ltfatarghelper does not complain
if ~isempty(candCellId) && any(candCellId)
    candCellIdLast = find(candCellId,1,'last');
    winCell = varargincell{candCellIdLast};
    varargincell(candCellId) = []; % But remove all
    varargincell{end+1} = winCell{1};
end