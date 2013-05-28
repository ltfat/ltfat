function p = blockpanel(params)
%BLOCKPANEL Control panel
%   Usage: blockpanel(params)
%
%   Input parameters:
%      params: Cell-array of parameters specifications.
%
%   Output parameters:
%      p : Control panel Java object
%
%   `blockpanel(params)` creates a Java object containing GUI for changing
%   parameters during the playback. `params` should be a cell-array, whose 
%   elements are another cell array of the followong format:
%
%      {'var','label',minVal,maxVal,defVal,valCount}
%
%   Example:
%
%   params = {
%               {'G','Gain',-20,20,0,21}
%            }
%

if nargin<1
    error('%s: Too few input parameters.',upper(mfilename));
end

if ~iscell(params)
    error('%s: Input should be cell array.',upper(mfilename));
end

if ~iscell(params{1})
    params = {params};
end

p = javaObject('net.sourceforge.ltfat.ContFrame');
paramList = javaObject('java.util.LinkedList');
    

for ii = 1:numel(params)
   param = params{ii};
   paramListEl = javaObject('java.util.LinkedList');
   for jj=1:numel(param)
        javaMethod('add',paramListEl,param{jj});
   end
   javaMethod('add',paramList,paramListEl);
    
   
end
 javaMethod('addControlElements',p,paramList);