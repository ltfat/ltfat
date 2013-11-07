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

if ~iscell(params) || isempty(params)
    error('%s: Input should be a nonempty cell array.',upper(mfilename));
end

if ~iscell(params{1})
    params = {params};
end

try
  p = javaObject('net.sourceforge.ltfat.ContFrame');
catch% err
   error(['%s: Could not load net.sourceforge.ltfat.ContFrame. It is not ',...
          'compiled or it is not in Matlab classpath. In the latter case, ',...
          'ltfatstart should do the trick.'],upper(mfilename));
end

% Using Java LinkedList class for passing the cell-array
% (since there is no way how to pass the cell-array directly)
paramList = javaObject('java.util.LinkedList');
    

for ii = 1:numel(params)
   param = params{ii};
   if numel(param)<6
      error('%s: Parameter %i is not in format {''var'',''label'',minVal,maxVal,defVal,valCount}.',upper(mfilename),ii);
   end
   if param{3}>=param{4}
      error('%s: In parameter %i: minVal cannot be greater or equal to maxVal.',upper(mfilename),ii);
   end
   
   if param{5}<param{3} || param{5}>param{4}
      error('%s: In parameter %i: defVal is not in range minVal-maxVal.',upper(mfilename),ii);
   end
   
   if param{6}<=1
      error('%s: In parameter %i: valCount has to be >=2.',upper(mfilename),ii);
   end
   
   % Each element of the linked list paramList is again a linked list
   paramListEl = javaObject('java.util.LinkedList');
   for jj=1:numel(param)
        javaMethod('add',paramListEl,param{jj});
   end
   javaMethod('add',paramList,paramListEl);
    
   
end

% Pass the data
javaMethod('addControlElements',p,paramList);
 
% Give the object time to inilialize properly.
pause(0.1);
 