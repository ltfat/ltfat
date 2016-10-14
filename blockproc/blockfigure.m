function p = blockfigure(varargin)
%BLOCKFIGURE Block figure object
%   Usage: p=blockfigure();
%          p=blockfigure('cm',cmaps);
%
%   Output parameters:
%         p     : JAVA object of the class net.sourceforge.ltfat.SpectFrame
%
%   `p=blockfigure()` initializes a JAVA object for the |blockplot|
%   routine. 
%
%   Optional key-value pairs:
%
%      `'cm',cmaps`   : Custom colormap (a $L \times 3$ matrix) or colormaps
%                       (cell array of matrixes).
%
%   The function takes in the additional optional arguments:
%
%       `'location',location`:   Window inital position. `location`
%                                has to be 2 element row vector `[x,y]`
%                                defining distance from the top-left
%                                corner of the screen. 

definput.keyvals.cm=[];
definput.keyvals.location=[50,120];
[flags,kv]=ltfatarghelper({},definput,varargin);

complainif_isjavaheadless('BLOCKFIGURE');
p = [];
try
   p = javaObject('net.sourceforge.ltfat.SpectFrame');
   javaMethod('setLocation',p,kv.location(1),kv.location(2));
catch% err
   error(['%s: Could not load net.sourceforge.ltfat.SpectFrame. It is not ',...
          'compiled or it is not in Matlab classpath. In the latter case, ',...
          'ltfatstart should do the trick.'],upper(mfilename));
end



if ~isvector(kv.location) || any(size(kv.location)~=[1,2]) ||...
    any(kv.location<0) || ~isreal(kv.location)
   error(['%s: Location has to be a 2 element row vector of ',...
         ' positive numbers.'],upper(mfilename));  
end


% Gives p time to be correctly initialized. Otherwise, calling methods
% of the object can fail.
pause(0.1);


if isempty(kv.cm)
   % Create a default colormap.
   kv.cm = {ltfat_inferno()};
elseif isnumeric(kv.cm)
   % Enclose the colormap to a single element cell array.
   kv.cm = {kv.cm};
elseif iscell(kv.cm)
   if numel(kv.cm)>1
      error('%s: TO DO: More than one colormap is not supported yet.',upper(mfilename));
   end
end

% Set the colormap.
if isoctave
   javaMethod('setColormap',p,kv.cm{1}(:),size(kv.cm{1},1),size(kv.cm{1},2));
else
   javaMethod('setColormap',p,kv.cm{1});
end

javaMethod('show',p);




