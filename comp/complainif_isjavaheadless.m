function complainif_isjavaheadless(callfun)
% COMPLAINIF_ISJAVAHEADLESS 
%
%   Prints warning if the available JRE ius in headless mode.

if nargin<1
    callfun=mfilename;
end

try
   ge = javaMethod('getLocalGraphicsEnvironment','java.awt.GraphicsEnvironment');
catch
    % No java support at all. Either we are running matlab -nojvm
    % or octave(<3.8.0) without octave-java package.
    % Both cases shoud have already been caught somewhere.   
    return;
end

if javaMethod('isHeadless',ge)
       error(['%s: JRE is available in headless mode only. ',...
              'Block processing GUI will not work. Consider ',...
              'installing full JRE.'],upper(callfun));
end
