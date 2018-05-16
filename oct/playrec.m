%PLAYREC Multichannel non-blocking audio recording and playback
%
% Run: 
%     playrec('help') 
%
% Tho get the actual help.
%
% PLAYREC is a MEX file which means it has to be compiled in order to work.
% Please see INSTALL-Matlab or INSTALL-Octave for details. 
% 
% Originally created by Robert Humphrey (see license_playrec.txt) and 
% currently hosted at https://github.com/PlayrecForMatlab/playrec
% For further details visit the webpage or run playrec('about') and
% playrec('overview')
%
% This is however a sightly modified version doing an on-the-fly
% sample rate conversion if the target sample rate is not supported by
% the device.  
% 
