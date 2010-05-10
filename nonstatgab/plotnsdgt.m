function [] = plotnsdgt(c,a_new,dynrange,sr)
%PLOTNSDGT Plot spectrogram from nonstationary Gabor coefficients
%   Usage:  plotnsdgt(c,a,dynrange,sr);
%
%   Input parameters:
%         c        : Cell array of coefficients.
%         a        : Vector of time positions of windows.
%         dynrange : Colorscale dynamic range in dB (default 60 dB).
%         sr       : signal sample rate in Hz (default 1 Hz).
%
%   PLOTNSDGT the spectrogram from coefficients computed with the 
%   function ndgt. For more details on the format of the variables c and a 
%   format, please read the ndgt function help.
%
%   PLOTNSDGT uses a dB colorscale, and the dynrange value can be used to
%   specify the dynamic of this colorscale, as the produced image uses a 
%   colormap in the interval [chigh-dynrange,chigh], where chigh is the 
%   highest value in the plot.
%
%   Limitation: PLOTNSDGT only works for coefficients c obtained from a
%   monochannel signal.
%
%   SEE ALSO:  NSDGT

%   AUTHOR : Florent Jaillet
%   TESTING: 
%   REFERENCE: 
%   Last changed 2009-05

% Todo :
% - Check the validity of the input
% - Handle the case of transform of multichannel signals

if nargin < 4
  % Default value for sampling frequency.
  sr=1; 
end

timepos=cumsum(a_new)-a_new(1);

if nargin < 3
  % Default value for colorscale dynamic.
  dynrange=60;
end 
  
% Compute time limit for the representation of each window.
tlim=diff(timepos)/2;
tlim=[timepos(1)-tlim(1);timepos(1:end-1)+tlim;timepos(end)+tlim(end)];
tlim=tlim/sr;

% Compute maximum of the representation for colorscale dynamic handling.
temp=cell2mat(c);
ma=20*log10(max(abs(temp(:))));

% Plot the representation: as the sampling grid in the time frequency plane
% is irregular, the representation by done by plotting many images next to 
% each other, with one image for each window
hold('on');
for ii=1:length(a_new)
  temp = 20*log10(abs(c{ii})+eps); % +eps is here to avoid log of 0
  % Octave cannot plot images that are only one point wide, so we use
  % images that are to points wide
  imagesc(adapt(tlim(ii:ii+1)),[0,1-1/length(c{ii})]*sr,[temp,temp],...
    [ma-dynrange,ma]);
end
hold('off');
axis('tight');

end

function [res]=adapt(lim) 
% we have to adapt the time values to fit the way the image function handle
% the x position
res=[(3*lim(1)+lim(2))/4,(lim(1)+3*lim(2))/4];
end
