function plottfjigsawsep(fplanes,cplanes,info,varargin)
%PLOTTFJIGSAWSEP Plots separated layers from tfjigsawsep
%   Usage: plottfjigsawsep(fplanes,cplanes,info);
%          plottfjigsawsep(fplanes,cplanes,info,fs); 
%       
%   `plottfjigsawsep(fplanes,cplanes,info)` shows the original and the three 
%   separated layers in both the time and the time-domains.
%   All parameters can be obtained from |tfjigsawsep|.
%
%   `plottfjigsawsep(fplanes,cplanes,info,fs)` works as above but scales the
%   axes labels according to the sampling rate *fs*. 
%   
%   The function calls |plotdgtreal| and forwards any arguments it understands.
%
%   Additional parameters
%   ---------------------
%
%      'equalyrange'   
%           The y-axis range is equal in all time domain plots.
%           By default, the individual y-axes are scaled automatically.
%
%      'showbuttons'
%           The figure will contain buttons for convenient playback of the
%           separated layers. Note that the sampling rate *fs* must be provided.
%           This is disabled by default and does not work on Octave.
%           
%   See also:   tfjigsawsep 
%
%   References: jato07

%AUTHOR: Daniel Haider, 2017 

complainif_notenoughargs(nargin,3,mfilename);

definput.import={'ltfattranslate','tfplot'};
definput.importdefaults = {'dynrange',90};
definput.flags.timeplot = {'noyscale','equalyrange'};
definput.flags.buttons = {'noshowbuttons','showbuttons'};
[flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);

if isempty(fs), fs = 1; end

clf; f = gcf;% set(f,'Visible','off');

xplot=(0:size(fplanes,1)-1)./fs;

orig = sum(fplanes,2);
ylims = [min(orig),max(orig)];
ton = fplanes(:,1);
trans = fplanes(:,2);
res = fplanes(:,3);

subplot(3,4,1)
plot(xplot,orig)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('original')
xlabel('Time (samples)')
subplot(3,4,2)
plot(xplot,ton)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('tonal')
xlabel('Time (samples)')
subplot(3,4,3)
plot(xplot,trans)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('transient')
xlabel('Time (samples)')
subplot(3,4,4)
plot(xplot,res)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('residual')
xlabel('Time (samples)')

if flags.do_equalyrange, set(gca,'ylim',ylims); end

subplot(2,4,5)
plotdgtreal(dgtreal(orig,info.g3,info.a3,info.M3),info.a3,info.M3,'argimport',flags,kv);
title('original')
subplot(2,4,6)
plotdgtreal(cplanes{1},info.a1,info.M1,'argimport',flags,kv);
title('tonal')
subplot(2,4,7)
plotdgtreal(cplanes{2},info.a2,info.M2,'argimport',flags,kv);
title('transient')
subplot(2,4,8)
plotdgtreal(cplanes{3},info.a3,info.M3,'argimport',flags,kv);
title('residual')

colormap(ltfat_inferno)

if flags.do_showbuttons && ~isoctave() && usejava('desktop')

if fs == 1
    warning(sprintf(['%s: Invalid sampling rate for audio playback.',...
        ' Pass ''fs'',fs to be used for audio playback.'],upper(mfilename)));
end

% The callbacks here used to be 'sound...' instead of '@(varargin)sound...'
% The difference is that the anonymous function version actually 
% stores fplanes.
try
btn1 = uicontrol('Style','pushbutton','String','Play Original Sound','Callback',@(varargin) sound(sum(fplanes,2),fs));
    
set(btn1,'Units','normalized');
set(btn1,'Position',[0.125 0.54 0.15 0.05]);
set(btn1,'BackgroundColor','m');

btn2 = uicontrol('Style','pushbutton','String','Play Tonal Part','Callback',@(varargin) sound(fplanes(:,1),fs));
set(btn2,'Units','normalized');
set(btn2,'OuterPosition',[0.335 0.54 0.15 0.05]);
set(btn2,'BackgroundColor','y');

btn3 = uicontrol('Style','pushbutton','String','Play Transient Part','Callback',@(varargin)sound(fplanes(:,2),fs));
set(btn3,'Units','normalized');
set(btn3,'OuterPosition',[0.545 0.54 0.15 0.05]);
set(btn3,'BackgroundColor','c');

btn4 = uicontrol('Style','pushbutton','String','Play Residual','Callback',@(varargin)sound(fplanes(:,3),fs));
set(btn4,'Units','normalized');
set(btn4,'OuterPosition',[0.752 0.54 0.15 0.05]);
set(btn4,'BackgroundColor','g');
catch
    warning('%s: Buttons are not supported',upper(mfilename));
end
end


%set(f,'Visible','on');
