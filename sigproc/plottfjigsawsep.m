function plottfjigsawsep(fplanes,cplanes,info,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

complainif_notenoughargs(nargin,3,mfilename);


definput.import={'ltfattranslate','tfplot'};
definput.importdefaults = {'dynrange',90};
definput.flags.timeplot = {'noyscale','equalyrange'};
[flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);

if isempty(fs)
    fs = 1;
end

f = figure('Visible','off');
    
xplot=(0:size(fplanes,1)-1)./fs;

orig = sum(fplanes,2);
ylims = [min(orig),max(orig)];
ton = fplanes(:,1);
trans = fplanes(:,2);
res = fplanes(:,3);
    
subplot(3,4,1)
p1 = plot(xplot,orig)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('original')
xlabel('Time (s)')
subplot(3,4,2)
p2 = plot(xplot,ton)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('tonal')
xlabel('Time (s)')
subplot(3,4,3)
p3 = plot(xplot,trans)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('transient')
xlabel('Time (s)')
subplot(3,4,4)
p4 = plot(xplot,res)
if flags.do_equalyrange, set(gca,'ylim',ylims); end
title('residual')
xlabel('Time (s)')

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

% if ~isoctave
% btn1 = uicontrol('Style','pushbutton','String','Play Original Sound','Callback','sound(orig,fs)');
% set(btn1,'Units','normalized');
% set(btn1,'Position',[0.125 0.54 0.15 0.05]);
% set(btn1,'BackgroundColor','m');
% 
% btn2 = uicontrol('Style','pushbutton','String','Play Tonal Part','Callback','sound(ton,fs)');
% set(btn2,'Units','normalized');
% set(btn2,'OuterPosition',[0.335 0.54 0.15 0.05]);
% set(btn2,'BackgroundColor','y');
% 
% btn3 = uicontrol('Style','pushbutton','String','Play Transient Part','Callback','sound(trans,fs)');
% set(btn3,'Units','normalized');
% set(btn3,'OuterPosition',[0.545 0.54 0.15 0.05]);
% set(btn3,'BackgroundColor','c');
% 
% btn4 = uicontrol('Style','pushbutton','String','Play Residual','Callback','sound(res,fs)');
% set(btn4,'Units','normalized');
% set(btn4,'OuterPosition',[0.752 0.54 0.15 0.05]);
% set(btn4,'BackgroundColor','g');
% end


set(f,'Visible','on');
