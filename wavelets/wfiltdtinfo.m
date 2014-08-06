function wfiltdtinfo(dw,varargin)
%WFILTDTINFO Plots dual-tree filters info
%   Usage: wfiltdtinfo(dw);
%
%   Input parameters:
%         dw     : Wavelet dual-tree filterbank
%
%   `wfiltdtinfo(w)` plots impulse responses, frequency responses and 
%   approximation of the scaling and of the wavelet function(s) associated
%   with the dual-tree wavelet filters defined by *w* in a single figure. 
%   Format of *dw* is the same as in |dtwfb|.
%
%   Optionally it is possible to define scaling of the y axis of the
%   frequency seponses. Supported are:
%
%   'db','lin'   
%       dB or linear scale respectivelly. By deault a dB scale is used.
%
%   Examples:
%   ---------
%   :::
%      wfiltdtinfo('qshift4');
%   

complainif_notenoughargs(nargin,1,'WFILTDTINFO');


definput.flags.freqzscale = {'db','lin'};
[flags]=ltfatarghelper({},definput,varargin);


[dwstruct,info] = dtwfbinit({'strict',{dw,6,'dwt'}});
dw = info.dw;

filtNo = size(dw.g,1);
grayLevel = [0.6,0.6,0.6];
clf;

colorAr = repmat('rbkcmg',1,filtNo);

for ii=1:2
    subplot(5,filtNo,filtNo*(ii-1)+1);
    title(sprintf('Scaling imp. response, tree %i',ii));
    loAna = dw.g{1,ii}.h;
    loShift = -dw.g{1,ii}.offset;
    xvals = -loShift + (0:length(loAna)-1);
    hold on;
    if ~isempty(loAna(loAna==0))
       stem(xvals(loAna==0),loAna(loAna==0),'Color',grayLevel);
    end

    loAnaNZ = find(loAna~=0);
    stem(xvals(loAnaNZ),loAna(loAnaNZ),colorAr(1));
    axis tight;
    hold off;
end

for ii=1:2
    for ff=2:filtNo
        subplot(5,filtNo,ff+filtNo*(ii-1));
        title(sprintf('Wavelet imp. response no: %i, tree %i',ff-1,ii));
        filtAna = dw.g{ff,ii}.h;
        filtShift = -dw.g{ff,ii}.offset;
        xvals = -filtShift + (0:length(filtAna)-1);
        filtNZ = find(filtAna~=0);
        hold on;

        if ~isempty(filtAna(filtAna==0))
           stem(xvals(filtAna==0),filtAna(filtAna==0),'Color',grayLevel);
        end

        stem(xvals(filtNZ),filtAna(filtNZ),colorAr(ff));
        axis tight;
        hold off;
    end
end

L = 1024;
Lc = wfbtclength(L,dwstruct,'per');
c = wavpack2cell(zeros(sum([Lc;Lc(end:-1:1)]),1),...
                 [Lc;Lc(end:-1:1)]);
c{1}(1) = 1;
sfn = flipud(idtwfb(c,dwstruct,L));


subplot(5,filtNo,[2*filtNo+1]);
xvals = ((-floor(L/2)+1:floor(L/2)).');
plot(xvals,fftshift([abs(sfn),real(sfn),imag(sfn)]));
axis tight;
title('Scaling function');
legend({'abs','real','imag'},'Location','west')

for ff=2:filtNo
   subplot(5,filtNo,[2*filtNo+ff]);
   
   c{ff-1}(1) = 0;
   c{ff}(1) = 1;
   wfn = flipud(idtwfb(c,dwstruct,L));
   
   plot(xvals,fftshift([abs(wfn),real(wfn),imag(wfn)]));
   axis tight;
   legend({'abs','real','imag'},'Location','west')
   title(sprintf('Wavelet function: %i',ff-1));
end

subplot(5,filtNo,3*filtNo + (1:filtNo) );
title('Magnitude frequency response');
maxLen=max(cellfun(@(gEl) numel(gEl.h),dw.g));
Ls = nextfastfft(max([maxLen,1024]));
H = filterbankfreqz(dw.g(:),[dw.a(:);dw.a(:)],Ls);



%[H] = wtfftfreqz(w.g);
if flags.do_db
   plotH = 20*log10(abs(H));
elseif flags.do_lin
   plotH = abs(H);   
else
    error('%s: Unknown parameter',upper(mfilaname));
end
xVals = linspace(0,1,numel(H(:,1)));
hold on;
for ff=1:filtNo*2
   plot(xVals,plotH(:,ff),colorAr(ff));
   axis tight;
end
if flags.do_db
    ylim([-30,max(plotH(:))])
end
ylabel('|\itH|[dB]');
xlabel('\omega [-]')
hold off;

subplot(4,filtNo,3*filtNo + (1:filtNo) );
title('Phase frequency response');
hold on;
for ff=1:filtNo*2
   plot(xVals,unwrap(angle((H(:,ff))))/pi,colorAr(ff));
   axis tight;
end
ylabel('arg H(\omega)[\pi rad]');
xlabel('\omega [-]')

axis tight;
% plot(unwrap(angle([H])));
% axis tight;
hold off;
