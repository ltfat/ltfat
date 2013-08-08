function wfiltinfo(w)
%WFILTINFO Plots filters info
%   Usage: wfiltinfo(w);
%
%   Input parameters:
%         w     : Wavelet filterbank
%
%   Plots impulse responses, frequency responses and approximation of
%   the scaling of the wavelet function(s) associated withthe wavelet filters
%   defined by *w* in a single figure. Format of *w* is the same as in |fwt|.
%
%   Examples:
%   ---------
%   
%   Details of the 'syn:spline8:8' wavelet filters (see |wfilt_spline|):::
%   
%      wfiltinfo('syn:spline8:8');
%   
%   Details of the 'ana:spline8:8' wavelet filters:::
%
%      wfiltinfo('ana:spline8:8');
%
%   See also: wfilt_db 


w = fwtinit({'strict',w});
filtNo = length(w.g);

grayLevel = [0.6,0.6,0.6];


colorAr = repmat('rbkcmg',1,filtNo);

subplot(4,filtNo,1);
title('Scaling imp. response');
loAna = w.g{1}.h;
loShift = -w.g{1}.offset;
xvals = -loShift + (1:length(loAna));
hold on;
if ~isempty(loAna(loAna==0))
   stem(xvals(loAna==0),loAna(loAna==0),'Color',grayLevel);
end

loAnaNZ = find(loAna~=0);
stem(xvals(loAnaNZ),loAna(loAnaNZ),colorAr(1));
axis tight;
hold off;

for ff=2:filtNo
    subplot(4,filtNo,ff);
    title(sprintf('Wavelet imp. response no: %i',ff-1));
    filtAna = w.g{ff}.h;
    filtShift = -w.g{ff}.offset;
    xvals = -filtShift + (1:length(filtAna));
    filtNZ = find(filtAna~=0);
    hold on;
    
    if ~isempty(filtAna(filtAna==0))
       stem(xvals(filtAna==0),filtAna(filtAna==0),'Color',grayLevel);
    end
    
    stem(xvals(filtNZ),filtAna(filtNZ),colorAr(ff));
    axis tight;
    hold off;
end

[wfn,sfn,xvals] = wavfun(w,'fft');
subplot(4,filtNo,[filtNo+1]);

plot(xvals(:,end),sfn,colorAr(1));
axis tight;
title('Scaling function');

for ff=2:filtNo
   subplot(4,filtNo,[filtNo+ff]);
   plot(xvals(:,ff-1),wfn(:,ff-1),colorAr(ff));
   axis tight;
   title(sprintf('Wavelet function: %i',ff-1));
end

subplot(4,filtNo,2*filtNo + (1:filtNo) );
title('Magnitude frequency response');
maxLen=max(cellfun(@(gEl) numel(gEl.h),w.g));
Ls = nextfastfft(max([maxLen,1024]));
H = zeros(Ls,filtNo);
for ii=1:filtNo
   H(:,ii) = comp_transferfunction(w.g{ii},Ls);
end

%[H] = wtfftfreqz(w.g);
plotH = 20*log10(abs(H));
xVals = linspace(0,1,numel(H(:,1)));
hold on;
for ff=1:filtNo
   plot(xVals,plotH(:,ff),colorAr(ff));
   axis tight;
end
ylim([-30,max(plotH(:))])
ylabel('|\itH|[dB]');
xlabel('\omega [-]')

subplot(4,filtNo,3*filtNo + (1:filtNo) );
title('Phase frequency response');
hold on;
for ff=1:filtNo
   plot(xVals,unwrap(angle((H(:,ff))))/pi,colorAr(ff));
   axis tight;
end
ylabel('arg H(\omega)[\pi rad]');
xlabel('\omega [-]')

 axis tight;
% plot(unwrap(angle([H])));
% axis tight;

