clc;
clear all;

fontname = 'Times New Roman';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
set(0,'defaulttextfontsize',12);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','none');

J = 3;
wavelet = 'db4';
[lo,hi,lo_r,hi_r] = wfilters(wavelet);
multi = makeMultirateIdentity( lo_r,hi_r, J);

figure(1);
clf;
freLen = 512;
fre = zeros(freLen,length(multi));

for i = 1:length(multi)
    [H,f] = freqz(multi{i},1,freLen);
    fre(:,i) = abs(H);
end



plot(f/pi,fre);
l = xlabel('xlabel');
set(l,'FontSize',12);
l=ylabel('ylabel');
set(l,'FontSize',12);
l = legend('g','h2g','h2h4g','h2h4hmmmmm');
set(l,'FontSize',16);

%savetoeps(gcf,'ch01multirateFreqz.eps',15);

figure(2);
clf;
hold on;
cc=lines(12);
handles = zeros(1,length(multi));
for i = 1:length(multi)
   handles(i) = stem(0:length(multi{i})-1,multi{i},'color',cc(i,:),'LineWidth',0.5);
   %set(handles(i),'LineWidth',1);
end

xlim([0 length(multi{end}) ]);
ylim([min(multi{1}) max(multi{1}) + 0.5]);
l = xlabel('xlabel');
set(l,'FontSize',12);
l=ylabel('ylabel');
set(l,'FontSize',12);
l = legend(handles',{'g','h2g','h2h4g','h2h4hmmmmmmm'});
set(l,'FontSize',16);
hold off;

%savetoeps(gcf,'ch01multirateImpz.eps',[15 10]);