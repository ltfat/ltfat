function test_failed = test_gabphasederiv(dflag,varargin)
test_failed = 0;

L = 1000;
l = 0:L-1; l = l(:);

f1 = exp(1i*pi*0.1*l);
f2 = pchirp(L,1);
f3 = expchirp(L,0.4,2.2,'fc',-1.2);
f4 = exp(1i*2*pi*0.1*L*(l./(L)).^2);
tfr = 1;
g = {'gauss',tfr};
%g = {'hann',64};
a = 1;
M = L/1;
minlvl = 1e-10;

f= f4;
%f = dft(f);



c = dgt(f,g,a,M);

dynrange = 90;
dynrangerat = 10^(dynrange/20)
mincoef = max(abs(c(:)))/dynrangerat;

phase = unwrap(angle(c),[],2);
phaselocked = (angle(phaselock(c,a)));
phase(abs(c)<mincoef) = 0;
phaselocked(abs(c)<mincoef) = 0;
%  figure(10);
%  plotdgt(phase/pi,a,'lin','clim',[-pi,pi]);
%   figure(11);
%   plotdgt(phaselocked/pi,a,'lin','clim',[-pi,pi]);

clim = [-1,1];

figure(1);
phasedt2 = gabphasederiv(dflag,'dgt',f,g,a,M,'minlvl',minlvl,varargin{:});
phasedt2(abs(c)<mincoef) = 0;
plotdgt(phasedt2,a,'lin','clim',clim,'tc');

figure(2);
phasedf2 = gabphasederiv(dflag,'abs',abs(c),g,a,varargin{:});
phasedf2(abs(c)<mincoef) = 0;
plotdgt(phasedf2,a,'lin','clim',clim);

figure(3);
phasedff2 = gabphasederiv(dflag,'phase',angle(c),a,varargin{:});
phasedff2(abs(c)<mincoef) = 0;
plotdgt(phasedff2,a,'lin','clim',clim);

figure(4);
phasedt2 = gabphasederiv(dflag,'cross',f,g,a,M,'minlvl',minlvl,varargin{:});
phasedt2(abs(c)<mincoef) = 0;
plotdgt(phasedt2,a,'lin','clim',clim);

%figure(4);
%plotdgt(abs(phasedt2-phasedff2),a,'lin','clim',[0,50]);

figure(5);
plotdgt(c,a,'dynrange',dynrange);
