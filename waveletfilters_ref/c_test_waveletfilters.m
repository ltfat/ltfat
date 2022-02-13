
%test waveletfilters
addpath('/run/media/clara/3143d7fe-0bef-4c9d-8983-732cfe02d2c9/ltfat/waveletfilters_ref');
pedantic = 0;%switch pedantic on to compare fb coefficients, but this may take a while
 [f, fs] = gspi;
 Ls = length(f);
 scales = linspace(10,0.1,100);
 alpha = 1-2/(1+sqrt(5));
 delays = @(n,a) a*(mod(n*alpha+.5,1)-.5);

[g,a,fc,L, info]=waveletfilters(Ls,scales, {'fbsp', 4, 3}, 'repeat', 'delay', delays);
%[g,a,fc,L,info] = waveletfilters(fs, fmin, fmax, M, Ls, 'redtar', 0.5);
[gtemp,atemp,fctemp,Ltemp, infotemp]=waveletfilterstemp(Ls,scales, {'fbsp', 4, 3}, 'repeat', 'delay', delays);
%filterbankfreqz(g,a,L,fs,'plot','linabs','posfreq');
if pedantic
    %check coefficients
    c = filterbank(f, g, a);
    ctemp = filterbank(f, gtemp, atemp);
    assert(abs(mean(c{2,1}-ctemp{2,1})) == 0)
    assert(sum(c{1,1}-ctemp{1,1}) == 0)
end
%abs(sum(sum(c{1,1}/ctemp{1,1})))/numel(c{1,1})
%check basic params
assert(isequal(a(:,1), atemp(:,1)))
assert(isequal(fc, fctemp))
assert(isequal(L, Ltemp))
%check info struct
assert(isequal(info.fc, infotemp.fc))
assert(isequal(info.basefc, infotemp.basefc))
assert(isequal(info.foff, infotemp.foff))
assert(isequal(info.fsupp, infotemp.fsupp))
assert(isequal(info.scale, infotemp.scale))
assert(isequal(info.dilation, infotemp.dilation))
assert(isequal(info.bw, infotemp.bw))
assert(isequal(info.tfr, infotemp.tfr))
assert(isequal(info.aprecise, infotemp.aprecise))
assert(isequal(info.a_natural, infotemp.a_natural))
assert(isequal(info.cauchyAlpha, infotemp.cauchyAlpha))
assert(isequal(info.lowpassstart, infotemp.startindex))
assert(isequal(g{1,95}.delay, gtemp{1,95}.delay))
assert(isequal(g{1,1}.delay, gtemp{1,1}.delay))
