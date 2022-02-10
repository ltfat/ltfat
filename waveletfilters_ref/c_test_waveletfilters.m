%test waveletfilters
addpath('/run/media/clara/3143d7fe-0bef-4c9d-8983-732cfe02d2c9/ltfat/waveletfilters_ref');
pedantic = 1;%switch pedantic on to compare fb coefficients, but this may take a while
[f, fs] = gspi;
Ls = length(f);
scales = linspace(10,0.1,100);

[g,a,fc,L, info]=waveletfilters(Ls,scales, {'fbsp', 4, 3}, 'redtar', 0.5);
[gtemp,atemp,fctemp,Ltemp, infotemp]=waveletfilterstemp(Ls,scales, {'fbsp', 4, 3}, 'redtar', 0.5);

if pedantic
    %check coefficients
    c = filterbank(f, g, a);
    ctemp = filterbank(f, gtemp, atemp);
    assert(abs(mean(c{1,1}-ctemp{1,1})) < 0.001)
end

%check basic params
assert(isequal(a, atemp))
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
