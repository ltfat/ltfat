function test_libphaseret_overlaynthframe
f = greasy;
a = 128;
M = 1024;
M2 = floor(M/2) + 1; 
gl = 1024;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gamma = 0.25645*gl^2;
gd = long2fir(gabdual(g,a,M),gl);
gshift = fftshift(g);
N = 10;
idx = 9;
g2 = repmat(g.*gd,1,N);
g1 = repmat(fftshift(g.*gd),1,N);
cout = zeros(gl,1);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','phaseret_overlaynthframe_d',g1,gl,N,a,idx,coutPtr);

[~,out2] = comp_overlayframes(g2,a,gl,idx);
[~,out3orig] = overlayframes(g2,a,gl,idx);
out1orig = M*coutPtr.Value;
out3orig = out3orig*M;

out1 = out1orig;
out3 = out3orig;
out1(out1orig==0) = 1;
out1(out1<1e-6) = 1e-6;
out3(out3orig==0) = 1;
out3(out3<1e-6) = 1e-6;

figure(1); plot([gshift./out1,gshift./out3])



function [partrec,frame] = overlayframes(cframes,a,M,n)

N = size(cframes,2);
bufLen = N*a - (a-1) + M-1;
partrec = zeros(bufLen,1);

startidx = ceil(M/2)-1;
idxrange = startidx + [0:floor(M/2),-ceil(M/2)+1:-1];
for ii=0:N-1
    idx = ii*a + idxrange + 1;
    partrec(idx) = partrec(idx) + cframes(:,ii+1);
end

frame = partrec(a*n+1:a*n + M);
