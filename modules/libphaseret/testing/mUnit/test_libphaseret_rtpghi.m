clear all;
f = cocktailparty;
a = 128;
M = 2048;
M2 = floor(M/2) + 1; 
gl = 2048;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;


corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);
cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);
gamma = gl^2*0.25645;

tic
calllib('libphaseret','phaseret_rtpghioffline_d',s,L,1,a,M,gamma,1e-6,1,coutPtr);
t = toc;
t/N*1000

coutsingle = zeros(2*M2,N,'single');
coutsinglePtr = libpointer('singlePtr',coutsingle);

ssingle = cast(s,'single');

tic
calllib('libphaseret','phaseret_rtpghioffline_s',ssingle,L,1,a,M,gamma,1e-6,1,coutsinglePtr);
t =toc;
t/N*1000

coutsingle2 = interleaved2complex(coutsinglePtr.Value);
cout2 = interleaved2complex(coutPtr.Value);


frec = idgtreal(coutsingle2,{'dual',{'hann',gl}},a,M,'timeinv');
s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)


cout2 = rtpghi(s,gamma,a,M,'timeinv','tol',1e-6,'causal');

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');
s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)




