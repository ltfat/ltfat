clear all;
f = gspi;
a = 256;
M = 1024;
M2 = floor(M/2) + 1; 
gl = 1024;
L = dgtlength(numel(f),a,M);
g = firwin('blackman',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;


corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);
cout = zeros(2*M2,N);
cout(1:numel(s)) = s(:);
coutPtr = libpointer('doublePtr',cout);
calllib('libphaseret','phaseret_spsi_d',coutPtr,L,1,a,M,libpointer(),coutPtr);

%cout2 = interleaved2complex(coutPtr.Value);
cout2 = spsi(s,a,M,'timeinv');

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');

s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)




