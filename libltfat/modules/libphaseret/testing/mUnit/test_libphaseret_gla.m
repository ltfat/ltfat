
f = greasy;
a = 128;
M = 1024;
M2 = floor(M/2) + 1; 
gl = M;
L = dgtlength(numel(f),a,M);
g = firwin('blackman',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;
maxit = 100;

corig = dgtreal(f,{'blackman',gl},a,M);
s = abs(corig) + 1i*zeros(size(corig));

cinPtr = libpointer('doublePtr',complex2interleaved(s));
cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','phaseret_gla_d',cinPtr,libpointer(),g,L,gl,1,a,M,maxit,coutPtr);

cout2 = interleaved2complex(coutPtr.Value);

frec = idgtreal(cout2,{'dual',{'blackman',gl}},a,M);

s2 = dgtreal(frec,{'blackman',gl},a,M);
magnitudeerrdb(s,s2)



