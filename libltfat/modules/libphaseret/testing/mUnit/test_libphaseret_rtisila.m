clear all;
f = greasy;
a = 128;
M = 1024;
M2 = floor(M/2) + 1; 
gl = 1024;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;
lookahead = 0;
maxit = 10;

corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);

cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','phaseret_rtisilaoffline_d',s,g,L,gl,1,a,M,lookahead,maxit,coutPtr);

cout2 = interleaved2complex(coutPtr.Value);

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');

s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)



c=rtisila(s,g,a,M,'lookahead',lookahead,'maxit',maxit,'timeinv');
frec = idgtreal(c,{'dual',{'hann',gl}},a,M,'timeinv');
magnitudeerrdb(s,dgtreal(frec,{'hann',gl},a,M,'timeinv'))

plotdgtreal(abs(c-cout2),a,M,'linabs')

