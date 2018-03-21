
f = greasy;
a = 128;
M = 1024;
M2 = floor(M/2) + 1; 
gl = 1024;
L = dgtlength(numel(f),a,M);
g = firwin('hann',gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;


corig = dgtreal(f,{'hann',gl},a,M,'timeinv');
s = abs(corig);
cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);
gamma = gl^2*0.25645;
calllib('libphaseret','phaseret_pghi_d',s,L,1,a,M,gamma,coutPtr);

cout2 = interleaved2complex(coutPtr.Value);
%cout2 = pghi(s,gamma,a,M,'timeinv');

frec = idgtreal(cout2,{'dual',{'hann',gl}},a,M,'timeinv');

s2 = dgtreal(frec,{'hann',gl},a,M,'timeinv');
magnitudeerrdb(s,s2)




