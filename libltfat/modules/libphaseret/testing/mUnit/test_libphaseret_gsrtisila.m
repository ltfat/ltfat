;
f = gspi;
a = 256;
M = 2048;
M2 = floor(M/2) + 1; 
gl = 2048;
L = dgtlength(numel(f),a,M);
win = 'hann';
g = firwin(win,gl);
gd = long2fir(gabdual(g,a,M),gl);
N = L/a;
lookahead = 0;
maxit = 16;

corig = dgtreal(f,{win,gl},a,M,'timeinv');
s = abs(corig);

cout = zeros(2*M2,N);
coutPtr = libpointer('doublePtr',cout);

calllib('libphaseret','phaseret_gsrtisilaoffline_d',s,g,L,gl,1,a,M,lookahead,maxit,coutPtr);

coutsingle = zeros(2*M2,N,'single');
coutsinglePtr = libpointer('singlePtr',coutsingle);
ssingle = cast(s,'single');
gsingle = cast(g,'single');
calllib('libphaseret','phaseret_gsrtisilaoffline_s',ssingle,gsingle,L,gl,1,a,M,lookahead,maxit,coutsinglePtr);


coutsingle2 = interleaved2complex(coutsinglePtr.Value);
cout2 = interleaved2complex(coutPtr.Value);


frec = idgtreal(coutsingle2,{'dual',{win,gl}},a,M,'timeinv');

s2 = dgtreal(frec,{win,gl},a,M,'timeinv');
magnitudeerrdb(s,s2)



c=gsrtisila(s,g,a,M,'lookahead',lookahead,'maxit',maxit,'timeinv');
frec = idgtreal(c,{'dual',{win,gl}},a,M,'timeinv');
magnitudeerrdb(s,dgtreal(frec,{win,gl},a,M,'timeinv'))


