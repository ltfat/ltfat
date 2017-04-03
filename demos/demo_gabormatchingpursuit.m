[f,fs] = gspi; a = 512; M = 2048; 
maxit = 40000;

%[c,relresmp] = dgtrealmp(f,a,M,'mp','tol',1e-14,'maxit',maxit,'relresstep',1000,'g1','hann');
[c,relreslocmp] = dgtrealmp(f,a,M,'locomp','tol',1e-14,'maxit',maxit,'relresstep',1000,'g1','hann');
% 
plot(20*log10([relresmp,relreslocmp]))
