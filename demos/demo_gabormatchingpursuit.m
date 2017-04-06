[f,fs] = cocktailparty; a = 1024; M = 8192; 
maxit = 100000;
Ls = numel(f);
L = dgtlength(Ls,a,M);
tfr = a*M/L;

g1 = pherm(L,0,tfr);
g2 = pherm(L,0,tfr/2);

g3 = pherm(L,0,tfr);
%  [c,relresmp1,fhat] = dgtrealmp(f,a,M,'localomp','tol',0,'maxit',maxit,...
%                   'g1',g1,...
%                   'relresstep',1000,'printdb','printstep',1000);
            
[c,relresmp2,fhat] = dgtrealmp(f,a,M,'mp','tol',1e-10,'maxit',maxit,...
                'g1',{'hann',M/2},...
                'printstep',10000,'relresstep',10000,'printdb');            
            
% [c,relresmp2,fhat] = dgtrealmp(f,1024,M,'mp','tol',1e-5,'maxit',maxit,...
%                 'g1',{'hann',4096},...
%                 'relresstep',1000);
%             
% [c,relresmp3,fhat] = dgtrealmp(f,32,M,'mp','tol',1e-5,'maxit',maxit,...
%                 'g1',{'hann',128},...
%                 'relresstep',1000);
            
%plot(20*log10([relresmp1,relresmp2]))
%[c,relreslocmp] = dgtrealmp(f,a,M,'locomp','tol',1e-4,'maxit',maxit,...
%                    'relresstep',1000,'g1','hann');
% 

