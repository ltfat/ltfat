for ii= 70:70

filename = sprintf('~/Desktop/SQAM/%02d.wav',ii);
disp(filename);

[f,fs] =wavload(filename); 
f = f(1:10*fs,1);

%f = f(1:fs);
[f,fs] =gspi;

a1 = 32; M1 = 128;
a2 = 512; M2 = 2048;
a3 = 512; M3 = 2048;
a4 = 1024; M4 = 4096;


Ls = numel(f);
L = dgtlength(Ls,a2,M2);
atoms = 0.5*L;
f = postpad(f,L);
f(:) = 1;

step = 10;
g = 'gauss';


[c1,fhat1,info1] = dgtrealmp(f,'mp','tol',1e-4,'atoms',atoms,...
                   'g2',{'blackman',512,2048},...
                   'atsellim',10000,...
                   'relresstep',step,'errappr','printdb',...
                   'relrestoldb',-60);
%                     'g3',{'gauss',tfr/2^2},...
%                    'g4',{'gauss',tfr/2^2},...
%                    'g5',{'gauss',tfr/2^2},...
%                    'g6',{'gauss',tfr/2^2},...              

 info1
 info1.exitmsg
 20*log10(info1.relres)
end
% a = 256;
% [c2,fhat2,info2] = dgtrealmp(f,a,M,'localomp','tol',1e-4,'atoms',atoms,...
%                     'g1',g,...
%                     'printstep',step,'relresstep',step,'printdb');   
%     
% a = 128;
% [c3,fhat3,info3] = dgtrealmp(f,a,M,'localomp','tol',1e-4,'atoms',atoms,...
%                      'g1',g,...
%                      'printstep',step,'relresstep',step,'printdb'); 
                
% [c4,fhat4,info4] = dgtrealmp(f,a,M,'compllocalomp','tol',1e-4,'atoms',atoms,...
%                     'g1',g,...
%                     'printstep',step,'relresstep',step,'printdb'); 
            
% [c,relresmp2,fhat] = dgtrealmp(f,1024,M,'mp','tol',1e-5,'maxit',maxit,...
%                 'g1',{'hann',4096},...
%                 'relresstep',1000);
%             
% [c,relresmp3,fhat] = dgtrealmp(f,32,M,'mp','tol',1e-5,'maxit',maxit,...
%                 'g1',{'hann',128},...
%                 'relresstep',1000);
            
%plot(20*log10([info1.relres,info2.relres,info3.relres]))
%[c,relreslocmp] = dgtrealmp(f,a,M,'locomp','tol',1e-4,'maxit',maxit,...
%                    'relresstep',1000,'g1','hann');
% 

