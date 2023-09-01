
% This test example is taken from demo_gabmulappr

% Setup parameters for the Gabor system and length of the signal
L=576; % Length of the signal
a=32;   % Time shift 
M=72;  % Number of modulations
N=L/a;
fs=44100; % assumed sampling rate
SNRtv=63; % signal to noise ratio of change rate of time-variant system

% construction of slowly time variant system
% take an initial vector and multiply by random vector close to one
A = [];
c1=(1:L/2); c2=(L/2:-1:1); c=[c1 c2].^(-1); % weight of decay x^(-1)
A(1,:)=(tester_rand(1,L)-0.5).*c;  % convolution kernel
Nlvl = exp(-SNRtv/10);
Slvl = 1-Nlvl;
for ii=2:L;
     A(ii,:)=(Slvl*circshift(A(ii-1,:),[0 1]))+(Nlvl*(tester_rand(1,L)-0.5)); 
end;
A = A/norm(A)*0.99; % normalize matrix

% perform best approximation by gabor multiplier
g=gabtight(a,M,L);
sym1=gabmulappr(A,g,a,M);


% Now do the same using the general frame algorithm.

F=frame('dgt',g,a,M);

sym2=framemulappr(F,F,A);


norm(sym1-reshape(sym2,M,N))

% Test for exactness

testsym=tester_crand(M,N);
FT=frsynmatrix(F,L);

T=FT*diag(testsym(:))*FT';


sym1b=gabmulappr(T,g,a,M);
sym2b=framemulappr(F,F,T);

norm(testsym-sym1b)
norm(testsym-reshape(sym2b,M,N))
