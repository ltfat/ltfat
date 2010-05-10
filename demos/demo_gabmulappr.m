%DEMO_GABMULAPPR Approximate a slowly time variant system by a Gabor multiplier
%   
%   This script construct a slowly time variant system and performs the 
%   best approximation by a Gabor multiplier with specified parameters
%   (a and L see below). Then it shows the action of the slowly time variant
%   system as well as of the Gabor multiplier on three different test
%   signals (sinusoids.)
%
% FIGURE 1 
%
%  shows three different input signals (sinusoids) as well as the output signals
%  after applying the slowly time variant system (A) and the best 
%  approximation by a Gabor multiplier (G) on the signals.   
%
% FIGURE 2
%  shows the spectograms of the signal in Figure 1
%
% FIGURE 3
%  shows the Wigner distribution of the signal in Figure 1
%
% FIGURE 4
%  shows the original matrix and the approximation by Gabor multipliers
%
% FIGURE 5
%  shows the time-frequency spread of the matrices in Figure 4.
%
% See also: gabmulappr, gabmul

%   AUTHOR : Nina Engelputzeder, Peter Balazs.

disp('Type "help demo_gabmulappr" to see a description of how this demo works.');

% Setup parameters for the Gabor system and length of the signal
L=144; % Length of the signal
a=4;   % Time shift 
M=36;  % Number of modulations

% construction of slowly time variant system
% take an initial vector and multiply by random vector close to one
A = [];
c1=(1:L/2); c2=(L/2:-1:1); c=[c1 c2].^(-1); % weight of decay x^(-1)
A(1,:)=rand(1,L).*c;  % convolution kernel
for ii=2:L;
    A(ii,:)=circshift(A(ii-1,:).*(0.85+0.3.*rand(1,L)),[0 1]); 
end;

% perform best approximation by gabor multiplier
sym=gabmulappr(A,a,M);

% creation of 3 different input signals (sinusoids)
x=[2*pi/L:2*pi/L:2*pi]; 
s1=sin(10.*x);
s2=sin(20.*x);
s3=sin(30.*x);

% application of the slowly time variant system
As1=A*s1';  
As2=A*s2';
As3=A*s3';

% application of the Gabor multiplier
Gs1=gabmul(s1,sym,a); 
Gs2=gabmul(s2,sym,a);
Gs3=gabmul(s3,sym,a);

% the Gabor multiplier as matrix
GM = tfmat('gabmul',sym,a);

% Plotting the results

%% -------------- figure 1 ----------------------------------------
figure(1);
subplot(3,3,1); plot(s1); title('Input Signal  f_1');
subplot(3,3,2); plot(s2); title('Input Signal  f_2');
subplot(3,3,3); plot(s3); title('Input Signal  f_3');
subplot(3,3,4); plot(real(As1)); title ('Output Signal  Af_1');
subplot(3,3,5); plot(real(As2)); title ('Output Signal  Af_2');
subplot(3,3,6); plot(real(As3)); title ('Output Signal  Af_3');
subplot(3,3,7); plot(real(Gs1)); title ('Approx. by Gabor Multiplier  Gf_1');
subplot(3,3,8); plot(real(Gs2)); title ('Approx. by Gabor Multiplier  Gf_2');
subplot(3,3,9); plot(real(Gs3)); title ('Approx. by Gabor Multiplier  Gf_3');

%% ------------- figure 2 ------------------------------------------

figure(2);
subplot(3,3,1); sgram(s1,'tfr',10,'clim',[-40,13]); title('Spectogram: Input Signal  f_1');
subplot(3,3,2); sgram(s2,'tfr',10,'clim',[-40,13]); title('Spectogram: Input Signal  f_2');
subplot(3,3,3); sgram(s3,'tfr',10,'clim',[-40,13]); title('Spectogram: Input Signal  f_3');
subplot(3,3,4); sgram(real(As1),'tfr',10,'clim',[-40,13]); 
title ('Spectogram: Output Signal  Af_1');

subplot(3,3,5); sgram(real(As2),'tfr',10,'clim',[-40,13]);
title ('Spectogram: Output Signal  Af_2');

subplot(3,3,6); sgram(real(As3),'tfr',10,'clim',[-40,13]);
title ('Spectogram: Output Signal  Af_3');

subplot(3,3,7); sgram(real(Gs1),'tfr',10,'clim',[-40,13]);
title ('Spectogram: Approx. by Gabor Multiplier  Gf_1');

subplot(3,3,8); sgram(real(Gs2),'tfr',10,'clim',[-40,13]);
title ('Spectogram: Approx. by Gabor Multiplier  Gf_2');

subplot(3,3,9); sgram(real(Gs3),'tfr',10,'clim',[-40,13]);
title ('Spectogram: Approx. by Gabor Multiplier  Gf_3');

%% ----------- figure 3 -------------------------------

figure(3)
subplot(2,2,1);
surf(abs(A));
title ('Original Matrix ');
shading('interp');

subplot(2,2,2);
surf(abs(GM));
title ('Approximation by Gabor Multiplier');
shading('interp');

subplot(2,2,3);
surf(log(abs(A)));
title ('Original Matrix (log.scale)');
shading('interp');

subplot(2,2,4);
surf(log(abs(GM)));
title ('Approximation by Gabor Multiplier (log.scale)');
shading('interp');
     

%% ----------- figure 4 -------------------------------

figure(4)
TA = spreadfun(A); % the T-F spreading representation of the operator 
TGM = spreadfun(GM);
TA=fftshift(TA);
TGM=fftshift(TGM);

subplot(2,2,1);
surf(abs(TA));
title ('Time-Frequency spread of Original Matrix ');
shading('interp');

subplot(2,2,2);
surf(abs(TGM));
title ('Time-Frequency spread of Approximation by Gabor Multiplier ');
shading('interp');

subplot(2,2,3);
surf(log(abs(TA)));
title ('Time-Frequency spread of Original Matrix (log.scale)');
shading('interp');

subplot(2,2,4);
surf(log(abs(TGM)));
title (['Time-Frequency spread of Approximation by Gabor Multiplier ' ...
        '(log.scale)']);
shading('interp');
