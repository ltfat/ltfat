%DEMO_DGT  How to call a DGT
%
%   This script sets up a Gabor frame with the specified
%   parameters (a, M, L defined below), computes a window
%   and the corresponding canonical dual  window and a test signal,
%   and plots the windows and the energy of the coefficients.
%
%   FIGURE 1 window+dual
%
%    This figure shows a Gaussian used as the window function
%    and its canonical dual. 
%
%   FIGURE 2 absolute value of coefficients
%
%    This figure shows a (color coded) image of the dgt coefficient
%    modulus. 
%
%   See also:  dgt, idgt, pgauss

disp('Type "help demo_dgt" to see a description of how this demo works.');

% Setup parameters and length of signal.
% Not that it must hold that L=M*b=N*a for some integers
% M and N, and that a*b <= L

L=480;  % Length of signal.
a=20;   % Time shift.
M=24;   % Number of modulations.

% From the parameters L, a and M we can calculate:
b=L/M;  % Length of frequency shift.
N=L/a;  % Number of time shifts.

% Get a good window and its canonical dual.
g=pgauss(L,a/b);
gamma=gabdual(g,a,M);

% Plot them:

% Standard note on plotting:
%
% - The windows are all centered around zero, but this
%   is not visually pleasing, so the window must be
%   shifted to the middle by an FFTSHIFT

figure(1);

subplot(2,1,1);
plot(fftshift(g));
title('Periodized gaussian.');
legend('off');

subplot(2,1,2);
plot(fftshift(gamma));
title('Canonical dual window.');
legend('off');

% Setup a complex test signal.
ftest=ctestfun(L);

% Calculate coefficients.
coef=dgt(ftest,gamma,a,M);

% Plot
figure(2);
imagesc(log(abs(coef)));
title('DGT magnitude logarithm');
if ~isoctave
    set(gca,'ydir','normal');
    xlabel('Time');
    ylabel('Frequency');
end

% Test reconstruction
ftest_r=idgt(coef,g,a);

% Print relative error of reconstruction.
disp('');
%disp('Relative error of reconstruction, should be close to zero.');
rec_err = norm(ftest-ftest_r)/norm(ftest);

fprintf('Relative error of reconstruction (should be close to zero.):    %e \n',rec_err);



