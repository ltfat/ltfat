%DEMO_DGT  Basic introduction to DGT analysis/synthesis
%
%   This demo shows how to compute Gabor coefficients of a signal.
%
%   FIGURE 1 Spectrogram of the 'bat' signal.
%
%     The figure shows a spectrogram of the 'bat' signal. The
%     coefficients is shown on a linear scale.
%
%   FIGURE 2 Gabor coefficient of the 'bat' signal.
%
%     The figure show a set of Gabor coefficients for the 'bat' signal,
%     computed using a DGT with a Gaussian window. The coefficients
%     contains all the information to reconstruct the signal, even though
%     there a far fewer coefficients than the spectrogram contains.
%
%   FIGURE 3 Real-valued Gabor analysis
%
%     This figure shows only the coefficients for the positive
%     frequencies. As the signal is real-value, these coefficients
%     contain all the necessary information. Compare to the shape of the
%     spectrogram shown on Figure 1.
%
%   FIGURE 4 DGT coefficients on a spectrogram
%
%     This figure shows how the coefficients from DGTREAL can be picked
%     from the coefficients computed by a full Short-time Fourier
%     transform, as visualized by a spectrogram.
%
%   See also: sgram, dgt, dgtreal

disp('Type "help demo_dgt" to see a description of how this demo works.');

% Load a test signal
f=bat;

% sampling rate of the test signal, only important for plotting
fs=143000;

% Length of signal
Ls=length(f);

disp(' ');
disp('------ Spectrogram analysis -----------------------------------');

disp('Figure 1 displays a spectrogram of the bat test signal.');
figure(1);
c_sgram=sgram(f,fs,'lin');


% Number of coefficients in the Spectrogram
no_sgram=numel(c_sgram);

disp(' ');
disp('The spectrogram is highly redundant.');
fprintf('No. of coefficients in the signal:       %i\n',Ls);
fprintf('No. of coefficients in the spectrogram:  %i\n',no_sgram);
fprintf('Redundacy of the spectrogram:            %f\n',no_sgram/Ls);

% WARNING: In the above code, the spectrogram routine SGRAM returns the
% coefficients use to plot the image. These coefficients are ONLY
% intended to be used by post-processing image tools, and in this
% example, the are only used to illustrate the redundancy of the
% spectogram. Numerical Gabor signal analysis and synthesis should ALWAYS
% be done using the DGT, IDGT, DGTREAL and IDGTREAL functions, see the
% following sections of this example.

disp(' ');
disp('---- Simple Gabor analysis using a standard Gaussian window. ----');

disp('Setup parameters for a Discrete Gabor Transform.')
disp('Time shift:')
a=20

disp('Number of frequency channels.');
M=40

disp(' ');
disp('Note that it must hold that L = M*b = N*a for some integers b, N and L,');
disp('and that a<M. L is the transform length, and the DGT will choose the');
disp('smallest possible value of L that is larger or equal to the length of the');
disp('signal. Choosing a<M makes the transform redundant, otherwise the');
disp('transform will be lossy, and reconstruction will not be possible.');

% Simple DGT using a standard Gaussian window.
c=dgt(f,'gauss',a,M);

disp('Number of time shifts in transform:')
N = size(c,2);

disp('Length of transform:')
L = N*a


disp('Figure 2 visualize the Gabor coefficients.');
figure(2);
imagesc(abs(c).^2);
axis('xy');
colorbar;

disp(' ');
disp(['The redundancy of the Gabor transform can be reduced without loosing ' ...
      'information.']);
fprintf('No. of coefficients in the signal:       %i\n',Ls);
fprintf('No. of output coefficients from the DGT: %i\n',numel(c));
fprintf('Redundacy of the DGT (in this case)      %f\n',numel(c)/Ls);

disp(' ');
disp('---- Real valued Gabor analysis. ----');

% Figure 1 and Figure 2 looks quite different, because Figure 2 also
% displays the coefficients for the n

% Simple real valued DGT using a standard Gaussian window.
c_real=dgtreal(f,'gauss',a,M);

disp('Figure 3 shows the positive-frequency DGT coefficients (DGTREAL).');
figure(3);
imagesc(abs(c_real).^2);
axis('xy');
colorbar;

disp('Figure 4 shows placement of the DGTREAL coefficients on the spectrogram.');
figure(4)
b=L/M;
[X,Y]=meshgrid(1:a:L+a,1:b:L/2+b);

hold on;
imagesc(c_sgram);
plot([X(:),X(:)]',[Y(:),Y(:)]','wo','Linewidth',1);
axis('xy','image');
hold off;

disp(' ');
disp('---- Perfect reconstruction. ----');

% Reconstruction from the full DGT coefficients
r      = idgt(c,'gaussdual',a);

% Reconstruction from the DGTREAL coefficients
% The parameter M cannot be deduced from the size of the coefficient
% array c_real, so it is an explicit input parameter.
r_real = idgtreal(c_real,'gaussdual',a,M);

fprintf('Reconstruction error using IDGT:      %e\n',norm(f-r));
fprintf('Reconstruction error using IDGTREAL:  %e\n',norm(f-r_real));


