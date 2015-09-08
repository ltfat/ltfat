% Testing script for audfilters

[f,fs]=greasy;  % Get the test signal
Ls=length(f);
fmin = 0;
fmax = fs/2;

%% Generate a fractional ERBlet FB with V= 1 and a Hann window
disp('--- Fractional non-uniform ERBlet FB, V=1, Hann window ---')
[g_erb,a_erb,fc_erb]=audfilters(fs,Ls,fmin,fmax,'erb','fractional');
c_erb=filterbank(f,{'realdual',g_erb},a_erb);
r_erb=2*real(ifilterbank(c_erb,g_erb,a_erb));
disp('Reconstruction error:')
norm(f-r_erb)

% Plot the FB response
figure(1);
R1=filterbankresponse(g_erb,a_erb,Ls,fs,'real','plot');
title('ERBlet FB response')
ylabel('Magnitude');

% Plot frequency responses of individual filters
gd_erb=filterbankrealdual(g_erb,a_erb,Ls);
figure(2);
subplot(2,1,1);
filterbankfreqz(g_erb,a_erb,Ls,fs,'plot','linabs','posfreq');
title('ERBlet analysis filters')
subplot(2,1,2);
filterbankfreqz(gd_erb,a_erb,Ls,fs,'plot','linabs','posfreq');
title('Synthesis, dual filters')

%% Generate a uniform Barklet FB with V= 3 and a cosine window
disp('--- Uniform Barklet FB, V=3, cosine window ---')
[g_bark,a_bark,fc_bark]=audfilters(fs,Ls,fmin,fmax,'bark','uniform','spacing',1/3,'cosine');
c_bark=filterbank(f,{'realdual',g_bark},a_bark);
r_bark=2*real(ifilterbank(c_bark,g_bark,a_bark));
disp('Reconstruction error:')
norm(f-r_bark)

% Plot the FB response
figure(3);
R2=filterbankresponse(g_bark,a_bark,Ls,fs,'real','plot');
title('Barklet FB response')
ylabel('Magnitude');

% Plot frequency responses of individual filters
gd_bark=filterbankrealdual(g_bark,a_bark,Ls);
figure(4);
subplot(2,1,1);
filterbankfreqz(g_bark,a_bark,Ls,fs,'plot','linabs','posfreq');
title('Barklet analysis filters')
subplot(2,1,2);
filterbankfreqz(gd_bark,a_bark,Ls,fs,'plot','linabs','posfreq');
title('Synthesis, dual filters')

%% Generate a fractional uniform Mel FB with 30 filters and a triangular window
disp('--- Fractional uniform Mel FB, M=30, triangular window ---')
[g_mel,a_mel,fc_mel,L_mel]=audfilters(fs,Ls,fmin,fmax,'mel','fractionaluniform','M',30,'tria');
if L_mel > Ls
    f = postpad(f,L_mel);
end
c_mel=filterbank(f,{'realdual',g_mel},a_mel);
r_mel=2*real(ifilterbank(c_mel,g_mel,a_mel));
disp('Reconstruction error:')
norm(f-r_mel)

% Plot the FB response
figure(5);
R3=filterbankresponse(g_mel,a_mel,Ls,fs,'real','plot');
title('Mel FB response')
ylabel('Magnitude');

% Plot frequency responses of individual filters
gd_mel=filterbankrealdual(g_mel,a_mel,L_mel);
figure(6);
subplot(2,1,1);
filterbankfreqz(g_mel,a_mel,Ls,fs,'plot','linabs','posfreq');
title('Mel analysis filters')
subplot(2,1,2);
filterbankfreqz(gd_mel,a_mel,L_mel,fs,'plot','linabs','posfreq');
title('Synthesis, dual filters')