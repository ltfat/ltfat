data_fb =load('longer_fb.log');
data_fac=load('longer_fac.log');

% Columns in data, fb : a M L W gl time
% Columns in data, fac: a M L W time
Ls=data_fac(:,3);
t_fb =data_fb(:,6);
t_fac=data_fac(:,5);
figure(1);

M=data_fac(1,2);
L=data_fac(:,3);
NL=nextfastfft(L/M);
Nt_fac=t_fac(NL-9);

plot(Ls,t_fb,Ls,Nt_fac);
legend('FB','FAC','Location','SouthEast');
xlabel('Signal length / samples');
ylabel('Running time / seconds');


