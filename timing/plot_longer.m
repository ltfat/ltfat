data_fb =load("longer_fb.log");
data_fac=load("longer_fac.log");

% Columns in data, fb : a M L W gl time
% Columns in data, fac: a M L W time
Ls=data_fb(:,3);
t_fb =data_fb(:,6);
t_fac=data_fac(:,5);
figure(1);

plot(Ls,t_fb,Ls,t_fac);
legend('FB','FAC');
xlabel('Signal length / samples');
ylabel('Running time / seconds');
