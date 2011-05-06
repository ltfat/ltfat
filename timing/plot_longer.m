lw=2;

data_fb =load('longer_fb.log');
data_fac=load('longer_fac.log');
data_fb_real =load('longer_fb_real.log');
data_fac_real=load('longer_fac_real.log');



% Columns in data, fb : a M L W gl time
% Columns in data, fac: a M L W time
Ls=data_fac(:,3);
t_fb =data_fb(:,6);
t_fac=data_fac(:,5);
t_fb_real =data_fb_real(:,6);
t_fac_real=data_fac_real(:,5);

if 0
  % Color legend
  l1='b';
  l2='b--';
  l3='r';
  l4='r--';
else
  % bw legend
  l1='b';
  l2='b--';
  l3='b-.';
  l4='b:';
end;

figure(1);

M=data_fac(1,2);
L=data_fac(:,3);

plot(Ls,t_fb,l1,...  
     Ls,t_fac,l2,...
     Ls,t_fb_real,l3,...  
     Ls,t_fac_real,l4,'LineWidth',lw);

legend('Portnoff','Fac','Portnoff, real','Fac, real',...
       'Location','NorthWest');

xlabel('Signal length / samples');
ylabel('Running time / seconds');
set(gca,'fontsize',16);

print -deps plot_longer_1.eps

