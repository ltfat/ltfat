lw=2;
fz=20;

data_fb_real =load('much_longer_fb_real.log');
data_fac_real=load('much_longer_fac_real.log');
data_ola_real=load('much_longer_ola_real.log');


% Columns in data, fb : a M L W gl time
% Columns in data, fac: a M L W time
% Columns in data, ola: a M L W gl bl time

Ls=data_fac_real(:,3);
t_fb_real =data_fb_real(:,6);
t_fac_real=data_fac_real(:,5);
t_ola_real=data_ola_real(:,7);

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

M=data_fac_real(1,2);
L=data_fac_real(:,3);

plot(Ls,t_fb_real,l1,...  
     Ls,t_fac_real,l2,...
     Ls,t_ola_real,l3,'LineWidth',lw);
set(gca,'Fontsize',fz);

h=legend('Portnoff, real','Fac, real','Fac-OLA, real',...
       'Location','NorthWest');

% Grow the box a little, otherwise the export to .eps is messed up.
q=get(h,'Position');
set(h,'Position',[q(1)*.9 q(2)*.95 q(3)*1.8 q(4)]);

xlabel('Signal length / samples','Fontsize',fz);
ylabel('Running time / seconds','Fontsize',fz);

print -deps plot_much_longer_1.eps



