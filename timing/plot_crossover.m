%PLOT_CROSSOVER

fontsize=16;
lw=2;

data=load('crossover.log');
dataref=load('crossover.ref');
data_real=load('crossover_real.log');
dataref_real=load('crossover_real.ref');

% Columns in data: a M L W gl time

% overlapping factor
%x=data(:,5)./data(:,1);

% gl
x=data(:,5);

t=data(:,6);
t_real=data_real(:,6);

tref=dataref(:,5);
tref_real=dataref_real(:,5);

% Make a horizontal line in the plot.
plotref=ones(size(x,1),1)*tref;
plotref_real=ones(size(x,1),1)*tref_real;

% Compute the flopcounts based on the setup from the data.
N=size(data,1);
flop_fac = zeros(N,1);
flop_fb  = zeros(N,1);
flop_fac_real = zeros(N,1);
flop_fb_real  = zeros(N,1);

for ii=1:N
  [fcfac,fcfb,fcfac_real,fcfb_real] = flopcounts(data(ii,1),data(ii,2),data(ii,3),data(ii,5));
  
  flop_fac(ii)=fcfac*data(ii,4);
  flop_fb(ii)=fcfb*data(ii,4);
  flop_fac_real(ii)=fcfac_real*data(ii,4);
  flop_fb_real(ii)=fcfb_real*data(ii,4);

end;

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
set(gca,'fontsize',fontsize);

plot(x,flop_fb,l1,...
     x,flop_fac,l2,...
     x,flop_fb_real,l3,...
     x,flop_fac_real,l4,'LineWidth',lw);
set(gca,'Fontsize',fz);

h=legend('Portnoff','Fac','Portnoff, real','Fac, real',...
       'Location','NorthWest');

% Grow the box a little, otherwise the export to .eps is messed up.
q=get(h,'Position');
set(h,'Position',[q(1)*0.90 q(2)*.95 q(3)*1.6 q(4)]);

%xlabel('Overlapping factor');
xlabel('Window length / samples','Fontsize',fz);
ylabel('Flop count / flop','Fontsize',fz);
title('Flop count comparison','Fontsize',fz);

print -deps plot_crossover_1.eps

figure(2);
set(gca,'fontsize',fontsize);

plot(x,t,l1,...
     x,plotref,l2,...
     x,t_real,l3,...
     x,plotref_real,l4,'LineWidth',lw);
set(gca,'Fontsize',fz);

h=legend('Portnoff','Fac','Portnoff, real','Fac, real',...
       'Location','NorthWest');

% Grow the box a little, otherwise the export to .eps is messed up.
q=get(h,'Position');
set(h,'Position',[q(1)*0.90 q(2)*.95 q(3)*1.6 q(4)]);

%xlabel('Overlapping factor');
xlabel('Window length / samples','Fontsize',fz);
ylabel('Execution time / seconds','Fontsize',fz);
title('Execution time comparison','Fontsize',fz);

print -deps plot_crossover_2.eps

