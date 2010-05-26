%PLOT_CROSSOVER

data=load('crossover.log');
dataref=load('crossover.ref');

% Columns in data: a M L W gl time
x=data(:,5)./data(:,1);
t=data(:,6);
tref=dataref(:,5);
% Make a horizontal line in the plot.
plotref=ones(size(x,1),1)*tref;

% Compute the flopcounts based on the setup from the data.
N=size(data,1);
flop_fac = zeros(N,1);
flop_fb  = zeros(N,1);
for ii=1:N
  [fcfac,fcfb] = flopcounts(data(ii,1),data(ii,2),data(ii,3),data(ii,5));
  
  flop_fac(ii)=fcfac*data(ii,4);
  flop_fb(ii)=fcfb*data(ii,4);
end;

l1='b';
l2='b--';


figure(1);
set(gca,'fontsize',16);
set(gca,'LineWidth',2);

plot(x,flop_fb,l1,...
     x,flop_fac,l2);
legend('Poisson','Fac');

xlabel('Overlapping factor');
ylabel('Flop count / flop');
title('Flop count comparison');

figure(2);
set(gca,'fontsize',16);
set(gca,'LineWidth',2);

plot(x,t,l1,...
     x,plotref,l2);

legend('Poisson','Fac');
xlabel('Overlapping factor');
ylabel('Execution time / seconds');
title('Execution time comparison');