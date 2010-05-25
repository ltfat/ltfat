data=load("crossover.log");
dataref=load("crossover.ref");

% Columns in data: a M L W gl time
x=data(:,5)./data(:,1);
t=data(:,6);
tref=dataref(:,5);
% Make a horizontal line in the plot.
plotref=ones(size(x,1),1)*tref;

figure(1);

plot(x,t,x,plotref);
legend('FB','FAC');
xlabel('Overlapping factor');
ylabel('Running time / seconds');
