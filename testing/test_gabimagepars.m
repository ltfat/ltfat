%TEST_GABIMAGEPARS
%
%   This will run a simple test of the gabimagepars routine over a range
%   of sizes, and make a plot of the efficiancy in the end.

x=800;
y=1200;

Ntests=10000;
res=zeros(Ntests,5);
offset=99;
for ii=1:Ntests;  
  Ls=ii+offset;    
  [res(ii,1),res(ii,2),res(ii,3),res(ii,4),res(ii,5)]=gabimagepars(Ls,x,y);  
end;

figure(1);
% res(:,3) is L
plot(res(:,3)./((1+offset:Ntests+offset)'))

