f = gspi;

J = 10;
w = fwtinit({'db',8});
c = fwt(f,w,J);
figure(1);
plotfwt(c);
figure(2);
plot(f);
