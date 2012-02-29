maxval=25;
exponent=1.5;

Lx=40;
Ly=101;

scale=log(maxval+1).^(1/exponent);

x=linspace(0,2*pi,Lx);
y=linspace(0,1,Ly);
yy=exp((scale*y).^exponent)-1;
dom=repmat(yy.',1,Lx).*exp(i*repmat(x,Ly,1));

figure(1);
image(x,y,domaincolor(dom));
pbaspect([.5,1,1])
xlabel('Phase (rad)');
ylabel('Magnitude');

figure(2);
plot(y,yy)

figure(3);

c=dgtreal(greasy,'gauss',10,588);

s=20*log10(abs(c));
s=s-max(s(:))+90;
s(s<0)=0;
s=s/90;

ss=exp((scale*s).^exponent)-1;
dom=ss.*exp(i*angle(c));

image(domaincolor(dom));