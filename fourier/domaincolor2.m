function c = domaincolor(w)
%DOMAINCOLOR  Convert complex values to color
%   Usage: c = domaincolor(x);
%
%   Example:
%   --------
%
%   The following example shows a patch of the complex plane:::
%
%     L=51;
%     R=2;
%     x=linspace(-R,R,L);
%     dom=repmat(x,L,1)+i*repmat(x.',1,L);
%     image(x,x,domaincolor(dom));
%     axis('xy');
%     xlabel('Real axis');
%     ylabel('Imaginary axis');
%
%   The following example create a `colorbar`-like representation for a
%   mapping of values $0\leq x\leq 1$ using the following mapping for the magnitude:
%
%   .. math:: y = e^{(2.2x)^{1.5}}-1
%
%   :::
%
%     Lx=40;
%     Ly=101;
%     x=linspace(0,2*pi,Lx);
%     y=linspace(0,1,Ly);
%     dom=(exp(repmat(2.2*y.',1,Lx).^1.5)-1).*exp(i*repmat(x,Ly,1));
%     image(x,y,domaincolor(dom));
%     pbaspect([.5,1,1])
%     xlabel('Phase (rad)');
%     ylabel('Magnitude');
%
  
remembersize=size(w);
N=numel(w);
c = zeros(N,3);

%c = hsv2rgb([angle(w(:))/(2*pi)+.5,1-abs(w(:)),abs(w(:))]);
c = hsv2rgb([angle(w(:))/(2*pi)+.5,ones(N,1),abs(w(:))]);

max(abs(w(:)))
c=reshape(c,[remembersize,3]);
