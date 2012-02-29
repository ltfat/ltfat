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

for ii = 1:N
  c(ii,:) = modify(argcolorone(w(ii)),w(ii));
end

% Handle NaNs
idx=isnan(w(:));
c(idx,:)=repmat([1,1,1],sum(idx),1);

c=reshape(c,[remembersize,3]);


function c2 = modify(c1,w)

wf = 2*abs(w)/(abs(w).^2+1);
if abs(w) < 1
    c2 = c1 * wf;
else
    c2 = c1 + (1-c1)*(1-wf); 
end



function col = argcolorone(w)

a = angle(w);
idx=a<0;
a(idx)=a(idx)+2*pi;

a = a/pi*3; % a in [0,6)

switch floor(a)
 case 0
  col = [1,a,0];
 case 1
  col = [2-a,1,0];
 case 2
  col = [0,1,a-2];
 case 3        
  col = [0,4-a,1];
 case 4
  col = [a-4,0,1];
 otherwise
  col = [1,0,6-a];
end;

%
