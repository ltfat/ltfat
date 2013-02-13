function [H] = wfreq_mexh(N,s);

if rem(N,2) == 0
    % N is even
    w = 2*pi*[0:N/2-1]/N;             % frequency axis: [0, pi)
else
    % N is odd
    w = 2*pi*[0:(N-1)/2]/N;           % frequency axis: [0, pi)
end

% frequency axis
%w=[0:2*pi/L:2*pi*(1-1/L)];
%w(1)=eps;
%w(L/2+1)=w(L/2+1)+1e-15;

w = linspace(-1,1,N);
sigma = 1;
Htmp = exp(-sigma^2*(s*w).^2).*(s*w).^2.*(-sqrt(8)*sigma^(5/2)*pi^(1/4)/sqrt(3));

H = Htmp;

% if rem(N,2) == 0
%     % N is even
%     H = [Htmp, 0, Htmp(end:-1:2)];
% else
%     % N is odd
%     H = [Htmp, Htmp(end:-1:2)];
% end

