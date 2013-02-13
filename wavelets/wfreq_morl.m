function [H] = wfreq_morl(s,N);



%H = cell(length(s),1);
if rem(N,2) == 0
    % N is even
    w = 2*pi*[0:N/2-1]/N;             % frequency axis: [0, pi)
else
    % N is odd
    w = 2*pi*[0:(N-1)/2]/N;           % frequency axis: [0, pi)
end

w = linspace(-8,8,N);
w0 = 6;



H = (1/sqrt(abs(s))).*sqrt(2)*pi^(1/4).*exp(-((w.*s-w0).^2)/2);
H = H./max(H);

% for ii=1:length(s)
%  
%   if rem(N,2) == 0
%     % N is even
%     H{ii} = [Htmp; 0; Htmp(end:-1:2)];
%   else
%     % N is odd
%     H{ii} = [Htmp; Htmp(end:-1:2)];
%   end
% end





