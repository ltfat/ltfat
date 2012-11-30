function [H, aa] = wfreq_radwt(p,q,s,N);
% Ilker Bayram, Ivan W. Selesnick: Frequency-domain design of overcomplete rational-dilation wavelet transforms.
% According to http://eeweb.poly.edu/iselesni/rational/RADWT_2009_preprint.pdf
% 
% Ivan Selesnick and Ilker Bayram, Polytechnic Institute, New
% York, November 2008
%
% Two channels: H(omega)->ups(p)->downs(q)
%               G(omega)->downs(s)
%
%  
%
%
% Multirate Identity possible! The lowpass frequency response is localized
% in abs(supp(H(omega))) <= pi*p/q;
%
% input length = lcm(q,s)


aa = [q/p;s;];

% Make H

wp = (s-1)*pi/s;
ws = pi*p/q;

if rem(N,2) == 0
    % N is even
    w = 2*pi*[0:N/2-1]/N;             % frequency axis: [0, pi)
else
    % N is odd
    w = 2*pi*[0:(N-1)/2]/N;           % frequency axis: [0, pi)
end

k_pass = (abs(w) <= wp);                 % pass-band
k_stop = (abs(w) >= ws);                 % stop-band
k_trans = (abs(w) > wp) & (abs(w) < ws);    % transition-band

a = (1-1/s)*pi;
b = p/q - (1 - 1/s);
w_scaled = (abs(w) - a)/b;

Htmp = zeros(size(w.'));
Htmp(k_pass) = 1;
Htmp(k_trans) = (1+cos(w_scaled(k_trans))) .* sqrt(2-cos(w_scaled(k_trans)))/2;


if rem(N,2) == 0
    % N is even
    H(:,1) = [Htmp; 0; Htmp(end:-1:2)];
else
    % N is odd
    H(:,1) = [Htmp; Htmp(end:-1:2)];
end


% Make G

ws = (s-1)*pi/s;
wp = p*pi/q;

if rem(N,2) == 0
    % N is even
    w = 2*pi*[0:N/2-1]/N;             % frequency axis: [0, pi)
else
    % N is odd
    w = 2*pi*[0:(N-1)/2]/N;           % frequency axis: [0, pi)
end

k_pass = (abs(w) >= wp);                 % pass-band
k_stop = (abs(w) <= ws);                 % stop-band
k_trans = (abs(w) > ws) & (abs(w) < wp);    % transition-band

a = (1-1/s)*pi;
b = p/q - (1 - 1/s);
w_scaled = (abs(w) - a)/b;

Htmp = zeros(size(w.'));
Htmp(k_pass) = 1;
Htmp(k_trans) = (1-cos(w_scaled(k_trans))) .* sqrt(2+cos(w_scaled(k_trans)))/2;

if rem(N,2) == 0
    % N is even
    H(:,2) = [Htmp; 1; Htmp(end:-1:2)];
else
    % N is odd
    H(:,2) = [Htmp; Htmp(end:-1:2)];
end

% max(abs(imag(ifft(G))))


% normalization

H(:,1) = H(:,1)*sqrt(p*q);
H(:,2) = H(:,2)*sqrt(s);

