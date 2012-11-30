function [H, G, a] = wfilt_db(N,varargin)
%DBFILT    Daubechies FIR filterbank
%   Usage:  [H,G] = dbfilt(N);
%           [H,G] = dbfilt(N,J);
%           [H,G,a] = dbfilt(...);
%
%   Input parameters:
%         N     : Order of Daubechies filters. 
%         J     : Number of filterbank iterations, default J=1.
%   Output parameters:
%         H     : J+1 cell array of analysing filters impulse reponses
%         G     : J+1 cell array of synthetizing filters impulse reponses
%         a     : J+1 array of subsampling (or hop) factors accociated with
%                 corresponding filters
%
%   `[H,G] = dbfilt(N)` computes two-channel Daubechies FIR filterbank from prototype
%   maximum-phase analysing lowpass filter obtained by spectral factorization of the Lagrange interpolator filter. 
%   *N* also denotes number of zeros at z=-1 of the lowpass filters of length 2*N. 
%   The prototype lowpass filter is in a form (all roots of R(z) are outside of the unit circle):
%                 
%   .. H_l(z)=(1+z^-1)^N*R(z),
%
%   .. math:: H_l(z)=\left(1+z^{-1}\right)R(z),
%   
%   where R(z) is a spectral factor of P(z)=2R(z)*R(z^-1) 
%   All subsequent filters of the two-channel filterbank are derived as
%   follows:
%
%   .. H_h(z)=H_l((-z)^-1)
%   .. G_l(z)=H_l(z^-1)
%   .. G_h(z)=-H_l(-z)
%
%    making them orthogonal causal perfect-reconstruction QMF.
%
%   `[H,G] = dbfilt(N,J)` computes one-iteration J+1 channel Daubechies FIR filterbank
%   equivalent to the J iterations of the basic two-channel filterbank
%   using multirate identity.

%   References: Daubechies, Kovacevic, MultirateId

if(nargin<1)
   error('%s: Too few input parameters.',upper(mfilename));
end

if(N<1)
    error('Minimum N is 1.');
end
if(N>20)
    warning('Instability may occur when N is too large.');
end

definput.keyvals.J = 1;
[flags,kv,J]=ltfatarghelper({'J'},definput,varargin);

H = cell(2,1);
flen = 2*N;

% Calculating Lagrange interpolator coefficients
sup = [-N+1:N];
a = zeros(1,N);
for n = 1:N
    non  = sup(sup ~= n);
    a(n) = prod(0.5-non)/prod(n-non);
end
P = zeros(1,2*N-1);
P(1:2:end) = a;
P = [P(end:-1:1),1,P];

R = roots(P);
R = R(abs(R)>1 & real(R)>0);

% roots of the 2*conv(lo_a,lo_r) filter
hroots = [R(:);-ones(N,1)];


% building synthetizing low-pass filter from roots
H{1}= real(poly(hroots));
% normalize to norm(H{1})=sqrt(2)
H{1}= sqrt(2)*H{1}/sum(H{1});
% QMF modulation low-pass -> highpass
H{2}= (-1).^(1:flen).*H{1}(end:-1:1);


if(nargout>1)
   % Building reconstruction filterbank
   G = cell(2,1); 
   % flip
   G{1} = H{1}(end:-1:1);
   % modulation
   G{2} = -(-1).^(1:flen).*H{1};

   if(nargout>2)
       a = [2;2];
   end
   if(J>1)
        % make J level multirate identity filterbank
        G = makeMultirateIdentity(G{1},G{2},J);
       if(nargout>2)
         a = zeros(J+1,1);
         a(1) = 2^J;
         for j=1:J a(end+1-j)= 2^j; end;
       end
   end
end

if(J>1)
    % make J level multirate identity filterbank
    H = makeMultirateIdentity(H{1},H{2},J);
end





