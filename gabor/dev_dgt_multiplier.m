%this script
%1. calculates a lattice lat based on L, a, and b
%2. calculates the Gabor matrix G [L x MN]
%3. calculates the STFT matrix [L x L] therefrom

clear all;
clc;
L=100;
symbols = ones(L, L);
sg = [0.01 0.1 0.5 1 1.5]; %TODO: find suitable values

%create (separable) lattice of length L
a=1;
b=1;
sep = 0; %"separability constant" when given in in lower triangular Hermite normal form
M=abs(L/b);
N=abs(L/a);

% make (separable) lattice lat
lat=zeros(2,M*N);
for n=0:N-1
  for m=0:M-1
    offset = mod(sep*n,b);
    % Determine gridpoint in rectangular coordinates.
    lat(1,m+n*M+1) = n*a;
    lat(2,m+n*M+1) = m*b + offset;
  end
end
lat=mod(lat,L);

% calculate Gabor matrix G
MN=size(lat,2);
G=zeros(L,MN);
jj=(0:L-1).';

%figure;
for ii = 1:numel(sg)
    %s = sqrt(sg(ii));
    s=sg(ii);
    g = pgauss(L, 'width', s*L); %calculates pgauss(L,pi*s^2/4L*log(2))

    for p=1:MN
        G(:,p)=exp(-2*pi*1i*lat(2,p)*jj/L).*circshift(g,lat(1,p));
    end

%derive STFT matrix via summation over G
Gs = reshape(G, L, L, ceil(size(G,2)/L));
V = sum(Gs, 3)./L;
V_inv = conj(V);

%calculate multiplier (which dims? - shorten to [M x N]?)
%mult = V_inv.* V.* symbols;

%    mult = idftmatrix.*premul(1:L, 1:L).*symbols;

%    plot(real(eig(mult)))
%    hold on
end
%imagesc(abs(mult))
