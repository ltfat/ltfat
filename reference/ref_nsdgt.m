function [c,Ls] = ref_nsdgt(f,g,a,M)
%NSDGT  Non-stationary Discrete Gabor transform
%   Usage:  c=nsdgt(f,g,a,M);
%           [c,Ls]=nsdgt(f,g,a,M);

timepos=cumsum(a)-a(1)+1;

Ls=length(f);
L=nsdgtlength(Ls,a);
N=numel(g);

F=zeros(L,sum(M));


% Construct the analysis operator matrix explicitly
Y = 1;
for n = 1:length(timepos)
  X = length(g{n});
  win_range = mod(timepos(n)+(-floor(X/2):ceil(X/2)-1)-1,L)+1;
  F(win_range,Y) = fftshift(g{n}); 
  for m = 1:M(n)-1
    F(win_range,Y+m) = F(win_range,Y).*exp(2*pi*i*m*(-floor(X/2):ceil(X/2)-1)/M(n)).';
  end
  Y=Y+M(n);
end


cmat=F'*f;

c=mat2cell(cmat,M);

