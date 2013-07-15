function [h,g,a]=wfilt_maxflat(N)
%WFILT_MAXFLAT Generates maximally flat FIR filters
%   Usage: [h,g,a] = wfilt_maxflat(N);
%
%   `[h,g,a]=wfilt_maxflat(N)` calculates half-band maximally flat FIR filters,
%   where $(N-1)$ is the degree of flatness at $w=0$ and $w=\pi$ radians. 
%
%
%   References: vaidy93mult

% Original copyright goes to:
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es



Npi=N;

poly=[];		% Calculate trigonometric polynomial
for i=1:N
	poly=[poly , 2*numcomb(Npi+i-2,i-1)];
end
poly=poly(length(poly):-1:1);
zerospoly=roots(poly);	% Calculate roots

% Transform roots

rootsz=[];

for i=1:length(zerospoly)
    rooty=zerospoly(i);
    rootz1=(1-2*rooty)-2*sqrt(rooty*(rooty-1));
    rootz2=(1-2*rooty)+2*sqrt(rooty*(rooty-1));
    rootsz=[rootsz,rootz1,rootz2];
end     

zeros=rootsz;
N=length(zeros);

% To construct rh for the minimum phase choice, we choose all the zeros 
% inside the unit circle. 

modulus=abs(zeros);

j=1;
for i=1:N
	if modulus(i)<1
		zerosinside(j)=zeros(i);
		j=j+1;
	end
end

An=poly(1);

realzeros=[];
imagzeros=[];
numrealzeros=0;
numimagzeros=0;


Ni=length(zerosinside);

for i=1:(Ni)
	if (imag(zerosinside(i))==0)
		numrealzeros=numrealzeros+1;
		realzeros(numrealzeros)=zerosinside(i);
	else
		numimagzeros=numimagzeros+1;
		imagzeros(numimagzeros)=zerosinside(i);	
		
	end
end

% Construction of rh from its zeros

rh=[1 1];

for i=2:N
	rh=conv(rh,[1 1]);
end

for i=1:numrealzeros
	rh=conv(rh,[1 -realzeros(i)]);
end

for i=1:2:numimagzeros
	rh=conv(rh,[1 -2*real(imagzeros(i)) abs(imagzeros(i))^2]);
end

% Normalization

rh=sqrt(2)/sum(rh)*rh;

g{1} = rh;
g{2} = -(-1).^(1:length(rh)).*g{1}(end:-1:1);
 
h = g;

a= [2;2];

function y=numcomb(n,k)

if n==k,
   y=1;
elseif k==0,
   y=1;
elseif k==1,
   y=n;
else 
   y=fact(n)/(fact(k)*fact(n-k));
end

function y=fact(x)

for j=1:length(x)
    if x(j)==0,
       y(j)=1;
    else
       y(j)=x(j)*fact(x(j)-1);
    end
end



      




