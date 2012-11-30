function [h,g,a]=wfilt_spline(m,n)

% WSPLINE  Generates spline wavelets.
%
%          [H,G,RH,RG]=WSPLINE(M,N) returns the analysis and
%          synthesis filters corresponding to a biortoghonal 
%          scheme with spline wavelets of compact support. 
%          H is the analysis lowpass filter, RH the synthesis 
%          lowpass filter, G the analysis highpass filter and
%          RG the synthesis highpass filter.
%          N and M specify the number of zeros in z=-1 required 
%          for H(z) and RH(z) respectively.
% 
%          N+M must be even. 
%
%          With these examples is possible to achieve arbitrarily 
%          high regularity. For large M, the analysis wavelet
%          will belong to C^k if N>4.165M+5.165(k+1).

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es


if(rem(m+n,2)~=0)
    error('M+N must be even.');
end

% Calculate rh coefficients, RH(z)=sqrt(2)*((1+z^-1)/2)^m;

rh=sqrt(2)*(1/2)^m*binewton(m);

% Calculate h coefficients, H(-z)=sqrt(2)*((1+z^-1)/2)^n*P(z)

% First calculate P(z) (pol)

if (rem(n,2)==0)
   N=n/2+m/2;
else
   N=(n+m-2)/2+1;
end

pol=trigpol(N);

% Now calculate ((1+z*-1)/2)^n;

r0=(1/2)^n*binewton(n);


hrev=sqrt(2)*conv(r0,pol);

l=length(hrev);
hh=hrev(l:-1:1);


[h{2}, g{2}]=calhpf(hh,rh);
h{1} = hh;
g{1} = rh;

if(length(h{1})>length(h{2}))
    if(rem(length(h{1}),2)~=1)
       r0 = (length(h{1})-length(h{2}))/2;
       l0 = r0;
    else
       r0 = (length(h{1})-length(h{2}))/2+1;
       l0 = (length(h{1})-length(h{2}))/2-1;
    end
      h{2} = [zeros(1,l0), h{2}, zeros(1,r0) ];
else
    if(rem(length(h{1}),2)~=1)
       r0 = (length(h{2})-length(h{1}))/2;
       l0 = r0;
    else
       r0 = (length(h{2})-length(h{1}))/2+1;
       l0 = (length(h{2})-length(h{1}))/2-1;
    end
      h{1} = [zeros(1,l0), h{1}, zeros(1,r0) ];
end

if(length(g{1})>length(g{2}))
    if(rem(length(g{1}),2)~=1)
       r0 = (length(g{1})-length(g{2}))/2;
       l0 = r0;
    else
       r0 = (length(g{1})-length(g{2}))/2+1;
       l0 = (length(g{1})-length(g{2}))/2-1;
    end
      g{2} = [zeros(1,l0), g{2}, zeros(1,r0) ];
else
    if(rem(length(g{1}),2)~=1)
       r0 = (length(g{2})-length(g{1}))/2;
       l0 = r0;
    else
       r0 = (length(g{2})-length(g{1}))/2+1;
       l0 = (length(g{2})-length(g{1}))/2-1;
    end
      g{1} = [zeros(1,l0), g{1}, zeros(1,r0) ];
end

a= [2;2];




