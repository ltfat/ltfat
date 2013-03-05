function c=comp_ufilterbank_td(f,g,a,skip,ext)  
%COMP_UFILTERBANK_TD   Uniform filterbank by conv2
%   Usage:  c=comp_ufilterbank_fft(f,g,a,skip,ext);
%
%   Input parameters:
%         f   : Input data.
%         g   : Filterbank filters. 
%         a   : Subsampling factor.
%         skip: Delay of the filters.
%         ext : Border exension technique.
%
%   Output parameters:
%         c  : Coefficients 
%


%input data length
L=size(f,1);
%input channel number
W=size(f,2);
%filter number
M=size(g,2);
%length of filters
filtLen = size(g,1);
% Allow filter delay only in the filter support range
if(all(skip>=filtLen) || all(skip<0))
  error('%s: The filter zero index position outside of the filter support.', upper(mfilename));  
end

if(numel(skip)==1)
    skip = skip*ones(M,1);
end

% Determine output length
% Lext -- length of the signal after convolution before subsampling
% N -- after subsampling
if(strcmp(ext,'per'))
   Lext = L;
   N = ceil(Lext/a);
else
   Lext = (L+filtLen-1);
   N = ceil((Lext-skip)/a); 
end
%The minimum input signal length which produces N output samples
Lreq = a*(N-1) + 1;

% Output memory allocation
c=zeros(N,M,W);

% Explicitly extend the input. length(fext) = length(f) + 2*(filtLen-1)
fext = comp_extBoundary(f,filtLen-1,ext);
% CONV2 does 2-D linear convolution. 
% length(fextconv2) = length(f) + 3*(filtLen-1)
for m=1:M
  c(:,m,:) = comp_downs(conv2(fext,g(:,m),'valid'),a,skip(m),Lreq); 
end;

 

