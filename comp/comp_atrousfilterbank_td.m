function c=comp_atrousfilterbank_td(f,g,a,offset)  
%COMP_ATROUSFILTERBANK_TD   Uniform filterbank by conv2
%   Usage:  c=comp_atrousfilterbank_fft(f,g,a,skip);
%
%   Input parameters:
%         f   : Input data - L*W array.
%         g   : Filterbank filters - filtLen*M array. 
%         a   : Filter upsampling factor - scalar.
%         offset: Delay of the filters - scalar or array of length M. 
%
%   Output parameters:
%         c  : L*M*W array of coefficients
%


%input data length
L=size(f,1);
%input channel number
W=size(f,2);
%filter number
M=size(g,2);
g = comp_ups(g,a,1);

%length of filters
filtLen = size(g,1);
skip = -offset;
% Allow filter delay only in the filter support range
if(all(skip>=filtLen) || all(skip<0))
  error('%s: The filter zero index position outside of the filter support.', upper(mfilename));  
end

if(numel(skip)==1)
    skip = skip*ones(M,1);
end

% Output memory allocation
c=zeros(L,M,W,assert_classname(f,g));

% Explicitly extend the input. length(fext) = length(f) + 2*(filtLen-1)
fext = comp_extBoundary(f,filtLen-1,'per','dim',1);
% CONV2 does 2-D linear convolution. 'valid' option crops tails
% length(fextconv2) = length(f) + (filtLen-1)
% length(c(:,m,:)) = N
% W channels done simultaneously by conv2
for m=1:M
  c(:,m,:) = comp_downs(conv2(fext,g(:,m),'valid'),1,skip(m),L); 
end;

 

