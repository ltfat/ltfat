function c=comp_filterbank_td(f,g,a,skip,ext)  
%COMP_UFILTERBANK_TD   Non-uniform filterbank by conv2
%   Usage:  c=comp_ufilterbank_fft(f,g,a,ext);
%
%   Input parameters:
%         f   : Input data.
%         g   : Filterbank filters
%         a   : Subsampling factors
%         skip: Delay of the filters.
%         ext : Border exension technique.
%
%   Output parameters:
%         c  : Coefficients

%input data length
L=size(f,1);
%input channel number
W=size(f,2);
%filter number
M=size(g,2);
%filter length
filtLen = size(g,1);
% Allow filter delay only in the filter support range
if(any(skip>=filtLen) || any(skip)<0)
  error('%s: The filter zero index position outside of the filter support.', upper(mfilename));  
end


% Determine output lengths
% Lext -- length of the signal after convolution before subsampling
% N -- after subsampling
if(strcmp(ext,'per'))
   Lext = L;
   N = ceil(Lext./a);
else
   Lext = (L+filtLen-1);
   N = ceil((Lext-skip)./a); 
end
Lreq = a.*(N-1) + 1;

% Output memory allocation
c=cell(M,1);
for m=1:M
  c{m}=zeros(N(m),W);
end;

% Explicitly extend the input. length(fext) = length(f) + 2*(filtLen-1)
fext = comp_extBoundary(f,filtLen-1,ext);
% CONV2 with 'valid' does 2-D linear convolution and crops (filtLen-1) samples from both ends.  
% length(fextconv2) = length(f) + (filtLen-1)

% W channels are done simultaneously
for m=1:M
   c{m} = comp_downs(conv2(fext,g(:,m),'valid'),a(m),skip(m),Lreq(m)); 
end


 

