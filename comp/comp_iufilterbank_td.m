function f=comp_iufilterbank_td(c,g,a,Ls,skip,ext)  
%COMP_UFILTERBANK_TD   Synthesis Uniform filterbank by conv2
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

%input length
N=size(c,1);
%input channel number
W=size(c,2);
%filter number
M=size(g,2);
%length of filters
filtLen = size(g,1);
% Allow filter delay only in the filter support range
if(skip>=filtLen || skip<0)
  error('%s: The filter zero index position outside of the filter support.', upper(mfilename));  
end

% Output memory allocation
f=zeros(Ls,W);

if(~strcmp(ext,'per'))
    ext = 'zero';
end

skipOut = a*(filtLen-1)+skip;

% W channels are done simultaneously
for m=1:M
   cext = comp_extBoundary(squeeze(c(:,m,:)),filtLen-1,ext,'dim',1); 
   ftmp = conv2(g(:,m),comp_ups(cext,a));
   f = f + ftmp(1+skipOut(m):Ls+skipOut(m),:); 
end


 

