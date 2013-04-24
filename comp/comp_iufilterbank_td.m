function f=comp_iufilterbank_td(c,g,a,Ls,skip,ext)  
%COMP_IUFILTERBANK_TD   Synthesis Uniform filterbank by conv2
%   Usage:  f=comp_iufilterbank_td(c,g,a,Ls,skip,ext);
%
%   Input parameters:
%         c    : N*M*W array of coefficients.
%         g    : Filterbank filters - filtLen*M array. 
%         a    : Upsampling factor - scalar.
%         Ls   : Output length.
%         skip : Delay of the filters - scalar or array of length M.
%         ext  : Border exension technique.
%
%   Output parameters:
%         f  : Output Ls*W array. 
%

%input channel number
W=size(c,3);
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

% Output memory allocation
f=zeros(Ls,W,assert_classname(c,g));

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


 

