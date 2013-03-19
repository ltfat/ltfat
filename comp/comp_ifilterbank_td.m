function f=comp_ifilterbank_td(c,g,a,Ls,skip,ext)  
%COMP_UFILTERBANK_TD   Synthesis Uniform filterbank by conv2
%   Usage:  f=comp_iufilterbank_fft(c,g,a,Ls,skip,ext);
%
%   Input parameters:
%         c    : Cell array of length M, each element is N(m)*W matrix.
%         g    : Filterbank filters - length M cell-array, each element is vector of length filtLen(m) 
%         a    : Upsampling factors - array of length M.
%         skip : Delay of the filters - scalar or array of length M.
%         Ls   : Output length.
%         ext  : Border exension technique.
%
%   Output parameters:
%         f  : Output Ls*W array. 
%

%input channel number
W=size(c{1},2);
%filter number
M=numel(g);
%length of filters
filtLen = cellfun(@(x) numel(x),g(:));


% Allow filter delay only in the filter support range
if(any(skip(:)>=filtLen) || any(skip)<0)
  error('%s: The filter zero index position outside of the filter support.', upper(mfilename));  
end

if(numel(skip)==1)
    skip = skip*ones(M,1);
end


% Output memory allocation
f=zeros(Ls,W);

if(~strcmp(ext,'per'))
    ext = 'zero';
end

skipOut = a.*(filtLen-1)+skip(:);

% W channels are done simultaneously
for m=1:M
   cext = comp_extBoundary(c{m},filtLen(m)-1,ext,'dim',1); 
   ftmp = conv2(g{m}(:),comp_ups(cext,a(m)));
   f = f + ftmp(1+skipOut(m):Ls+skipOut(m),:); 
end
