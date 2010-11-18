function h=pfilt(f,g,a,dim)
%PFILT  Apply filter with periodic boundary conditions
%   Usage:  h=pfilt(f,g);
%           h=pfilt(f,g,a,dim);
%
%   PFILT(f,g) applies the filter g to the input f. If f is a matrix, the
%   filter is applied along each column.
%
%   PFILT(f,g,a) does the same, but downsamples the output keeping only
%   every a'th sample (starting with the first one).
%
%   PFILT(f,g,a,dim) filters along dimension dim.
%
%   See also: pconv
  
% Assert correct input.
error(nargchk(2,4,nargin));

if nargin<4
  dim=1;
end;

if nargin<3
  a=1;
end;

L=[];

[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,'PFILT');

[g,info] = comp_fourierwindow(g,L,'PFILT');

g=fir2long(g,L);

% Force FFT along dimension 1, since we have permuted the dimensions
% manually
h=ifft(fft(f,L,1).*repmat(fft(g,L,1),1,W),L,1);

h=h(1:a:end,:);

permutedsize(1)=size(h,1);

h=assert_sigreshape_post(h,dim,permutedsize,order);
