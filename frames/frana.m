function outsig=frana(F,insig);
%FRANA  Frame analysis operator
%   Usage: c=frana(F,f);
%
%   `c=frana(F,f)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |newframe|_.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   The output coefficients are stored as columns. This is usually
%   **not** the same format as the 'native' format of the frame. As an
%   examples, the output from |frana|_ for a gabor frame cannot be
%   passed to |idgt|_ without a reshape.
%
%   See also: newframe, frsyn, plotframe
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
 case 'gen'
  outsig=F.ga'*insig;
 
 case 'dgt'
  outsig=dgt(insig,F.ga,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 
 case 'dgtreal'
  outsig=dgtreal(insig,F.ga,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 
 case 'dwilt'
  outsig=dwilt(insig,F.ga,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 
 case 'wmdct'
  outsig=wmdct(insig,F.ga,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 
 case 'ufilterbank'
  outsig=filterbank(insig,F.ga,F.a);
  [N,M,W]=size(outsig);
  outsig=reshape(outsig,[N*M,W]);
  
 case {'dft','fft','fftreal',...
       'dcti','dctii','dctiii','dctiv',...
       'dsti','dstii','dstiii','dstiv'}
  outsig=feval(F.type,insig);
end;

  