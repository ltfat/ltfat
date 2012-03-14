function outsig=frsyn(F,insig);
%FRSYN  Frame synthesis operator
%   Usage: f=frsyn(F,c);
%
%   `f=frsyn(F,c)` constructs a signal *f* from the frame coefficients *c*
%   using the frame *F*. The frame object *F* must have been created using
%   |newframe|_.
%
%   See also: newframe, frana, plotframe
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;

switch(F.type)
 case 'gen'
  outsig=F.gs*insig;
  
 case 'dgt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);
  outsig=idgt(insig,F.gs,F.a);
  
 case 'dgtreal'
  [MN,W]=size(insig);
  M2=floor(F.M/2)+1;
  N=MN/M2;
  insig=reshape(insig,[M2,N,W]);
  outsig=idgtreal(insig,F.gs,F.a,F.M);
  
 case 'dwilt'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[2*F.M,N/2,W]);
  outsig=idwilt(insig,F.gs);
  
 case 'wmdct'
  [MN,W]=size(insig);
  N=MN/F.M;
  insig=reshape(insig,[F.M,N,W]);
  outsig=iwmdct(insig,F.gs);
  
 case 'ufilterbank'
  [MN,W]=size(insig);
  M=numel(F.gs);
  N=MN/M;
  insig=reshape(insig,[N,M,W]);
  outsig=ifilterbank(insig,F.gs,F.a);   
  
 case {'dcti','dctiv','dsti','dstiv'}
  outsig=feval(F.type,insig);
  
 case 'dctii'
  outsig=dctiii(insig);  
 case 'dctiii'
  outsig=dctii(insig);
 case 'dstii'
  outsig=dstiii(insig);
 case 'dstiii'
  outsig=dstii(insig);
 case 'dft'
  outsig=idft(insig);
 case 'fft'
  outsig=ifft(insig);
 case 'fftreal'
  outsig=ifftreal(insig,F.L);  
end;

  