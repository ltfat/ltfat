function outsig=frsynadj(F,insig);
%FRSYNADJ  Adjoint synthesis operator
%   Usage: c=frsynadj(F,f);
%
%   `c=frsynadj(F,f)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |newframe|_.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   The output coefficients are stored as columns. This is usually
%   **not** the same format as the 'native' format of the frame. As an
%   examples, the output from |frsynadj|_ for a gabor frame cannot be
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
  case 'identity'
    outsig=insig;    
  case 'gen'
    outsig=F.gs'*insig; 
  case 'dgt'
    outsig=framenative2coef(F,dgt(insig,F.gs,F.a,F.M,F.vars{:}));
  case 'dgtreal'
    outsig=framenative2coef(F,dgtreal(insig,F.gs,F.a,F.M));
  case 'dwilt'
    outsig=framenative2coef(F,dwilt(insig,F.gs,F.M));
  case 'wmdct'
    outsig=framenative2coef(F,wmdct(insig,F.gs,F.M));
    
  case 'filterbank'
    outsig=framenative2coef(F,filterbank(insig,F.gs,F.a));
  case 'filterbankreal'
    outsig=framenative2coef(F,filterbank(insig,F.gs,F.a));
  case 'ufilterbank'
    outsig=framenative2coef(F,ufilterbank(insig,F.gs,F.a));
  case 'ufilterbankreal'
    outsig=framenative2coef(F,ufilterbank(insig,F.gs,F.a));

  case 'nsdgt'
    outsig=framenative2coef(F,nsdgt(insig,F.gs,F.a));
  case 'unsdgt'
    outsig=framenative2coef(F,unsdgt(insig,F.gs,F.a));
  case 'nsdgtreal'
    outsig=framenative2coef(F,nsdgtreal(insig,F.gs,F.a));
  case 'unsdgtreal'
    outsig=framenative2coef(F,unsdgtreal(insig,F.gs,F.a));
    
  case {'dft',...
        'dcti','dctii','dctiii','dctiv',...
        'dsti','dstii','dstiii','dstiv'}
    outsig=feval(F.type,insig);
  case 'fft'
    L=size(insig,1);
    outsig=fft(insig)/L;
  case 'fftreal'
    outsig=fftreal(insig,F.L)/F.L;  
  case 'fusion'
    % All frames must use the same length signal.
    L=framelength(F,size(insig,1));
    insig=postpad(insig,L);

    coefs = cell(F.Nframes,1);
    for ii=1:F.Nframes
        coefs(ii)=frsynadj(F.frames{ii},insig)/(F.Nframes*F.w(ii));
    end;
    outsig=cell2mat(coefs);
end;

  
