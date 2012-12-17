function outsig=frana(F,insig);
%FRANA  Frame analysis operator
%   Usage: c=frana(F,f);
%
%   `c=frana(F,f)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |frame|_ or |framepair|_.
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
%   See also: frame, framepair, frsyn, plotframe
  
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
    outsig=F.g'*insig;  
  case 'dgt'
    outsig=framenative2coef(F,dgt(insig,F.g,F.a,F.M,F.vars{:}));
  case 'dgtreal'
    outsig=framenative2coef(F,dgtreal(insig,F.g,F.a,F.M,F.vars{:}));
  case 'dwilt'
    outsig=framenative2coef(F,dwilt(insig,F.g,F.M));
  case 'wmdct'
    outsig=framenative2coef(F,wmdct(insig,F.g,F.M));
    
  case 'filterbank'
    outsig=framenative2coef(F,filterbank(insig,F.g,F.a));
  case 'filterbankreal'
    outsig=framenative2coef(F,filterbank(insig,F.g,F.a));
 case 'ufilterbank'
   outsig=framenative2coef(F,ufilterbank(insig,F.g,F.a));
  case 'ufilterbankreal'
    outsig=framenative2coef(F,ufilterbank(insig,F.g,F.a));
    
  case 'nsdgt'
    outsig=framenative2coef(F,nsdgt(insig,F.g,F.a,F.M));
  case 'unsdgt'
    outsig=framenative2coef(F,unsdgt(insig,F.g,F.a,F.M));
  case 'nsdgtreal'
    outsig=framenative2coef(F,nsdgtreal(insig,F.g,F.a,F.M));
  case 'unsdgtreal'
    outsig=framenative2coef(F,unsdgtreal(insig,F.g,F.a,F.M));
    
  case {'dft',...
        'dcti','dctii','dctiii','dctiv',...
        'dsti','dstii','dstiii','dstiv'}
    outsig=feval(F.type,insig);
  case 'fusion'
    % All frames must use the same length signal.
    L=framelength(F,size(insig,1));
    insig=postpad(insig,L);
    
    coefs = cell(F.Nframes,1);
    for ii=1:F.Nframes
        coefs(ii)={F.w(ii)*frana(F.frames{ii},insig)};
    end;
    outsig=cell2mat(coefs);
end;

  
