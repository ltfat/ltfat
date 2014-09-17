function f = ifwt2(c,w,J,varargin)
%IFWT2   Inverse Fast Wavelet Transform 
%   Usage:  f = ifwt2(c,w,J)
%           f = ifwt2(c,w,J,Ls,...)
%
%   Input parameters:
%         c     : Coefficients stored in a matrix.
%         w     : Wavelet filters definition.
%         J     : Number of filterbank iterations.
%         Ls    : Size of the reconstructed signal.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   `f = ifwt2(c,w,J)` reconstructs signal *f* from the wavelet coefficients
%   *c* using a *J*-iteration synthesis filterbank build from the basic 
%   synthesis filterbank defined by *w*. *f* is a matrix with 
%   `size(f)==size(c)`.
%
%   `f = ifwt2(c,w,J,Ls)` works as above but the result *f* is cut or
%   extended to size *Ls* if *Ls* is a two-element vector or to `[Ls,Ls]`
%   if *Ls* is a scalar.
%
%   This function takes the same optional parameters as |fwt2|. Please see
%   the help on |fwt2| for a description of the parameters.
%   
%   See also:  fwt2, fwtinit
%
%   Demos: demo_imagecompression
%
%   References: ma98

% AUTHOR: Zdenek Prusa

if nargin<3
   error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(c)
  error('%s: Unrecognized coefficient format.',upper(mfilename));
end;

% Initialize the wavelet filters structure
w = fwtinit(w);


%% PARSE INPUT
definput.keyvals.Ls=[];    
definput.import = {'fwt','fwt2'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

if (isempty(Ls))
   Ls = size(c);
end

if (numel(Ls)==1)
  Ls = [Ls,Ls];
end

Lcrows = fwtclength(Ls(1),w,J,'per');
Lccols = fwtclength(Ls(2),w,J,'per');
nFilts = numel(w.g);

if flags.do_standard
  Jstep = 1;
  for jj=1:J-1
    LcIdx =  jj*(nFilts-1)+2;
    colRange = 1:Lcrows(LcIdx);
    rowRange = 1:Lccols(LcIdx);
    c(colRange,rowRange) = ifwt(c(colRange,rowRange),w,Jstep,Lcrows(LcIdx),'dim',1,'per');
    c(colRange,rowRange) = ifwt(c(colRange,rowRange),w,Jstep,Lccols(LcIdx),'dim',2,'per');
  end

  c = ifwt(c,w,Jstep,Ls(1),'dim',1,'per');
  f = ifwt(c,w,Jstep,Ls(2),'dim',2,'per');
  
end;

if flags.do_tensor
  f = ifwt(c,w,J,Ls(1),'dim',1,'per');
  f = ifwt(f,w,J,Ls(2),'dim',2,'per');
end;

