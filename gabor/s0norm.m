function y = s0norm(f,varargin)
%S0NORM S0-norm of signal
%   Usage: y = s0norm(f);
%          y = s0norm(f,...);
%
%   S0NORM(f) computes the S0-norm of a vector.
%
%   If the input is a matrix or ND-array, the RMS is computed along the
%   first (non-singleton) dimension, and a vector of values is returned.
%
%   WARNING: The S0-norm is computed by computing a full Short-time
%   Fourier transform of a signal, which can be quite time-consuming. Use
%   this function with care for long signals.
%
%   S0NORM takes the following flags at the end of the line of input
%   parameters:
%
%-     'dim',d  : Work along specified dimension. The default value of []
%                 means to work along the first non-singleton one.
%

%   AUTHOR : Peter L. Soendergaard
  
%% ------ Checking of input parameters ---------

if ~isnumeric(f) 
  error('%s: Input must be numerical.',upper(mfilename));
end;

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.dim=[];
[flags,kv]=ltfatarghelper({},definput,varargin);

%% ------ Computation --------------------------
 
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim, ...
                                                  upper(mfilename));
permutedsize(1)=1;
y=zeros(permutedsize);

g=pgauss(L);

for ii=1:W  
  c=dgt(f(:,ii),g,1,L);    
  y(1,ii)=sum(abs(c(:)))/L;
end;
  
y=assert_sigreshape_post(y,kv.dim,permutedsize,order);
