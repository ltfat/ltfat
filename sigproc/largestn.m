function [xo,Nout]=largestn(xi,N,varargin)
%LARGESTN   Keep N largest coefficients
%   Usage:  xo=largestn(x,N);
%           xo=largestn(x,N,mtype);
%
%   `largestn(x,N)` returns an array of the same size as *x* keeping
%   the *N* largest coefficients.
%
%   `largestn` takes the following flags at the end of the line of input
%   arguments:
%
%     'hard'    Perform hard thresholding. This is the default.
%
%     'wiener'  Perform empirical Wiener shrinkage. This is in between
%               soft and hard thresholding.
%
%     'soft'    Perform soft thresholding.  
%
%     'full'    Returns the output as a full matrix. This is the default.
%
%     'sparse'  Returns the output as a sparse matrix.   
%
%   If the coefficients represents a signal expanded in an orthonormal
%   basis then this will be the best N-term approximation.
%
%   **Note:** If soft- or Wiener thresholding is selected, only $N-1$
%   coefficients will actually be returned. This is caused by the *N*'th
%   coefficient being set to zero.
%
%   See also:  largestr
%
%   References: ma98

%   AUTHOR : Peter L. Søndergaard and Bruno Torresani.  
%   TESTING: OK
%   REFERENCE: OK

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'thresh'};
[flags,keyvals]=ltfatarghelper({},definput,varargin);

if (prod(size(N))~=1 || ~isnumeric(N))
  error('N must be a scalar.');
end;

if flags.do_sparse
  if ndims(xi)>2
    error('Sparse output is only supported for 1D/2D input. This is a limitation of Matlab/Octave.');
  end;
end;

% Determine the size of the array.
ss=numel(xi);

% Sort the absolute values of the coefficients.
sxi=sort(abs(xi(:)));


% Find the coeffiecient sitting at position N through the array,
% and use this as a threshing value. 
if N<=0
    % Choose a thresh value higher than max
    lambda=sxi(end)+1;
else
    lambda=sxi(ss-N+1);
end;

[xo,Nout]=thresh(xi,lambda,flags.outclass,flags.iofun);

