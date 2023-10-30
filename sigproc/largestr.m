function [xo,Nout]=largestr(xi,p,varargin)
%LARGESTR   Keep fixed ratio of largest coefficients
%   Usage:  xo=largestr(x,p);
%           xo=largestr(x,p,mtype);  
%           [xo,N]=largestr(...);
%
%   `largestr(x,p)` returns an array of the same size as *x* keeping
%   the fraction *p* of the coefficients. The coefficients with the largest
%   magnitude are kept.
%
%   `[xo,n]=largestr(xi,p)` additionally returns the number of coefficients
%   kept.
% 
%   **Note:** If the function is used on coefficients coming from a
%   redundant transform or from a transform where the input signal was
%   padded, the coefficient array will be larger than the original input
%   signal. Therefore, the number of coefficients kept might be higher than
%   expected.
%
%   `largestr` takes the following flags at the end of the line of input
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
%   **Note:** If soft- or Wiener thresholding is selected, one less
%   coefficient will actually be returned. This is caused by that
%   coefficient being set to zero.
%
%   See also:  largestn
%
%   References: ma98

%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.import={'thresh'};
[flags,keyvals]=ltfatarghelper({},definput,varargin);

if (prod(size(p))~=1 || ~isnumeric(p))
  error('p must be a scalar.');
end;

wascell=iscell(xi);

if wascell
  %[xi,shape]=cell2vec(xi);
  xi = cell2mat(xi);
end;

% Determine the size of the array.
ss=numel(xi);

% Determine if p was given as a fraction or as an integer
if floor(p)==p
    N = p;
else
    N=round(ss*p);
end
  
[xo,Nout]=largestn(xi,N,flags.outclass,flags.iofun);

if wascell
  %xo=vec2cell(xo,shape);
  xo = mat2cell(xo,size(xo,1), size(xo,2));
end;
