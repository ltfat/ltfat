function cout = gabprojkern(cin,g,a,varargin)
%PROJKERN  Projection onto generating kernel space
%   Usage:  cout=projkern(cin,g,a);
%           cout=projkern(cin,g,a,varargin);
%
%   Input parameters:
%         cin   : Input coefficients
%         g     : analysis/synthesis window
%         a     : Length of time shift.
%   Output parameters:
%         cout  : Output coefficients
%
%   `cout=projkern(cin,g,a)` projects a set of Gabor coefficients *cin* onto the
%   space of possible Gabor coefficients, specified by the window g, time-shift
%   a and the same number of channels as cin. This means that *cin* and *cout*
%   synthesize to the same signal.
%
%   `cout=projkern(cin,g,a,M,'real')` has to be used if cin comes from
%   dgtreal with M channels
% 
%
%   Additionally, the function accepts the following flag:
%
%       'tight' A tight window generated from g will be used for
%               analysis and synthesis.
%
%       'same' The window g is used for analysis and synthesis.
%
%   The rationale for this function is a follows: After modification of 
%   the Gabor coefficient one does in general not stay in the range of the 
%   Gabor frame analysis operator. This function projects the modified
%   coefficients *cin*, but you are in reality working with *cout*.
%
%   See also: dgt, idgt

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.ga = g;
definput.keyvals.M = [];
definput.flags.proj = {'dual','tight','same'};
definput.flags.real = {'complex','real'};
[flags,kv] = ltfatarghelper({'M'},definput,varargin);

if isempty(kv.M)
    kv.M=size(cin,1);
end

M = kv.M;
N = size(cin,2);
L = N*a;

ga = gabwin(kv.ga,a,M);

assert_squarelat(a,M,1,'PROJKERN');

switch flags.proj
    case 'dual'
        gs=gabdual(ga,a,M,L);
    case 'tight'
        gs=gabtight(ga,a,M,L);
        ga=gs;
    case 'same'
        gs=ga;
end

switch flags.real
    case 'complex'
        cout=dgt(idgt(cin,gs,a,M),ga,a,M,L);
    case 'real'
        cout=dgtreal(idgtreal(cin,gs,a,M),ga,a,M);
end



