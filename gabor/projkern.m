function c=projkern(c,g,a,varargin)
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
%   `cout=projkern(cin,g,a,M)` does the same but projects onto Gabor
%   coefficients with M channels
% 
%   `cout=projkern(cin,g,a,M,'ga',ga)` uses the analysis window ga
%
%   Additionally, the function accepts the following flag:
%
%       'tight' A tight window generated from g will be used for
%               analysis and synthesis.
%
%   The rationale for this function is a follows: Because the coefficient
%   space of a Gabor frame is larger than the signal space (since the frame
%   is redundant) then there are many coefficients that correspond to the
%   same signal. Therefore, you might desire to work with the coefficients
%   *cin*, but you are in reality working with *cout*.
%
%   See also: dgt, idgt

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.M = [];
definput.keyvals.ga = g;
definput.flags.proj = {'same','tight'};
[flags,kv] = ltfatarghelper({'M'},definput,varargin);

M = kv.M;
ga = kv.ga;

if isempty(kv.M); M=size(c,1); end

assert_squarelat(a,M,1,'PROJKERN');

switch flags.proj
    case 'same'
        gs=g;
        c=dgt(idgt(c,gs,a),ga,a,M);
    case 'tight'
        gs=gabtight(g,a,M);
        ga=gs;
        c=dgt(idgt(c,gs,a),ga,a,M);
end



