function c = gga(f,indvec)
%GGA Generalized Goertzel algorithm
%   Usage:  y = gga(x,indvec)
%
%   Input parameters:
%         x      : Input data.
%         indvec : Indices to calculate. 
%
%   Output parameters:
%         c      : Coefficient vector.
%
% `gga(f,indvec)` computes DTFT of one-dimensional signal *f* at 'indices'
% contained in `indvec`, using the generalized second-order Goertzel algorithm.
% Thanks to the generalization, the 'indices' can be non-integer valued
% in the range 0 to N-1, where N is the length of vector X.
% (Index 0 corresponds to the DC component.)
% Integers in INDVEC result in the classical DFT coefficients.
%
% The output *c* is a column complex vector of length LENGTH(INDVEC) containing
% the desired DTFT values.
%
% If the *f* is a matrix, the Goertzel algorithm is applied to each of *W*
% columns. 
%
% Remark:
% Besides the generalization the algorithm is also shortened by one
% iteration compared to the conventional Goertzel.
%
% References: syra2012goertzel
       
% The original copyright goes to
% 2013 Pavel Rajmic, Brno University of Technology, Czech Rep.


%% Check the input arguments
if nargin < 2
    error('%s: Not enough input arguments',upper(mfilename))
end

if isempty(f)
    error('%s: X must be a nonempty vector or a matrix',upper(mfilename))
end

if ~isvector(indvec) || isempty(indvec)
    error('%s: INDVEC must be a nonempty vector',upper(mfilename))
end

if ~isreal(indvec)
    error('%s: INDVEC must contain real numbers',upper(mfilename))
end

% if isinteger(indvec)
%     disp('Warning: The traditional Goertzel algorithm is a bit more effective in case of INDVEC being integer-valued')
% end

[f,L] =comp_sigreshape_pre(f,upper(mfilename),0);

% if any(indvec>L)
%     error('%s: INDVEC must contain real numbers less than %i.',upper(mfilename),L)
% end

c = comp_gga(f,indvec);

