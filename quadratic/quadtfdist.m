function p = quadtfdist(f, q)
%QUADTFDIST Quadratic time-frequency distribution
%   Usage p = quadtfdist(f, q);
%
%   Input parameters:
%         f  : Input vector:w
%         q  : Kernel
%
%   Output parameters:
%         p  : Quadratic time-frequency distribution
% 
%   For an input vector of length L, the kernel should be a L x L matrix.
%   `quadtfdist(f, q);` computes a discrete quadratic time-frequency 
%   distribution. 
%

% AUTHOR: Jordy van Velthoven

complainif_notenoughargs(nargin, 2, 'QUADTFDIST');

[M,N] = size(q);

if ~all(M==N)
  error('%s: The kernel should be a square matrix.', upper(mfilename));
end

[f,~,Ls,W,~,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

p = comp_quadtfdist(f, q);

