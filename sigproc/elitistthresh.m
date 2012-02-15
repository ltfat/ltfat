function [xo]=elitistthresh(xi,lambda,varargin)
%ELITISTTHRESH   elitist (hard/soft) thresholding
%   Usage:  xo=elitistthresh(xi,lambda);
%
%   `elitistthresh(xi,lambda)` performs hard elitist thresholding on *xi*,
%   with threshold *lambda*. The input *xi* must be a two-dimensional array,
%   the first dimension labelling groups, and the second one labelling
%   members.  All coefficients within a given group are shrunk according to
%   the value of the $l^1$ norm of the group in comparison to the threshold
%   value *lambda*.
%
%   `elitistthresh(x,lambda,'soft')` will do the same using soft
%   thresholding.
%
%   `elitistthresh` accepts the following flags at the end of the line of input
%   arguments:
%
%     'hard'    Perform hard thresholding. This is the default.
%
%     'soft'    Perform soft thresholding.  
%
%     'full'    Return the output as a full matrix. This is the default.
%
%     'sparse'  Return the output as a sparse matrix.
%  
%   See also:  groupthresh
%
%   Demos:  demo_audioshrink
%
%   References: Kowalski08sparsity kowalski2009mixed

%   AUTHOR : Bruno Torresani.  
 
if nargin<2
  error('Too few input parameters.');
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.flags.iofun={'hard','soft'};
definput.flags.outclass={'full','sparse'};

[flags,keyvals]=ltfatarghelper({},definput,varargin,mfilename);

NbGroups = size(xi,1);
NbMembers = size(xi,2);

if flags.do_sparse
  xo = sparse(size(xi));
else
  xo = zeros(size(xi));
end;

for g=1:NbGroups,
    y = sort(abs(xi(g,:)),'descend');
    rhs = cumsum(y);
    rhs = rhs .* lambda ./ (1 + lambda * (1:NbMembers));
    M_g = find(diff(sign(y-rhs)));
    if (M_g~=0)
        tau_g = lambda * norm(y(1:M_g),1)/(1+lambda*M_g);
    else
        tau_g = 0;
    end
    xo(g,:) = thresh(xi(g,:),tau_g,flags.iofun,flags.outclass);
end
