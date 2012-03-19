function [xo]=groupthresh(xi,lambda,varargin)
%GROUPTHRESH   Group thresholding
%   Usage:  xo=groupthresh(xi,lambda);
%
%   `groupthresh(x,lambda)` performs group thresholding on *x*, with
%   threshold *lambda*.  *x* must be a two-dimensional array, the first
%   dimension labelling groups, and the second one labelling
%   members. Several types of grouping behaviour are available:
%
%   * `groupthresh(x,lambda,'group')` shrinks all coefficients within a given
%     group according to the value of the $l^2$ norm of the group in
%     comparison to the threshold *lambda*. This is the default.
%
%   * `groupthresh(x,lambda,'elite')` shrinks all coefficients within a
%     given group according to the value of the $l^1$ norm of the
%     group in comparison to the threshold value *lambda*.
%
%   `groupthresh` accepts all the flags of |thresh|_ to choose the
%   thresholding type within each group and the output type (full / sparse
%   matrix). Please see the help of |thresh|_ for the available
%   options. Default is to use soft thresholding and full matrix output.
%  
%   See also:  thresh
%
%   Demos:  demo_audioshrink
%
%   References: Kowalski08sparsity kowalski2009mixed

%   AUTHOR : Bruno Torresani.  
%   REFERENCE: OK
 
if nargin<2
  error('Too few input parameters.');k
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.import={'thresh'};
definput.importdefaults={'soft'};
definput.flags.gtype={'group','elite'};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

NbGroups = size(xi,1);
NbMembers = size(xi,2);

if flags.do_sparse
  xo = sparse(size(xi));
else
  xo = zeros(size(xi));
end;


if flags.do_group
  for g=1:NbGroups,
    threshold = norm(xi(g,:));
    mask = (1-lambda/threshold);
    mask = mask * (mask>0);
    if flags.do_hard
      mask = (mask>0);
    else      
      mask = mask * (mask>0);
    end;
    xo(g,:) = xi(g,:) * mask;
  end
end;

if flags.do_elite
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
end;
