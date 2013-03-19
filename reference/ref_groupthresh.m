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
%   By default, `groupthresh` chooses the groups to be the row vectors of
%   the input (the vectors along the 2nd dimension). This can be changed by
%   calling `groupthresh(x,lambda,dim)`, where *dim* is the dimension to
%   **not** group along.
%
%   `groupthresh` accepts all the flags of |thresh| to choose the
%   thresholding type within each group and the output type (full / sparse
%   matrix). Please see the help of |thresh| for the available
%   options. Default is to use soft thresholding and full matrix output.
%  
%   See also:  thresh
%
%   Demos:  demo_audioshrink
%
%   References: Kowalski08sparsity kowalski2009mixed yu2008audio

%   AUTHOR : Bruno Torresani.  
%   REFERENCE: OK
 
if nargin<2
  error('Too few input parameters.');k
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.import={'thresh','groupthresh'};
definput.importdefaults={'soft'};
definput.keyvals.dim=1;

[flags,keyvals,dim]=ltfatarghelper({'dim'},definput,varargin);

[xi,L,NbGroups,NbMembers,dim,permutedsize,order]=assert_sigreshape_pre(xi,[],dim,'GROUPTHRESH');

if flags.do_sparse
  xo = sparse(size(xi));
else
  xo = zeros(size(xi));
end;

if flags.do_group
  
  groupnorm = sqrt(sum(abs(xi).^2, 2));
  w = thresh(groupnorm, lambda, flags.iofun,flags.outclass)./groupnorm;
  
  % Clean w for NaN. NaN appears if the input has a group with norm
  % exactly 0.
  w(isnan(w)) = 0;
  
  xo = xi.*repmat(w,1,NbMembers);

end

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
    
    % FIXME: The following line does not work for sparse matrices.
    xo(g,:) = thresh(xi(g,:),tau_g,flags.iofun,flags.outclass);
  end
end;

xo=assert_sigreshape_post(xo,dim,permutedsize,order);

