function L=filterbanklengthcoef(coef,a)
%FILTERBANKLENGTHCOEF  Filterbank length from coefficients
%   Usage: L=filterbanklengthcoef(coef,a);
%
%   `filterbanklengthcoef(coef,a)` returns the length of a filterbank with
%   time-shifts *a*, such that the filterbank is long enough to expand the
%   coefficients *coef*.
%
%   If instead a signal is given, call |filterbanklength|.
%
%   See also: filterbank, filterbanklength

complainif_notenoughargs(nargin,2,upper(mfilename));


if iscell(coef)
  cl=cellfun(@(x) size(x,1),coef);
else
  Mcoef=size(coef,2);
  cl=ones(1,Mcoef)*size(coef,1);
end;

cl=cl(:);

% Make 'a' have the length of '
if isvector(a)
    a=bsxfun(@times,a,ones(numel(cl),1));
    a=a(:);

    L=a.*cl;
else
    L=a(:,1).*cl./a(:,2);
end;


if var(L)>0
  error(['%s: Invalid set of coefficients. The product of the no. of ' ...
         'coefficients and the channel time shift must be the same for ' ...
         'all channels.'],upper(mfilename));
end;

L=L(1);


