function L=filterbanklengthcoef(coef,a);
%FILTERBANKLENGTHCOEF  Filterbank length from coefficients
%   Usage: L=filterbanklengthcoef(coef,a);
%
%   `filterbanklengthcoef(coef,a)` returns the length of a filterbank with
%   time-shifts *a*, such that the filterbank is long enough to expand the
%   coefficients *coef*.
%
%   If instead a signal is given, call |filterbanklengthsignal|_.
%
%   See also: filterbank, filterbanklengthsignal
  
if iscell(coef)
  Mcoef=numel(coef);
  cl=cellfun(@numel,coef);
else
  Mcoef=size(coef,2);
  cl=ones(1,Mcoef)*size(coef,1);    
end;

a=a(:);
cl=cl(:);

[a,~]=scalardistribute(a,cl);

L=a.*cl;

if var(L)>0
  error(['%s: Invalid set of coefficients. The product of the no. of ' ...
         'coefficients and the channel time shift must be the same for ' ...
         'all channels.'],upper(mfilename));
end;

L=L(1);


