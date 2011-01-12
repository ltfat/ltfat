function order=assert_groworder(order)
%ASSERT_GROWORDER  Grow the order parameter
%
%   ASSERT_GROWORDER is meant to be used in conjunction with
%   ASSERT_SIGRESHAPE_PRE and ASSERT_SIGRESHAPE_POST. It is used to
%   modify the order parameter in between call in order to expand the
%   processed dimension by 1, i.e. for use in a routine that creates 2D
%   output from 1D input, for instance in DGT or FILTERBANK.
  
if numel(order)>1
  % We only need to handle the non-trivial order, where dim>1
  
  p=order(1);
  
  % Shift orders higher that the working dimension by 1, to make room for
  % the new dimension, but leave lower dimensions untouched.
  order(order>p)=order(order>p)+1;
  
  order=[p,p+1,order(2:end)];
end;
