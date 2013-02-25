function order=assert_groworder(order)
%ASSERT_GROWORDER  Grow the order parameter
%
%   `assert_groworder` is meant to be used in conjunction with
%   `assert_sigreshape_pre` and `assert_sigreshape_post`. It is used to
%   modify the *order* parameter in between calls in order to expand the
%   processed dimension by 1, i.e. for use in a routine that creates 2D
%   output from 1D input, for instance in `dgt` or `filterbank`.
  
if numel(order)>1
  % We only need to handle the non-trivial order, where dim>1
  
  p=order(1);
  
  % Shift orders higher that the working dimension by 1, to make room for
  % the new dimension, but leave lower dimensions untouched.
  order(order>p)=order(order>p)+1;
  
  order=[p,p+1,order(2:end)];
end;


