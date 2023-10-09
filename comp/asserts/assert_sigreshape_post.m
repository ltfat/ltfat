function f=assert_sigreshape_post(f,dim,permutedsize,order)
%ASSERT_SIGRESHAPE_POST  Restore dimension input.
%
%   Input parameters:
%      f            : Input signal as matrix
%      dim          : Verified dim
%      permutedsize : pass to assert_sigreshape_post
%      order        : pass to assert_sigreshape_post
%   Output parameters:
%      f            : signal, possibly ND-array
%
%   `assert_sigreshape_post` works in conjunction with
%   `assert_sigreshape_pre` and restores the original
%   dimensions of the input array


% Restore the original, permuted shape.
f=reshape(f,permutedsize);

if dim>1
  % Undo the permutation.
  f=ipermute(f,order);
end;