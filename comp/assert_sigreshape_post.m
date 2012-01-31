function f=assert_sigreshape_post(f,dim,permutedsize,order)

% Restore the original, permuted shape.
f=reshape(f,permutedsize);

if dim>1
  % Undo the permutation.
  f=ipermute(f,order);
end;




