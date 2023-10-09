function [f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,L,dim,callfun)
%ASSERT_SIGRESHAPE_PRE  Preprocess and handle dimension input.
%
%   Input parameters:
%      f            : signal, possibly ND-array
%      L            : L parameter
%      dim          : dim parameter
%      callfun      : Name of calling function
%   Output parameters:
%      f            : Input signal as matrix
%      L            : Verified L
%      Ls           : Length of signal along dimension to be processed
%      W            : Number of transforms to do.
%      dim          : Verified dim
%      permutedsize : pass to assert_sigreshape_post
%      order        : pass to assert_sigreshape_post
  
  
if ~isnumeric(f)
  error('%s: The input must be numeric.',callfun);
end;

D=ndims(f);

% Dummy assignment.
order=1;

if isempty(dim)
  dim=1;

  if sum(size(f)>1)==1
    % We have a vector, find the dimension where it lives.
    dim=find(size(f)>1);
  end;

else
  if (numel(dim)~=1 || ~isnumeric(dim))
    error('%s: dim must be a scalar.',callfun);
  end;
  if rem(dim,1)~=0
    error('%s: dim must be an integer.',callfun);
  end;
  if (dim<1) || (dim>D)
    error('%s: dim must be in the range from 1 to %d.',callfun,D);
  end;  

end;

if (numel(L)>1 || ~isnumeric(L))
  error('%s: L must be a scalar.',callfun);
end;
if (~isempty(L) && rem(L,1)~=0)
  error('%s: L must be an integer.',callfun);
end;


if dim>1
  order=[dim, 1:dim-1,dim+1:D];

  % Put the desired dimension first.
  f=permute(f,order);

end;

Ls=size(f,1);

% If L is empty it is set to be the length of the transform.
if isempty(L)
  L=Ls;
end;  

% Remember the exact size for later and modify it for the new length
permutedsize=size(f);
permutedsize(1)=L;
  
% Reshape f to a matrix.
if ~isempty(f)
  f=reshape(f,size(f,1),numel(f)/size(f,1));
end;
W=size(f,2);





