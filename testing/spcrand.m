function f=spcrand(n1,n2,p);
%SPCRAND   Sparse Random complex numbers for testing.
%   Usage: f=sptester_crand(n1,n2,p);

% Make a random real valued matrix, extract the indices, put complex
% numbers in and recollect.
f=sprand(n1,n2,p);

[row,col,val]=find(f);

L=numel(val);
val=rand(L,1)-.5+i*(rand(L,1)-.5);

f=sparse(row,col,val,n1,n2);

