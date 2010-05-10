function cout=comp_col2diag(cin);
%COMP_COL2DIAG  Compute matrix represenation from spreading symbol
%
%  This function is its own inverse.

%   AUTHOR : Peter Soendergaard.
%   TESTING: OK
%   REFERENCE: OK

L=size(cin,1);
cout=zeros(L);

jj=(0:L-1).';
for ii=0:L-1
  cout(ii+1,:)=cin(ii+1,mod(ii-jj,L)+1);
end;

