function cout=comp_col2diag(cin);
%COMP_COL2DIAG  transforms columns to diagonals (in a special way)
%
%  This function transforms the first column to the main diagonal. The
%  second column to the first side-diagonal below the main diagonal and so
%  on. 
% 
%  This way fits well the connection of matrix and spreading function, see
%  spreadfun.
%
%  This function is its own inverse.

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

L=size(cin,1);
cout=zeros(L,assert_classname(cin));

jj=(0:L-1).';
for ii=0:L-1
  cout(ii+1,:)=cin(ii+1,mod(ii-jj,L)+1);
end;



