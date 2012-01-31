function cout=ref_col2diag(cin);
%REF_COL2DIAG  Compute matrix represenation from spreading symbol
%
%  This function is its own inverse.
  
L=size(cin,1);
cout=zeros(L);

for ii=0:L-1
  for jj=0:L-1
    cout(ii+1,jj+1)=cin(ii+1,mod(ii-jj,L)+1);
  end;
end;


