function cout=ref_col2diag_1(cin);
%REF_COL2DIAG_1  Compute matrix represenation from spreading symbol
%
%  This function is its own inverse.
  
L=size(cin,1);
cout=zeros(L);

for jj=0:L-1
  for ii=0:jj-1
    cout(ii+1,jj+1)=cin(ii+1,ii-jj+L+1);
  end;
    for ii=jj:L-1
    cout(ii+1,jj+1)=cin(ii+1,ii-jj+1);
  end;
end;

% The second code also works.
if 0

  for ii=0:L-1
    for jj=0:ii
      cout(ii+1,jj+1)=cin(ii+1,ii-jj+1);
    end;
    for jj=ii+1:L-1
      cout(ii+1,jj+1)=cin(ii+1,ii-jj+L+1);
    end;  
  end;
end;


