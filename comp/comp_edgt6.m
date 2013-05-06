function cout=comp_edgt6(cin,a)
%COMP_EDGT6   Compute Even DGT type 6
%

M=size(cin,1);
N=size(cin,2)/2;
W=size(cin,3);

cout=zeros(M,N,W,assert_classname(cin));

cout=cin(:,1:N,:);

cout=reshape(cout,M*N,W);

