function [cout]=comp_iedgt6(cin,a,M)
%COMP_IEDGT6   Compute inverse even DGT type 6
%

N=size(cin,1)/M;
W=size(cin,2);
L=N*a;

cin=reshape(cin,M,N,W);

cout=zeros(M,2*N,W,assert_classname(cin));
cout(:,1:N,:)=cin;

% Copy the non modulated coefficients.
cout(1,N+1:2*N,:)=cin(1,N:-1:1,:);

% Copy the modulated coefficients.
cout(2:M,N+1:2*N,:)=-cin(M:-1:2,N:-1:1,:);



