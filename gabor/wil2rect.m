function cout=wil2rect(cin);
%WIL2RECT  Arrange Wilson coefficients in a rectangular layout
%   Usage:  c=wil2rect(c);
%
%   `wil2rect(c)` rearranges the coefficients *c* in a rectangular shape. The
%   coefficients must have been obtained from |dwilt|. After rearrangement
%   the coefficients are placed correctly on the time/frequency-plane.
%
%   The rearranged array is larger than the input array: it contains
%   zeros on the spots where the Wilson transform is missing a DC or
%   Nyquest component.
%
%   See also: rect2wil, dwilt, wmdct
  
complainif_argnonotinrange(nargin,1,1,mfilename);
  
M=size(cin,1)/2;
N=size(cin,2)*2;
W=size(cin,3);

cout=zeros(M+1,N,W,assert_classname(cin));

if rem(M,2)==0
  for ii=0:N/2-1
    cout(1:M+1,2*ii+1,:)=cin(1:M+1  ,ii+1,:);
    cout(2:M,2*ii+2,:)  =cin(M+2:2*M,ii+1,:);
  end;
else
  for ii=0:N/2-1
    cout(1:M,2*ii+1,:)  =cin(1:M    ,ii+1,:);
    cout(2:M,2*ii+2,:)  =cin(M+2:2*M,ii+1,:);
    cout(M+1,2*ii+2,:)  =cin(M+1    ,ii+1,:);
  end;  
end;

