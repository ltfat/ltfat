function cout=rect2wil(cin);
%RECT2WIL  Inverse of WIL2RECT
%   Usage:  c=rect2wil(c);
%   
%   `rect2wil(c)` takes Wilson coefficients processed by |wil2rect| and
%   puts them back in the compact form used by |dwilt| and |idwilt|. The
%   coefficients can then be processed by |idwilt|.
%
%   See also: wil2rect, dwilt, idwilt
  
%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

complainif_argnonotinrange(nargin,1,1,mfilename);
  
M=size(cin,1)-1;
N=size(cin,2);
W=size(cin,3);

cout=zeros(2*M,N/2,W,assert_classname(cin));

if rem(M,2)==0
  for ii=0:N/2-1
    cout(1:M+1  ,ii+1,:)=cin(1:M+1,2*ii+1,:);
    cout(M+2:2*M,ii+1,:)=cin(2:M,2*ii+2,:);
  end;
else
  for ii=0:N/2-1
    cout(1:M    ,ii+1,:)=cin(1:M,2*ii+1,:);
    cout(M+2:2*M,ii+1,:)=cin(2:M,2*ii+2,:);
    cout(M+1    ,ii+1,:)=cin(M+1,2*ii+2,:);
  end;  
end;

