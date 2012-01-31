function f=ref_iedgt_1(c,g,a,M)
%REF_IEDGT_1   Reference Inverse Even DGT by IDGTII
%   Usage  c=ref_edgt(f,g,a,M);
%
%   The input window must be odd-centered of length 2L.
%
%   a must be divisable by 2.

L=size(g,1)/2;
W=size(c,2);

N=L/a;
M=L/b;

clong=zeros(M,2*N,W);

% Copy the first M/2 coefficients of the first/last time shift.
clong(1:M/2,1,:)=c(1:M/2,:);
clong(1:M/2,N+1,:)=c(1+(N-1)*M+M/2:N*M,:);

% Scale the first coefficients correctly
clong(1,1,:)=clong(1,1,:)*sqrt(2);
clong(1,N+1,:)=clong(1,N+1,:)*sqrt(2);

% Copy the remaining coefficients to the first/last time shift, such
% that we get an odd symmetry.
clong(M:-1:M/2+2,1,:)   = -clong(2:M/2,1,:);
clong(M:-1:M/2+2,N+1,:) = -clong(2:M/2,N+1,:);

% Copy the body of the coefficients
clong(:,2:N,:)=c(M/2+1:M/2+M*(N-1),:);

% Copy the unmodulated coefficients for the second half
clong(1,N+2:2*N,:)=clong(1,N:-1:2,:);

% Copy the modulated coefficients for the second half
clong(2:M,N+2:2*N,:)=-clong(M:-1:2,N:-1:2,:);

clong=reshape(clong,2*M*N,W);

fdouble=ref_idgtii(clong,g,a,2*b);

f=fdouble(1:L,:);


