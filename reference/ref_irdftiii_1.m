function f=ref_irdftiii_1(c)
%REF_RDFTIII_1  Reference IRDFTIII by IDFTIII
%   Usage:  f=ref_irdftiii(c);
%
%   Only works for real functions.
%
%   Compute IRDFT by doubling the signal and doing an IDFTIII
%   to obtain the reconstructed signal

L=size(c,1);
Lhalf=floor(L/2);
Lend=Lhalf*2;


cc=zeros(size(c));

cc(1:Lhalf,:)=(c(1:2:Lend,:)- i*c(2:2:Lend,:))/sqrt(2);
cc(L-Lhalf+1:end,:)= (c(Lend-1:-2:1,:)  +i*c(Lend:-2:2,:))/sqrt(2);

% If f has an even length, we must also copy the Nyquest-wave
% (it is real)
if mod(L,2)==1
  cc((L+1)/2,:)=c(L,:);
end;

f=real(ref_idftiii(cc));








