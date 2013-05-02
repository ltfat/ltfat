function cout=comp_irdgt(cin,a,M)
%COMP_IRDGT  Compute inverse real DGT.

N=size(cin,1)/M;
W=size(cin,2);
L=N*a;

cin=reshape(cin,M,N,W);

Mhalf=ceil(M/2);
Mend=Mhalf*2-1;

cout=zeros(M,N,W,assert_classname(cin));

% Copy the first coefficient, it is real
cout(1,:,:)=cin(1,:,:);

cout(2:Mhalf,:,:)=(cin(2:2:Mend,:,:)- i*cin(3:2:Mend,:,:))/sqrt(2);
cout(M-Mhalf+2:M,:,:)= (cin(Mend-1:-2:2,:,:)  +i*cin(Mend:-2:3,:,:))/sqrt(2);

% If f has an even length, we must also copy the Nyquest-wave
% (it is real)
if mod(M,2)==0
  cout(M/2+1,:,:)=cin(M,:,:);
end;




