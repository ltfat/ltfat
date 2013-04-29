% This script test the quality of the Hermite implementation.

% There seem to be some numerical problems for ii>66

L=200;
H=pherm(L,0:159,'fast','qr');


H1=H(:,1:4:end);
H2=H(:,2:4:end);
H3=H(:,3:4:end);
H4=H(:,4:4:end);

norm(H1'*H1)
norm(H2'*H2)
norm(H3'*H3)
norm(H4'*H4)

norm(H1'*H2)
norm(H1'*H3)
norm(H1'*H4)

norm(H2'*H3)
norm(H2'*H4)

norm(H3'*H4)

norm(abs(H./dft(H))-1)

