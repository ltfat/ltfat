function [A,B]=ref_wilbounds(g,a,M);

L=length(g);

N=L/a;

F=ref_idwilt(eye(M*N),g,a,M);

S=F'*F;

d=eig(S);

A=min(d);
B=max(d);

