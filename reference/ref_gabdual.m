function gd=ref_gabdual(g,a,M)
%REF_GABDUAL   Reference GABDUAL
%   Usage:  gd=ref_gabdual(g,a,M);
%
%   Calculate the canonical dual window by simple linear algebra
  
G=tfmat('dgt',g,a,M);

gd=(G*G')\g;