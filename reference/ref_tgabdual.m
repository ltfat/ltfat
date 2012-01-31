function gd=ref_tgabdual(ttype,g,L,a,M)
%REF_WIN   Compute appropriate dual window for transform
%

if nargin<5
  error('Too few input parameters.');
end;

info=ref_transforminfo(ttype,L,a,M);

gd=info.winscale*gabdual(g,info.a,info.M);


