outsig = comp_frana_tensor(F,insig)

outsig=frsyn(F.frames{1},insig);
perm=circshift((1:F.Nframes).',-1);
for ii=2:F.Nframes
    outsig=permute(outsig,perm);
    outsig=frsyn(F.frames{ii},outsig);
end;
outsig=permute(outsig,perm);
