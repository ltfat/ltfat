function outsig=comp_frsyn_fusion(F,insig)

W=size(insig,2);
L=size(insig,1)/framered(F);

outsig=zeros(L,W);

idx=0;    
for ii=1:F.Nframes
    coeflen=L*framered(F.frames{ii});
    outsig=outsig+frsyn(F.frames{ii},insig(idx+1:idx+coeflen,:))*F.w(ii);
    idx=idx+coeflen;
end;
