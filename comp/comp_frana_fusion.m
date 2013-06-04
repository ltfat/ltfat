function outsig=comp_frana_fusion(F,insig)

% All frames must use the same length signal.
L=F.length(size(insig,1));
insig=postpad(insig,L);

coefs = cell(F.Nframes,1);
for ii=1:F.Nframes
    coefs(ii)={F.w(ii)*frana(F.frames{ii},insig)};
end;
outsig=cell2mat(coefs);

