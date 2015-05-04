s=ltfattext;

figure(1);
sig_iter = isgramreal(s,'gauss',8,800);
C=sgram(sig_iter,'dynrange',100);
C=C-min(C(:));
C=C/max(C(:))*255;
C=flipud(C);
imwrite(C,jet(256),'frontpage.png','png','Transparency',0);
