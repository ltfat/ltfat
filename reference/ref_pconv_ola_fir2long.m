function h=ref_pconv_ola_fir2long(f,g,Lb)
%
%  This function implements periodic convolution using overlap-add. The
%  window g is supposed to be extended by fir2iir.
  
L=length(f);
Lg=length(g);
  
% Number of blocks
Nb=L/Lb;

% Length of extended block and padded g
Lext=Lb+Lg;
gpad=fir2long(g,Lext);
b2=Lg/2;

h=zeros(L,1);
for ii=0:Nb-1
  block=pconv(postpad(f(ii*Lb+1:(ii+1)*Lb),Lext),gpad);
  h(ii*Lb+1:(ii+1)*Lb)+=block(1:Lb);  % Large block

  s_ii=mod(ii+1,Nb);
  h(s_ii*Lb+1:s_ii*Lb+b2)+=block(Lb+1:Lb+b2);  % Small block +

  s_ii=mod(ii-1,Nb)+1;
  h(s_ii*Lb-b2+1:s_ii*Lb)+=block(Lb+b2+1:Lext);  % Small block -

end;

