function gf=filterbankresponse(g,a,L)
  
  gf=zeros(L,1);
  M=numel(g);
  
  for m=1:M
    gf=gf+abs(fft(middlepad(g{m},L))).^2;
  end;
  
  