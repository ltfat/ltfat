% This script test the quality of the Hermite implementation.

% There seem to be some numerical problems for ii>66

L=200;

rsym=zeros(L,1);
  
if 0
  rnan=zeros(L,1);
  
  for ii=0:L-1
    g=pherm(L,ii);
    rnan(ii+1)=any(isnan(g));
    if mod(ii,2)==0
      rsym(ii+1)=norm(imag(fft(g)));
    else
      rsym(ii+1)=norm(real(fft(g)));
    end;
  end;
  
else

  [V,D]=pherm3(L);

  for ii=0:floor(L/2-1)
    if mod(ii,2)==0
      rsym(ii+1)=norm(imag(fft(V(:,ii+1))));
    else
      rsym(ii+1)=norm(real(fft(V(:,ii+1))));
    end;
  end;

  for ii=floor(L/2-1)+1:L-1
    if mod(ii,2)==1
      rsym(ii+1)=norm(imag(fft(V(:,ii+1))));
    else
      rsym(ii+1)=norm(real(fft(V(:,ii+1))));
    end;
  end;

  
end;

semilogy(rsym);
