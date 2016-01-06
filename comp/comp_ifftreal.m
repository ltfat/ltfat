function f=comp_ifftreal(c,N)
%COMP_IFFTREAL  Compute an IFFTREAL
  
% Force IFFT along dimension 1, since we have permuted the dimensions
% manually
if rem(N,2)==0
  f=[c;...
     flipud(conj(c(2:end-1,:)))];
else
  f=[c;...
     flipud(conj(c(2:end,:)))];
end;

f=real(ifft(f,N,1));


