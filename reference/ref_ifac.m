function f=ref_ifac(ff,W,c,d,p,q,permutation);
%REF_IFAC  Reference inverse factorization.
%

% Output array
f=zeros(c*d*p*q,W);

work=zeros(p*q*d,W);
if d>1
   
  for ko=0:c-1
    
    work(:) = ifft(reshape(ff(:,1+ko*d:(ko+1)*d),q*W*p,d).');

    % Permute again.
    f(permutation+ko,:)=work;

  end;

else

  for ko=0:c-1
    
    % Multiply
    work(:)=ff(:,ko+1);
    
    % Permute again.
    f(permutation+ko,:)=work;          

  end;
  
  
end;


