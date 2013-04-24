function test_failed=test_zak

disp(' ===============  TEST_ZAK ================');

Lr=[1,12,12,12,12,12,15,15,15];
Kr=[1, 1, 2, 3, 4, 6, 3, 5,15];

ref_failed=0;
inv_failed=0;

W=1;

for ii=1:length(Lr)
  
  L=Lr(ii);
  K=Kr(ii);

  f=tester_rand(L,W)+i*rand(L,W)-.5-i*.5;
      
  ccref=ref_zak(f,K);
  cc=zak(f,K);
  r=izak(cc);

  res=ccref-cc;
  
  nres=norm(res(:));
  ninv=norm(f-r);

  s=sprintf('REF L:%3i K:%3i %0.5g',L,K,nres);
  disp(s)

  [ref_failed,fail]=ltfatdiditfail(nres,ref_failed);
  % if nres>10e-10
    % disp('FAILED');
    % ref_failed=ref_failed+1;
  % end;

  s=sprintf('INV L:%3i K:%3i %0.5g',L,K,ninv);
  disp(s)

  [inv_failed,fail]=ltfatdiditfail(ninv,inv_failed);
  % if nres>10e-10
    % disp('FAILED');
    % inv_failed=inv_failed+1;
  % end;

end;

test_failed=ref_failed+inv_failed;

