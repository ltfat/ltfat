function test_failed=test_fwt(verbose)
%
%  Check perfect reconstruction of some parameter combinations
%

f = randn(17897,1);
J = 7;
% Examples
if(nargin>0)
  verbose=1;
else
  verbose=0;
end
test_failed = 0;

w = fwtinit({'db',7});
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
err=norm(f-fhat);
if(err>1e-6)
 test_failed = 1; 
 if verbose, stem([f,fhat]); end;
 return;
end

[c] = fwt(f,{'db',7},J);
fhat = ifwt(c,{'db',7},J,length(f));
err=norm(f-fhat);
if(err>1e-6)
 test_failed = 1; 
 if verbose, stem([f,fhat]); end;
 return;
end


w = fwtinit({'algmband',1});
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
err=norm(f-fhat);
if(err>1e-6)
 test_failed = 1; 
 if verbose, stem([f,fhat]); end;
 return;
end
%%
w = fwtinit({'db',8});
[c] = fwt(f,w,J);
fhat = ifwt(c,w,J,length(f));
err=norm(f-fhat);
if(err>1e-6)
 test_failed = 1; 
 if verbose, stem([f,fhat]); end;
 return;
end
%%

[h,g] = wfilt_db(8);
[c] = fwt(f,h,J);
fhat = ifwt(c,g,J,length(f));
err=norm(f-fhat);
if(err>1e-6)
 test_failed = 1; 
 if verbose, stem([f,fhat]); end;
 return;
end
%%
% TODO: put somewhere else
% w = fwtinit({'db',8},'undec');
% [c] = fwt(f,w,J);
% fhat = ifwt(c,w,J,length(f));
% err=norm(f-fhat);
% if(err>1e-6)
%  test_failed = 1; 
%  if verbose, stem([f,fhat]); end;
%  return;
% end
% 
% w = waveletfb({'db',8},'undec','sym');
% [c] = fwt(f,w,J);
% fhat = ifwt(c,w,J,length(f));
% err=norm(f-fhat);
% if(err>1e-6)
%  test_failed = 1; 
%  if verbose, stem([f,fhat]); end;
%  return;
% end
% 
% 
% w = waveletfb({'db',8},'undec','sym');
% [c] = fwt(f,w,J,'per');
% fhat = ifwt(c,w,J,length(f),'per');
% err=norm(f-fhat);
% if(err>1e-6)
%  test_failed = 1; 
%  if verbose, stem([f,fhat]); end;
%  return;
% end
% 
% w = waveletfb({'db',8},'undec','sym');
% [c] = fwt(f,w,J);
% fhat = ifwt(c,w,J,length(f));
% err=norm(f-fhat);
% if(err>1e-6)
%  test_failed = 1; 
%  if verbose, stem([f,fhat]); end;
%  return;
% end
% 
% 
% w = waveletfb({'db',8});
% [c] = fwt(f,w,J,'undec','sym');
% fhat = ifwt(c,w,J,length(f),'undec','sym');
% err=norm(f-fhat);
% if(err>1e-6)
%  test_failed = 1; 
%  if verbose, stem([f,fhat]); end;
%  return;
% end

