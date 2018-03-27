function test_failed = test_libltfat_dgtreal2dgt(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr  = [350 350   9   1];
glarr = [ 20  10   9   1];
aarr  = [ 10  10   9   1];
Marr  = [ 35  35   3   1];
Warr  = [  1   3   3   1];

for idx = 1:numel(Larr)
    L = Larr(idx);
    W = Warr(idx);
    a = aarr(idx);
    M = Marr(idx);
    M2 = floor(M/2) + 1;
    gl = glarr(idx);
    
    N = L/a;

    g = randn(gl,1,flags.complexity);
    c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
    cout = complex2interleaved(c);
    f = randn(L,W,flags.complexity);

    cin = complex2interleaved(dgtreal(f,g,a,M));
    truecfull = dgt(f,g,a,M);
    
    cinPtr = libpointer(dataPtr,cin);
    coutPtr = libpointer(dataPtr,cout);
    N = size(truecfull,2);

    funname = makelibraryname('dgtreal2dgt',flags.complexity,1);
    status = calllib('libltfat',funname, cinPtr, M,N*W, coutPtr);

    res = norm(reshape(truecfull,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['dgtreal2dgt OP    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,gl,W,a,M,flags.complexity,ltfatstatusstring(status),fail);
  
    cintmp = dgtreal(f,g,a,M);
    cin = complex2interleaved(postpad(cintmp(:),M*N*W,1));
    cinPtr = libpointer(dataPtr,cin);
    
    funname = makelibraryname('dgtreal2dgt',flags.complexity,1);
    status = calllib('libltfat',funname, cinPtr, M,N*W, cinPtr);
    
    res = norm(reshape(truecfull,M,N*W) - reshape(interleaved2complex(cinPtr.Value),M,N*W),'fro');
    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['dgtreal2dgt IP    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,gl,W,a,M,flags.complexity,ltfatstatusstring(status),fail);  
    
    
    %%
    c = cast(randn(M2,N,W)+1i*randn(M2,N,W),flags.complexity);
    cout = complex2interleaved(c);
    
    cin = complex2interleaved(dgt(f,g,a,M));
    truecreal = dgtreal(f,g,a,M);
    
    cinPtr = libpointer(dataPtr,cin);
    coutPtr = libpointer(dataPtr,cout);
    N = size(truecreal,2);

    funname = makelibraryname('dgt2dgtreal',flags.complexity,1);
    status = calllib('libltfat',funname, cinPtr, M,N*W, coutPtr);

    res = norm(reshape(truecreal,M2,N*W) - interleaved2complex(coutPtr.Value),'fro');
    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['dgt2dgtreal OP    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,gl,W,a,M,flags.complexity,ltfatstatusstring(status),fail);
  
    cintmp = dgt(f,g,a,M);
    cin = complex2interleaved(cintmp);
    cinPtr = libpointer(dataPtr,cin);
    
    funname = makelibraryname('dgt2dgtreal',flags.complexity,1);
    status = calllib('libltfat',funname, cinPtr, M,N*W, cinPtr);
    
    cinvalue = interleaved2complex(cinPtr.Value);
    res = norm(reshape(truecreal,M2,N*W) - reshape(postpad(cinvalue(:),M2*N*W),M2,N*W),'fro');
    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['dgtreal2dgt IP    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,gl,W,a,M,flags.complexity,ltfatstatusstring(status),fail);
    
end


