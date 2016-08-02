function test_failed = test_libltfat_phaselock(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr = [350 350 10  9   1];
Warr = [1     3  1  3   1];
aarr = [10   10 10  9   1];
Marr = [35   35  2  3   1];


for idx = 1:numel(Larr)
    L = Larr(idx);
    W = Warr(idx);
    a = aarr(idx);
    M = Marr(idx);
    
    N = L/a;

    z = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
    zi = complex2interleaved(z);
    zout = randn(size(zi),flags.complexity);

    ziPtr = libpointer(dataPtr,zi);
    zoutPtr = libpointer(dataPtr,zout);

    trueres = phaselock(z,a);

    funname = makelibraryname('dgt_phaselock',flags.complexity,1);
    status = calllib('libltfat',funname,ziPtr,L,W,a,M,zoutPtr);

    res = norm(reshape(trueres,M,N*W) - interleaved2complex(zoutPtr.Value),'fro');

    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['PHASELOCK OP L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,W,a,M,flags.complexity,ltfatstatusstring(status),fail);

    status = calllib('libltfat',funname,ziPtr,L,W,a,M,ziPtr);

    res = norm(reshape(trueres,M,N*W) - interleaved2complex(ziPtr.Value),'fro');

    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['PHASELOCK IP L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,W,a,M,flags.complexity,ltfatstatusstring(status),fail);
    
    
    M2 = floor(M/2) + 1;
    
    z = cast(randn(M2,N,W)+1i*randn(M2,N,W),flags.complexity);
    zi = complex2interleaved(z);
    zout = randn(size(zi),flags.complexity);

    ziPtr = libpointer(dataPtr,zi);
    zoutPtr = libpointer(dataPtr,zout);

    trueres = phaselockreal(z,a,M);

    funname = makelibraryname('dgtreal_phaselock',flags.complexity,1);
    status = calllib('libltfat',funname,ziPtr,L,W,a,M,zoutPtr);

    res = norm(reshape(trueres,M2,N*W) - interleaved2complex(zoutPtr.Value),'fro');

    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['PHASELOCKREAL OP L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,W,a,M,flags.complexity,ltfatstatusstring(status),fail);

    status = calllib('libltfat',funname,ziPtr,L,W,a,M,ziPtr);

    res = norm(reshape(trueres,M2,N*W) - interleaved2complex(ziPtr.Value),'fro');

    [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
    fprintf(['PHASELOCKREAL IP L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],L,W,a,M,flags.complexity,ltfatstatusstring(status),fail);   
 
    
end


