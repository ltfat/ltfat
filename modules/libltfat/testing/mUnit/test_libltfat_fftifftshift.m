function test_failed = test_libltfat_fftifftshift(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr = [1,9,11,110,111];

for L = Larr

        z = cast((1:L)',flags.complexity);
        zi = complex2interleaved(fft(z));
        zout = randn(size(zi),flags.complexity);

        ziPtr = libpointer(dataPtr,zi);
        zoutPtr = libpointer(dataPtr,zout);

        trueres = ifftshift(z);

        funname = makelibraryname('fftifftshift',flags.complexity,1);
        status = calllib('libltfat',funname,ziPtr,L,zoutPtr);

        res = norm(trueres - real(ifft(interleaved2complex(zoutPtr.Value))));

        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTIFFTSHIFT OP L:%3i, %s %s %s\n'],L,flags.complexity,ltfatstatusstring(status),fail);

        status = calllib('libltfat', funname, ziPtr,L,ziPtr);

        res = norm(trueres - real(ifft(interleaved2complex(ziPtr.Value))));

        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTIFFTSHIFT IP L:%3i, %s %s %s\n'],L,flags.complexity,ltfatstatusstring(status),fail);
        
        z = cast((1:L)',flags.complexity);
        zi = complex2interleaved(fftreal(z));
        zout = randn(size(zi),flags.complexity);

        ziPtr = libpointer(dataPtr,zi);
        zoutPtr = libpointer(dataPtr,zout);

        trueres = ifftshift(z);

        funname = makelibraryname('fftrealifftshift',flags.complexity,1);
        status = calllib('libltfat',funname,ziPtr,L,zoutPtr);

        res = norm(trueres - ifftreal(interleaved2complex(zoutPtr.Value),L));

        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTREALIFFTSHIFT OP L:%3i, %s %s %s\n'],L,flags.complexity,ltfatstatusstring(status),fail);

        status = calllib('libltfat', funname, ziPtr,L,ziPtr);

        res = norm(trueres - ifftreal(interleaved2complex(ziPtr.Value),L));

        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTREALIFFTSHIFT IP L:%3i, %s %s %s\n'],L,flags.complexity,ltfatstatusstring(status),fail);

end


