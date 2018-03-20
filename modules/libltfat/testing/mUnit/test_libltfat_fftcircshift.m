function test_failed = test_libltfat_fftcircshift(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr = [1,9,11,110,111];
shiftarr = [-10, 100, 0, -1, 1, 3, 52];

for L = Larr
    for shift = shiftarr

        z = cast((1:L)' + 1i*(L:-1:1)',flags.complexity);
        zi = complex2interleaved(fft(z));
        zout = randn(size(zi),flags.complexity);

        ziPtr = libpointer(dataPtr,zi);
        zoutPtr = libpointer(dataPtr,zout);

        trueres = circshift(z,shift);

        status = calllib('libltfat',makelibraryname('fftcircshift',flags.complexity,1),...
            ziPtr,L,shift,zoutPtr);

        res = norm(trueres - ifft(interleaved2complex(zoutPtr.Value)));


        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTCIRCSHIFT OP L:%3i, shift:%3i %s %s %s\n'],L,shift,flags.complexity,ltfatstatusstring(status),fail);

        status = calllib('libltfat',makelibraryname('fftcircshift',flags.complexity,1),...
            ziPtr,L,shift,ziPtr);

        res = norm(trueres - ifft(interleaved2complex(ziPtr.Value)));


        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTCIRCSHIFT IP L:%3i, shift:%3i %s %s %s\n'],L,shift,flags.complexity,ltfatstatusstring(status),fail);
        
        
        z = cast((1:L)',flags.complexity);
        zi = complex2interleaved(fftreal(z));
        zout = randn(size(zi),flags.complexity);

        ziPtr = libpointer(dataPtr,zi);
        zoutPtr = libpointer(dataPtr,zout);

        trueres = circshift(z,shift);

        status = calllib('libltfat',makelibraryname('fftrealcircshift',flags.complexity,1),...
            ziPtr,L,shift,zoutPtr);

        res = norm(trueres - real(ifftreal(interleaved2complex(zoutPtr.Value),L)));

        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTREALCIRCSHIFT OP L:%3i, shift:%3i %s %s %s\n'],L,shift,flags.complexity,ltfatstatusstring(status),fail);

        status = calllib('libltfat',makelibraryname('fftrealcircshift',flags.complexity,1),...
            ziPtr,L,shift,ziPtr);

        res = norm(trueres - real(ifftreal(interleaved2complex(ziPtr.Value),L)));


        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['FFTREALCIRCSHIFT IP L:%3i, shift:%3i %s %s %s\n'],L,shift,flags.complexity,ltfatstatusstring(status),fail);
    end
end


