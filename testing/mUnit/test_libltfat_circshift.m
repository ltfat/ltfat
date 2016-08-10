function test_failed = test_libltfat_circshift(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr = [1,9,11,110];
shiftarr = [-10, 100, 0, -1, 1, 3];

for do_complex = 0:1
    complexstring = '';
    if do_complex, complexstring = 'complex'; end
    for L = Larr
        for shift = shiftarr
            
            if do_complex
                z = cast((1:L)'+1i*(L:-1:1)',flags.complexity);
                zi = complex2interleaved(z);
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            else
                z = cast((1:L)',flags.complexity);
                zi = z;
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            end
            
            trueres = circshift(z,shift);
            
            
            status = calllib('libltfat',makelibraryname('circshift',flags.complexity,do_complex),...
                ziPtr,L,shift,zoutPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(zoutPtr.Value));
            else
                res = norm(trueres - zoutPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed,0);
            fprintf(['CIRCSHIFT OP L:%3i, shift:%3i %s %s %s %s\n'],L,shift,flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
            status = calllib('libltfat',makelibraryname('circshift',flags.complexity,do_complex),...
                ziPtr,L,shift,ziPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(ziPtr.Value));
            else
                res = norm(trueres - ziPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed,0);
            fprintf(['CIRCSHIFT IP L:%3i, shift:%3i %s %s %s %s\n'],L,shift,flags.complexity,complexstring,ltfatstatusstring(status),fail);
        end
    end
end

