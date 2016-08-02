function test_failed = test_libltfat_gabdual(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

Larr = [12  15  120  120];
aarr = [3   3    10   10];
Marr = [6   5    20   20];

for do_complex = 0:1
    complexstring = '';
    if do_complex, complexstring = 'complex'; end
    for idx = 1:numel(Larr)
        L = Larr(idx);
        a = aarr(idx);
        M = Marr(idx);
        
            if do_complex
                z = cast(randn(L,1)+1i*randn(L,1),flags.complexity);
                zi = complex2interleaved(z);
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            else
                z = cast(randn(L,1),flags.complexity);
                zi = z;
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            end
            
            trueres = gabdual(z,a,M);
            
            funname = makelibraryname('gabdual_long',flags.complexity,do_complex);
            status = calllib('libltfat',funname, ziPtr,L,a,M,zoutPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(zoutPtr.Value));
            else
                res = norm(trueres - zoutPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
            fprintf(['GABDUAL OP L:%3i, %s %s %s %s\n'],L,flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
            status = calllib('libltfat',funname,ziPtr,L,a,M,ziPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(ziPtr.Value));
            else
                res = norm(trueres - ziPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
            fprintf(['GABDUAL IP L:%3i, %s %s %s %s\n'],L,flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
            
            
            if do_complex
                z = cast(randn(M,1)+1i*randn(M,1),flags.complexity);
                zi = complex2interleaved(z);
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            else
                z = cast(randn(M,1),flags.complexity);
                zi = z;
                zout = randn(size(zi),flags.complexity);
                
                ziPtr = libpointer(dataPtr,zi);
                zoutPtr = libpointer(dataPtr,zout);
            end
            
            trueres = gabdual(z,a,M);
            
            funname = makelibraryname('gabdual_painless',flags.complexity,do_complex);
            status = calllib('libltfat',funname, ziPtr,M,a,M,zoutPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(zoutPtr.Value));
            else
                res = norm(trueres - zoutPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
            fprintf(['GABDUAL PAINLESS OP L:%3i, %s %s %s %s\n'],L,flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
            status = calllib('libltfat',funname,ziPtr,M,a,M,ziPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(ziPtr.Value));
            else
                res = norm(trueres - ziPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
            fprintf(['GABDUAL PANLESS IP L:%3i, %s %s %s %s\n'],L,flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
    end
end

