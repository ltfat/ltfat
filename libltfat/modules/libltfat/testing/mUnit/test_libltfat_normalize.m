function test_failed = test_libltfat_normalize(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

[~,~,enuminfo]=libltfatprotofile;
Cenumnorms = enuminfo.ltfat_norm_t;

Larr = [1,9,11,110];
normpairs = {{'null',Cenumnorms.LTFAT_NORM_NULL},...
             {'1', Cenumnorms.LTFAT_NORM_1},...
             {'2', Cenumnorms.LTFAT_NORM_2}};

for do_complex = 0:1
    complexstring = '';
    if do_complex, complexstring = 'complex'; end
    
    funname = makelibraryname('normalize',flags.complexity,do_complex);
    
    for L = Larr
        for npId = 1:numel(normpairs)
            normpair = normpairs{npId};
            
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
            
            trueres = normalize(z,normpair{1});
            
            
            status = calllib('libltfat',funname, ziPtr,L,normpair{2},zoutPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(zoutPtr.Value));
            else
                res = norm(trueres - zoutPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res-status,test_failed);
            fprintf(['NORMALIZE OP L:%3i, norm: %s %s %s %s %s\n'],L,normpair{1},flags.complexity,complexstring,ltfatstatusstring(status),fail);
            
            status = calllib('libltfat',funname, ziPtr,L,normpair{2},ziPtr);
            
            if do_complex
                res = norm(trueres - interleaved2complex(ziPtr.Value));
            else
                res = norm(trueres - ziPtr.Value);
            end
            
            [test_failed,fail]=ltfatdiditfail(res-status,test_failed);
            fprintf(['NORMALIZE IP L:%3i, norm: %s %s %s %s %s\n'],L,normpair{1},flags.complexity,complexstring,ltfatstatusstring(status),fail);
        end
    end
end

