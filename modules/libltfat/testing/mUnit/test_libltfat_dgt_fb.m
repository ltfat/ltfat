function test_failed = test_libltfat_dgt_fb(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

[~,~,enuminfo]=libltfatprotofile;
phaseconv = enuminfo.ltfat_phaseconvention;

fftwflags = struct('FFTW_MEASURE',0,'FFTW_ESTIMATE',64,'FFTW_PATIENT',32,'FFTW_DESTROY_INPUT',1,...
    'FFTW_UNALIGNED',2,'FFTW_EXHAUSTIVE',8,'FFTW_PRESERVE_INPUT',16);

Larr  = [350 350   9   1];
glarr = [ 20  10   9   1];
aarr  = [ 10  10   9   1];
Marr  = [ 35  35   3   1];
Warr  = [  1   3   3   1];

for do_complex = 0:1
    complexstring = '';
    if do_complex, complexstring = 'complex'; end

    for idx = 1:numel(Larr)
        L = Larr(idx);
        W = Warr(idx);
        a = aarr(idx);
        M = Marr(idx);
        gl = glarr(idx);

        N = L/a;

        if do_complex
            g = randn(gl,1,flags.complexity) + 1i*randn(gl,1,flags.complexity);
            f = randn(L,W,flags.complexity) + 1i*randn(L,W,flags.complexity);   
            gin = complex2interleaved(g);
            fin = complex2interleaved(f);
        else
            g = randn(gl,1,flags.complexity);
            f = randn(L,W,flags.complexity);
            gin = g; fin = f;
        end
        
        fPtr = libpointer(dataPtr,fin);
        gPtr = libpointer(dataPtr,gin);
        c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
        cout = complex2interleaved(c);
        coutPtr = libpointer(dataPtr,cout);

        truec = dgt(f,g,a,M);

        funname = makelibraryname('dgt_fb',flags.complexity,do_complex);
        status = calllib('libltfat',funname,fPtr,gPtr,L,gl,W,a,M,phaseconv.LTFAT_FREQINV,coutPtr);

        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['DGT FREQINV    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],L,gl,W,a,M,complexstring,flags.complexity,ltfatstatusstring(status),fail);

        truec = dgt(f,g,a,M,'timeinv');
        status = calllib('libltfat',funname,fPtr,gPtr,L,gl,W,a,M,phaseconv.LTFAT_TIMEINV,coutPtr);

        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
        [test_failed,fail]=ltfatdiditfail(res+status,test_failed);
        fprintf(['DGT TIMEINV    L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],L,gl,W,a,M,complexstring,flags.complexity,ltfatstatusstring(status),fail);

        % With plan
        c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
        cout = complex2interleaved(c);
        coutPtr = libpointer(dataPtr,cout);

        plan = libpointer();
        funname = makelibraryname('dgt_fb_init',flags.complexity,do_complex);
        statusInit = calllib('libltfat',funname,gPtr,gl,a,M,phaseconv.LTFAT_FREQINV,fftwflags.FFTW_MEASURE,plan);

        funname = makelibraryname('dgt_fb_execute',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan, fPtr,L,W,coutPtr);

        funname = makelibraryname('dgt_fb_done',flags.complexity,do_complex);
        statusDone = calllib('libltfat',funname,plan);

        truec = dgt(f,g,a,M);
        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
        [test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
        fprintf(['DGT FREQINV WP L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],L,gl,W,a,M,complexstring,flags.complexity,ltfatstatusstring(status),fail);

        %%%%%%
        c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
        cout = complex2interleaved(c);
        coutPtr = libpointer(dataPtr,cout);

        plan = libpointer();
        funname = makelibraryname('dgt_fb_init',flags.complexity,do_complex);
        statusInit = calllib('libltfat',funname,gPtr,gl,a,M,phaseconv.LTFAT_TIMEINV,fftwflags.FFTW_MEASURE,plan);

        funname = makelibraryname('dgt_fb_execute',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan, fPtr,L,W,coutPtr);

        funname = makelibraryname('dgt_fb_done',flags.complexity,do_complex);
        statusDone = calllib('libltfat',funname,plan);

        truec = dgt(f,g,a,M,'timeinv');
        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
        [test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
        fprintf(['DGT TIMEINV WP L:%3i, gl:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],L,gl,W,a,M,complexstring,flags.complexity,ltfatstatusstring(status),fail);
    end
end

