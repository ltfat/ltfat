function test_failed = test_libltfat_dgtwrapper(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

[~,~,enuminfo]=libltfatprotofile;
phaseconv = enuminfo.ltfat_phaseconvention;
hintstruct = enuminfo.ltfat_dgt_hint;

fftwflags = struct('FFTW_MEASURE',0,'FFTW_ESTIMATE',64,'FFTW_PATIENT',32,'FFTW_DESTROY_INPUT',1,...
    'FFTW_UNALIGNED',2,'FFTW_EXHAUSTIVE',8,'FFTW_PRESERVE_INPUT',16);

Larr  = [350 360   9   2];
glarr = [ 20  10   9   1];
aarr  = [ 10  10   3   1];
Marr  = [ 35  36   3   2];
Warr  = [  1   3   3   1];

for do_complex = [1,0]
    complexstring = '';
    if do_complex, complexstring = 'complex'; end

for anasyn = {'ana','syn'}
dirstr = anasyn{1};
for hint = fieldnames(hintstruct).'
    hint = hint{1};
for idx = 1:numel(Larr)
    L = Larr(idx);
    W = Warr(idx);
    a = aarr(idx);
    M = Marr(idx);
    gl = glarr(idx);
    
    N = L/a;
    g = randn(gl,1,flags.complexity);
    if do_complex
        g = g + 1i*randn(gl,1,flags.complexity);
        gPtr = libpointer(dataPtr,complex2interleaved(g)); 
    else
       gPtr = libpointer(dataPtr,g); 
    end
    
    f = cast(randn(L,W)+1i*randn(L,W),flags.complexity);
    fPtr = libpointer(dataPtr,complex2interleaved(f));
    if ~do_complex, f = reshape(postpad(fPtr.value(:),L*W),L,W); end

    c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
    cout = complex2interleaved(c);
    coutPtr = libpointer(dataPtr,cout);
    

    plan = libpointer();
    funname = makelibraryname('dgt_init',flags.complexity,do_complex);
    statusInit = calllib('libltfat',funname,gPtr,gl,L,W,a,M,fPtr,coutPtr,libpointer(),plan);
    
    if strcmp(dirstr,'ana') 
        funname = makelibraryname('dgt_execute_ana',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan);
    
        truec = dgt(f,g,a,M);
        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
    elseif strcmp(dirstr,'syn')
        funname = makelibraryname('dgt_execute_syn',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan);
    
        truef = idgt(c,{'dual',g},a);
        res = norm(truef - interleaved2complex(fPtr.Value),'fro');      
    end
    
    funname = makelibraryname('dgt_done',flags.complexity,do_complex);
    statusDone = calllib('libltfat',funname,plan);
    
    [test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    fprintf(['dgt FREQINV WP auto %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],dirstr,L,W,a,M,flags.complexity,complexstring,ltfatstatusstring(statusExecute),fail);

    
    %%%%%%
    f = cast(randn(L,W)+1i*randn(L,W),flags.complexity);
    fPtr = libpointer(dataPtr,complex2interleaved(f));
    if ~do_complex, f = reshape(postpad(fPtr.value(:),L*W),L,W); end
       
    c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
    cout = complex2interleaved(c);
    coutPtr = libpointer(dataPtr,cout);
    
    params = calllib('libltfat','ltfat_dgt_params_allocdef');
    calllib('libltfat','ltfat_dgt_setpar_phaseconv',params,phaseconv.LTFAT_TIMEINV);
    calllib('libltfat','ltfat_dgt_setpar_hint',params,hintstruct.(hint));
    calllib('libltfat','ltfat_dgt_setpar_fftwflags',params,fftwflags.FFTW_MEASURE);
    
    plan = libpointer();
    funname = makelibraryname('dgt_init',flags.complexity,do_complex);
    statusInit = calllib('libltfat',funname,gPtr,gl,L,W,a,M,fPtr,coutPtr,params,plan);
    
    if strcmp(dirstr,'ana') 
        funname = makelibraryname('dgt_execute_ana',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan);
    
        truec = dgt(f,g,a,M,'timeinv');
        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
    elseif strcmp(dirstr,'syn')
        funname = makelibraryname('dgt_execute_syn',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan);
    
        truef = idgt(c,{'dual',g},a,'timeinv');
        res = norm(truef - interleaved2complex(fPtr.Value),'fro');      
    end
    
    funname = makelibraryname('dgt_done',flags.complexity,do_complex);
    statusDone = calllib('libltfat',funname,plan);
    
    calllib('libltfat','ltfat_dgt_params_free',params);
    
    [test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    fprintf(['dgt TIMEINV WP %-4s %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],hint(15:end),dirstr,L,W,a,M,flags.complexity,complexstring,ltfatstatusstring(statusExecute),fail);
    
    %%%%%%%
    f = cast(randn(L,W)+1i*randn(L,W),flags.complexity);
    fPtr = libpointer(dataPtr,complex2interleaved(f));
    if ~do_complex, f = reshape(postpad(fPtr.value(:),L*W),L,W); end
   
    c = cast(randn(M,N,W)+1i*randn(M,N,W),flags.complexity);
    cout = complex2interleaved(c);
    coutPtr = libpointer(dataPtr,cout);
    
    params = calllib('libltfat','ltfat_dgt_params_allocdef');
    calllib('libltfat','ltfat_dgt_setpar_phaseconv',params,phaseconv.LTFAT_TIMEINV);
    calllib('libltfat','ltfat_dgt_setpar_hint',params,hintstruct.(hint));
    calllib('libltfat','ltfat_dgt_setpar_fftwflags',params,fftwflags.FFTW_ESTIMATE);
    
    plan = libpointer();
    funname = makelibraryname('dgt_init',flags.complexity,do_complex);
    statusInit = calllib('libltfat',funname,gPtr,gl,L,W,a,M,libpointer(),libpointer(),params,plan);
    
    if strcmp(dirstr,'ana') 
        funname = makelibraryname('dgt_execute_ana_newarray',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan,fPtr,coutPtr);
    
        truec = dgt(f,g,a,M,'timeinv');
        res = norm(reshape(truec,M,N*W) - interleaved2complex(coutPtr.Value),'fro');
    elseif strcmp(dirstr,'syn')
        funname = makelibraryname('dgt_execute_syn_newarray',flags.complexity,do_complex);
        statusExecute = calllib('libltfat',funname,plan,coutPtr,fPtr);
    
        truef = idgt(c,{'dual',g},a,'timeinv');
        res = norm(truef - interleaved2complex(fPtr.Value),'fro');      
    end
    
    funname = makelibraryname('dgt_done',flags.complexity,do_complex);
    statusDone = calllib('libltfat',funname,plan);
    
    calllib('libltfat','ltfat_dgt_params_free',params);
    
    [test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    fprintf(['dgt TIMEINV WP %-4s %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s %s\n'],hint(15:end),dirstr,L,W,a,M,flags.complexity,complexstring,ltfatstatusstring(statusExecute),fail);    
    
    
end
end
end
end


