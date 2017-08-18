function test_failed = test_libltfat_dgtrealmp(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];
dataPtrPtr = [flags.complexity, 'PtrPtr'];

[~,~,enuminfo]=libltfatprotofile;
phaseconv = enuminfo.ltfat_phaseconvention;
hintstruct = enuminfo.ltfat_dgtreal_hint;

hintmpstruct = enuminfo.ltfat_dgtrealmp_hint;
algmpstruct = enuminfo.ltfat_dgtrealmp_alg;

fftwflags = struct('FFTW_MEASURE',0,'FFTW_ESTIMATE',64,'FFTW_PATIENT',32,'FFTW_DESTROY_INPUT',1,...
    'FFTW_UNALIGNED',2,'FFTW_EXHAUSTIVE',8,'FFTW_PRESERVE_INPUT',16);

base = 2048;
Larr  = [128* base   9   2];
glarr = [ base  10   9   1];
aarr  = [   base/4   10   3   1];
Marr  = [ base  36   3   2];
Warr  = [  1   3   3   1];

for idx = 1:1%numel(Larr)
    
%     L = Larr(idx);
%     W = Warr(idx);
%     a = aarr(idx);
%     M = Marr(idx);
%     M2 = floor(M/2) + 1;
%     gl = glarr(idx);
% f = randn(L,W,flags.complexity);

[f,fs] = gspi;
%f = postpad(f,fs);
[Ls,W] = size(f);
a = 64;
M = 2048
M2 = floor(M/2) + 1;
gl = 2048;
L = dgtlength(Ls,a,M);
f = postpad(f,L);


    
    N = L/a;
    %g = randn(gl,1,flags.complexity);
    g = firwin('blackman',gl,'2');
    gPtr = libpointer(dataPtrPtr,g);
    
    
    fPtr = libpointer(dataPtr,f);
    
    fout = randn(L,W,flags.complexity);
    foutPtr = libpointer(dataPtr,fout);

    c = cast(randn(M2,N,W)+1i*randn(M2,N,W),flags.complexity);
    cout = complex2interleaved(c);
    coutPtr = libpointer(dataPtr,cout);
    
    ctrue = dgt(f,g,a,M);
    atoms = 90000;
    
    
    plan = libpointer();
    funname = makelibraryname('dgtrealmp_init',flags.complexity,0);
    statusInit = calllib('libltfat',funname,gPtr,libpointer('int64Ptr',gl),...
        L,1,libpointer('int64Ptr',a),libpointer('int64Ptr',M),libpointer(),plan);
    
    funname = makelibraryname('dgtrealmp_reset',flags.complexity,0);
    statusReset = calllib('libltfat',funname,plan,fPtr);
    
  
    funname = makelibraryname('dgtrealmp_set_iterstep',flags.complexity,0);
    sttatusSet =  calllib('libltfat',funname, plan, 10000);
    
    funname = makelibraryname('dgtrealmp_set_maxatoms',flags.complexity,0);
    sttatusSet =  calllib('libltfat',funname, plan, atoms);        

     tic
     funname = makelibraryname('dgtrealmp_execute',flags.complexity,0);
     statusExecute = calllib('libltfat',funname,plan,fPtr,coutPtr,foutPtr)
     t =toc;
     cout2 = interleaved2complex(coutPtr.value);
     atoms = numel(find(abs(cout2(:))));
     fprintf('%i atoms, %i atoms/s\n',atoms,atoms/t);
     
     
     
     plot([fPtr.value, foutPtr.value])
     
     20*log10(norm(fPtr.value -foutPtr.value)/norm(fPtr.value))
        
    
    funname = makelibraryname('dgtrealmp_done',flags.complexity,0);
    statusDone = calllib('libltfat',funname,plan);
    
    %[test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    %fprintf(['DGTREAL FREQINV WP auto %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],dirstr,L,W,a,M,flags.complexity,ltfatstatusstring(statusExecute),fail);
   
   
end



