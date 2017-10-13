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


algmpstruct = enuminfo.ltfat_dgtrealmp_alg;
statusenum = enuminfo.ltfat_dgtrealmp_status;

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



for ii= [57]

filename = sprintf('~/Desktop/SQAM/%02d.wav',ii);
disp(filename);
    
[f,fs] = wavload(filename);
[f,fs] = gspi;

%f = postpad(f,fs);
f = cast(f,flags.complexity);

[Ls,W] = size(f(:,1));
Ls = min([Ls,10*fs]);
f = postpad(f(:,1),Ls);
%f = [zeros(100,1);f(6*fs+1:7*fs)];
%f = [zeros(1000,1);f;zeros(1000,1);];
%f = pconv(f,fir2long(firwin('hann',10),numel(f)));
%f = postpad(f,Ls+fs);
%f = [zeros(fs,1);postpad(f(fs:end,1),Ls);zeros(fs,1)];
W = 1;
f = f(:,1);
Ls = numel(f);
a  = [ 512,  64,  64,  64,  64];
M  = [2048, 1024,8192,8192,8192];
gl = [2048, 512, 512,1024,2048];
M2 = floor(M/2) + 1;
P = [1];
Psize = numel(P);
L = dgtlength(Ls,max(a),max(M));
f = postpad(f,L);
cphaseconv = phaseconv.LTFAT_TIMEINV;
mphaseconv = 'timeinv';
if cphaseconv == phaseconv.LTFAT_FREQINV
    mphaseconv = 'freqinv';
end
%f(:) = linspace(0,1,L); 
%f(:) = 1;
    
N = L./a;
%g = randn(gl,1,flags.complexity);

gCell = cell(Psize,1);
for p=1:Psize
    gCell{p} = cast(firwin('blackman',gl(P(p)),'2'),flags.complexity);
end
g = cell2mat(gCell);
 
 
gPtr = libpointer(dataPtr,g);
glPtr = libpointer('int64Ptr',cast(gl(P),'int64'));
aPtr  = libpointer('int64Ptr',cast(a(P),'int64'));
MPtr  = libpointer('int64Ptr',cast(M(P),'int64'));


fPtr = libpointer(dataPtr,f);

fout = randn(L,W,flags.complexity);
foutPtr = libpointer(dataPtr,fout);

sizeaccum = 0;
for p=1:Psize
    sizeaccum = sizeaccum + M2(P(p))*N(P(p))*W;
end

cout = complex2interleaved(...
cast(zeros(sizeaccum,1)+...
         1i*zeros(sizeaccum,1),flags.complexity));
coutPtr = libpointer(dataPtr,cout);


%ctrue = dgt(f,g(1:gl(1)),a(1),M(1));
atoms = 0.8*L;
%atoms = 16650;

params = calllib('libltfat','ltfat_dgtrealmp_params_allocdef');
calllib('libltfat','ltfat_dgtrealmp_setpar_maxatoms',params,atoms);
%calllib('libltfat','ltfat_dgtrealmp_setpar_maxit',params,atoms);
calllib('libltfat','ltfat_dgtrealmp_setpar_errtoldb',params,-40);
calllib('libltfat','ltfat_dgtrealmp_setpar_kernrelthr',params,1e-4);
calllib('libltfat','ltfat_dgtrealmp_setpar_phaseconv',params,cphaseconv);
%calllib('libltfat','ltfat_dgtrealmp_setpar_alg',params,algmpstruct.ltfat_dgtrealmp_alg_LocCyclicMP);
calllib('libltfat','ltfat_dgtrealmp_setpar_iterstep',params,1e6);
calllib('libltfat','ltfat_dgtrealmp_setpar_cycles',params,1);


plan = libpointer();
funname = makelibraryname('dgtrealmp_init_compact',flags.complexity,0);
statusInit = calllib('libltfat',funname,gPtr,glPtr,...
    L,Psize,aPtr,MPtr,params,plan);

calllib('libltfat','ltfat_dgtrealmp_params_free',params);

tic
funname = makelibraryname('dgtrealmp_reset',flags.complexity,0);
statusReset = calllib('libltfat',funname,plan,fPtr);
t1 = toc;

cres1 = complex2interleaved(...
cast(randn(sizeaccum,1)+...
         1i*randn(sizeaccum,1),flags.complexity));
cresPtr = libpointer(dataPtr,cres1);

funname = makelibraryname('dgtrealmp_getresidualcoef_compact',flags.complexity,0);
calllib('libltfat',funname,plan,cresPtr);
cres2 = reshape(interleaved2complex(cresPtr.value),M2(P(1)),N(P(1)));
figure(2); plotdgtreal(cres2,1,100,'clim',[-90,10]);




funname = makelibraryname('dgtrealmp_getresidualcoef_compact',flags.complexity,0);
calllib('libltfat',funname,plan,cresPtr);
cres1 = reshape(interleaved2complex(cresPtr.value),M2(P(1)),N(P(1)));
%figure(2); plotdgtreal(cres1,1,100,'clim',[-90,10]);ylim([0,0.005]);




tic
 %funname = makelibraryname('dgtrealmp_execute_compact',flags.complexity,0);
 %statusExecute = calllib('libltfat',funname,plan,fPtr,coutPtr,foutPtr);
 funname = makelibraryname('dgtrealmp_execute_niters',flags.complexity,0);
 statusExecute = calllib('libltfat',funname,plan,atoms,coutPtr);
t2 =toc;

%%%%%%%%%%%%%%
errdb = libpointer('doublePtr',[1]);
funname = makelibraryname('dgtrealmp_get_errdb',flags.complexity,0);
calllib('libltfat',funname,plan,errdb);
err1 = errdb.value;
%%%%%%%%%%%%%%

cout2 = interleaved2complex(coutPtr.value);
figure(5); plotdgtreal(reshape(cout2,M2(P(1)),N(P(1))),1,100,'clim',[-90,10]);

% funname = makelibraryname('dgtrealmp_revert',flags.complexity,0);
% calllib('libltfat',funname,plan,coutPtr);

%%%%%%%%%%%%%%
errdb = libpointer('doublePtr',[1]);
funname = makelibraryname('dgtrealmp_get_errdb',flags.complexity,0);
calllib('libltfat',funname,plan,errdb);
err2 = errdb.value;
%%%%%%%%%%%%%%

cres2 = complex2interleaved(...
cast(randn(sizeaccum,1)+...
         1i*randn(sizeaccum,1),flags.complexity));
cresPtr = libpointer(dataPtr,cres2);
 
funname = makelibraryname('dgtrealmp_getresidualcoef_compact',flags.complexity,0);
calllib('libltfat',funname,plan,cresPtr);
cres2 = reshape(interleaved2complex(cresPtr.value),M2(P(1)),N(P(1)));
figure(3); plotdgtreal(cres2,1,100,'clim',[-90,10]);
%figure(4); plotdgtreal(cres1-cres2,1,100,'dynrange',90);


fprintf('Init %.3f, execute %.3f, both %.3f seconds, status %s.\n',t1,t2,t1+t2,dgtrealmpstring(statusExecute));



clear coutPtr cout
atoms = numel(find(abs(cout2(:))));

coutCell = cell(Psize,1);
sizeaccum = 0;
fout(:) = 0;
for p=1:Psize
    coutCell{p} = cout2(sizeaccum +1: sizeaccum + M2(P(p))*N(P(p))*W);
    coutCell{p} = reshape(coutCell{p},M2(P(p)),N(P(p)),W);
    sizeaccum = sizeaccum + M2(P(p))*N(P(p))*W;
    fout = fout + idgtreal(coutCell{p},gCell{P(p)},a(P(p)),M(P(p)),mphaseconv);
end

clear cout2
figure(1);plot((0:L-1)/fs,[fPtr.value, fout]);

errdb = 20*log10(norm(fPtr.value -fout)/norm(fPtr.value));

 fprintf('%i atoms, sparsity %.3f, %i atoms/s, Err: True: %.2f dB, En: %.2f dB, Rev: %.2f dB\n',atoms,atoms/L,atoms/t2, errdb,err1,err2);

shg; 

funname = makelibraryname('dgtrealmp_done',flags.complexity,0);
statusDone = calllib('libltfat',funname,plan);


  
    %[test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    %fprintf(['DGTREAL FREQINV WP auto %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],dirstr,L,W,a,M,flags.complexity,ltfatstatusstring(statusExecute),fail);
   
   
end
end


function sstring=dgtrealmpstring(status)

[~,~,enuminfo]=libltfatprotofile;

map = structfun(@(a) a==status ,enuminfo.ltfat_dgtrealmp_status);
names = fieldnames(enuminfo.ltfat_dgtrealmp_status);
sstring = names{map};


