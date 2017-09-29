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



for ii= 1:1

filename = sprintf('~/Desktop/SQAM/%02d.wav',ii);
disp(filename);
    
[f,fs] = wavload(filename);
[f,fs] = gspi;

%f = postpad(f,fs);
f = cast(f,flags.complexity);
%f = [zeros(10000,1);f];
[Ls,W] = size(f(:,1));
Ls = min([Ls,10*fs]);
W = 1;
f = f(:,1);
a  = [  512,  64,  64,  64,  64];
M  = [2048, 4*1024,8192,8192,8192];
M2 = floor(M/2) + 1;
gl = [2048, 256, 512,1024,2048];
P = [1];
Psize = numel(P);
L = dgtlength(Ls,max(a),max(M));
f = postpad(f,L);
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
cast(randn(sizeaccum,1)+...
         1i*randn(sizeaccum,1),flags.complexity));
coutPtr = libpointer(dataPtr,cout);

%ctrue = dgt(f,g(1:gl(1)),a(1),M(1));
atoms = 0.8*L;
atoms = 19000;

params = calllib('libltfat','ltfat_dgtrealmp_params_allocdef');
calllib('libltfat','ltfat_dgtrealmp_setpar_maxatoms',params,atoms);
calllib('libltfat','ltfat_dgtrealmp_setpar_errtoldb',params,-40);
calllib('libltfat','ltfat_dgtrealmp_setpar_kernrelthr',params,1e-4);
calllib('libltfat','ltfat_dgtrealmp_setpar_phaseconv',params,phaseconv.LTFAT_TIMEINV);
calllib('libltfat','ltfat_dgtrealmp_setpar_alg',params,algmpstruct.ltfat_dgtrealmp_alg_locomp);
calllib('libltfat','ltfat_dgtrealmp_setpar_iterstep',params,1e6);

plan = libpointer();
funname = makelibraryname('dgtrealmp_init_compact',flags.complexity,0);
statusInit = calllib('libltfat',funname,gPtr,glPtr,...
    L,Psize,aPtr,MPtr,params,plan);

tic
funname = makelibraryname('dgtrealmp_reset',flags.complexity,0);
statusReset = calllib('libltfat',funname,plan,fPtr);
t1 = toc;

tic
 funname = makelibraryname('dgtrealmp_execute_compact',flags.complexity,0);
 statusExecute = calllib('libltfat',funname,plan,fPtr,coutPtr,foutPtr);

t2 =toc;

fprintf('Init %.3f, execute %.3f, both %.3f seconds, status %d.\n',t1,t2,t1+t2,statusExecute);


cout2 = interleaved2complex(coutPtr.value);
clear coutPtr cout
atoms = numel(find(abs(cout2(:))));

coutCell = cell(Psize,1);
sizeaccum = 0;
for p=1:Psize
    coutCell{p} = cout2(sizeaccum +1: sizeaccum + M2(P(p))*N(P(p))*W);
    coutCell{p} = reshape(coutCell{p},M2(P(p)),N(P(p)),W);
    sizeaccum = sizeaccum + M2(P(p))*N(P(p))*W;
end
clear cout2
plot((0:L-1)/fs,[fPtr.value, foutPtr.value]);

errdb = 20*log10(norm(fPtr.value -foutPtr.value)/norm(fPtr.value));

 fprintf('%i atoms, sparsity %.3f, %i atoms/s, err: %.2f dB\n',atoms,atoms/L,atoms/t2, errdb);

shg; 

funname = makelibraryname('dgtrealmp_done',flags.complexity,0);
statusDone = calllib('libltfat',funname,plan);


  
    %[test_failed,fail]=ltfatdiditfail(res+statusInit,test_failed);
    %fprintf(['DGTREAL FREQINV WP auto %s L:%3i, W:%3i, a:%3i, M:%3i %s %s %s\n'],dirstr,L,W,a,M,flags.complexity,ltfatstatusstring(statusExecute),fail);
   
   
end
end



