function test_failed = test_libltfat_slidgtrealmp(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

definput.flags.complexity={'double','single'};
[flags]=ltfatarghelper({},definput,varargin);
dataPtr = [flags.complexity, 'Ptr'];

[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

glarr =     [ 4*2048];
Warr =      [ 1];

bufLenInit = 100;
bufLenMax = 1000;

for initId = 0:1
    
for ii = 1:numel(glarr)
    gl = glarr(ii);
    W = Warr(ii);
    sliplan = libpointer();
    parbuf = libpointer();
    mpstate = libpointer();
    slimpstate = libpointer();

    
    taperLen = 1024;%round(floor(gl/100)/2)*2;
    zpadLen = 1024;
    procdelay = gl - zpadLen - 1;
%     funname = makelibraryname('slicing_processor_init',flags.complexity,0);
%     calllib('libltfat',funname, gl, taperLen, zpadLen, W,bufLenMax,sliplan);

    funname = makelibraryname('dgtrealmp_parbuf_init',flags.complexity,0);
    calllib('libltfat',funname, parbuf);
    
    funname = makelibraryname('dgtrealmp_parbuf_add_firwin',flags.complexity,0);
    calllib('libltfat',funname, parbuf, LTFAT_FIRWIN.LTFAT_BLACKMAN, 2048,  512, 2048);
    calllib('libltfat',funname, parbuf, LTFAT_FIRWIN.LTFAT_BLACKMAN, 512,  128, 512);
    
    funname = makelibraryname('dgtrealmp_setparbuf_maxit',flags.complexity,0);
    calllib('libltfat',funname, parbuf, gl);
    
        funname = makelibraryname('dgtrealmp_setparbuf_iterstep',flags.complexity,0);
    calllib('libltfat',funname, parbuf, gl);
    
    funname = makelibraryname('dgtrealmp_setparbuf_snrdb',flags.complexity,0);
    calllib('libltfat',funname, parbuf, 40);

%     funname = makelibraryname('dgtrealmp_init',flags.complexity,0);
%     calllib('libltfat',funname, parbuf, gl, mpstate);
    
    funname = makelibraryname('slidgtrealmp_init',flags.complexity,0);
    calllib('libltfat',funname, parbuf, gl, W, bufLenMax, slimpstate);
    
    funname = makelibraryname('slidgtrealmp_getprocdelay',flags.complexity,0);
    procdelay = calllib('libltfat',funname,slimpstate);    
    initstr = 'INIT WIN';


    [bufIn,fs] = gspi;
    bufIn = cast(bufIn,flags.complexity);
    bufIn = bsxfun(@times, repmat(bufIn,1,W), [1, rand(1,W-1,flags.complexity) + 1]);

    bufOut = 1000*ones(size(bufIn),flags.complexity);
    L = size(bufIn,1);
    status = 0;
    startIdx = 1;
    bufLen = bufLenInit;
    while startIdx <= L
        stopIdx = min([startIdx + bufLen - 1,L]);
        slice = startIdx : stopIdx;
        buf = bufIn(slice,:);
        bufInPtr = libpointer(dataPtr,buf);
        bufOutPtr = libpointer(dataPtr,randn(size(buf),flags.complexity));

        % Matlab automatically converts Ptr to PtrPtr
        funname = makelibraryname('slidgtrealmp_execute',flags.complexity,0);
        status = calllib('libltfat',funname,slimpstate,bufInPtr,numel(slice),W,bufOutPtr);
        if status
            break;
        end

        bufOut(slice,:) = bufOutPtr.Value;
        startIdx = stopIdx + 1;
        bufLen = randi(bufLenMax);
    end

    inshift = circshift(bufIn,(procdelay));
    inshift(1:(procdelay),:) = 0;
    plotthat = [bufOut - inshift];
    plotthat(end-(procdelay):end,:) = 0;

    [test_failed,fail]=ltfatdiditfail(20*log10(norm(inshift)/norm(plotthat)) < 35 + any(bufOut(:)>10),test_failed);
    fprintf(['DGTREAL_PROCESSOR OP %s gl:%3i, W:%3i, %s %s %s\n'],initstr,gl,W,flags.complexity,ltfatstatusstring(status),fail);

  

end
end






 