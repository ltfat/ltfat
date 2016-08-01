function test_libltfat_firwin
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;


gl = 18;
g = zeros(gl,1);
gPtr = libpointer('doublePtr',g);

calllib('libltfat','ltfat_firwin_d',LTFAT_FIRWIN.LTFAT_HAMMING,gl,gPtr);

gtrue = firwin('hamming',gl);

res = norm(gtrue - gPtr.Value)

if res>1e-14
    gtrue - gPtr.Value
end