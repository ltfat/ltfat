function test_libltfat_gabframediag
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

a = 14;
gl = 34;
M = 64;
g = zeros(gl,1);
d = zeros(a,1);
gPtr = libpointer('doublePtr',g);
dPtr = libpointer('doublePtr',d);

calllib('libltfat','ltfat_firwin_d',LTFAT_FIRWIN.LTFAT_HANN,gl,gPtr);


calllib('libltfat','ltfat_gabframediag_d',gPtr,gl,a,M,a,dPtr);


d =gabframediag(gPtr.Value,a,M,lcm(a,M));
d(1:a)-dPtr.Value



