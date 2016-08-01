function test_libltfat_gabdual
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

a = 16;
gl = 32;
M = 45;
g = zeros(gl,1);
gd = zeros(gl,1);
gPtr = libpointer('doublePtr',g);
gdPtr = libpointer('doublePtr',gd);
calllib('libltfat','ltfat_firwin_d',LTFAT_FIRWIN.LTFAT_HANN,gl,gPtr);
gdtrue = gabdual(gPtr.Value,a,M);

prd=calllib('libltfat','ltfat_gabdual_painless_d',gPtr,gl,a,M,gdPtr);
prd

norm(long2fir(gdtrue,gl) - gdPtr.Value)


prd=calllib('libltfat','ltfat_gabdual_painless_d',gPtr,gl,a,M,gPtr);
prd

norm(long2fir(gdtrue,gl) - gPtr.Value)

