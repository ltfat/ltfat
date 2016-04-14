function test_libltfat_gabdual
[~,~,enuminfo]=libltfatprotofile;
LTFAT_FIRWIN = enuminfo.LTFAT_FIRWIN;

a = 15;
gl = 33;
M = 33;
g = zeros(gl,1);
gd = zeros(gl,1);
gPtr = libpointer('doublePtr',g);
gdPtr = libpointer('doublePtr',g);
calllib('libltfat','firwin_d',LTFAT_FIRWIN.HANN,gl,gPtr);

prd=calllib('libltfat','gabdual_painless_d',gPtr,gl,a,M,gdPtr);
prd

gdtrue = gabdual(gPtr.Value,a,M);

norm(long2fir(gdtrue,gl) - gdPtr.Value)



