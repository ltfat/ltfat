function test_libltfat_circshift

L = 100000;
z = (1:L)+i*(1:L);
zi = complex2interleaved(z);
zout = zeros(size(zi));

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);

tic;
calllib('libltfat','circshift_cd',ziPtr,L,-50000,ziPtr);
toc;

tic;
calllib('libltfat','circshift_cd',ziPtr,10,50000,zoutPtr);
toc;


%interleaved2complex(zoutPtr.Value)

