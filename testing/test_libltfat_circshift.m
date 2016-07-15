function test_libltfat_circshift

L = 9;
z = (1:L)+i*(L:-1:1);
zi = complex2interleaved(z);
zout = zeros(size(zi));

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);

shift = -0;

calllib('libltfat','ltfat_circshift_dc',ziPtr,L,shift,zoutPtr);



calllib('libltfat','ltfat_circshift_dc',zoutPtr,L,-shift,zoutPtr);


norm(z - interleaved2complex(zoutPtr.Value))
