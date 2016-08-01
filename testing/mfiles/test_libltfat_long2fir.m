function test_libltfat_long2fir

L = 7;
Llong = 10;
zi = 1:Llong;

zout = zeros(1,L);

fir2longtrue = long2fir(zi,L);

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);


ret = calllib('libltfat','ltfat_long2fir_d',ziPtr,Llong,L,zoutPtr);

zoutPtr.Value


ret = calllib('libltfat','ltfat_long2fir_d',ziPtr,Llong,L,ziPtr)

ziPtr.Value(1:L)



errOutOfPlace = norm(fir2longtrue - zoutPtr.Value)
errInPlace = norm(fir2longtrue - ziPtr.Value(1:L))


%interleaved2complex(zoutPtr.Value)

