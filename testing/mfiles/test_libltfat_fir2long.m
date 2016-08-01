function test_libltfat_fir2long

L = 6;
Llong = 9;
ziorig = 1:L;
zi = postpad(ziorig,Llong);
zout = zeros(1,Llong);

fir2longtrue = fir2long(ziorig,Llong);

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);


ret = calllib('libltfat','ltfat_fir2long_d',ziPtr,L,Llong,zoutPtr)

zoutPtr.Value


ret = calllib('libltfat','ltfat_fir2long_d',ziPtr,L,Llong,ziPtr)

ziPtr.Value



errOutOfPlace = norm(fir2longtrue - zoutPtr.Value)
errInPlace = norm(fir2longtrue - ziPtr.Value)


%interleaved2complex(zoutPtr.Value)

