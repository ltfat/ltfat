function test_libltfat_pgauss


z = randn(100) + 1i*randn(100);

zi = complex2interleaved(z);
ziout = complex2interleaved(zeros(size(z)));

ziPtr = libpointer('doublePtr',zi);
zioutPtr = libpointer('doublePtr',ziout);

calllib('libltfat','col2diag_cd',ziPtr,100,zioutPtr);


res = interleaved2complex(zioutPtr.Value);

restrue = col2diag(z);

norm(res - restrue)