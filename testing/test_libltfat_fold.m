function test_libltfat_fold

L = 10;
Lfold = 10;
zi = (1:L) ;
zout = zeros(1,Lfold);
shift = 102;

periods = ceil(L/Lfold);
fext = postpad(zi,periods*Lfold);
ffoldtrue = circshift(sum(reshape(fext, Lfold, periods),2)',[0,shift]);
ffold2ndtrue = sum(reshape(circshift(fext,[0,shift]), Lfold, periods),2)';

errTrues = norm(ffoldtrue - ffold2ndtrue)

ziPtr = libpointer('doublePtr',zi);
zoutPtr = libpointer('doublePtr',zout);


ret = calllib('libltfat','fold_array_d',ziPtr,L,shift,Lfold,zoutPtr)

zoutPtr.Value;


ret = calllib('libltfat','fold_array_d',ziPtr,L,shift,Lfold,ziPtr)

ziPtr.Value(1:min(Lfold,L));



errOutOfPlace = norm(ffoldtrue - zoutPtr.Value)
errInPlace = norm(ffoldtrue - ziPtr.Value(1:(min(Lfold,L))))


%interleaved2complex(zoutPtr.Value)

