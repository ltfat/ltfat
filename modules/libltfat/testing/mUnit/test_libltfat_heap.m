function test_failed = test_libltfat_heap(varargin)
test_failed = 0;

fprintf(' ===============  %s ================ \n',upper(mfilename));

L = 1025;
f = rand(L,1);
fPtr = libpointer('doublePtr',f);
idxArr = randperm(L);

heap = calllib('libltfat','ltfat_heap_init_d', L, fPtr);

for w = idxArr
    calllib('libltfat','ltfat_heap_insert_d', heap, w-1);
end

sortedIdxArr = zeros(L,1);

for w = 1:L
    sortedIdxArr(w) = calllib('libltfat','ltfat_heap_delete_d', heap);
end

calllib('libltfat','ltfat_heap_done_d',heap);

res = norm(f(sortedIdxArr+1) - sort(f,'descend'));

plot([f,f(sortedIdxArr+1),sort(f,'descend')])

if res>0
    test_failed=1;
end





%plot(sort(f,'descend'))

