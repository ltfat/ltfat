function lname = makelibraryname(fname,complexity,is_complex)

suffix = '';

if strcmp(complexity,'double')
    suffix = [suffix,'_d'];
elseif strcmp(complexity,'single')
    suffix = [suffix,'_s'];
end

if is_complex
    suffix = [suffix,'c'];
end

lname = ['ltfat_',fname,suffix];