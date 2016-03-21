function fi = complex2interleaved(fc)

if ~isnumeric(fc)
    error('%s: Input must be numeric.',upper(mfilename))
end

fsize = size(fc);
fsize(1)=fsize(1)*2;

fi = zeros(fsize,class(fc));

fi(1:2:end,:,:) = real(fc);

if ~isreal(fc)
    fi(2:2:end,:,:) = imag(fc);
end