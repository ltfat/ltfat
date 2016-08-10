function fc = interleaved2complex(fi)

if ~isreal(fi)
    error('%s: Input must be real.',upper(mfilename))
end

fsize = size(fi);
fsize(1)=fsize(1)/2;

if rem(fsize(1),1)~=0
     error('%s: Wrong dimension.',upper(mfilename))
end


fc=fi(1:2:end,:,:) + 1i*fi(2:2:end,:,:);

