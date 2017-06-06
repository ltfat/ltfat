function out = comp_circshiftdim1(in,shift)

[L,W] = size(in);

p = mod(L-shift,L);
if p < 0
    p = p + L;
end

out = zeros(size(in));

rest = L-p;
out(1:rest,:) = in(p+1:end,:);
out(rest+1:end,:) = in(1:p,:);
