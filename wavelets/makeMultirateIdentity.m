function out = makeMultirateIdentity( lo, hi, J)

out=cell(J+1,1);

out{end}=hi(:);
 h = lo(:);
 g= hi(:);
for j=2:J
  out{end+1-j}=conv(h,ups(g,2^(j-1),1));
  h= conv(h,ups(lo,2^(j-1),1));
end

out{1} = h;