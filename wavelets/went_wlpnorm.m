function E = went_wlpnorm(x,p)

%u = f./norm(f);
v = sort(abs(x(:)),'descend');

k = 1:length(v);
E = max((k.'.^(1/p)).*v);