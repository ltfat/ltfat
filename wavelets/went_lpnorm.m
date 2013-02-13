function E = went_lpnorm(x,p)

assert(0<p&&p<2,'%s: Parameter p have to be in the range ]0,2[',upper(mfilename));
%u = x./norm(x);

%E = sum(abs(x(:)).^p)^(1/p);
E = sum(abs(x(:)).^p);