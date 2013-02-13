function E = went_shannon(x)
%u = x./norm(x);
u = x(x~=0).^2;
E = -sum(u.*log(u));